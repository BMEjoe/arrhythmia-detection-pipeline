"""
Phase 10 runner (resumable; one JSON line per task; a task is skipped if its id is already present).

    python -m experiments.phase10_final.run10 --phase dev  --part real   --workers 4
    python -m experiments.phase10_final.run10 --phase conf --part all    --workers 4   # guarded
    python -m experiments.phase10_final.run10 --phase conf --download                  # guarded

Parts
  real    per window (512 intervals): K1 raw, K1 edited, K3, K4, K3 with RR-rule labels, titration
          raw / masked / edited; window metadata (ectopy burden, clock hour, HRV)
  robust  windows j in ROBUST_J: K1 and K3 with 499 surrogates; K1 and K3 at 256 and 1024 intervals
  spike   Question 1: one base window per subject x family x parameter -> additive spike-ins at every
          f in F_GRID and the replacement design; K1 and K3 on each
  seg     Question 2c/2d: 12-min segments, titration raw / masked / edited / Wu-NN
  synth   Question 2e: titration raw / masked / edited on Phase 5-6 nulls and non-chaotic coupled vdP
          with ectopy, 512 and 800 intervals

The CONFIRMATORY phase (and the confirmatory download) refuses to start unless PREREGISTRATION.md is
committed, unmodified and pushed (the Phase 3 guard, experiments.phase3_lle.run_phase3).
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import time
import traceback
from multiprocessing import Pool

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402

from experiments.phase10_final import data10 as D  # noqa: E402
from experiments.phase10_final import harness10 as H  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
PREREG = HERE / "PREREGISTRATION.md"
METHODS_LISTED = ("K1", "K3", "K4", "K3RR", "TIT")
SYNTH_SEEDS = {"dev": range(10000, 10005), "conf": range(10100, 10200)}
SYNTH_LENGTHS = (512, 800)
DEV_SPIKE_F = (0.1, 0.3, 0.5, 0.9)     # development used a reduced f grid (runtime)


def guard():
    from experiments.phase3_lle import run_phase3 as R3
    R3.PREREG = PREREG
    R3.check_preregistration(list(METHODS_LISTED))


def dbs(phase):
    if phase == "dev":
        return D.DEV_DIR, ("nsrdb", "chfdb", "mitdb")
    return D.CONF_DIR, D.CONF_DBS


GROUP = {"nsrdb": "NSR", "nsr2db": "NSR", "chfdb": "CHF", "chf2db": "CHF", "mitdb": "MIT"}


def subject_of(db, rec):
    if db == "mitdb" and rec == "202":
        return "mitdb:201"
    return f"{db}:{rec}"


def n_per(db):
    return 3 if db == "mitdb" else H.N_PER


_CACHE = {}


def series(root, db, rec):
    k = (str(root), db, rec)
    if k not in _CACHE:
        _CACHE.clear()
        t, lab, sym, fs = D.beat_series(db, rec, root)
        _CACHE[k] = (t, lab, sym, fs, D.header_info(db, rec, root))
    return _CACHE[k]


def excluded(phase):
    """Subjects excluded after the overlap check (confirmatory only; written by overlap10.py)."""
    p = HERE / "results" / "conf" / "exclusions.json"
    if phase == "conf" and p.exists():
        return set(json.load(open(p))["excluded_subjects"])
    return set()


# --------------------------------------------------------------------------------------- tasks
def tasks(phase, part, ids=None):
    root, names = dbs(phase)
    ex = excluded(phase)
    T = []
    for db in names:
        recs = D.records(db, root)
        for si, rec in enumerate(recs):
            if subject_of(db, rec) in ex:
                continue
            if part in ("real", "robust", "spike"):
                t = series(root, db, rec)[0]
                ws = H.window_starts(t, n_per=n_per(db))
                if part == "real":
                    T += [{"id": f"real|{db}|{rec}|{j}", "part": "real", "db": db, "rec": rec, "j": j, "i0": i0}
                          for j, i0 in ws.items()]
                elif part == "robust" and db != "mitdb":
                    T += [{"id": f"robust|{db}|{rec}|{j}", "part": "robust", "db": db, "rec": rec, "j": j,
                           "i0": ws[j]} for j in H.ROBUST_J if j in ws]
                elif part == "spike" and db != "mitdb" and ws:
                    jb = base_window(si, ws)
                    for fam, lst in H.FAMILIES.items():
                        for p in range(len(lst)):
                            T.append({"id": f"spike|{db}|{rec}|{jb}|{fam}|{p}", "part": "spike", "db": db,
                                      "rec": rec, "j": jb, "i0": ws[jb], "family": fam, "pidx": p, "si": si,
                                      "f_grid": list(DEV_SPIKE_F if phase == "dev" else H.F_GRID)})
            elif part == "seg" and db != "mitdb":
                T.append({"id": f"seg|{db}|{rec}", "part": "seg", "db": db, "rec": rec})
    if part == "synth":
        from experiments.phase8_cardiac import nulls as NL8
        conds = [("null", c) for c in NL8.PASS_NULLS] + [("vdp", (r, e)) for r in H.VDP_NONCHAOTIC for e in H.ECT]
        for kind, c in conds:
            for n in SYNTH_LENGTHS:
                for s in SYNTH_SEEDS[phase]:
                    cid = c if kind == "null" else f"vdp{c[0]}_{c[1]}"
                    T.append({"id": f"synth|{cid}|{n}|{s}", "part": "synth", "kind": kind, "cond": c, "n": n,
                              "seed": s})
    if ids is not None:
        T = [t for t in T if t["id"] in ids]
    return T


def base_window(si, ws):
    """Question 1 base window of subject si (index in RECORDS): j = 1 + (si mod 12), else the next
    available j cyclically (spreads base windows over the time of day)."""
    for k in range(H.N_PER):
        j = 1 + (si + k) % H.N_PER
        if j in ws:
            return j
    return None


# --------------------------------------------------------------------------------------- execution
def _window(root, task, L):
    t, lab, sym, fs, info = series(root, task["db"], task["rec"])
    i0 = task["i0"]
    if i0 + L >= len(t):
        return None
    tt = t[i0:i0 + L + 1]
    return {"rr": np.diff(tt), "lab": lab[i0:i0 + L + 1], "sym": sym[i0:i0 + L + 1], "fs": fs, "t0": float(tt[0]),
            "info": info, "clean": bool(np.all((np.diff(tt) >= H.LO) & (np.diff(tt) <= H.HI)))}


def run_real(root, task):
    w = _window(root, task, H.L_PRIMARY)
    rr, lab, sym = w["rr"], w["lab"], w["sym"]
    m = H.k3_mask(lab)
    out = {"subject": subject_of(task["db"], task["rec"]), "group": GROUP[task["db"]], "j": task["j"],
           "t_start": w["t0"], "clock_hour": H.clock_hour(w["info"]["base_time"], w["t0"] - 0.0),
           "burden": int(np.sum(~np.isin(sym, list(D.NORMAL)))), "n_masked": int(m.sum()),
           "nn_fraction": float(1 - m.mean()), "mean_rr": float(rr.mean()), "sdnn": float(rr.std()),
           "rmssd": float(np.sqrt(np.mean(np.diff(rr) ** 2)))}
    out["K1"] = H.k1(rr)
    xe, nnf = H.edit(rr, lab)
    if not m.any():
        out["K1_edited"] = {"same_as_raw": True, **out["K1"]}
    elif xe is not None and nnf >= H.NN_EDIT_MIN:
        out["K1_edited"] = {"same_as_raw": False, **H.k1(xe)}
    else:
        out["K1_edited"] = {"eligible": False, "nn_fraction": nnf}
    out["K3"] = H.k3(rr, lab)
    out["K4"] = H.k4(rr)
    out["K3RR"] = H.k3(rr, H.rr_rule_labels(rr))
    for arm in ("raw", "masked", "edited"):
        out[f"TIT_{arm}"] = H.titration(rr, lab, arm)
    return out


def run_robust(root, task):
    out = {"subject": subject_of(task["db"], task["rec"]), "group": GROUP[task["db"]], "j": task["j"]}
    w = _window(root, task, H.L_PRIMARY)
    out["K1_499"] = H.k1(w["rr"], n_sur=499)
    out["K3_499"] = H.k3(w["rr"], w["lab"], n_sur=499)
    t = series(root, task["db"], task["rec"])[0]
    for L in (256, 1024):
        sL = H.window_starts(t, n_per=n_per(task["db"]), L=L).get(task["j"])
        wl = None if sL is None else _window(root, {**task, "i0": sL}, L)
        if wl is None or not wl["clean"]:
            out[f"L{L}"] = {"available": False}
            continue
        out[f"L{L}"] = {"available": True, "burden": int(np.sum(~np.isin(wl["sym"], list(D.NORMAL)))),
                        "K1": H.k1(wl["rr"]), "K3": H.k3(wl["rr"], wl["lab"])}
    return out


def run_spike(root, task, phase):
    w = _window(root, task, H.L_PRIMARY)
    rr, lab, sym, fs = w["rr"], w["lab"], w["sym"], w["fs"]
    fam, p = task["family"], task["pidx"]
    dbc = {"nsrdb": 1, "chfdb": 2, "nsr2db": 3, "chf2db": 4}[task["db"]]
    key = [dbc, task["si"], task["j"]]
    try:
        c = H.chaos_component(fam, p, len(rr), key)
    except Exception as exc:                                    # noqa: BLE001
        return {"error": f"chaos component: {type(exc).__name__}: {exc}"}
    out = {"subject": subject_of(task["db"], task["rec"]), "group": GROUP[task["db"]], "j": task["j"],
           "family": fam, "param": H.FAMILIES[fam][p][0], "lam": H.FAMILIES[fam][p][2],
           "lam_unit": H.FAMILIES[fam][p][3], "fs": fs, "burden": int(np.sum(~np.isin(sym, list(D.NORMAL)))),
           "designs": {}}
    for f in task["f_grid"]:
        q, meta = H.spike_additive(rr, lab, c, f, fs)
        out["designs"][f"add_{f}"] = {"f": f, **meta, "K1": H.k1(q), "K3": H.k3(q, lab)}
    rng = np.random.default_rng(np.random.SeedSequence([H.ENTROPY, 2, *key, H.FAM_CODE[fam], p]))
    q, meta = H.spike_replacement(rr, lab, sym, c, fs, rng)
    out["designs"]["replace"] = {**meta, "K1": H.k1(q), "K3": H.k3(q, lab)}
    return out


def run_seg(root, task):
    t, lab, sym, fs, info = series(root, task["db"], task["rec"])
    out = {"subject": subject_of(task["db"], task["rec"]), "group": GROUP[task["db"]], "segments": []}
    for s in H.segments(t, lab):
        rec = {"seg": s["seg"], "n": s["n"], "n_removed": s["n_removed"], "cover": s["cover"], "gap": s["gap"],
               "t_start": float(t[s["i0"]]), "clock_hour": H.clock_hour(info["base_time"], float(t[s["i0"]]))}
        if not s["gap"]:
            rr_all = np.diff(t[s["i0"]:s["i1"] + 1])
            lb = lab[s["i0"]:s["i1"] + 1]
            sy = sym[s["i0"]:s["i1"] + 1]
            keep = (rr_all >= H.LO) & (rr_all <= H.HI)
            rr, m = rr_all[keep], H.k3_mask(lb)[keep]
            rec["burden"] = int(np.sum(~np.isin(sy, list(D.NORMAL))))
            rec["n_masked"] = int(m.sum())
            for arm in ("raw", "masked", "edited"):
                rec[f"TIT_{arm}"] = H.titration(rr, None, arm, imask=m)
            nn = H.wu_nn(rr, m)
            rec["TIT_wu"] = (H.titration(nn, None, "raw") if len(nn) >= H.SEG_MIN_N
                             else {"eligible": False, "n": int(len(nn)), "positive": False, "NL": 0.0})
        out["segments"].append(rec)
    return out


def run_synth(task):
    rr, types = H.synthetic_window(task["kind"], tuple(task["cond"]) if task["kind"] == "vdp" else task["cond"],
                                   task["seed"], task["n"])
    out = {"kind": task["kind"], "cond": task["cond"], "n": task["n"], "seed": task["seed"],
           "n_ectopic": int(np.sum(np.asarray(types) != "N"))}
    for arm in ("raw", "masked", "edited"):
        out[f"TIT_{arm}"] = H.titration(rr, types, arm)
    out["TIT_raw_q128"] = H.titration(H.quantize_times(0.0, rr, 128.0), None, "raw")
    return out


def execute(args):
    phase, task = args
    root = dbs(phase)[0]
    t0 = time.time()
    out = {"id": task["id"], "part": task["part"]}
    try:
        if task["part"] == "real":
            out.update(run_real(root, task))
        elif task["part"] == "robust":
            out.update(run_robust(root, task))
        elif task["part"] == "spike":
            out.update(run_spike(root, task, phase))
        elif task["part"] == "seg":
            out.update(run_seg(root, task))
        elif task["part"] == "synth":
            out.update(run_synth(task))
    except Exception as exc:                                    # noqa: BLE001
        out["error"] = f"{type(exc).__name__}: {exc}"
        out["traceback"] = traceback.format_exc(limit=4)
    out["t_total"] = time.time() - t0
    return out


def _done(path):
    ids = set()
    if path.exists():
        for line in open(path):
            try:
                ids.add(json.loads(line)["id"])
            except Exception:                                   # noqa: BLE001
                pass
    return ids


def _default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        v = float(o)
        return None if not np.isfinite(v) else v
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", choices=["dev", "conf"], required=True)
    ap.add_argument("--part", default="all")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--download", action="store_true")
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args(argv)
    if a.phase == "conf":
        guard()
        if a.download:
            for db in D.CONF_DBS:
                recs = D.download(db, D.CONF_DIR)
                print(db, len(recs), "records; SHA-256 OK")
            return
    parts = ("real", "seg", "synth", "spike", "robust") if a.part == "all" else tuple(a.part.split(","))
    out_dir = HERE / "results" / a.phase
    out_dir.mkdir(parents=True, exist_ok=True)
    for part in parts:
        path = out_dir / f"{part}.jsonl"
        done = _done(path)
        T = [t for t in tasks(a.phase, part) if t["id"] not in done]
        if a.limit:
            T = T[:a.limit]
        print(f"{part}: {len(T)} tasks to run ({len(done)} done)", flush=True)
        if not T:
            continue
        t0 = time.time()
        with Pool(a.workers) as pool, open(path, "a") as fh:
            for i, r in enumerate(pool.imap_unordered(execute, [(a.phase, t) for t in T], chunksize=1)):
                fh.write(json.dumps(r, default=_default) + "\n")
                fh.flush()
                if (i + 1) % 25 == 0 or i + 1 == len(T):
                    el = time.time() - t0
                    print(f"  {part} {i + 1}/{len(T)} {el / 60:.1f} min; eta {el / (i + 1) * (len(T) - i - 1) / 60:.1f} min",
                          flush=True)


if __name__ == "__main__":
    main()
