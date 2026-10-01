"""
Phase 7 EXPLORATORY spike-in positive controls (chosen after the preregistered
results were seen).  For every primary raw-RR window (305), the frozen DETECTOR is
run on synthetic chaotic RR built to that window's own scale and ectopy:

  (a) Henon x (a = 1.4, b = 0.3) and (b) logistic r = 4 (Phase 2E generators,
      experiments.phase2e.systems), standardized and rescaled to the real window's
      mean and SD, then quantized to 1/360 s (Phase 5 quantize).  If the window's SD
      would push the minimum interval below CFG.rr_min_seconds (0.25 s; windows with
      long pauses), the SD is reduced to the largest value keeping min >= 0.25 s
      (counted: sd_capped).
  (c) the same chaotic series (paired: same realization) with the real window's
      abnormal beats inserted at their real positions with the Phase 6 Part B
      ectopy model, applied to the unquantized series before quantization:
        - a maximal run of k consecutive abnormal ending beats starting at interval j
          is one group; its type is the first beat's symbol: supraventricular
          (A a J S) -> atrial model, otherwise (V E F f Q) -> ventricular model;
        - ventricular (full compensatory pause, sinus node not reset):
          premature interval j = c * x_j, c ~ U(0.6, 0.8); k = 2: second premature
          c' * x_j; k >= 3 (run): later premature intervals an absolute cycle
          U(0.40, 0.55) s; compensatory interval j+k = sum(x_j..x_j+k) - sum(premature)
          (floored at 0.25 s, counted);
        - atrial (sinus node reset, Phase 6 E5): premature c * x_j (later premature
          beats of a run c' * x_{j+i}), post-ectopic interval j+k = d * x_{j+k},
          d ~ U(1.0, 1.1);
        - a group ending at the window's last interval has no compensatory interval.
  (d) (c) plus R-peak timing jitter at every ventricular-type ectopic beat (V E F):
      an offset (samples) drawn from the empirical distribution of
      (detected - annotated) sample offsets of all V/E/F beats matched at 150 ms in
      the 44 non-paced records (QC); the interval ending at the beat gets +offset,
      the next interval -offset (stays on the 1/360 s grid).

    python -m experiments.phase7_mitbih.spike_in [--workers 4]      # resumable
    python -m experiments.phase7_mitbih.spike_in --summary
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import time
from multiprocessing import Pool

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import final_pipeline as fp  # noqa: E402
from experiments.phase2e import systems as S2  # noqa: E402
from experiments.phase5_rr import systems as S5  # noqa: E402
from experiments.phase6_robust import config as C6  # noqa: E402
from experiments.phase7_mitbih import analysis as A  # noqa: E402
from experiments.phase7_mitbih import data as D  # noqa: E402
from experiments.phase7_mitbih import detector as DT  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "results" / "spike_in"
ENTROPY = 20261002
SYSTEMS = {"henon": (1, S2.henon), "logistic": (2, S2.logistic)}
VARIANTS = ("a_clean", "c_ectopy", "d_ectopy_jitter")
S_TYPE = frozenset("A a J S".split())
V_TYPE = frozenset("V E F".split())
RR_FLOOR = fp.CFG.rr_min_seconds
FS = 360.0


def v_offsets():
    """Empirical (detected - annotated) offsets, in samples, of V/E/F beats matched at 150 ms."""
    pk = np.load(HERE / "results" / "qc" / "detected_peaks.npz")
    offs = []
    for r in D.RECORDS:
        if r in D.PACED_RECORDS:
            continue
        _, fs, s, y = D.load(r)
        b = np.isin(y, list(D.BEAT_SYMBOLS))
        rs, ry = s[b], y[b]
        m = D.greedy_match(pk[r], rs, fs, D.TOL_PRIMARY_MS)
        for i in np.flatnonzero(m >= 0):
            if ry[m[i]] in V_TYPE:
                offs.append(int(pk[r][i] - rs[m[i]]))
    return np.asarray(offs, dtype=int)


def window_specs():
    """Per primary window: mean, SD, per-interval abnormal mask and symbols."""
    raw = A.load(HERE / "results" / "run", "raw")
    specs = []
    for r in D.RECORDS:
        g = raw[raw.record == r]
        if not len(g):
            continue
        ser = D.record_series(r, cache_dir=os.path.join(D.DATA_DIR, "_phase7_peaks"))
        b = ser["is_beat"]
        rs, ry = ser["ann_s"][b], ser["ann_y"][b]
        m = D.greedy_match(ser["beats"], rs, ser["fs"], D.TOL_PRIMARY_MS)
        sym = np.where(m >= 0, ry[np.maximum(m, 0)], "")
        for _, st in g.iterrows():
            s0 = int(st.start_rr)
            x = ser["rr"][s0:s0 + 256]
            lab = ser["labels"][s0 + 1:s0 + 257]
            specs.append({"record": r, "subject": st.subject, "window": int(st.window), "label": int(st.label),
                          "beat_type": st.beat_type, "abnormal_fraction": float(st.abnormal_fraction),
                          "mean": float(np.mean(x)), "sd": float(np.std(x)),
                          "abnormal": (lab == 1).tolist(), "symbols": sym[s0 + 1:s0 + 257].tolist()})
    return specs


def groups(abnormal, symbols):
    out, j, n = [], 0, len(abnormal)
    while j < n:
        if abnormal[j]:
            k = 1
            while j + k < n and abnormal[j + k]:
                k += 1
            out.append((j, k, "atrial" if symbols[j] in S_TYPE else "ventricular"))
            j += k
        else:
            j += 1
    return out


def insert_ectopy(x, grp, rng):
    x = x.copy()
    n_floor = 0
    for j, k, kind in grp:
        base = x[j:j + k + 1].copy()
        c = rng.uniform(*C6.PREMATURE_COUPLING_RANGE)
        prem = [c * base[0]]
        for i in range(1, k):
            if kind == "ventricular" and k >= 3:
                prem.append(rng.uniform(*C6.RUN_CYCLE_RANGE_S))           # run: absolute cycle
            elif kind == "ventricular":
                prem.append(rng.uniform(*C6.PREMATURE_COUPLING_RANGE) * base[0])   # couplet: c' * x_j
            else:
                prem.append(rng.uniform(*C6.PREMATURE_COUPLING_RANGE) * base[i])   # atrial run: c' * x_{j+i}
        x[j:j + k] = prem
        if j + k < len(x):
            if kind == "ventricular":
                comp = float(base.sum() - np.sum(prem))
                if comp < RR_FLOOR:
                    comp, n_floor = RR_FLOOR, n_floor + 1
                x[j + k] = comp
            else:
                x[j + k] = rng.uniform(*C6.ATRIAL_POST_FACTOR_RANGE) * base[k]
    return x, n_floor


def build(spec, system, offsets):
    code, gen = SYSTEMS[system]
    key = [ENTROPY, code, int(spec["record"]), spec["window"]]
    z = S5._standardize(gen(256, np.random.default_rng(np.random.SeedSequence(key))))
    sd = spec["sd"]
    capped = False
    if spec["mean"] + sd * z.min() < RR_FLOOR:
        sd, capped = (spec["mean"] - RR_FLOOR) / (-z.min()), True
    x = spec["mean"] + sd * z
    grp = groups(spec["abnormal"], spec["symbols"])
    xe, n_floor = insert_ectopy(x, grp, np.random.default_rng(np.random.SeedSequence(key + [10])))
    a, c = S5.quantize(x), S5.quantize(xe)
    d = c.copy()
    jr = np.random.default_rng(np.random.SeedSequence(key + [20]))
    n_jit = 0
    for j, sym in enumerate(spec["symbols"]):
        if spec["abnormal"][j] and sym in V_TYPE:
            off = int(jr.choice(offsets)) / FS
            d[j] += off
            if j + 1 < len(d):
                d[j + 1] -= off
            n_jit += 1
    meta = {"sd_used": float(sd), "sd_capped": capped, "n_groups": len(grp), "n_comp_floored": n_floor,
            "n_jittered": n_jit, "min_rr_d": float(d.min())}
    return {"a_clean": a, "c_ectopy": c, "d_ectopy_jitter": d}, meta


def execute(task):
    tid, spec, system, variant, rr, meta = task
    rec = DT.evaluate(rr)
    keep = ("and_detected", "lle_detected", "upo_detected", "lle_p", "lle_z", "upo_score", "upo_status",
            "error", "detrend_applied", "runtime_s")
    return {"task_id": tid, "system": system, "variant": variant,
            **{k: spec[k] for k in ("record", "subject", "window", "label", "beat_type", "abnormal_fraction")},
            **meta, **{k: rec.get(k) for k in keep}}


def run(workers):
    OUT.mkdir(parents=True, exist_ok=True)
    offsets = v_offsets()
    np.save(OUT / "v_offsets_samples.npy", offsets)
    specs = window_specs()
    tasks = []
    for spec in specs:
        for system in SYSTEMS:
            series, meta = build(spec, system, offsets)
            for v in VARIANTS:
                tid = f"{system}:{v}:{spec['record']}:{spec['window']}"
                tasks.append((tid, spec, system, v, series[v], meta))
    path = OUT / "spike_in.jsonl"
    done = set()
    if path.exists():
        for line in open(path):
            try:
                done.add(json.loads(line)["task_id"])
            except (json.JSONDecodeError, KeyError):
                pass
    todo = [t for t in tasks if t[0] not in done]
    print(f"{len(specs)} windows, {len(tasks)} tasks, {len(todo)} to run; "
          f"{len(offsets)} V offsets, median {np.median(offsets):.0f} samples", flush=True)
    t0 = time.time()
    with open(path, "a") as fh, Pool(workers, maxtasksperchild=50) as pool:
        for i, res in enumerate(pool.imap_unordered(execute, todo, chunksize=1), 1):
            fh.write(json.dumps(res, default=lambda o: o.item() if hasattr(o, "item") else str(o)) + "\n")
            fh.flush()
            if i % 50 == 0 or i == len(todo):
                el = time.time() - t0
                print(f"{i}/{len(todo)} {el:.0f}s eta {el / i * (len(todo) - i):.0f}s", flush=True)


def summary():
    df = pd.DataFrame([json.loads(l) for l in open(OUT / "spike_in.jsonl")]).drop_duplicates("task_id")
    for c in ("and_detected", "lle_detected", "upo_detected"):
        df[c] = df[c].fillna(False).astype(bool)
    off = np.load(OUT / "v_offsets_samples.npy")
    md = ["# Phase 7 EXPLORATORY spike-in positive controls (generated by spike_in.py --summary)", "",
          f"Tasks {len(df)} (errors {int(df.error.notna().sum())}). V/E/F timing offsets used for (d): "
          f"n = {len(off)}, median {np.median(off) / FS * 1000:.1f} ms, "
          f"{100 * np.mean(off / FS > 0.075):.1f} % > +75 ms, {100 * np.mean(off / FS < -0.075):.1f} % < -75 ms, "
          f"IQR [{np.percentile(off, 25) / FS * 1000:.1f}, {np.percentile(off, 75) / FS * 1000:.1f}] ms.", ""]
    first = df[df.variant == "a_clean"]
    md += [f"Windows with SD capped (pauses): {int(first.sd_capped.sum() / 2)} of {int(len(first) / 2)}; "
           f"compensatory intervals floored: {int(first.n_comp_floored.sum())} (both systems).", ""]
    rows = []
    for (system, variant, label), g in df.groupby(["system", "variant", "label"]):
        row = {"system": system, "variant": variant, "label": "abnormal" if label else "normal", "n": len(g)}
        for c, name in (("and_detected", "AND"), ("lle_detected", "LLE alone"), ("upo_detected", "UPO alone")):
            k = int(g[c].sum())
            lo, hi = A.wilson(k, len(g))
            row[name] = f"{k}/{len(g)} ({100 * k / len(g):.1f} %; {100 * lo:.1f}–{100 * hi:.1f})"
        row["errors"] = int(g.error.notna().sum())
        rows.append(row)
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "rates.csv", index=False)
    md += ["## Detection rates (Wilson 95 %)", "", tab.to_markdown(index=False), ""]
    ab = df[df.label == 1].copy()
    ab["abn_bin"] = pd.cut(ab.abnormal_fraction, [0.0999, 0.2, 0.3, 0.5, 1.0])
    rows = []
    for (system, variant, b), g in ab.groupby(["system", "variant", "abn_bin"], observed=True):
        rows.append({"system": system, "variant": variant, "abnormal fraction": str(b), "n": len(g),
                     "AND": f"{int(g.and_detected.sum())}/{len(g)}", "UPO": f"{int(g.upo_detected.sum())}/{len(g)}",
                     "LLE": f"{int(g.lle_detected.sum())}/{len(g)}"})
    t2 = pd.DataFrame(rows)
    t2.to_csv(OUT / "abnormal_by_fraction.csv", index=False)
    md += ["## Abnormal windows by abnormal fraction", "", t2.to_markdown(index=False), ""]
    rows = []
    for (system, variant, bt), g in ab.groupby(["system", "variant", "beat_type"]):
        rows.append({"system": system, "variant": variant, "beat type": bt, "n": len(g),
                     "AND": f"{int(g.and_detected.sum())}/{len(g)}", "UPO": f"{int(g.upo_detected.sum())}/{len(g)}"})
    t3 = pd.DataFrame(rows)
    md += ["## Abnormal windows by beat type", "", t3.to_markdown(index=False), ""]
    nm = df[df.label == 0].copy()
    nm["has_ectopy"] = nm.n_groups > 0
    rows = []
    for (system, variant, h), g in nm.groupby(["system", "variant", "has_ectopy"]):
        rows.append({"system": system, "variant": variant, "normal window with >= 1 ectopic beat": bool(h),
                     "n": len(g), "AND": f"{int(g.and_detected.sum())}/{len(g)}",
                     "UPO": f"{int(g.upo_detected.sum())}/{len(g)}"})
    md += ["## Normal windows with / without any inserted ectopic beat", "", pd.DataFrame(rows).to_markdown(index=False), ""]
    open(OUT / "spike_in.md", "w").write("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--summary", action="store_true")
    a = ap.parse_args()
    summary() if a.summary else run(a.workers)
