"""
Phase 9 Part E: threshold calibration of the candidates on DEVELOPMENT material only.

    python -m experiments.phase9_waveform.calibrate9 --run --workers 4    # raw statistics -> results/dev/calib_raw.jsonl
    python -m experiments.phase9_waveform.calibrate9 --calibrate          # -> results/dev/calibration9.json

Material (512 intervals; seeds in the DEV range 9000-9499 only):
  NULL conditions: the 17 PASS nulls, each at
     'input' (true beats, no recording noise: RR at 1/360 s resolution, generator beat labels) and
     'real'  (ECGSYN -> nstdb mix 12 dB, DEV signal 0 -> MIT-BIH quantization -> Pan-Tompkins beats ->
              Phase 7 V jitter; labels = annotations of the matched beats, 'X' for unmatched),
     seeds 9000-9099.
  NON-CHAOTIC conditions: Mackey-Glass tau = 14, 15, 16 and the six circle-map controls, each at
     'input' (ectopy none) and 'real' with ectopy none, S2_5, E1, E3_10; seeds 9100-9149.
  CHAOTIC (power only, never used for thresholds): Mackey-Glass chaotic regimes at the same five
     variants, seeds 9200-9219; Phase 5 positives P1, P2, G1, G2, G3 (input, all beats 'N'), seeds
     9200-9219.
Rules (predeclared, Phase 8 recipe):
  Z_DET (Z_DET_RR): the smallest z such that EVERY null condition has <= 2 % of its windows with
     zmax >= z (annotation mask; RR-rule mask), plus 0.5.
  G_MIN (G_MIN_RR): the smallest g such that EVERY null and non-chaotic condition has <= 2 % of its
     windows with (zmax >= Z_DET and G >= g), plus 0.05.
  F_MIN: the same for the FSLE gate value (fsle10), plus 0.05.
A window that is not analysable never counts as a detection.
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

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "results" / "dev"
N = 512


def window(t):
    from experiments.phase8_cardiac import nulls as NL8
    from experiments.phase9_waveform import library as L
    if t["source"] == "positive":
        x = NL8.generate(t["name"], t["seed"], N)
        return np.asarray(x, float), np.array(["N"] * (len(x) + 1)), {}
    ect = "none" if t["source"] == "null" else t["ectopy"]
    if t["variant"] == "input":
        w = L.build(t["source"], t["name"], ect, N, t["seed"], "none", "DEV", "true", False)
    else:
        w = L.build(t["source"], t["name"], ect, N, t["seed"], "mix12", "DEV", "detected", True)
    if "error" in w:
        raise RuntimeError(w["error"])
    info = {k: w[k] for k in ("n_missed", "n_extra") if k in w}
    return np.asarray(w["rr"], float), np.asarray(w["beat_types"]), info


def one(t):
    from experiments.phase9_waveform import candidates9 as C
    t0 = time.time()
    out = dict(t)
    try:
        x, types, info = window(t)
        out.update(info)
        out["n_nonN"] = int(np.sum(types != "N"))
        out.update(C.raw_stats(x, types))
    except Exception as exc:                               # noqa: BLE001
        out["error"] = f"{type(exc).__name__}: {exc}"
    out["t"] = time.time() - t0
    return out


def tasks():
    from experiments.phase8_cardiac import nulls as NL8
    from experiments.phase9_waveform import devmaps as DM
    from experiments.phase9_waveform import families as F
    T = []
    for c in NL8.PASS_NULLS:
        for s in range(9000, 9100):
            for v in ("input", "real"):
                T.append({"group": "null", "source": "null", "name": c, "variant": v, "ectopy": "native", "seed": s})
    regs = [r for r in F.load_regimes() if r["family"] == "mackey_glass"]
    non = [("model", r["name"]) for r in regs if r["label"] == "NON-CHAOTIC"] + [("devmap", r["name"]) for r in DM.REGIMES]
    cha = [("model", r["name"]) for r in regs if r["label"] == "CHAOTIC"]
    for grp, lst, seeds in (("noncha", non, range(9100, 9150)), ("chaos", cha, range(9200, 9220))):
        for src, nm in lst:
            for s in seeds:
                T.append({"group": grp, "source": src, "name": nm, "variant": "input", "ectopy": "none", "seed": s})
                for e in ("none", "S2_5", "E1", "E3_10"):
                    T.append({"group": grp, "source": src, "name": nm, "variant": "real", "ectopy": e, "seed": s})
    for c in ("P1_henon_rr", "P2_logistic_rr", "G1_lorenz_maxima", "G2_rossler_flow", "G3_mackey_glass"):
        for s in range(9200, 9220):
            T.append({"group": "chaos", "source": "positive", "name": c, "variant": "input", "ectopy": "none", "seed": s})
    for t in T:
        t["tid"] = f"{t['group']}|{t['name']}|{t['variant']}|{t['ectopy']}|{t['seed']}"
    return T


def cond(r):
    return f"{r['name']}|{r['variant']}|{r['ectopy']}"


def smallest_threshold(groups, frac=0.02):
    """groups: list of arrays of statistic values (nan = cannot exceed).  Smallest thr with every
    group's fraction of values >= thr at most frac."""
    thr = -np.inf
    for v in groups:
        v = np.asarray(v, float)
        n = len(v)
        allowed = int(np.floor(frac * n))
        fin = np.sort(v[np.isfinite(v)])[::-1]
        if len(fin) > allowed:
            # need thr > fin[allowed]
            thr = max(thr, np.nextafter(fin[allowed], np.inf))
    return float(thr)


def calibrate():
    import pandas as pd
    R = [json.loads(line) for line in open(OUT / "calib_raw.jsonl")]
    R = [r for r in R if "error" not in r]
    res = {"n_windows": len(R)}
    for key, zname, gname in (("ann", "Z_DET", "G_MIN"), ("rr", "Z_DET_RR", "G_MIN_RR")):
        nulls = {}
        for r in R:
            if r["group"] == "null":
                a = r[key]
                nulls.setdefault(cond(r), []).append(a["zmax"] if a.get("analysable") else np.nan)
        z = smallest_threshold(list(nulls.values())) + 0.5
        res[zname] = z
        res[zname + "_binding"] = max(nulls, key=lambda c: np.nanpercentile(np.r_[nulls[c], -1e9], 98))
        gs = {}
        fs = {}
        for r in R:
            if r["group"] in ("null", "noncha"):
                a = r[key]
                ok = a.get("analysable") and a["zmax"] >= z
                gs.setdefault(cond(r), []).append(a["G"] if ok else np.nan)
                if key == "ann":
                    fs.setdefault(cond(r), []).append(r["fsle10"] if ok else np.nan)
        res[gname] = smallest_threshold(list(gs.values())) + 0.05
        if key == "ann":
            res["F_MIN"] = smallest_threshold(list(fs.values())) + 0.05
    # development power / specificity table with the calibrated thresholds
    rows = []
    for r in R:
        a, b = r["ann"], r["rr"]
        det_a = bool(a.get("analysable") and a["zmax"] >= res["Z_DET"])
        det_b = bool(b.get("analysable") and b["zmax"] >= res["Z_DET_RR"])
        rows.append({"group": r["group"], "cond": cond(r),
                     "k2": det_a,
                     "k3": det_a and a["G"] >= res["G_MIN"],
                     "k4": det_b and b["G"] >= res["G_MIN_RR"],
                     "k5": det_a and np.isfinite(r["fsle10"]) and r["fsle10"] >= res["F_MIN"]})
    D = pd.DataFrame(rows)
    tab = D.groupby(["group", "cond"])[["k2", "k3", "k4", "k5"]].agg(["sum", "count"])
    res["by_condition"] = {f"{g}|{c}": {k: f"{int(tab.loc[(g, c), (k, 'sum')])}/{int(tab.loc[(g, c), (k, 'count')])}"
                                         for k in ("k2", "k3", "k4", "k5")} for g, c in tab.index}
    json.dump(res, open(OUT / "calibration9.json", "w"), indent=1, default=float)
    print(json.dumps({k: v for k, v in res.items() if k != "by_condition"}, indent=1, default=float))
    return res


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--calibrate", action="store_true")
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args(argv)
    if a.run:
        OUT.mkdir(parents=True, exist_ok=True)
        path = OUT / "calib_raw.jsonl"
        done = {json.loads(line)["tid"] for line in open(path)} if path.exists() else set()
        T = [t for t in tasks() if t["tid"] not in done]
        np.random.default_rng(2).shuffle(T)
        print(len(T), "to run", flush=True)
        with open(path, "a") as fh, Pool(a.workers, maxtasksperchild=50) as pool:
            for i, r in enumerate(pool.imap_unordered(one, T, chunksize=2), 1):
                fh.write(json.dumps(r, default=float) + "\n")
                fh.flush()
                if i % 200 == 0:
                    print(i, "done", flush=True)
    if a.calibrate:
        calibrate()


if __name__ == "__main__":
    main()
