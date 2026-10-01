"""
Phase 8 Part A: ground-truth largest Lyapunov exponents from the model equations.

    python -m experiments.phase8_cardiac.ground_truth --scan        # predeclared grids, 3 ICs each
    python -m experiments.phase8_cardiac.ground_truth --select      # label + pick regimes

Labelling rule (predeclared here, applied mechanically to the scan):
  per grid point, K = 3 independent initial conditions / seeds; for each the LE estimate
  lam_k with batch-means SE_k and a convergence curve.
  combined lam = mean_k lam_k; SE = sqrt(mean_k SE_k^2 / K + var_k(lam_k) / K)
  CHAOTIC      iff  lam - 2.576 SE > 0  and  min_k (lam_k - 2.576 SE_k) > 0  and the estimate has
               converged (|lam(full) - lam(half run)| <= 0.25 |lam(full)| for every k)
  NON-CHAOTIC  iff  lam + 2.576 SE <= tol, tol = 0 for maps; for flows / delay systems the
               largest exponent of a periodic or quasi-periodic attractor is 0, so
               tol = 2.576 * SE (i.e. "at or below zero within uncertainty") with the extra
               condition lam <= 0.1 * (smallest CHAOTIC lam of the same family) to keep weak
               chaos out
  otherwise    AMBIGUOUS (discarded)
Units: per iteration (maps: per stimulus / per beat) or per model time unit (flows).
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
from multiprocessing import Pool

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "results" / "ground_truth"
K = 3

GRIDS = {
    "phase_reset": [{"aggregate": "A", "tau": round(t, 3)} for t in np.arange(0.30, 2.001, 0.02)],
    "av_node": [{"H": h} for h in (36.0, 40.0, 45.0, 50.0, 52.0, 55.0, 56.0, 58.0, 60.0, 70.0, 90.0)],
    "mackey_glass": [{"tau": t} for t in (12.0, 14.0, 15.0, 16.0, 16.5, 17.0, 18.0, 20.0, 23.0, 26.0, 30.0)],
    "coupled_vdp": [{"rho": r, "omega": w} for r in (0.0, 2.0, 4.0, 6.0, 8.0, 10.0)
                    for w in (1.5, 2.1, 2.7, 3.3, 4.0, 5.6)
                    if not (r == 0.0 and w != 1.5)] + [{"rho": 5.45, "omega": 5.6}, {"rho": 8.625, "omega": 2.1},
                                                         {"rho": 9.6, "omega": 2.1}],
}
IS_MAP = {"phase_reset": True, "av_node": True, "mackey_glass": False, "coupled_vdp": False}


def one(task):
    family, params, k = task
    if family == "phase_reset":
        from experiments.phase8_cardiac.models import phase_reset as M
        p = M.prc_params(params["aggregate"])
        lam, se, conv = M.lyapunov(params["tau"], p, n=200_000, phi0=[0.123, 0.517, 0.871][k])
    elif family == "av_node":
        from experiments.phase8_cardiac.models import av_node as M
        lam, se, conv = M.lyapunov(params["H"], n=100_000, A0=[100.0, 120.0, 140.0][k])
    elif family == "mackey_glass":
        from experiments.phase8_cardiac.models import mackey_glass as M
        lam, se, conv = M.lyapunov(params["tau"], total_time=2.0e5, dt=0.05, seed=k)
    else:
        from experiments.phase8_cardiac.models import coupled_vdp as M
        lam, se, conv, _ = M.lyapunov({"rho": params["rho"], "omega": params["omega"]}, model_time=8000.0,
                                      transient=500.0 + 100.0 * k, seed=k)
    keys = sorted(conv)
    half = conv[keys[-2]]
    return {"family": family, **{f"p_{a}": b for a, b in params.items()}, "k": k, "lam": lam, "se": se,
            "lam_half": half, "converged": bool(abs(lam - half) <= 0.25 * abs(lam)) if lam != 0 else True,
            "conv": json.dumps(conv)}


def scan(workers):
    OUT.mkdir(parents=True, exist_ok=True)
    tasks = [(f, p, k) for f, grid in GRIDS.items() for p in grid for k in range(K)]
    rows = []
    with Pool(workers) as pool:
        for i, r in enumerate(pool.imap_unordered(one, tasks), 1):
            rows.append(r)
            if i % 20 == 0:
                print(i, "/", len(tasks), flush=True)
    pd.DataFrame(rows).to_csv(OUT / "scan_raw.csv", index=False)


def select():
    raw = pd.read_csv(OUT / "scan_raw.csv")
    pcols = [c for c in raw.columns if c.startswith("p_")]
    rows = []
    for (fam, *pv), g in raw.groupby(["family"] + pcols, dropna=False):
        lam = g.lam.mean()
        se = float(np.sqrt((g.se ** 2).mean() / len(g) + g.lam.var(ddof=1) / len(g)))
        lo_each = (g.lam - 2.576 * g.se).min()
        rows.append({"family": fam, **{c: v for c, v in zip(pcols, pv)}, "lam": lam, "se": se,
                     "lam_min_seed": g.lam.min(), "lam_max_seed": g.lam.max(), "lo_each": lo_each,
                     "converged_all": bool(g.converged.all()), "n": len(g)})
    df = pd.DataFrame(rows)
    labels = []
    for fam, g in df.groupby("family"):
        ch = g[(g.lam - 2.576 * g.se > 0) & (g.lo_each > 0) & g.converged_all]
        min_ch = ch.lam.min() if len(ch) else np.inf
        for i, r in g.iterrows():
            if i in ch.index:
                lab = "CHAOTIC"
            elif IS_MAP[fam] and r.lam + 2.576 * r.se <= 0:
                lab = "NON-CHAOTIC"
            elif (not IS_MAP[fam]) and r.lam - 2.576 * r.se <= 0 and r.lam <= 0.1 * min_ch:
                lab = "NON-CHAOTIC"
            else:
                lab = "AMBIGUOUS"
            labels.append((i, lab))
    for i, lab in labels:
        df.loc[i, "label"] = lab
    df.to_csv(OUT / "scan_labels.csv", index=False)
    print(df.groupby(["family", "label"]).size())
    return df


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan", action="store_true")
    ap.add_argument("--select", action="store_true")
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    if a.scan:
        scan(a.workers)
    if a.select:
        select()
