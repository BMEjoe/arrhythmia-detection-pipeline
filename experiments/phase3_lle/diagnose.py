"""
Phase 3, B2 -- diagnosis of the Rosenstein LLE bias (DEVELOPMENT seeds only).

For each development window, compute the mean log-divergence curve of the
frozen fp.rosenstein_lle under four embeddings:
    V0  production: TDMI tau*, Cao m* at tau*          (current pipeline)
    V1  tau = 1,  m = m*                               (delay changed only)
    V2  tau = 1,  m = Cao at lag 1                     (map-compatible delay, Cao dimension)
    V3  tau = 1,  m = true map dimension (1 logistic, 2 Henon; 2 for other systems)
and several fit regions on the same curve:
    prod       the production early fit (first 6 valid steps, fp.rosenstein_lle)
    k1_6       steps 1..6 (drops the k = 0 point)
    sat        steps 1..k_sat, k_sat = first step where the curve reaches
               half of its rise from y[1] to the plateau (mean of y[25:40]); >= 3 points
Saturation index and the k = 0 -> 1 jump are recorded to measure the noise floor.

    python -m experiments.phase3_lle.diagnose            # writes results/dev/diagnosis.jsonl
    python -m experiments.phase3_lle.diagnose --summary  # tables + plots
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import json
import multiprocessing as mp
import pathlib
import sys
import warnings

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import final_pipeline as fp  # noqa: E402
from experiments.phase2e import systems as S  # noqa: E402
from experiments.phase3_lle import config as C  # noqa: E402

CFG = fp.CFG
HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "results" / "dev" / "diagnosis.jsonl"
N_DIAG_SEEDS = 30          # first 30 development seeds of every condition


def _f(v):
    v = float(v)
    return v if np.isfinite(v) else None


def curve(x, tau, m, theiler):
    emb = fp._embed_backward(x, int(tau), int(m))
    lle, d = fp.rosenstein_lle(emb, theiler=theiler, max_iter=CFG.lle_max_iter,
                               min_fit_points=CFG.lle_min_fit_points,
                               max_fit_fraction=CFG.lle_max_fit_fraction, min_r2=CFG.lle_min_r2,
                               require_valid_fit=False)
    y = d.get("mean_log_divergence")
    return lle, d, (None if y is None else np.asarray(y, float)), len(emb)


def _slope(y, k0, k1):
    k = np.arange(k0, k1 + 1)
    yy = y[k]
    if len(k) < 2 or not np.all(np.isfinite(yy)):
        return None
    return _f(np.polyfit(k, yy, 1)[0])


def fits(y):
    if y is None or not np.all(np.isfinite(y[:8])):
        return {}
    plateau = float(np.nanmean(y[25:40]))
    rise = plateau - y[1]
    k_sat = None
    if rise > 0:
        above = np.flatnonzero(y[1:] >= y[1] + 0.5 * rise)
        k_sat = int(above[0] + 1) if len(above) else None
    k_end = max(3, k_sat) if k_sat is not None else 3
    k95 = None
    if plateau - y[0] > 0:
        a = np.flatnonzero(y >= y[0] + 0.95 * (plateau - y[0]))
        k95 = int(a[0]) if len(a) else None
    return {"k1_6": _slope(y, 1, 6), "k0_5": _slope(y, 0, 5), "sat": _slope(y, 1, k_end),
            "k_half": k_sat, "k95": k95, "plateau": _f(plateau), "y0": _f(y[0]), "y1": _f(y[1]),
            "jump01": _f(y[1] - y[0]), "jump12": _f(y[2] - y[1])}


def diagnose(task):
    system, snr, n, seed = task
    x = S.generate(system, n, seed, snr)
    th = fp._mean_period_beats(x)
    _, tau_p, m_p, *_ = fp.takens_embed(x, config=CFG)
    m1 = fp.cao_method(x, 1, CFG.cao_max_dim, CFG.cao_tol, CFG.cao_theiler)[0]
    m_true = C.TRUE_DIMENSION.get(system, 2)
    rec = {"system": system, "snr_db": snr, "window_length": n, "seed": seed, "theiler": th,
           "tau_prod": int(tau_p), "m_prod": int(m_p), "m_cao_lag1": int(m1), "m_true": m_true,
           "variants": {}}
    for name, (t, m) in {"V0": (tau_p, m_p), "V1": (1, m_p), "V2": (1, min(m1, CFG.cao_max_dim)),
                         "V3": (1, m_true)}.items():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            lle, d, y, ne = curve(x, t, m, th)
        fi = d.get("fit_indices")
        rec["variants"][name] = {
            "tau": int(t), "m": int(m), "n_embedded": ne, "prod": _f(lle) if lle is not None else None,
            "prod_fit_k": None if fi is None else [int(fi[0]), int(fi[-1])],
            "r2": _f(d["r2"]) if "r2" in d else None, "reason": d.get("reason"),
            "curve": None if y is None else [_f(v) for v in y], **fits(y)}
    return rec


def all_tasks():
    out = []
    for n in C.WINDOW_LENGTHS:
        for system, snr in C.conditions():
            for seed in C.seeds("dev", system, n, snr)[:N_DIAG_SEEDS]:
                out.append((system, snr, n, seed))
    return out


def run(workers):
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with mp.get_context("fork").Pool(workers) as pool, open(OUT, "w") as fh:
        for rec in pool.imap(diagnose, all_tasks(), chunksize=2):
            fh.write(json.dumps(rec) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--summary", action="store_true")
    a = ap.parse_args(argv)
    if a.summary:
        from experiments.phase3_lle import diagnose_summary
        diagnose_summary.main()
    else:
        run(a.workers)


if __name__ == "__main__":
    main()
