"""
Phase 3, B3 -- candidate tuning on DEVELOPMENT seeds only (256 samples).

For every development window and every grid combination, compute the divergence
curve of the window and of N_SUR surrogates, then several fit rules on the same
curves.  Writes results/dev/tuning.jsonl; summarize with --summary.

Grid
  m          2, 3                     (tau = 1 always; B2: TDMI tau and Cao m bias the LLE)
  nbr        nn1   Rosenstein single nearest neighbour
             eps   single nearest neighbour beyond eps = EPS_FACTOR x median NN distance
             knn5  mean distance of the 5 nearest neighbours (Kantz-style)
  surrogate  aaft (fp.aaft_surrogate), iaaft
  fit        f05   steps 0..5 (production-style)
             f15   steps 1..5
             s50   steps 1..k_end, k_end = half-rise to plateau (3..12 points)
             s30   same with 30 % of the rise
Also the 0-1 test K (Gottwald-Melbourne) with the same surrogates.
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

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import final_pipeline as fp  # noqa: E402
from experiments.phase2e import systems as S  # noqa: E402
from experiments.phase3_lle import config as C  # noqa: E402
from experiments.phase3_lle import estimators as E  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "results" / "dev" / "tuning.jsonl"
N_SUR = 39
EPS_FACTOR = 1.0
MS = (2, 3)
NBRS = ("nn1", "eps", "knn5")
SURS = ("aaft", "iaaft")
FITS = {"f05": dict(k_start=0, fixed=5), "f15": dict(k_start=1, fixed=5),
        "s50": dict(k_start=1, frac=0.5), "s30": dict(k_start=1, frac=0.3)}


def fit(y, rule):
    if "fixed" in rule:
        ks = np.arange(rule["k_start"], rule["fixed"] + 1)
        yy = y[ks]
        if not np.all(np.isfinite(yy)):
            return np.nan
        return float(np.polyfit(ks, yy, 1)[0])
    return E.saturation_fit(y, k_start=rule["k_start"], frac=rule["frac"])[0]


def nn_median(X, theiler):
    D = np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(-1))
    i = np.arange(len(X))
    D[np.abs(i[:, None] - i[None, :]) <= theiler] = np.inf
    d = D.min(1)
    d = d[np.isfinite(d) & (d > 1e-15)]
    return float(np.median(d)) if len(d) else 0.0


def stats(x, m, nbr, theiler):
    X = fp._embed_backward(x, 1, m)
    if nbr == "nn1":
        y, _ = E.divergence_curve(X, theiler)
    elif nbr == "eps":
        y, _ = E.divergence_curve(X, theiler, min_dist=EPS_FACTOR * nn_median(X, theiler))
    else:
        y, _ = E.divergence_curve(X, theiler, n_neighbors=5)
    return {f: fit(y, r) for f, r in FITS.items()}


def one(task):
    system, snr, n, seed = task
    x = S.generate(system, n, seed, snr)
    th = fp._mean_period_beats(x)
    rec = {"system": system, "snr_db": snr, "window_length": n, "seed": seed, "theiler": th, "res": {}}
    with E.warnings_ignored():
        for kind in SURS:
            rng = E.window_rng(x, stream={"aaft": 1, "iaaft": 2}[kind])
            surs = E.surrogates(x, kind, N_SUR, rng)
            rec["res"][f"K01|{kind}"] = {"obs": E.zero_one_K(x, np.random.default_rng(0)),
                                        "sur": [E.zero_one_K(s, np.random.default_rng(0)) for s in surs]}
            for m in MS:
                for nbr in NBRS:
                    obs = stats(x, m, nbr, th)
                    sur = [stats(s, m, nbr, th) for s in surs]
                    for f in FITS:
                        rec["res"][f"m{m}|{nbr}|{kind}|{f}"] = {"obs": obs[f], "sur": [d[f] for d in sur]}
    return rec


def all_tasks(n=256):
    out = []
    for system, snr in C.conditions():
        for seed in C.seeds("dev", system, n, snr):
            out.append((system, snr, n, seed))
    return out


def pval(obs, sur):
    sur = np.asarray([s for s in sur if s is not None and np.isfinite(s)])
    if obs is None or not np.isfinite(obs) or len(sur) == 0:
        return np.nan
    return (1 + np.sum(sur >= obs)) / (len(sur) + 1)


def summary():
    R = [json.loads(line) for line in open(OUT)]
    combos = list(R[0]["res"].keys())
    labs = [C.condition_label(s, snr) for s, snr in C.conditions()]
    rows = []
    for cb in combos:
        row = {"combo": cb}
        for lab in labs:
            g = [r for r in R if C.condition_label(r["system"], r["snr_db"]) == lab]
            p = [pval(r["res"][cb]["obs"], r["res"][cb]["sur"]) for r in g]
            row[lab] = np.mean([(v <= 0.05) if np.isfinite(v) else False for v in p])
            ref = C.LLE_REFERENCE.get(lab.split("@")[0])
            if ref and not cb.startswith("K01"):
                o = [r["res"][cb]["obs"] for r in g]
                o = [v for v in o if v is not None and np.isfinite(v)]
                row["bias:" + lab] = (np.median(o) - ref) if o else np.nan
        rows.append(row)
    rows.sort(key=lambda r: -(r["logistic@20dB"] + r["henon@20dB"]))
    head = ["combo", "white_noise", "ar1", "sinusoid", "logistic_p4", "logistic", "henon",
            "logistic@30dB", "henon@30dB", "logistic@20dB", "henon@20dB", "logistic@10dB", "henon@10dB",
            "bias:logistic", "bias:henon", "bias:logistic@20dB", "bias:henon@20dB"]
    print(f"rates = fraction with p <= 0.05 ({N_SUR} surrogates, nominal 2/40 = 0.05); dev seeds, 256")
    print(" | ".join(h.replace("white_noise", "WN") for h in head))
    for r in rows:
        print(" | ".join([r["combo"]] + [f"{r.get(h, float('nan')):.2f}" for h in head[1:]]))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--summary", action="store_true")
    a = ap.parse_args(argv)
    if a.summary:
        summary()
        return
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with mp.get_context("fork").Pool(a.workers) as pool, open(OUT, "w") as fh:
        for rec in pool.imap(one, all_tasks(), chunksize=1):
            fh.write(json.dumps(rec, default=float) + "\n")
            fh.flush()


if __name__ == "__main__":
    main()
