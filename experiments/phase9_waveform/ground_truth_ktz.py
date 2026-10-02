"""
Phase 9 Part A2: ground-truth Lyapunov exponents of the KTz repolarization (APD) dynamics FOR
THE EXACT RR INPUT USED (families.ktz_input: paced at 0.8 s, optionally with Phase 8 ectopy).

    python -m experiments.phase9_waveform.ground_truth_ktz --workers 4

Candidate regimes (predeclared from the constant-pacing scan in results/verification/ktz_scan.jsonl):
  chaotic candidates     P_nom = 92, 100, 146, 230, 232, 260, 278 ts
  non-chaotic candidates P_nom = 120, 150, 200, 250, 300 ts
Inputs: 'none' (constant pacing), 'S2_5', 'E1', 'E3_10'.
Per (P_nom, input): K = 3 realizations (ground-truth seeds 990000-990002: initial state and
ectopy draw), 3,000 transient beats discarded, then 8,000 beats with the largest LE of the
3-D map accumulated by QR at every time step (ktz.run, want_le=True).
Per realization: lambda_k (per time step) and its batch-means SE_k (20 blocks of beats, LE per
block = summed log growth / steps in the block); convergence: |lambda(all) - lambda(first half)|
<= 0.25 |lambda(all)|.
Labelling rule (as Phase 8 for maps):
  combined lambda = mean_k lambda_k; SE = sqrt(mean_k SE_k^2 / K + var_k(lambda_k) / K)
  CHAOTIC      iff lambda - 2.576 SE > 0 and min_k(lambda_k - 2.576 SE_k) > 0 and all converged
  NON-CHAOTIC  iff lambda + 2.576 SE <= 0
  otherwise    AMBIGUOUS (discarded)
A (P_nom, input) pair is a REGIME; only CHAOTIC / NON-CHAOTIC pairs are kept.  APD_ref (mean APD
over the 8,000 beats of all realizations, time steps) is stored for the T-wave mapping.
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
from multiprocessing import Pool

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402

from experiments.phase9_waveform import families as F  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "results" / "ground_truth"
CHAOTIC_CAND = (92, 100, 146, 230, 232, 260, 278)
NONCHAOTIC_CAND = (120, 150, 200, 250, 300)
INPUTS = ("none", "S2_5", "E1", "E3_10")
GT_SEEDS = (990000, 990001, 990002)
TRANSIENT, N_LE, N_BLOCKS = 3000, 8000, 20


def one(task):
    from experiments.phase9_waveform.morph import ktz
    P, inp, k = task
    rng = np.random.default_rng(F.ss9(F.FAMILY9["ktz"], 0, 0, F.ECT_CODE[inp], GT_SEEDS[k]))
    per, types, step = F.ktz_input(P, inp, TRANSIENT + N_LE, rng)
    x0, y0, z0 = rng.uniform(-0.1, 0.1, size=3)
    _, _, _, _, x, y, z, _ = ktz.run(per[:TRANSIENT], x0, y0, z0)
    apd, _, le, ns, *_, le_beat = ktz.run(per[TRANSIENT:], x, y, z, want_le=True)
    p = per[TRANSIENT:].astype(float)
    lam = le / ns
    blocks = np.array_split(np.arange(N_LE), N_BLOCKS)
    bl = np.array([le_beat[b].sum() / p[b].sum() for b in blocks])
    half = N_LE // 2
    lam_half = le_beat[:half].sum() / p[:half].sum()
    a = apd[~np.isnan(apd)]
    return {"P_nom": P, "input": inp, "k": k, "lam": float(lam), "se": float(bl.std(ddof=1) / np.sqrt(N_BLOCKS)),
            "lam_half": float(lam_half), "converged": bool(abs(lam - lam_half) <= 0.25 * abs(lam)) if lam != 0 else True,
            "apd_mean_ts": float(a.mean()), "apd_sd_ts": float(a.std()), "n_apd_nan": int(np.isnan(apd).sum()),
            "n_distinct_apd": int(len(np.unique(a))), "step_s": step}


def label(rows):
    lam = np.array([r["lam"] for r in rows])
    se = np.array([r["se"] for r in rows])
    K = len(rows)
    L = float(lam.mean())
    SE = float(np.sqrt(np.mean(se ** 2) / K + np.var(lam, ddof=1) / K))
    if L - 2.576 * SE > 0 and np.min(lam - 2.576 * se) > 0 and all(r["converged"] for r in rows):
        lab = "CHAOTIC"
    elif L + 2.576 * SE <= 0:
        lab = "NON-CHAOTIC"
    else:
        lab = "AMBIGUOUS"
    return lab, L, SE


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    tasks = [(P, inp, k) for P in CHAOTIC_CAND + NONCHAOTIC_CAND for inp in INPUTS for k in range(3)]
    with Pool(a.workers) as pool:
        rows = list(pool.imap_unordered(one, tasks))
    with open(OUT / "ktz_scan_raw.jsonl", "w") as fh:
        for r in sorted(rows, key=lambda r: (r["P_nom"], r["input"], r["k"])):
            fh.write(json.dumps(r) + "\n")
    regimes, table = [], []
    for P in CHAOTIC_CAND + NONCHAOTIC_CAND:
        for inp in INPUTS:
            rr = [r for r in rows if r["P_nom"] == P and r["input"] == inp]
            lab, L, SE = label(rr)
            apd_ref = float(np.mean([r["apd_mean_ts"] for r in rr]))
            table.append({"P_nom": P, "input": inp, "label": lab, "lam_per_ts": L, "se": SE,
                          "lam_per_beat": L * P, "apd_ref_ts": apd_ref,
                          "apd_sd_ts": float(np.mean([r["apd_sd_ts"] for r in rr])),
                          "n_distinct_apd": int(np.median([r["n_distinct_apd"] for r in rr])),
                          "lam_k": [r["lam"] for r in rr]})
    json.dump(table, open(OUT / "ktz_labels.json", "w"), indent=1)
    for t in table:
        print(f"P={t['P_nom']:3d} {t['input']:6s} {t['label']:12s} lam/ts={t['lam_per_ts']:+.5f} +/- {t['se']:.5f} "
              f"lam/beat={t['lam_per_beat']:+.3f} APDref={t['apd_ref_ts']:.1f} sd={t['apd_sd_ts']:.1f} distinct={t['n_distinct_apd']}")


if __name__ == "__main__":
    main()
