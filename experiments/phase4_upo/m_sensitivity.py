"""
Phase 4, item 6 -- sensitivity of the Phase 3 winner (fp.lle_chaos_test,
preregistered m = 2) to the embedding dimension.  REPORT ONLY; nothing is
changed or re-selected.

Windows: the Phase 3 TEST seeds 1000-1099 at 256 samples (Phase 2E generators),
conditions white_noise, ar1, sinusoid, logistic, henon, logistic@20dB,
henon@20dB.  fp.lle_chaos_test is run with lle_chaos_test_m = 2, 3, 4 (all
other preregistered parameters unchanged).  The m = 2 arm must reproduce the
stored Phase 3 C1 results exactly (checked in --summary).

    python -m experiments.phase4_upo.m_sensitivity            # results/m_sensitivity.jsonl
    python -m experiments.phase4_upo.m_sensitivity --summary  # results/tables/m_sensitivity.md
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
from dataclasses import replace

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import final_pipeline as fp  # noqa: E402
from experiments.phase2e import systems as S2  # noqa: E402
from experiments.phase3_lle import analysis as A3  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "results" / "m_sensitivity.jsonl"
TABLE = HERE / "results" / "tables" / "m_sensitivity.md"
CONDS = [("white_noise", None), ("ar1", None), ("sinusoid", None), ("logistic", None), ("henon", None),
         ("logistic", 20.0), ("henon", 20.0)]
MS = (2, 3, 4)
SEEDS = range(1000, 1100)
N = 256


def one(task):
    system, snr, seed = task
    x = S2.generate(system, N, seed, snr)
    out = {"system": system, "snr_db": snr, "seed": seed, "m": {}}
    for m in MS:
        r = fp.lle_chaos_test(x, config=replace(fp.CFG, lle_chaos_test=True, lle_chaos_test_m=m))
        out["m"][str(m)] = {"detected": bool(r["detected"]),
                            "lle": float(r["lle"]) if np.isfinite(r["lle"]) else None,
                            "p": float(r["p"]) if np.isfinite(r["p"]) else None}
    return out


def lab(system, snr):
    return system if snr is None else f"{system}@{int(snr)}dB"


def summary():
    R = [json.loads(line) for line in open(OUT)]
    p3 = {}
    for line in open(ROOT / "experiments/phase3_lle/results/test/c1_rosenstein_m2_iaaft.jsonl"):
        r = json.loads(line)
        t = r["task"]
        if t["window_length"] == N:
            p3[(t["system"], t["snr_db"], t["seed"])] = r["result"]
    mism = sum(1 for r in R if (p3[(r["system"], r["snr_db"], r["seed"])]["detected"] != r["m"]["2"]["detected"]
                                or p3[(r["system"], r["snr_db"], r["seed"])]["lle"] != r["m"]["2"]["lle"]))
    ref = {"logistic": np.log(2.0), "henon": 0.4192}
    md = ["# lle_chaos_test m-sensitivity (Phase 3 test seeds 1000-1099, 256 samples)", "",
          f"m = 2 reproduces the stored Phase 3 C1 results: {len(R) - mism}/{len(R)} windows identical "
          "(decision and LLE).", "",
          "Detections / 100 with 95 % Wilson CI; median LLE (bias vs reference for chaotic systems).", "",
          "| condition | m = 2 | m = 3 | m = 4 | median LLE m = 2 / 3 / 4 |", "|---|---|---|---|---|"]
    for s, snr in CONDS:
        g = [r for r in R if r["system"] == s and r["snr_db"] == snr]
        cells = [A3.frac(sum(r["m"][str(m)]["detected"] for r in g), len(g)) for m in MS]
        meds = []
        for m in MS:
            v = [r["m"][str(m)]["lle"] for r in g if r["m"][str(m)]["lle"] is not None]
            med = float(np.median(v)) if v else float("nan")
            meds.append(f"{med:.3f}" + (f" ({med - ref[s]:+.3f})" if s in ref else ""))
        md.append(f"| {lab(s, snr)} | " + " | ".join(cells) + " | " + " / ".join(meds) + " |")
    TABLE.parent.mkdir(parents=True, exist_ok=True)
    TABLE.write_text("\n".join(md) + "\n")
    print("\n".join(md))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--summary", action="store_true")
    a = ap.parse_args(argv)
    if a.summary:
        summary()
        return
    OUT.parent.mkdir(parents=True, exist_ok=True)
    todo = [(s, snr, seed) for s, snr in CONDS for seed in SEEDS]
    with mp.get_context("fork").Pool(a.workers) as pool, open(OUT, "w") as fh:
        for rec in pool.imap(one, todo, chunksize=2):
            fh.write(json.dumps(rec) + "\n")
            fh.flush()


if __name__ == "__main__":
    main()
