"""
Phase 5 Step 3 (report only): the Phase 4 C2 UPO configuration
(fp.phase4_upo_config()) with upo_fixed_dimension = 3 and 4 on the Phase 4
test seeds 2000-2149 at 256 samples, for white_noise, ar1, sinusoid, two_tone,
henon, henon@30dB, henon@20dB and logistic@20dB.  Detection =
run_upo_analysis(...)["instability_gate_detected"] (rng as in analyze_segment).
The stored m = 2 results (experiments/phase4_upo/results/test/test.jsonl, C2
rule on run m2_M15_S50) are shown beside them.

    python -m experiments.phase5_rr.phase4_m_check [--workers 4]
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
import time
import warnings
from dataclasses import replace

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import final_pipeline as fp  # noqa: E402
from experiments.phase4_upo import methods as MM  # noqa: E402
from experiments.phase4_upo import systems as S4  # noqa: E402
from experiments.phase5_rr.analysis import wilson  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "results" / "phase4_m_check.jsonl"
STORED = ROOT / "experiments/phase4_upo/results/test/test.jsonl"
CONDS = [("white_noise", None), ("ar1", None), ("sinusoid", None), ("two_tone", None), ("henon", None),
         ("henon", 30.0), ("henon", 20.0), ("logistic", 20.0)]
SEEDS = range(2000, 2150)
N = 256
M_VALUES = (3, 4)


def label(s, snr):
    return s if snr is None else f"{s}@{int(snr)}dB"


def execute(task):
    system, snr, seed, m = task
    cfg = replace(fp.phase4_upo_config(), upo_fixed_dimension=m)
    x = S4.generate(system, N, seed, snr)
    t = time.perf_counter()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            out = fp.run_upo_analysis(x, config=cfg, rng=np.random.default_rng(cfg.random_seed))
            det, status = bool(out.get("instability_gate_detected", False)), out["status"]
        except np.linalg.LinAlgError as exc:
            det, status = False, f"analysis_error {exc!r}"
    return {"condition": label(system, snr), "seed": seed, "m": m, "detected": det, "status": status,
            "runtime_s": time.perf_counter() - t}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--table-only", action="store_true")
    a = ap.parse_args(argv)
    done = set()
    if OUT.exists():
        for line in open(OUT):
            r = json.loads(line)
            done.add((r["condition"], r["seed"], r["m"]))
    todo = [(s, snr, seed, m) for m in M_VALUES for s, snr in CONDS for seed in SEEDS
            if (label(s, snr), seed, m) not in done]
    if todo and not a.table_only:
        print(f"{len(todo)} runs", flush=True)
        t0 = time.time()
        with mp.get_context("fork").Pool(a.workers) as pool, open(OUT, "a") as fh:
            for k, rec in enumerate(pool.imap_unordered(execute, todo, chunksize=1), 1):
                fh.write(json.dumps(rec) + "\n")
                fh.flush()
                if k % 50 == 0 or k == len(todo):
                    print(f"{k}/{len(todo)} {time.time() - t0:.0f} s", flush=True)
    new = {}
    for line in open(OUT):
        r = json.loads(line)
        new.setdefault((r["condition"], r["m"]), []).append(r["detected"])
    stored = {}
    for line in open(STORED):
        r = json.loads(line)
        t = r["task"]
        if t["window_length"] == N and t["seed"] in SEEDS:
            lab = label(t["system"], t["snr_db"])
            if lab in {label(*c) for c in CONDS}:
                stored.setdefault(lab, []).append(bool(MM.METHODS["c2_m2M15_mediangate"]["detect"](r["runs"]["m2_M15_S50"])))
    md = ["# Step 3: Phase 4 C2 UPO detector at upo_fixed_dimension = 2 / 3 / 4 (report only)", "",
          "Phase 4 test seeds 2000-2149, 256 samples; C2 otherwise unchanged (M = 15, 50 AAFT surrogates, "
          "median instability gate delta 0.2). m = 2: stored Phase 4 results.  Detections / windows "
          "(95 % Wilson interval).", "",
          "| condition | m = 2 (stored) | m = 3 | m = 4 |", "|---|---|---|---|"]
    for s, snr in CONDS:
        lab = label(s, snr)
        cells = []
        for v in [stored.get(lab, [])] + [new.get((lab, m), []) for m in M_VALUES]:
            k, n = sum(v), len(v)
            lo, hi = wilson(k, n)
            cells.append(f"{k}/{n} ({100 * lo:.0f}–{100 * hi:.0f} %)" if n else "–")
        md.append(f"| {lab} | " + " | ".join(cells) + " |")
    (HERE / "results" / "phase4_m_check.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
