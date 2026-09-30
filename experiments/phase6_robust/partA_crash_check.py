"""
Phase 6 Part A (report only): BASELINE-K = combined_chaos_config() with
keep_upo_on_short_lle_embedding = True, on
  - the three Phase 5 test windows where analyze_segment raised
    (P1_henon_rr 3008, 3154; P4_henon_rr_30dB 3014), and
  - N1-N6 Phase 5 test seeds 3000-3299 (1,800 windows).
Compares LLE / UPO / AND with the stored Phase 5 m = 2 decisions
(experiments/phase5_rr/results/test/test.jsonl).  Writes
results/partA/partA.{jsonl,md}.

    python -m experiments.phase6_robust.partA_crash_check [--workers 4]
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

ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import final_pipeline as fp  # noqa: E402
from experiments.phase5_rr import analysis as A5  # noqa: E402
from experiments.phase5_rr import systems as S5  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "results" / "partA" / "partA.jsonl"
BASELINE_K = replace(fp.combined_chaos_config(), keep_upo_on_short_lle_embedding=True)
CRASHED = [("P1_henon_rr", 3008), ("P1_henon_rr", 3154), ("P4_henon_rr_30dB", 3014)]
NULLS = ["N1_linear_rr", "N2_power_law", "N3_linear_rr_trend", "N4_linear_rr_step",
         "N5_linear_rr_warped", "N6_noisy_rsa"]


def one(task):
    cond, seed = task
    t = time.perf_counter()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            out = fp.analyze_segment(S5.generate(cond, seed), BASELINE_K)
        except Exception as exc:                                  # reported
            return {"condition": cond, "seed": seed, "error": f"{type(exc).__name__}: {exc}"}
    return {"condition": cond, "seed": seed, "error": None,
            "lle": bool(out["lle_chaos_test"]["detected"]), "lle_p": out["lle_chaos_test"]["p"],
            "upo": bool(out["upo"]["instability_gate_detected"]), "and": fp.combined_chaos_detected(out),
            "lle_embedding_too_short": bool(out.get("lle_embedding_too_short")),
            "runtime_s": time.perf_counter() - t}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args(argv)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if OUT.exists():
        done = {(r["condition"], r["seed"]) for r in map(json.loads, open(OUT))}
    tasks = [t for t in CRASHED + [(c, s) for c in NULLS for s in range(3000, 3300)] if t not in done]
    if tasks:
        with mp.get_context("fork").Pool(a.workers) as pool, open(OUT, "a") as fh:
            for k, rec in enumerate(pool.imap_unordered(one, tasks, chunksize=1), 1):
                fh.write(json.dumps(rec) + "\n")
                fh.flush()
                if k % 100 == 0:
                    print(f"{k}/{len(tasks)}", flush=True)
    new = {(r["condition"], r["seed"]): r for r in map(json.loads, open(OUT))}
    old = {(r["task"]["condition"], r["task"]["seed"]): r["m2"] for r in A5.load("test", "test")}
    md = ["# Phase 6 Part A: BASELINE-K (keep_upo_on_short_lle_embedding = True) vs Phase 5", "",
          "## The three Phase 5 windows that raised", "",
          "| window | Phase 5 | BASELINE-K error | LLE embedding too short | LLE | UPO | AND |", "|---|---|---|---|---|---|---|"]
    for c, s in CRASHED:
        n, o = new[(c, s)], old[(c, s)]
        md.append(f"| {c} {s} | raised ({o['error'].split(':')[0]}); AND counted 0 | {n['error'] or 'none'} | "
                  f"{n.get('lle_embedding_too_short')} | {n.get('lle')} | {n.get('upo')} | {n.get('and')} |")
    md += ["", "## N1-N6, Phase 5 test seeds 3000-3299", "",
           "| condition | AND Phase 5 | AND BASELINE-K | LLE P5 / K | UPO P5 / K | windows with any difference | errors (K) |",
           "|---|---|---|---|---|---|---|"]
    total_diff = 0
    for c in NULLS:
        keys = [(c, s) for s in range(3000, 3300)]
        diff = sum(any(new[k].get(f) != old[k][f] for f in ("lle", "upo")) or
                   (new[k]["error"] is None) != (old[k].get("error") is None) for k in keys)
        total_diff += diff
        ak = sum(bool(new[k].get("and")) for k in keys)
        ao = sum(old[k]["lle"] and old[k]["upo"] for k in keys)
        md.append(f"| {c} | {ao}/300 | {ak}/300 | {sum(old[k]['lle'] for k in keys)} / "
                  f"{sum(bool(new[k].get('lle')) for k in keys)} | {sum(old[k]['upo'] for k in keys)} / "
                  f"{sum(bool(new[k].get('upo')) for k in keys)} | {diff} | "
                  f"{sum(new[k]['error'] is not None for k in keys)} |")
    md += ["", f"Windows with any decision difference over N1-N6: **{total_diff}** / 1800."]
    (OUT.parent / "partA.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
