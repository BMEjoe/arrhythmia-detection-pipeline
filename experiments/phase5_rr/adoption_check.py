"""
Recomputes stored test-seed decisions through the adopted helpers:
fp.combined_chaos_detected(fp.analyze_segment(rr, fp.combined_chaos_config()))
for seeds 3000-3004 of every condition (120 windows); compares LLE, UPO and AND
with results/test/test.jsonl (m2).  Writes results/test/adoption_check.json.

    python -m experiments.phase5_rr.adoption_check
"""
from __future__ import annotations

import json
import multiprocessing as mp
import pathlib
import sys
import warnings

ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import final_pipeline as fp  # noqa: E402
from experiments.phase5_rr import analysis as A  # noqa: E402
from experiments.phase5_rr import config as C  # noqa: E402
from experiments.phase5_rr import systems as S  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
SEEDS = range(3000, 3005)


def one(task):
    cond, seed = task
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            out = fp.analyze_segment(S.generate(cond, seed), fp.combined_chaos_config())
        except ValueError as exc:
            return {"condition": cond, "seed": seed, "error": str(exc)}
    return {"condition": cond, "seed": seed, "error": None, "lle": bool(out["lle_chaos_test"]["detected"]),
            "upo": bool(out["upo"]["instability_gate_detected"]), "and": fp.combined_chaos_detected(out)}


def main():
    stored = {(r["task"]["condition"], r["task"]["seed"]): r["m2"] for r in A.load("test", "test")}
    tasks = [(c, s) for c in C.CONDITIONS for s in SEEDS]
    with mp.get_context("fork").Pool(4) as pool:
        rows = pool.map(one, tasks)
    n_diff = 0
    for r in rows:
        ref = stored[(r["condition"], r["seed"])]
        if r["error"] is not None or ref.get("error") is not None:
            r["same"] = (r["error"] is not None) == (ref.get("error") is not None)
        else:
            r["same"] = (r["lle"] == ref["lle"] and r["upo"] == ref["upo"] and r["and"] == (ref["lle"] and ref["upo"]))
        n_diff += not r["same"]
    (HERE / "results" / "test" / "adoption_check.json").write_text(
        json.dumps({"n": len(rows), "n_differences": n_diff, "rows": rows}, indent=1))
    print(f"{len(rows)} windows, {n_diff} differences")


if __name__ == "__main__":
    main()
