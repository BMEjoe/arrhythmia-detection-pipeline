"""
Recomputes stored test decisions of the Phase 6 winner through the adopted helper,
fp.combined_chaos_detected(fp.analyze_segment(rr, fp.combined_chaos_config(
replace(CFG, keep_upo_on_short_lle_embedding=True), detrend=True))), for test seeds
4000-4002 of every condition (111 windows), and compares LLE / UPO / AND with
results/test/test.jsonl (d2_linear_g07).  Writes results/test/adoption_check.json.

    python -m experiments.phase6_robust.adoption_check
"""
from __future__ import annotations

import json
import multiprocessing as mp
import pathlib
import sys
import warnings
from dataclasses import replace

ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import final_pipeline as fp  # noqa: E402
from experiments.phase6_robust import analysis as A  # noqa: E402
from experiments.phase6_robust import config as C  # noqa: E402
from experiments.phase6_robust import systems as S  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
CFG = fp.combined_chaos_config(replace(fp.CFG, keep_upo_on_short_lle_embedding=True), detrend=True)


def one(task):
    cond, seed = task
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out = fp.analyze_segment(S.generate(cond, seed), CFG)
    return {"condition": cond, "seed": seed, "lle": bool(out["lle_chaos_test"]["detected"]),
            "upo": bool(out["upo"].get("instability_gate_detected", False)), "and": fp.combined_chaos_detected(out)}


def main():
    stored = {(r["task"]["condition"], r["task"]["seed"]): r["methods"]["d2_linear_g07"]
              for r in A.load("test", "test")}
    tasks = [(c, s) for c in C.CONDITIONS for s in range(4000, 4003)]
    with mp.get_context("fork").Pool(4) as pool:
        rows = pool.map(one, tasks)
    for r in rows:
        ref = stored[(r["condition"], r["seed"])]
        r["same"] = all(r[k] == ref[k] for k in ("lle", "upo", "and"))
    n = sum(not r["same"] for r in rows)
    (HERE / "results" / "test" / "adoption_check.json").write_text(
        json.dumps({"n": len(rows), "n_differences": n, "rows": rows}, indent=1))
    print(f"{len(rows)} windows, {n} differences")


if __name__ == "__main__":
    main()
