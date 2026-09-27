"""
Same-environment before/after replicability check for final_pipeline.py changes.

The Phase 2E replicability check (run_phase2e --replicate) compares reruns with
the records stored in experiments/phase2e/results/, which were produced on a
different machine.  In this container that check gives 13/44 bitwise-identical
tasks even with the unmodified pipeline (floating-point differences at the
1e-15 level, amplified by ill-conditioned quantities; see HANDOFF.md, A2).

This script therefore
  --save PATH     reruns the same 44 tasks (first 4 of every Phase 2E experiment)
                  with the CURRENT pipeline and stores the full records;
  --compare PATH  reruns them again and requires every row and detail to be
                  bitwise identical to the saved records.
Use --save before changing final_pipeline.py and --compare after.

    python -m experiments.phase3_lle.replicate_check --save /tmp/rep_before.jsonl
    python -m experiments.phase3_lle.replicate_check --compare /tmp/rep_before.jsonl
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import json
import multiprocessing as mp
import sys

from experiments.phase2e import config as C
from experiments.phase2e import run_phase2e as R


def _chosen():
    chosen, seen = [], {}
    for t in R.all_tasks(list(R.EXPERIMENTS)):
        k = seen.get(t["experiment"], 0)
        if k < C.REPLICATE_PER_EXPERIMENT:
            chosen.append(t)
            seen[t["experiment"]] = k + 1
    return chosen


def _run(workers):
    with mp.get_context("fork").Pool(workers) as pool:
        return pool.map(R.execute, _chosen(), chunksize=1)


def _key(rec):
    strip = {k: v for k, v in rec["row"].items() if k != "n_runtime_warnings"}
    return R.dumps(strip), R.dumps(rec["detail"])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--save")
    g.add_argument("--compare")
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args(argv)
    recs = _run(a.workers)
    if a.save:
        with open(a.save, "w") as fh:
            for r in recs:
                fh.write(R.dumps(r) + "\n")
        print(f"saved {len(recs)} records to {a.save}")
        return 0
    before = {}
    with open(a.compare) as fh:
        for line in fh:
            r = json.loads(line)
            before[r["task"]["task_id"]] = r
    bad = []
    for r in recs:
        old = before.get(r["task"]["task_id"])
        if old is None or _key(json.loads(R.dumps(old))) != _key(json.loads(R.dumps(r))):
            bad.append(r["task"]["task_id"])
    print(f"same-environment replicability: {len(recs) - len(bad)}/{len(recs)} bitwise identical")
    if bad:
        print("differing:", bad)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
