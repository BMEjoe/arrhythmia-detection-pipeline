"""
Phase 5 runner: one JSONL line per (phase, condition, seed) with the frozen
detector's decisions at m = 2 (primary; fp.analyze_segment) and at m = 3 and
m = 4 (report only; both components' m set together; detector.decide).

    python -m experiments.phase5_rr.run_phase5 --phase dev --seeds 0-9 --tag runtime
    python -m experiments.phase5_rr.run_phase5 --phase test

Checkpointed: finished task ids are skipped, so an interrupted run is resumed
with the same command.  --phase test reuses the Phase 3 guard
(experiments/phase3_lle/run_phase3.check_preregistration) pointed at this
directory's PREREGISTRATION.md: committed, unmodified, pushed and listing the
method `and_detector`.
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import hashlib
import json
import multiprocessing as mp
import pathlib
import platform
import sys
import time
import warnings

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.phase3_lle import run_phase3 as R3  # noqa: E402
from experiments.phase5_rr import config as C  # noqa: E402
from experiments.phase5_rr import detector as D  # noqa: E402
from experiments.phase5_rr import systems as S  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
RESULTS = HERE / "results"
PREREG = HERE / "PREREGISTRATION.md"
M_VALUES = (2, 3, 4)
METHODS = ("and_detector",)


def tasks(phase, tag, conds, seeds=None):
    out = []
    for cond in conds:
        for seed in (seeds if seeds is not None else C.seeds(phase, cond)):
            out.append({"task_id": f"{phase}:{tag}:{cond}:{seed}", "phase": phase,
                        "condition": cond, "seed": int(seed)})
    return out


def execute(args):
    task, m_values = args
    t0 = time.perf_counter()
    rr = S.generate(task["condition"], task["seed"])
    rec = {"task": task, "data_seed": S.data_seed(task["condition"], task["seed"]), "runtime_s": {}}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for m in m_values:
            t = time.perf_counter()
            rec[f"m{m}"] = D.decide(rr, m)
            rec["runtime_s"][f"m{m}"] = time.perf_counter() - t
    rec["runtime_s"]["total"] = time.perf_counter() - t0
    return rec


def _done(path):
    ids = set()
    if path.exists():
        for line in open(path):
            try:
                ids.add(json.loads(line)["task"]["task_id"])
            except (json.JSONDecodeError, KeyError):
                pass
    return ids


def _seed_range(spec):
    if spec is None:
        return None
    lo, _, hi = spec.partition("-")
    return list(range(int(lo), int(hi or lo) + 1))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--phase", choices=("dev", "test"), required=True)
    ap.add_argument("--tag", default=None)
    ap.add_argument("--conditions", default=None, help="comma list (default: all)")
    ap.add_argument("--seeds", default=None, help="dev only: a-b range (default: config dev seeds)")
    ap.add_argument("--m", default="2,3,4")
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args(argv)
    m_values = tuple(int(v) for v in a.m.split(","))
    if a.phase == "test":
        if a.seeds or m_values != M_VALUES:
            raise SystemExit("--phase test runs the preregistered seeds and m values only")
        R3.PREREG = PREREG
        R3.check_preregistration(list(METHODS))
        tag = a.tag or "test"
    else:
        tag = a.tag or "dev"
    conds = list(C.CONDITIONS) if not a.conditions else a.conditions.split(",")
    outdir = RESULTS / a.phase
    outdir.mkdir(parents=True, exist_ok=True)
    out_path = outdir / f"{tag}.jsonl"
    sha = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
           for p in [ROOT / "final_pipeline.py", *sorted(HERE.glob("*.py")), *([PREREG] if PREREG.exists() else [])]}
    with open(outdir / "manifest.jsonl", "a") as fh:
        fh.write(json.dumps({"started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                             "phase": a.phase, "tag": tag, "m_values": list(m_values),
                             "git_head": R3._git("rev-parse", "HEAD").stdout.strip(),
                             "git_status": R3._git("status", "--porcelain").stdout, "sha256": sha,
                             "python": sys.version, "platform": platform.platform(),
                             "numpy": np.__version__, "workers": a.workers}) + "\n")
    done = _done(out_path)
    todo = [t for t in tasks(a.phase, tag, conds, _seed_range(a.seeds)) if t["task_id"] not in done]
    print(f"{len(todo)} windows to run; m = {m_values}", flush=True)
    t0 = time.time()
    with mp.get_context("fork").Pool(a.workers) as pool, open(out_path, "a") as fh:
        for k, rec in enumerate(pool.imap_unordered(execute, [(t, m_values) for t in todo], chunksize=1), 1):
            fh.write(json.dumps(rec, allow_nan=False, default=str) + "\n")
            fh.flush()
            if k % 25 == 0 or k == len(todo):
                print(f"{k}/{len(todo)} done, {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
