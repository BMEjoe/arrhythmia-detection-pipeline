"""
Phase 6 runner: one JSONL line per (phase, condition, seed) with the decisions
of every requested method (methods.py) on that window.

    python -m experiments.phase6_robust.run_phase6 --phase dev --methods baseline_k,linear --conditions P3_henon_rr_trend --seeds 0-29 --tag tuning
    python -m experiments.phase6_robust.run_phase6 --phase test

Checkpointed: a window is skipped when its line exists with every requested
method; missing methods are NOT merged into old lines (use a new --tag).
--phase test reuses the Phase 3 guard (committed, unmodified, pushed
PREREGISTRATION.md listing every preregistered method) and runs exactly the
preregistered methods, conditions and seeds.
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

import final_pipeline as fp  # noqa: E402
from experiments.phase3_lle import run_phase3 as R3  # noqa: E402
from experiments.phase6_robust import config as C  # noqa: E402
from experiments.phase6_robust import detector as D  # noqa: E402
from experiments.phase6_robust import methods as MM  # noqa: E402
from experiments.phase6_robust import systems as S  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
RESULTS = HERE / "results"
PREREG = HERE / "PREREGISTRATION.md"
TEST_CONDITIONS = None     # fixed with PREREGISTRATION.md (list of condition names)


def method_config(name):
    if name in MM.PREREGISTERED:
        return MM.PREREGISTERED[name]["config"]
    return MM.TUNING[name]


def execute(args):
    """Runs every method on one window.  Unless shortcut is False, methods with the same
    effective configuration on this window (methods.effective_key: a gated detrend below its
    gate is BASELINE-K exactly; two linear detrends above their gates are identical) share one
    analyze_segment call; rec["computed_as"] records the sharing."""
    task, methods, shortcut = args
    t0 = time.perf_counter()
    rr = S.generate(task["condition"], task["seed"])
    ratio = fp.linear_trend_ratio(rr)
    rec = {"task": task, "data_seed": S.data_seed(task["condition"], task["seed"]), "trend_ratio": ratio,
           "methods": {}, "runtime_s": {}, "computed_as": {}}
    cache = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for name in methods:
            cfg = method_config(name)
            key = MM.effective_key(cfg, ratio) if shortcut else (name,)
            t = time.perf_counter()
            if key not in cache:
                cache[key] = (name, D.decide(rr, cfg))
            rec["methods"][name] = cache[key][1]
            rec["computed_as"][name] = cache[key][0]
            rec["runtime_s"][name] = time.perf_counter() - t
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
    ap.add_argument("--methods", default=None)
    ap.add_argument("--conditions", default=None)
    ap.add_argument("--seeds", default=None, help="dev only: a-b (default: config dev seeds)")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--no-shortcut", action="store_true", help="dev only: run every method in full")
    a = ap.parse_args(argv)
    if a.phase == "test":
        if a.seeds or a.methods or a.conditions or a.no_shortcut:
            raise SystemExit("--phase test runs the preregistered methods, conditions and seeds only")
        methods = list(MM.PREREGISTERED)
        R3.PREREG = PREREG
        R3.check_preregistration(methods)
        conds = list(TEST_CONDITIONS)
        tag = a.tag or "test"
    else:
        methods = a.methods.split(",") if a.methods else ["baseline_k"]
        conds = a.conditions.split(",") if a.conditions else list(C.CONDITIONS)
        tag = a.tag or "dev"
    outdir = RESULTS / a.phase
    outdir.mkdir(parents=True, exist_ok=True)
    out_path = outdir / f"{tag}.jsonl"
    sha = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
           for p in [ROOT / "final_pipeline.py", *sorted(HERE.glob("*.py")), *([PREREG] if PREREG.exists() else [])]}
    with open(outdir / "manifest.jsonl", "a") as fh:
        fh.write(json.dumps({"started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                             "phase": a.phase, "tag": tag, "methods": methods, "conditions": conds,
                             "git_head": R3._git("rev-parse", "HEAD").stdout.strip(),
                             "git_status": R3._git("status", "--porcelain").stdout, "sha256": sha,
                             "python": sys.version, "platform": platform.platform(),
                             "numpy": np.__version__, "workers": a.workers}) + "\n")
    done = _done(out_path)
    seeds = _seed_range(a.seeds)
    todo = []
    for cond in conds:
        for seed in (seeds if seeds is not None else C.seeds(a.phase, cond)):
            tid = f"{a.phase}:{tag}:{cond}:{seed}"
            if tid not in done:
                todo.append({"task_id": tid, "phase": a.phase, "condition": cond, "seed": int(seed)})
    print(f"{len(todo)} windows to run; methods = {methods}", flush=True)
    t0 = time.time()
    with mp.get_context("fork").Pool(a.workers) as pool, open(out_path, "a") as fh:
        for k, rec in enumerate(pool.imap_unordered(execute, [(t, methods, not a.no_shortcut) for t in todo],
                                                          chunksize=1), 1):
            fh.write(json.dumps(rec, allow_nan=False, default=str) + "\n")
            fh.flush()
            if k % 25 == 0 or k == len(todo):
                print(f"{k}/{len(todo)} done, {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
