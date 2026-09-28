"""
Phase 4 runner: one JSONL line per (phase, window).  Each line holds the output
of every requested detector RUN configuration (detector.run) plus
fp.lle_chaos_test on the same window (for the combined detector).  Methods are
decision rules applied to these records (methods.py), so they share runs.

    python -m experiments.phase4_upo.run_phase4 --phase dev --runs cao:7:50,2:7:50 --tag explore
    python -m experiments.phase4_upo.run_phase4 --phase test --methods all

Checkpointed (finished task ids are skipped).  --phase test reuses the Phase 3
guard (experiments/phase3_lle/run_phase3.check_preregistration) pointed at this
directory's PREREGISTRATION.md: committed, unmodified, pushed, and listing every
method; the runs are then exactly those the listed methods need.
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
from dataclasses import replace

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import final_pipeline as fp  # noqa: E402
from experiments.phase3_lle import run_phase3 as R3  # noqa: E402
from experiments.phase4_upo import config as C  # noqa: E402
from experiments.phase4_upo import detector as D  # noqa: E402
from experiments.phase4_upo import systems as S  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
RESULTS = HERE / "results"
PREREG = HERE / "PREREGISTRATION.md"
LLE_CFG = replace(fp.CFG, lle_chaos_test=True)


def parse_runs(spec):
    out = []
    for item in spec.split(","):
        m, M, s = item.split(":")
        out.append((m if m == "cao" else int(m), int(M), int(s)))
    return out


def tasks(phase, windows, tag, conds=None):
    out = []
    for n in windows:
        for system, snr in (conds or C.conditions()):
            for seed in C.seeds(phase, system, n, snr):
                lab = C.condition_label(system, snr)
                out.append({"task_id": f"{phase}:{tag}:{lab}:{n}:{seed}", "phase": phase, "system": system,
                            "snr_db": snr, "window_length": int(n), "seed": int(seed)})
    return out


def execute(args):
    task, runs = args
    t0 = time.perf_counter()
    x = S.generate(task["system"], task["window_length"], task["seed"], task["snr_db"])
    rec = {"task": task, "data_seed": S.data_seed(task["system"], task["window_length"], task["seed"]),
           "runs": {}, "runtime_s": {}}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for (m_rule, M, n_sur) in runs:
            key = D.run_key(m_rule, M, n_sur)
            t = time.perf_counter()
            try:
                rec["runs"][key] = D.run(x, m_rule, M, n_sur)
            except Exception as exc:                     # recorded; counts as "not detected"
                rec["runs"][key] = {"status": f"raised {type(exc).__name__}: {exc}", "peaks": []}
            rec["runtime_s"][key] = time.perf_counter() - t
        t = time.perf_counter()
        lc = fp.lle_chaos_test(x, config=LLE_CFG)
        rec["runtime_s"]["lle_chaos_test"] = time.perf_counter() - t
    rec["lle_chaos_test"] = {"detected": bool(lc["detected"]), "p": D._f(lc["p"]), "lle": D._f(lc["lle"]),
                             "status": lc["status"]}
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


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--phase", choices=("dev", "test"), required=True)
    ap.add_argument("--runs", default=None, help="comma list m:M:surrogates (dev only)")
    ap.add_argument("--methods", default=None, help="test: 'all' preregistered methods")
    ap.add_argument("--tag", default=None)
    ap.add_argument("--windows", default=",".join(map(str, C.WINDOW_LENGTHS)))
    ap.add_argument("--systems", default=None, help="comma list of condition labels to restrict to")
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args(argv)
    if a.phase == "test":
        from experiments.phase4_upo import methods as MM
        names = list(MM.METHODS) if a.methods in (None, "all") else a.methods.split(",")
        R3.PREREG = PREREG
        R3.check_preregistration(names)
        runs = sorted({MM.METHODS[n]["run"] for n in names}, key=str)
        tag = a.tag or "test"
    else:
        runs = parse_runs(a.runs)
        tag = a.tag or "dev"
    windows = tuple(int(v) for v in a.windows.split(","))
    conds = None
    if a.systems:
        want = set(a.systems.split(","))
        conds = [c for c in C.conditions() if C.condition_label(*c) in want]
    outdir = RESULTS / a.phase
    outdir.mkdir(parents=True, exist_ok=True)
    out_path = outdir / f"{tag}.jsonl"
    sha = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
           for p in [ROOT / "final_pipeline.py", *sorted(HERE.glob("*.py")), *([PREREG] if PREREG.exists() else [])]}
    with open(outdir / "manifest.jsonl", "a") as fh:
        fh.write(json.dumps({"started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                             "phase": a.phase, "tag": tag, "runs": [list(map(str, r)) for r in runs],
                             "git_head": R3._git("rev-parse", "HEAD").stdout.strip(),
                             "git_status": R3._git("status", "--porcelain").stdout, "sha256": sha,
                             "python": sys.version, "platform": platform.platform(),
                             "numpy": np.__version__}) + "\n")
    done = _done(out_path)
    todo = [t for t in tasks(a.phase, windows, tag, conds) if t["task_id"] not in done]
    todo.sort(key=lambda t: -t["window_length"])
    print(f"{len(todo)} windows to run; runs = {[D.run_key(*r) for r in runs]}", flush=True)
    t0 = time.time()
    with mp.get_context("fork").Pool(a.workers) as pool, open(out_path, "a") as fh:
        for k, rec in enumerate(pool.imap_unordered(execute, [(t, runs) for t in todo], chunksize=1), 1):
            fh.write(json.dumps(rec, allow_nan=False, default=str) + "\n")
            fh.flush()
            if k % 25 == 0 or k == len(todo):
                print(f"{k}/{len(todo)} done, {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
