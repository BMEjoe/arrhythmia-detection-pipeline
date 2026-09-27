"""
Phase 3 runner: one JSONL line per (phase, method, system, snr, window, seed).

    python -m experiments.phase3_lle.run_phase3 --phase dev --methods baseline
    python -m experiments.phase3_lle.run_phase3 --phase test --methods all

Results go to experiments/phase3_lle/results/<phase>/<method>.jsonl and are
checkpointed: finished task ids are skipped on restart.  BLAS is pinned to one
thread.  --phase test refuses to run unless PREREGISTRATION.md is committed,
unmodified in the working tree, pushed to the remote branch, and lists the
method (line "method: <name>").
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
import re
import subprocess
import sys
import time
import warnings

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import final_pipeline as fp  # noqa: E402
from experiments.phase2e import systems as S  # noqa: E402
from experiments.phase3_lle import config as C  # noqa: E402
from experiments.phase3_lle import estimators as E  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
RESULTS = HERE / "results"
PREREG = HERE / "PREREGISTRATION.md"


def _git(*a):
    return subprocess.run(["git", *a], cwd=ROOT, capture_output=True, text=True)


def check_preregistration(methods):
    """Test seeds may only be touched after PREREGISTRATION.md is committed and pushed."""
    if not PREREG.exists():
        raise SystemExit("refusing --phase test: PREREGISTRATION.md does not exist")
    rel = str(PREREG.relative_to(ROOT))
    if _git("ls-files", "--error-unmatch", rel).returncode != 0:
        raise SystemExit("refusing --phase test: PREREGISTRATION.md is not committed")
    if _git("diff", "--quiet", "HEAD", "--", rel).returncode != 0:
        raise SystemExit("refusing --phase test: PREREGISTRATION.md has uncommitted changes")
    branch = _git("rev-parse", "--abbrev-ref", "HEAD").stdout.strip()
    _git("fetch", "-q", "origin", branch)
    if not _git("log", "-1", "--format=%H", f"origin/{branch}", "--", rel).stdout.strip():
        raise SystemExit("refusing --phase test: PREREGISTRATION.md is not pushed")
    if _git("diff", "--quiet", f"origin/{branch}", "--", rel).returncode != 0:
        raise SystemExit("refusing --phase test: local PREREGISTRATION.md differs from the pushed one")
    text = PREREG.read_text()
    listed = set(re.findall(r"^\s*(?:[-*]\s+)?method:\s*`?([A-Za-z0-9_]+)`?", text, flags=re.M))
    missing = [m for m in methods if m not in listed]
    if missing:
        raise SystemExit(f"refusing --phase test: methods not preregistered: {missing}")


def tasks(phase, methods, windows=C.WINDOW_LENGTHS, conds=None):
    out = []
    for method in methods:
        for n in windows:
            for system, snr in (conds or C.conditions()):
                for seed in C.seeds(phase, system, n, snr):
                    label = C.condition_label(system, snr)
                    out.append({"task_id": f"{phase}:{method}:{label}:{n}:{seed}", "phase": phase,
                                "method": method, "system": system, "snr_db": snr,
                                "window_length": int(n), "seed": int(seed)})
    return out


def execute(task):
    t0 = time.perf_counter()
    x = S.generate(task["system"], task["window_length"], task["seed"], task["snr_db"])
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        try:
            res = E.METHODS[task["method"]](x)
        except Exception as exc:          # recorded, counted as "not detected"
            res = {"lle": None, "p": None, "detected": False,
                   "status": f"raised {type(exc).__name__}: {exc}"}
    res["n_runtime_warnings"] = len(w)
    res["data_seed"] = S.data_seed(task["system"], task["window_length"], task["seed"])
    return {"task": task, "result": res, "runtime_s": time.perf_counter() - t0}


def _done(path):
    if not path.exists():
        return set()
    ids = set()
    with open(path) as fh:
        for line in fh:
            try:
                ids.add(json.loads(line)["task"]["task_id"])
            except (json.JSONDecodeError, KeyError):
                pass            # a torn last line from an interrupted run
    return ids


def manifest(phase, methods):
    sha = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
           for p in [ROOT / "final_pipeline.py", *sorted(HERE.glob("*.py")), *([PREREG] if PREREG.exists() else [])]}
    return {"started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "phase": phase,
            "methods": methods, "git_head": _git("rev-parse", "HEAD").stdout.strip(),
            "git_status": _git("status", "--porcelain").stdout, "sha256": sha,
            "python": sys.version, "platform": platform.platform(), "numpy": np.__version__}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--phase", choices=("dev", "test"), required=True)
    ap.add_argument("--methods", default="baseline", help="comma list or 'all'")
    ap.add_argument("--windows", default=",".join(map(str, C.WINDOW_LENGTHS)))
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args(argv)
    methods = list(E.METHODS) if a.methods == "all" else a.methods.split(",")
    for m in methods:
        if m not in E.METHODS:
            raise SystemExit(f"unknown method {m}")
    if a.phase == "test":
        check_preregistration(methods)
    windows = tuple(int(v) for v in a.windows.split(","))
    outdir = RESULTS / a.phase
    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "manifest.jsonl", "a") as fh:
        fh.write(json.dumps(manifest(a.phase, methods)) + "\n")
    todo = []
    for m in methods:
        done = _done(outdir / f"{m}.jsonl")
        todo += [t for t in tasks(a.phase, [m], windows) if t["task_id"] not in done]
    # longest first for load balance
    todo.sort(key=lambda t: -t["window_length"])
    print(f"{len(todo)} tasks to run", flush=True)
    handles = {m: open(outdir / f"{m}.jsonl", "a") for m in methods}
    t0, k = time.time(), 0
    with mp.get_context("fork").Pool(a.workers) as pool:
        for rec in pool.imap_unordered(execute, todo, chunksize=1):
            fh = handles[rec["task"]["method"]]
            fh.write(json.dumps(rec, allow_nan=False, default=str) + "\n")
            fh.flush()
            k += 1
            if k % 50 == 0 or k == len(todo):
                print(f"{k}/{len(todo)} done, {time.time() - t0:.0f} s", flush=True)
    for fh in handles.values():
        fh.close()


if __name__ == "__main__":
    main()
