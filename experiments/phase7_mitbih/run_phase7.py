"""
Phase 7 runner: the frozen detector (detector.DETECTOR) on every preregistered
MIT-BIH window.

    python -m experiments.phase7_mitbih.run_phase7 --phase test [--workers 4]
    python -m experiments.phase7_mitbih.run_phase7 --phase synthetic   # code check, no MIT-BIH

--phase test reuses the Phase 3 guard (experiments.phase3_lle.run_phase3.
check_preregistration): it refuses to start unless PREREGISTRATION.md is
committed, unmodified, pushed and lists `method: phase6_recommended_and`.

Arms (PREREGISTRATION.md B):
  raw         PRIMARY: raw RR of the Pan-Tompkins peaks, 150 ms labels (305 windows)
  edited      SECONDARY: Phase 6 annotation edit of the raw window, NN fraction >= 0.80
  annotation  SENSITIVITY: RR of the annotation beat times (368 windows)
The window lists are rebuilt from data.py and must equal results/qc/*.csv.
Results: results/run/{raw,edited,annotation}.jsonl (one line per window;
resumable: finished task ids are skipped) and results/run/manifest.jsonl.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import pathlib
import platform
import subprocess
import sys
import time
from multiprocessing import Pool

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from experiments.phase3_lle import run_phase3 as R3  # noqa: E402
from experiments.phase7_mitbih import data as D  # noqa: E402
from experiments.phase7_mitbih import detector as DT  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
PREREG = HERE / "PREREGISTRATION.md"
QC = HERE / "results" / "qc"
ARMS = ("raw", "edited", "annotation")
META = ("record", "subject", "window", "start_rr", "t0", "t1", "label", "abnormal_fraction", "n_abnormal",
        "n_v_type", "n_s_type", "n_q", "beat_type", "nn_fraction", "edited_eligible",
        "n_missed_beats_in_span", "start_beat_unusable", "af_fraction")


def _clean(v):
    if isinstance(v, (np.bool_,)):
        return bool(v)
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating, float)):
        return None if not np.isfinite(v) else float(v)
    return v


def build_tasks(arms=ARMS):
    """Every preregistered window as (task_id, arm, meta, rr).  Cross-checked against QC."""
    qc = {"detected": pd.read_csv(QC / "windows.csv", dtype={"record": str}),
          "annotation": pd.read_csv(QC / "windows_annotation_times.csv", dtype={"record": str})}
    tasks = []
    for source in ("detected", "annotation"):
        want = [a for a in arms if (a == "annotation") == (source == "annotation")]
        if not want:
            continue
        ref = qc[source].set_index(["record", "window"])
        n = 0
        for r in D.RECORDS:
            series = D.record_series(r, cache_dir=os.path.join(D.DATA_DIR, "_phase7_peaks"), source=source)
            for w in D.windows(series):
                q = ref.loc[(r, w["window"])]
                assert (int(q.t0), int(q.t1), int(q.label)) == (w["t0"], w["t1"], w["label"]), (r, w["window"])
                n += 1
                meta = {k: _clean(q[k]) if k in q.index else _clean(w.get(k)) for k in META
                        if k not in ("record", "window")}
                meta.update({"record": r, "window": int(w["window"]),
                             "retained_tol75": _clean(q["retained_tol75"]) if "retained_tol75" in q.index else None,
                             "label_tol75": _clean(q["label_tol75"]) if "label_tol75" in q.index else None})
                for arm in want:
                    if arm == "edited":
                        if not w["nn_fraction"] >= 0.80:
                            continue
                        rr = w["rr_edited"]
                    else:
                        rr = w["rr"]
                    tasks.append((f"{arm}:{r}:{w['window']}", arm, meta, np.asarray(rr, dtype=float)))
        assert n == len(qc[source]), f"{source}: {n} windows rebuilt, {len(qc[source])} in QC"
    return tasks


def synthetic_tasks():
    """Code check without MIT-BIH: Phase 5 synthetic windows with fake labels and subjects."""
    from experiments.phase5_rr import systems as S5
    tasks = []
    for i, (cond, lab) in enumerate([("N1_linear_rr", 0), ("N2_power_law", 0), ("S2_ectopic_10pct", 1),
                                     ("S2_ectopic_5pct", 1), ("P1_henon_rr", 1), ("N6_noisy_rsa", 0)]):
        for seed in range(4):
            meta = {"record": f"syn{i}", "subject": f"SYN{seed}", "window": seed, "label": lab,
                    "abnormal_fraction": 0.2 * lab, "beat_type": "ventricular" if lab else "normal",
                    "nn_fraction": 1.0 - 0.3 * lab, "af_fraction": 0.0, "n_missed_beats_in_span": 0,
                    "retained_tol75": True, "label_tol75": lab}
            for arm in ARMS:
                tasks.append((f"{arm}:syn{i}:{seed}", arm, meta, S5.generate(cond, seed)))
    return tasks


def execute(task):
    tid, arm, meta, rr = task
    rec = DT.evaluate(rr)
    return {"task_id": tid, "arm": arm, **meta, **D.hrv(rr), "n_rr": int(len(rr)),
            "rr_sha256": hashlib.sha256(np.ascontiguousarray(rr).tobytes()).hexdigest(), **rec}


def _done(path):
    ids = set()
    if path.exists():
        with open(path) as fh:
            for line in fh:
                try:
                    ids.add(json.loads(line)["task_id"])
                except (json.JSONDecodeError, KeyError):
                    pass
    return ids


def _json(o):
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, float) and not np.isfinite(o):
        return None
    raise TypeError(type(o))


def _finite(o):
    """Replace non-finite floats (JSON has no NaN/inf) by None, recursively."""
    if isinstance(o, dict):
        return {k: _finite(v) for k, v in o.items()}
    if isinstance(o, list):
        return [_finite(v) for v in o]
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


def manifest(phase):
    sha = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
           for p in [ROOT / "final_pipeline.py", *sorted(HERE.glob("*.py")), *([PREREG] if PREREG.exists() else [])]}
    git = lambda *a: subprocess.run(["git", *a], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    import scipy
    import sklearn
    import wfdb
    return {"started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "phase": phase,
            "method": DT.METHOD, "git_head": git("rev-parse", "HEAD"), "git_status": git("status", "--porcelain"),
            "sha256": sha, "python": sys.version, "platform": platform.platform(), "numpy": np.__version__,
            "scipy": scipy.__version__, "sklearn": sklearn.__version__, "wfdb": wfdb.__version__,
            "omp_num_threads": os.environ.get("OMP_NUM_THREADS")}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--phase", choices=("test", "synthetic"), required=True)
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args(argv)
    if a.phase == "test":
        R3.PREREG = PREREG
        R3.check_preregistration([DT.METHOD])
        tasks = build_tasks()
        outdir = HERE / "results" / "run"
    else:
        tasks = synthetic_tasks()
        outdir = HERE / "results" / "synthetic_check"
    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "manifest.jsonl", "a") as fh:
        fh.write(json.dumps({**manifest(a.phase), "n_tasks": len(tasks)}) + "\n")
    done = set().union(*[_done(outdir / f"{arm}.jsonl") for arm in ARMS])
    todo = [t for t in tasks if t[0] not in done]
    # longest-first is unknown; interleave arms so every file grows
    print(f"{len(tasks)} tasks, {len(done)} done, {len(todo)} to run", flush=True)
    t0 = time.time()
    fhs = {arm: open(outdir / f"{arm}.jsonl", "a") for arm in ARMS}
    try:
        with Pool(a.workers, maxtasksperchild=50) as pool:
            for i, res in enumerate(pool.imap_unordered(execute, todo, chunksize=1), 1):
                fhs[res["arm"]].write(json.dumps(_finite(res), default=_json) + "\n")
                fhs[res["arm"]].flush()
                if i % 25 == 0 or i == len(todo):
                    el = time.time() - t0
                    print(f"{i}/{len(todo)} {el:.0f}s eta {el / i * (len(todo) - i):.0f}s", flush=True)
    finally:
        for fh in fhs.values():
            fh.close()


if __name__ == "__main__":
    main()
