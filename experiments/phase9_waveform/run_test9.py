"""
Phase 9 Part E: the preregistered TEST run.

    python -m experiments.phase9_waveform.run_test9 --workers 4

Refuses to start unless PREREGISTRATION.md is committed, unmodified, pushed and lists every
candidate (Phase 3 guard), and unless every threshold in candidates9.py is set.  TEST seeds start at
9500 (never used before); TEST windows use nstdb signal 1.  Each window is generated once and
evaluated by every candidate.  Output: results/test/test9.jsonl (one line per window; resumable).

Groups (PREREGISTRATION.md Section 3), all windows 512 intervals:
  null_input    17 PASS nulls, true beats, no recording noise (labels = generator beat types),
                seeds 9500-9599
  null_real     17 PASS nulls, ECGSYN + nstdb mix 12 dB + Pan-Tompkins beats + Phase 7 V jitter
                (labels = annotations of matched beats, 'X' unmatched), seeds 9500-9599
  noncha_input  every Phase 8 NON-CHAOTIC regime (15 TEST + 3 DEV), input level, no ectopy,
                seeds 9500-9599
  noncha_real   every Phase 8 NON-CHAOTIC regime, realistic, ectopy cycling none / S2_5 / E1 / E3_10
                by seed (seed mod 4), seeds 9500-9599; every KTz NON-CHAOTIC (P, input) regime,
                realistic (its input is part of the regime), seeds 9500-9599
  chaos_primary TEST CHAOTIC regimes at the most realistic variant: Phase 8 TEST CHAOTIC x
                (S2_5, E1, E3_10), and KTz CHAOTIC (P, input) regimes, realistic; seeds 9500-9529
  chaos_second  TEST CHAOTIC regimes at the input level without ectopy (Phase 8) / with their input
                (KTz), seeds 9500-9529; DEV Mackey-Glass CHAOTIC x (S2_5, E1, E3_10) realistic,
                seeds 9500-9509
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import sys
import time
import zlib
from multiprocessing import Pool

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from experiments.phase3_lle import run_phase3 as R3  # noqa: E402
from experiments.phase9_waveform import candidates9 as C  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
PREREG = HERE / "PREREGISTRATION.md"
N = 512
S100 = range(9500, 9600)
S30 = range(9500, 9530)
S10 = range(9500, 9510)
ECT_CYCLE = ("none", "S2_5", "E1", "E3_10")


def tasks():
    from experiments.phase9_waveform import regimes9 as RG
    T = []
    src = RG.sources()
    for s in src:
        if s["source"] == "null":
            for seed in S100:
                T.append({"group": "null_input", "source": "null", "name": s["name"], "label": "NULL",
                          "ectopy": "none", "level": "input", "seed": seed})
                T.append({"group": "null_real", "source": "null", "name": s["name"], "label": "NULL",
                          "ectopy": "none", "level": "real", "seed": seed})
        elif s["source"] == "model":
            if s["label"] == "NON-CHAOTIC":
                for seed in S100:
                    T.append({"group": "noncha_input", "source": "model", "name": s["name"], "label": s["label"],
                              "family": s["family"], "ectopy": "none", "level": "input", "seed": seed})
                    T.append({"group": "noncha_real", "source": "model", "name": s["name"], "label": s["label"],
                              "family": s["family"], "ectopy": ECT_CYCLE[seed % 4], "level": "real", "seed": seed})
            elif s["split"] == "TEST":
                for e in ("S2_5", "E1", "E3_10"):
                    for seed in S30:
                        T.append({"group": "chaos_primary", "source": "model", "name": s["name"], "label": s["label"],
                                  "family": s["family"], "ectopy": e, "level": "real", "seed": seed})
                for seed in S30:
                    T.append({"group": "chaos_second", "source": "model", "name": s["name"], "label": s["label"],
                              "family": s["family"], "ectopy": "none", "level": "input", "seed": seed})
            else:   # DEV Mackey-Glass CHAOTIC (secondary)
                for e in ("S2_5", "E1", "E3_10"):
                    for seed in S10:
                        T.append({"group": "chaos_second", "source": "model", "name": s["name"], "label": s["label"],
                                  "family": s["family"], "ectopy": e, "level": "real", "seed": seed})
        elif s["source"] == "ktz":
            base = {"source": "ktz", "name": s["name"], "label": s["label"], "family": "ktz",
                    "ectopy": s["ectopy_fixed"]}
            if s["label"] == "NON-CHAOTIC":
                T += [{**base, "group": "noncha_real", "level": "real", "seed": seed} for seed in S100]
            else:
                T += [{**base, "group": "chaos_primary", "level": "real", "seed": seed} for seed in S30]
                T += [{**base, "group": "chaos_second", "level": "input", "seed": seed} for seed in S30]
    for t in T:
        t["task_id"] = f"{t['group']}|{t['name']}|{t['ectopy']}|{t['level']}|{t['seed']}"
    return T


def window(t):
    from experiments.phase9_waveform import library as L
    if t["level"] == "input":
        w = L.build(t["source"], t["name"], t["ectopy"], N, t["seed"], "none", "TEST", "true", False)
    else:
        w = L.build(t["source"], t["name"], t["ectopy"], N, t["seed"], "mix12", "TEST", "detected", True)
    if "error" in w:
        raise RuntimeError(w["error"])
    return w


def execute(t):
    t0 = time.time()
    out = dict(t)
    try:
        w = window(t)
        x = np.asarray(w["rr"], float)
        types = np.asarray(w["beat_types"])
        out.update({k: w[k] for k in ("n_missed", "n_extra") if k in w})
        out["n_nonN"] = int(np.sum(types != "N"))
        out["rr_mean"], out["rr_sd"] = float(np.mean(x)), float(np.std(x))
        out["gen_s"] = time.time() - t0
        res = C.run_all(x, types)
        for r in res.values():
            r.pop("traceback", None)
        out["results"] = res
    except Exception as exc:                                    # noqa: BLE001
        out["error"] = f"{type(exc).__name__}: {exc}"
    out["t_total"] = time.time() - t0
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args(argv)
    R3.PREREG = PREREG
    R3.check_preregistration(list(C.CANDIDATES))
    assert None not in (C.Z_DET, C.Z_DET_RR, C.G_MIN, C.G_MIN_RR, C.F_MIN), "thresholds not set"
    out = HERE / "results" / "test"
    out.mkdir(parents=True, exist_ok=True)
    path = out / "test9.jsonl"
    done = set()
    if path.exists():
        for line in open(path):
            try:
                done.add(json.loads(line)["task_id"])
            except Exception:                                   # noqa: BLE001
                pass
    T = [t for t in tasks() if t["task_id"] not in done]
    T.sort(key=lambda t: zlib.crc32(t["task_id"].encode()))      # interleave groups
    print(f"{len(T)} to run ({len(done)} done)", flush=True)
    t0 = time.time()
    with open(path, "a") as fh, Pool(a.workers, maxtasksperchild=50) as pool:
        for i, r in enumerate(pool.imap_unordered(execute, T, chunksize=1), 1):
            fh.write(json.dumps(r, default=lambda o: o.item() if hasattr(o, "item") else str(o)) + "\n")
            fh.flush()
            if i % 200 == 0 or i == len(T):
                el = time.time() - t0
                print(f"{i}/{len(T)} {el:.0f}s eta {el / i * (len(T) - i):.0f}s", flush=True)


if __name__ == "__main__":
    main()
