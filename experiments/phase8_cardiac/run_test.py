"""
Phase 8 Part D: the preregistered TEST run.

    python -m experiments.phase8_cardiac.run_test --workers 4

Reuses the Phase 3 guard (experiments.phase3_lle.run_phase3.check_preregistration): refuses
to start unless PREREGISTRATION.md is committed, unmodified, pushed, and lists every method
(`method: <name>` lines).  TEST seeds 5000-5099 (never used before).  Each window is
generated once per length and evaluated by every method of that length.

Groups (PREREGISTRATION.md Section 3):
  null          17 Phase 5-6 PASS nulls x seeds 5000-5099 (native, quantized)
  noncha        every NON-CHAOTIC regime (TEST and DEV), variant iv_q, seeds 5000-5099
  chaos_primary TEST CHAOTIC regimes x (iv_q, vi_S2_5, vi_E1, vi_E3_10) x seeds 5000-5029
  chaos_second  TEST CHAOTIC regimes x (i_clean, ii_dyn_lo, ii_dyn_hi, iii_meas30, iii_meas20,
                v_S2_5, v_E1, v_E3_10) x seeds 5000-5009
  dev_chaos     DEV CHAOTIC regimes x (iv_q, vi_S2_5, vi_E1, vi_E3_10) x seeds 5000-5009
  noncha_ect    every NON-CHAOTIC regime x (vi_S2_5, vi_E1, vi_E3_10) x seeds 5000-5009
  flows         Phase 5 G2_rossler_flow, G3_mackey_glass x seeds 5000-5099
Results: results/test/test_<length>.jsonl (one line per window, all methods of that length).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import pathlib
import sys
import time
from multiprocessing import Pool

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from experiments.phase3_lle import run_phase3 as R3  # noqa: E402
from experiments.phase8_cardiac import methods as MM  # noqa: E402
from experiments.phase8_cardiac import nulls as NL  # noqa: E402
from experiments.phase8_cardiac import series as SR  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
PREREG = HERE / "PREREGISTRATION.md"
METHODS_BY_LENGTH = {256: ("c3_nlp_ep256", "c4_titration256", "baseline_frozen256"),
                     512: ("c1_frozen512", "c2_nlp_iaaft512")}
S100 = range(5000, 5100)
S30 = range(5000, 5030)
S10 = range(5000, 5010)


def tasks(n):
    regs = SR.load_regimes()
    T = []
    for c in NL.PASS_NULLS:
        for s in S100:
            T.append({"group": "null", "regime": c, "label": "NULL", "variant": "native", "seed": s})
    for c in ("G2_rossler_flow", "G3_mackey_glass"):
        for s in S100:
            T.append({"group": "flows", "regime": c, "label": "SECONDARY", "variant": "native", "seed": s})
    for r in regs:
        base = {"regime": r["name"], "label": r["label"], "family": r["family"], "split": r["split"],
                "regime_idx": r["regime_idx"], "params": r["params"]}
        if r["label"] == "NON-CHAOTIC":
            T += [{**base, "group": "noncha", "variant": "iv_q", "seed": s} for s in S100]
            T += [{**base, "group": "noncha_ect", "variant": v, "seed": s}
                  for v in ("vi_S2_5", "vi_E1", "vi_E3_10") for s in S10]
        elif r["split"] == "TEST":
            T += [{**base, "group": "chaos_primary", "variant": v, "seed": s}
                  for v in ("iv_q", "vi_S2_5", "vi_E1", "vi_E3_10") for s in S30]
            T += [{**base, "group": "chaos_second", "variant": v, "seed": s}
                  for v in ("i_clean", "ii_dyn_lo", "ii_dyn_hi", "iii_meas30", "iii_meas20", "v_S2_5", "v_E1",
                            "v_E3_10") for s in S10]
        else:
            T += [{**base, "group": "dev_chaos", "variant": v, "seed": s}
                  for v in ("iv_q", "vi_S2_5", "vi_E1", "vi_E3_10") for s in S10]
    for t in T:
        t["n"] = n
        t["task_id"] = f"{t['group']}|{t['regime']}|{t['variant']}|{n}|{t['seed']}"
    return T


def execute(t):
    t0 = time.time()
    if t["group"] in ("null", "flows"):
        rr = NL.generate(t["regime"], t["seed"], t["n"])
    else:
        rr = SR.window(t["family"], t["regime_idx"], t["params"], t["n"], t["variant"], t["seed"])
    gen_s = time.time() - t0
    res = {}
    for m in METHODS_BY_LENGTH[t["n"]]:
        t1 = time.time()
        r = MM.run(m, rr)
        r.pop("traceback", None)
        r["runtime_s"] = time.time() - t1
        res[m] = r
    return {**{k: v for k, v in t.items() if k != "params"}, "gen_s": gen_s,
            "rr_sha": hashlib.sha256(np.ascontiguousarray(rr).tobytes()).hexdigest()[:16],
            "rr_mean": float(np.mean(rr)), "rr_sd": float(np.std(rr)), "results": res}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--lengths", default="256,512")
    a = ap.parse_args(argv)
    R3.PREREG = PREREG
    R3.check_preregistration(list(METHODS_BY_LENGTH[256]) + list(METHODS_BY_LENGTH[512]))
    assert MM.Z_C2 is not None and MM.Z_C3 is not None
    outdir = HERE / "results" / "test"
    outdir.mkdir(parents=True, exist_ok=True)
    for n in [int(x) for x in a.lengths.split(",")]:
        path = outdir / f"test_{n}.jsonl"
        done = set()
        if path.exists():
            for line in open(path):
                try:
                    done.add(json.loads(line)["task_id"])
                except Exception:
                    pass
        T = [t for t in tasks(n) if t["task_id"] not in done]
        # interleave heavy (coupled vdP) tasks
        T.sort(key=lambda t: (hash(t["task_id"]) % 997))
        print(f"length {n}: {len(T)} to run ({len(done)} done)", flush=True)
        t0 = time.time()
        with open(path, "a") as fh, Pool(a.workers, maxtasksperchild=50) as pool:
            for i, r in enumerate(pool.imap_unordered(execute, T, chunksize=1), 1):
                fh.write(json.dumps(r, default=lambda o: o.item() if hasattr(o, "item") else str(o)) + "\n")
                fh.flush()
                if i % 100 == 0 or i == len(T):
                    el = time.time() - t0
                    print(f"{n}: {i}/{len(T)} {el:.0f}s eta {el / i * (len(T) - i):.0f}s", flush=True)


if __name__ == "__main__":
    main()
