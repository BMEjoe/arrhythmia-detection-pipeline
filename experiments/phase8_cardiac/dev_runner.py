"""
Phase 8 Part B development runner: any detector method on DEVELOPMENT material only.

DEVELOPMENT material = DEV model regimes (results/ground_truth/regimes.json, split DEV),
the Phase 5-6 null and positive generators (nulls.py), development seeds 0-999.
It refuses TEST-split regimes.

    python -m experiments.phase8_cardiac.dev_runner --method frozen --lengths 256 --seeds 0-4 --tag base
Results: results/dev/<tag>.jsonl (resumable).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import pathlib
import time
from multiprocessing import Pool

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

from experiments.phase8_cardiac import methods as MM  # noqa: E402
from experiments.phase8_cardiac import nulls as NL  # noqa: E402
from experiments.phase8_cardiac import series as SR  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
DEV_VARIANTS = ("i_clean", "iii_meas20", "iv_q", "v_S2_5", "v_E1", "v_E3_10", "vi_S2_5", "vi_E1", "vi_E3_10")


def tasks(methods, lengths, seeds, variants, conditions, regimes_filter=None):
    regs = [r for r in SR.load_regimes() if r["split"] == "DEV"]
    out = []
    for n in lengths:
        for s in seeds:
            for r in regs:
                if regimes_filter and r["name"] not in regimes_filter:
                    continue
                for v in variants:
                    for m in methods:
                        out.append({"kind": "model", "family": r["family"], "regime_idx": r["regime_idx"],
                                    "regime": r["name"], "label": r["label"], "params": r["params"], "n": n,
                                    "variant": v, "seed": s, "method": m})
            for c in conditions:
                for m in methods:
                    out.append({"kind": "null" if c in NL.PASS_NULLS else "synthetic", "regime": c,
                                "label": "NON-CHAOTIC" if c in NL.PASS_NULLS else "secondary", "n": n,
                                "variant": "native", "seed": s, "method": m})
    for t in out:
        t["task_id"] = f"{t['method']}|{t['regime']}|{t['variant']}|{t['n']}|{t['seed']}"
    return out


def execute(t):
    if t["kind"] == "model":
        assert t["regime_idx"] in {r["regime_idx"] for r in SR.load_regimes() if r["split"] == "DEV"}
        rr = SR.window(t["family"], t["regime_idx"], t["params"], t["n"], t["variant"], t["seed"])
    else:
        rr = NL.generate(t["regime"], t["seed"], t["n"])
    t0 = time.time()
    rec = MM.run(t["method"], rr)
    return {**{k: v for k, v in t.items() if k != "params"}, **rec, "runtime_s": time.time() - t0,
            "rr_sha": hashlib.sha256(np.ascontiguousarray(rr).tobytes()).hexdigest()[:16]}


def _range(s):
    a, b = s.split("-")
    return range(int(a), int(b) + 1)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--method", required=True)
    ap.add_argument("--lengths", default="256")
    ap.add_argument("--seeds", default="0-4")
    ap.add_argument("--variants", default=",".join(DEV_VARIANTS))
    ap.add_argument("--conditions", default=",".join(NL.PASS_NULLS + ("G2_rossler_flow", "G3_mackey_glass", "P1_henon_rr")))
    ap.add_argument("--regimes", default="")
    ap.add_argument("--tag", required=True)
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args(argv)
    seeds = _range(a.seeds)
    assert max(seeds) < 1000, "development seeds are 0-999"
    T = tasks(a.method.split(","), [int(x) for x in a.lengths.split(",")], seeds,
              [v for v in a.variants.split(",") if v], [c for c in a.conditions.split(",") if c],
              set(a.regimes.split(",")) if a.regimes else None)
    out = HERE / "results" / "dev" / f"{a.tag}.jsonl"
    out.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if out.exists():
        for line in open(out):
            try:
                done.add(json.loads(line)["task_id"])
            except Exception:
                pass
    todo = [t for t in T if t["task_id"] not in done]
    print(len(T), "tasks", len(todo), "to run", flush=True)
    t0 = time.time()
    with open(out, "a") as fh, Pool(a.workers, maxtasksperchild=40) as pool:
        for i, r in enumerate(pool.imap_unordered(execute, todo, chunksize=1), 1):
            fh.write(json.dumps(r, default=lambda o: o.item() if hasattr(o, "item") else str(o)) + "\n")
            fh.flush()
            if i % 50 == 0 or i == len(todo):
                print(f"{i}/{len(todo)} {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
