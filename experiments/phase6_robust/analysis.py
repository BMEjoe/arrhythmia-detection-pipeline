"""
Phase 6 tables from runner JSONL files.

    python -m experiments.phase6_robust.analysis --phase test --file test

Per condition and method: AND (primary), LLE alone, UPO alone and OR with 95 %
Wilson intervals, errors (counted as not detected) and median runtime.
"""
from __future__ import annotations

import argparse
import json
import math
import pathlib

import numpy as np

from experiments.phase6_robust import config as C

HERE = pathlib.Path(__file__).resolve().parent
RULES = ("and", "lle", "upo", "or")


def load(phase, file):
    recs, seen = [], set()
    for line in open(HERE / "results" / phase / f"{file}.jsonl"):
        r = json.loads(line)
        key = (r["task"]["condition"], r["task"]["seed"])
        if key in seen:
            raise ValueError(f"duplicate window {key}")
        seen.add(key)
        recs.append(r)
    return recs


def wilson(k, n, z=1.959964):
    if n == 0:
        return (math.nan, math.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def rule(d, name):
    if d.get("error") is not None:
        return False
    return {"and": d["lle"] and d["upo"], "lle": d["lle"], "upo": d["upo"], "or": d["lle"] or d["upo"]}[name]


def summarize(recs, methods=None):
    methods = methods or list(dict.fromkeys(m for r in recs for m in r["methods"]))
    rows = []
    for cond in C.CONDITIONS:
        rs = [r for r in recs if r["task"]["condition"] == cond]
        if not rs:
            continue
        for m in methods:
            ds = [r["methods"][m] for r in rs if m in r["methods"]]
            if not ds:
                continue
            row = {"condition": cond, "group": C.group(cond), "method": m, "N": len(ds),
                   "seeds": sorted(r["task"]["seed"] for r in rs if m in r["methods"]),
                   "errors": sum(d.get("error") is not None for d in ds),
                   "runtime_median_s": float(np.median([r["runtime_s"][m] for r in rs if m in r["runtime_s"]]))}
            for name in RULES:
                k = sum(bool(rule(d, name)) for d in ds)
                row[f"{name}_k"] = k
                row[f"{name}_ci"] = wilson(k, len(ds))
            rows.append(row)
    return rows


def fmt(k, n, ci):
    return f"{k}/{n} ({100 * ci[0]:.0f}–{100 * ci[1]:.0f} %)"


def tables(rows, title):
    md = [f"# {title}", "", "Detections / windows (95 % Wilson interval). Errors: analyze_segment raised "
          "(counted as not detected).", ""]
    for m in dict.fromkeys(r["method"] for r in rows):
        md += [f"## {m}", "", "| condition | group | AND | LLE alone | UPO alone | OR | errors | median s |",
               "|---|---|---|---|---|---|---|---|"]
        for r in [r for r in rows if r["method"] == m]:
            md.append(f"| {r['condition']} | {r['group']} | " +
                      " | ".join(fmt(r[f'{n}_k'], r['N'], r[f'{n}_ci']) for n in RULES) +
                      f" | {r['errors']} | {r['runtime_median_s']:.1f} |")
        md.append("")
    methods = list(dict.fromkeys(r["method"] for r in rows))
    md += ["## AND detections by method", "", "| condition | N | " + " | ".join(methods) + " |",
           "|---|---|" + "---|" * len(methods)]
    for cond in dict.fromkeys(r["condition"] for r in rows):
        cell = {r["method"]: r for r in rows if r["condition"] == cond}
        n = next(iter(cell.values()))["N"]
        md.append(f"| {cond} | {n} | " + " | ".join(str(cell[m]["and_k"]) if m in cell else "–" for m in methods) + " |")
    return "\n".join(md) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", default="test", choices=("dev", "test"))
    ap.add_argument("--file", default="test")
    a = ap.parse_args(argv)
    rows = summarize(load(a.phase, a.file))
    tdir = HERE / "results" / a.phase / "tables"
    tdir.mkdir(parents=True, exist_ok=True)
    md = tables(rows, f"Phase 6 results ({a.phase} seeds, file {a.file})")
    (tdir / f"{a.file}_summary.md").write_text(md)
    keys = ["condition", "group", "method", "N", "errors", "and_k", "lle_k", "upo_k", "or_k", "runtime_median_s"]
    with open(tdir / f"{a.file}_summary.csv", "w") as fh:
        fh.write(",".join(keys) + "\n")
        for r in rows:
            fh.write(",".join(str(r[k]) for k in keys) + "\n")
    print(md)


if __name__ == "__main__":
    main()
