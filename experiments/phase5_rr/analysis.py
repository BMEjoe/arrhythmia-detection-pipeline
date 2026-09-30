"""
Phase 5 tables from runner JSONL files.

    python -m experiments.phase5_rr.analysis --phase test --file test

Per condition and m (2 = primary path through analyze_segment; 3, 4 = report
only): detections by AND (primary), the LLE component, the UPO component and
OR, with 95 % Wilson intervals; analyze_segment errors (counted as not
detected by every rule; the component decisions the direct calls would have
made are counted separately); median runtime.
"""
from __future__ import annotations

import argparse
import json
import math
import pathlib

import numpy as np

from experiments.phase5_rr import config as C

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
    return {"and": d["lle"] and d["upo"], "lle": d["lle"], "upo": d["upo"], "or": d["lle"] or d["upo"]}[name]


def summarize(recs, m_values=(2, 3, 4)):
    rows = []
    for cond in C.CONDITIONS:
        rs = [r for r in recs if r["task"]["condition"] == cond]
        if not rs:
            continue
        for m in m_values:
            ds = [r[f"m{m}"] for r in rs if f"m{m}" in r]
            if not ds:
                continue
            row = {"condition": cond, "group": C.group(cond), "m": m, "N": len(ds),
                   "errors": sum(d.get("error") is not None for d in ds),
                   "runtime_median_s": float(np.median([r["runtime_s"][f"m{m}"] for r in rs]))}
            for name in RULES:
                k = sum(bool(rule(d, name)) for d in ds)
                row[f"{name}_k"] = k
                row[f"{name}_ci"] = wilson(k, len(ds))
            errs = [d for d in ds if d.get("error") is not None]
            row["and_k_if_no_error"] = row["and_k"] + sum(d["lle_if_no_error"] and d["upo_if_no_error"] for d in errs)
            row["upo_status"] = {}
            for d in ds:
                row["upo_status"][d.get("upo_status")] = row["upo_status"].get(d.get("upo_status"), 0) + 1
            rows.append(row)
    return rows


def fmt(k, n, ci):
    return f"{k}/{n} ({100 * ci[0]:.0f}–{100 * ci[1]:.0f} %)"


def tables(rows, title):
    md = [f"# {title}", "",
          "Detections / windows (95 % Wilson interval). m = 2 is the frozen detector through "
          "`analyze_segment` (PRIMARY: AND). m = 3 and m = 4 set both `lle_chaos_test_m` and "
          "`upo_fixed_dimension` (report only). Errors: `analyze_segment` raised (counted as not detected).", ""]
    for m in sorted({r["m"] for r in rows}):
        md += [f"## m = {m}", "",
               "| condition | group | AND | LLE alone | UPO alone | OR | errors | median s |",
               "|---|---|---|---|---|---|---|---|"]
        for r in [r for r in rows if r["m"] == m]:
            md.append(f"| {r['condition']} | {r['group']} | " +
                      " | ".join(fmt(r[f'{n}_k'], r['N'], r[f'{n}_ci']) for n in RULES) +
                      f" | {r['errors']} | {r['runtime_median_s']:.1f} |")
        md.append("")
    md += ["## m-sensitivity (AND / LLE / UPO counts at m = 2, 3, 4)", "",
           "| condition | N | AND m2/m3/m4 | LLE m2/m3/m4 | UPO m2/m3/m4 |", "|---|---|---|---|---|"]
    for cond in dict.fromkeys(r["condition"] for r in rows):
        rs = {r["m"]: r for r in rows if r["condition"] == cond}
        if len(rs) < 2:
            continue
        cell = lambda n: "/".join(str(rs[m][f"{n}_k"]) if m in rs else "–" for m in (2, 3, 4))  # noqa: E731
        md.append(f"| {cond} | {rs[min(rs)]['N']} | {cell('and')} | {cell('lle')} | {cell('upo')} |")
    return "\n".join(md) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", default="test", choices=("dev", "test"))
    ap.add_argument("--file", default="test")
    a = ap.parse_args(argv)
    rows = summarize(load(a.phase, a.file))
    tdir = HERE / "results" / a.phase / "tables"
    tdir.mkdir(parents=True, exist_ok=True)
    md = tables(rows, f"Phase 5 results ({a.phase} seeds, file {a.file})")
    (tdir / f"{a.file}_summary.md").write_text(md)
    keys = ["condition", "group", "m", "N", "errors", "and_k", "lle_k", "upo_k", "or_k", "and_k_if_no_error",
            "runtime_median_s"]
    with open(tdir / f"{a.file}_summary.csv", "w") as fh:
        fh.write(",".join(keys) + ",upo_status\n")
        for r in rows:
            fh.write(",".join(str(r[k]) for k in keys) + ",\"" + json.dumps(r["upo_status"]).replace('"', "'") + "\"\n")
    print(md)


if __name__ == "__main__":
    main()
