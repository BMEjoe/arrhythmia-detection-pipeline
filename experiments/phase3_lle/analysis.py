"""
Phase 3 analysis: per method x window x condition tables.

    python -m experiments.phase3_lle.analysis --phase dev
    python -m experiments.phase3_lle.analysis --phase test

Rate = fraction of ALL windows with detected = True (failures count as "not
detected").  For null and periodic systems the rate is the false-positive
rate; for chaotic systems it is power.  95 % Wilson intervals.  Bias = median
LLE - reference (logistic ln 2, Henon 0.4192), over windows with a finite LLE.
Writes results/<phase>/tables/summary.{md,csv} and results/<phase>/tables/primary.md.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import pathlib

import numpy as np

from experiments.phase3_lle import config as C

HERE = pathlib.Path(__file__).resolve().parent
RESULTS = HERE / "results"


def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def frac(k, n):
    lo, hi = wilson(k, n)
    return f"{k}/{n} = {k / n:.3f} [{lo:.3f}, {hi:.3f}]" if n else "NA"


def load(phase, method):
    path = RESULTS / phase / f"{method}.jsonl"
    if not path.exists():
        return []
    return [json.loads(line) for line in open(path)]


def summarize(recs):
    """Group records by (window, condition) and compute the summary cells."""
    out = {}
    for r in recs:
        t, res = r["task"], r["result"]
        key = (t["window_length"], C.condition_label(t["system"], t["snr_db"]))
        out.setdefault(key, []).append(res)
    rows = []
    order = [C.condition_label(s, snr) for s, snr in C.conditions()]
    for (n, lab) in sorted(out, key=lambda k: (k[0], order.index(k[1]))):
        g = out[(n, lab)]
        N = len(g)
        k = sum(bool(x.get("detected")) for x in g)
        lle = [x["lle"] for x in g if x.get("lle") is not None and np.isfinite(x["lle"])]
        fail = sum(1 for x in g if x.get("status") != "ok")
        ref = C.LLE_REFERENCE.get(lab.split("@")[0])
        med = float(np.median(lle)) if lle else float("nan")
        q1, q3 = (np.percentile(lle, [25, 75]) if lle else (float("nan"), float("nan")))
        rows.append({
            "window": n, "condition": lab, "N": N, "detected_k": k, "rate": k / N if N else float("nan"),
            "rate_ci": frac(k, N), "failures": fail, "lle_n": len(lle),
            "lle_median": med, "lle_q1": float(q1), "lle_q3": float(q3),
            "bias": (med - ref) if ref is not None and lle else float("nan"),
            "abs_bias": abs(med - ref) if ref is not None and lle else float("nan"),
            "role": ("null" if lab in C.NULL_SYSTEMS else "periodic" if lab in C.PERIODIC_SYSTEMS
                     else "chaotic"),
        })
    return rows


def write(phase, methods):
    tdir = RESULTS / phase / "tables"
    tdir.mkdir(parents=True, exist_ok=True)
    allrows = []
    md = [f"# Phase 3 summary -- {phase} seeds", "",
          "Rate = detected / all windows (FPR for null/periodic, power for chaotic), 95 % Wilson CI. "
          "LLE median [Q1, Q3] over finite estimates; bias = median - reference.", ""]
    for m in methods:
        rows = summarize(load(phase, m))
        if not rows:
            continue
        for r in rows:
            r["method"] = m
        allrows += rows
        md += [f"## {m}", "", "| window | condition | role | rate | failures | LLE median [Q1, Q3] | bias |",
               "|---:|---|---|---|---:|---|---:|"]
        for r in rows:
            b = "" if not np.isfinite(r["bias"]) else f"{r['bias']:+.3f}"
            md.append(f"| {r['window']} | {r['condition']} | {r['role']} | {r['rate_ci']} | {r['failures']} | "
                      f"{r['lle_median']:.3f} [{r['lle_q1']:.3f}, {r['lle_q3']:.3f}] ({r['lle_n']}) | {b} |")
        md.append("")
    (tdir / "summary.md").write_text("\n".join(md))
    if allrows:
        with open(tdir / "summary.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(allrows[0].keys()))
            w.writeheader()
            w.writerows(allrows)
    return allrows


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", choices=("dev", "test"), required=True)
    ap.add_argument("--methods", default=None)
    a = ap.parse_args(argv)
    methods = (a.methods.split(",") if a.methods else
               sorted(p.stem for p in (RESULTS / a.phase).glob("*.jsonl")
                      if p.stem not in ("manifest", "diagnosis")))
    write(a.phase, methods)
    print((RESULTS / a.phase / "tables" / "summary.md").read_text())


if __name__ == "__main__":
    main()
