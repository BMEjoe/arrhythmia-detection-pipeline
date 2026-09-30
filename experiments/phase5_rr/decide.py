"""
Phase 5 PREREGISTERED primary rule (PREREGISTRATION.md Section 3), written
before any test window was generated.

    python -m experiments.phase5_rr.decide --phase test --file test

PASS  iff the AND detector (lle_chaos_test detected AND UPO
      instability_gate_detected, m = 2 through analyze_segment; a window where
      analyze_segment raised counts as not detected) fires on at most
      floor(0.07 N + 1e-9) windows (7/100) on EACH of N1-N6 (quantized, 256).
Power (AND on P1-P4 and G1-G3) is reported, not pass/fail.
"""
from __future__ import annotations

import argparse
import math

from experiments.phase5_rr import analysis as A
from experiments.phase5_rr import config as C

FP_FRACTION = 0.07


def evaluate(phase, file):
    rows = [r for r in A.summarize(A.load(phase, file), (2,)) if r["m"] == 2]
    cell = {r["condition"]: r for r in rows}
    checks = []
    for cond in C.PASS_CONDITIONS:
        r = cell[cond]
        limit = math.floor(FP_FRACTION * r["N"] + 1e-9)
        checks.append({"condition": cond, "k": r["and_k"], "N": r["N"], "limit": limit,
                       "ok": r["and_k"] <= limit, "errors": r["errors"],
                       "k_if_no_error": r["and_k_if_no_error"]})
    return checks, all(c["ok"] for c in checks), cell


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", default="test", choices=("dev", "test"))
    ap.add_argument("--file", default="test")
    a = ap.parse_args(argv)
    checks, passed, cell = evaluate(a.phase, a.file)
    md = [f"# Phase 5 primary decision ({a.phase} seeds, file {a.file})", "",
          "PASS iff the AND detector (m = 2, analyze_segment) fires on <= floor(0.07 N) windows on EACH of "
          "N1-N6 (quantized, 256 RR intervals). Errors count as not detected.", "",
          "| condition | AND detections | limit | within limit | analyze_segment errors | AND if errors were not counted |",
          "|---|---|---|---|---|---|"]
    for c in checks:
        md.append(f"| {c['condition']} | {c['k']}/{c['N']} | {c['limit']} | {'yes' if c['ok'] else '**no**'} | "
                  f"{c['errors']} | {c['k_if_no_error']} |")
    md += ["", f"**{'PASS' if passed else 'FAIL'}**", "",
           "Power (AND, reported, not pass/fail):", "", "| condition | AND | LLE alone | UPO alone | OR |",
           "|---|---|---|---|---|"]
    for cond, r in cell.items():
        if C.group(cond) in ("positive", "generalization"):
            md.append(f"| {cond} | {r['and_k']}/{r['N']} | {r['lle_k']} | {r['upo_k']} | {r['or_k']} |")
    tdir = A.HERE / "results" / a.phase / "tables"
    tdir.mkdir(parents=True, exist_ok=True)
    (tdir / f"{a.file}_primary.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
