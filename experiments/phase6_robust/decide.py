"""
Phase 6 Part C PREREGISTERED rule (PREREGISTRATION.md Section 4), written before
any test window was generated.

    python -m experiments.phase6_robust.decide --phase test --file test

A candidate PASSES only if, on the test windows,
  (1) its AND detector fires on at most floor(0.07 N + 1e-9) windows on EACH of
      N1-N6 and on EACH S2 rate (2, 5, 10 %), and
  (2) its AND detection on EACH of P1, P2, P4 (Henon / logistic at 30 / 20 dB)
      and G1 is at most floor(0.03 N + 1e-9) (= 6/200) windows below
      BASELINE-K on the same windows.
WINNER: among passing candidates, the highest pooled AND detection count on P3
  (Henon + trend, logistic + trend).
Tie-break: fewer pooled AND detections on the Part B ectopy conditions (all
  E* raw and annotation-edited conditions), then the lower candidate index.
A window where analyze_segment raised counts as not detected.  BASELINE-K is
evaluated on the same windows and is not eligible.
"""
from __future__ import annotations

import argparse
import math

from experiments.phase6_robust import analysis as A
from experiments.phase6_robust import config as C
from experiments.phase6_robust import methods as MM

FP_FRACTION = 0.07
POWER_LOSS_FRACTION = 0.03


def evaluate(phase, file):
    recs = A.load(phase, file)
    rows = A.summarize(recs, list(MM.PREREGISTERED))
    cell = {(r["method"], r["condition"]): r for r in rows}
    base = "baseline_k"
    out = []
    for m, spec in MM.PREREGISTERED.items():
        spec_checks = []
        for cond in list(C.PASS_NULLS) + list(C.S2_RATES):
            r = cell[(m, cond)]
            lim = math.floor(FP_FRACTION * r["N"] + 1e-9)
            spec_checks.append({"condition": cond, "k": r["and_k"], "N": r["N"], "limit": lim, "ok": r["and_k"] <= lim})
        power_checks = []
        for cond in C.POWER_CHECK:
            r, b = cell[(m, cond)], cell[(base, cond)]
            if r["seeds"] != b["seeds"]:
                raise ValueError(f"{m} and {base} not evaluated on the same windows for {cond}")
            lim = math.floor(POWER_LOSS_FRACTION * r["N"] + 1e-9)
            power_checks.append({"condition": cond, "k": r["and_k"], "base": b["and_k"], "N": r["N"],
                                 "max_loss": lim, "ok": b["and_k"] - r["and_k"] <= lim})
        p3 = sum(cell[(m, c)]["and_k"] for c in C.TREND_POSITIVES)
        p3n = sum(cell[(m, c)]["N"] for c in C.TREND_POSITIVES)
        ect = sum(cell[(m, c)]["and_k"] for c in C.PARTB_CONDITIONS)
        ectn = sum(cell[(m, c)]["N"] for c in C.PARTB_CONDITIONS)
        out.append({"method": m, "index": spec["index"], "eligible": spec["eligible"],
                    "spec": spec_checks, "power": power_checks,
                    "passes": all(c["ok"] for c in spec_checks + power_checks),
                    "p3": p3, "p3_n": p3n, "ectopy": ect, "ectopy_n": ectn})
    ranked = sorted([r for r in out if r["eligible"] and r["passes"]],
                    key=lambda r: (-r["p3"], r["ectopy"], r["index"]))
    return out, ranked, (ranked[0]["method"] if ranked else None)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", default="test", choices=("dev", "test"))
    ap.add_argument("--file", default="test")
    a = ap.parse_args(argv)
    out, ranked, winner = evaluate(a.phase, a.file)
    spec_c = list(C.PASS_NULLS) + list(C.S2_RATES)
    md = [f"# Phase 6 Part C primary decision ({a.phase} seeds, file {a.file})", "",
          "PASS: AND <= floor(0.07 N) on EACH of N1-N6 and each S2 rate, AND power on each of P1, P2, P4 (4) "
          "and G1 at most floor(0.03 N) below BASELINE-K on the same windows. WINNER: most pooled AND on P3; "
          "tie-break fewer pooled AND on the Part B ectopy conditions, then candidate index.", "",
          "## Specificity (AND detections; limit)", "",
          "| method | " + " | ".join(spec_c) + " |", "|---|" + "---|" * len(spec_c)]
    for r in out:
        md.append(f"| {r['method']} | " + " | ".join(
            f"{c['k']}/{c['N']} (≤{c['limit']}){'' if c['ok'] else ' **FAIL**'}" for c in r["spec"]) + " |")
    md += ["", "## Power on untrended positives (AND; BASELINE-K; max loss)", "",
           "| method | " + " | ".join(C.POWER_CHECK) + " |", "|---|" + "---|" * len(C.POWER_CHECK)]
    for r in out:
        md.append(f"| {r['method']} | " + " | ".join(
            f"{c['k']} (base {c['base']}, ≥{c['base'] - c['max_loss']}){'' if c['ok'] else ' **FAIL**'}"
            for c in r["power"]) + " |")
    md += ["", "## Summary", "", "| method | PASS | pooled P3 AND | pooled Part B ectopy AND |", "|---|---|---|---|"]
    for r in out:
        md.append(f"| {r['method']} | {'yes' if r['passes'] else 'no'}{'' if r['eligible'] else ' (not eligible)'} | "
                  f"{r['p3']}/{r['p3_n']} | {r['ectopy']}/{r['ectopy_n']} |")
    md += ["", "Ranking of passing candidates: " + (", ".join(r["method"] for r in ranked) or "none"), "",
           f"**WINNER: {winner}**" if winner else "**No candidate passes; there is no winner.**"]
    tdir = A.HERE / "results" / a.phase / "tables"
    tdir.mkdir(parents=True, exist_ok=True)
    (tdir / f"{a.file}_primary.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
