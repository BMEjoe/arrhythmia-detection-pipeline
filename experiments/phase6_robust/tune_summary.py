"""
Phase 6 C2: development tuning summary (results/dev/tuning.jsonl, dev seeds
only).  For every tuning configuration: AND on P3 (pooled), on the untrended
positives (each vs BASELINE-K), and on N1-N6 / S2 (each vs the 7 % limit at
the development N); the rule of PREREGISTRATION.md applied to development
windows, as a screen for choosing <= 3 candidates.

    python -m experiments.phase6_robust.tune_summary
"""
from __future__ import annotations

import math

import final_pipeline as fp
from experiments.phase6_robust import analysis as A
from experiments.phase6_robust import config as C
from experiments.phase6_robust import systems as S

GATES = (0.3, 0.5, 0.7)
GATED_FROM = ("linear", "sp_300", "movmed_101")


def add_gated(recs):
    """Gated variant '<m>_g<k>': the window is detrended only if fp.linear_trend_ratio(rr) >= k,
    so its decision equals method m's when the ratio passes and BASELINE-K's otherwise
    (exact: detrend_rr returns the input unchanged below the gate)."""
    for r in recs:
        ratio = fp.linear_trend_ratio(S.generate(r["task"]["condition"], r["task"]["seed"]))
        r["trend_ratio"] = ratio
        for m in GATED_FROM:
            if m not in r["methods"]:
                continue
            for k in GATES:
                src = m if ratio >= k else "baseline_k"
                r["methods"][f"{m}_g{k}"] = r["methods"][src]
                r["runtime_s"][f"{m}_g{k}"] = r["runtime_s"][src]
    return recs


def main():
    rows = A.summarize(add_gated(A.load("dev", "tuning")))
    cell = {(r["method"], r["condition"]): r for r in rows}
    methods = list(dict.fromkeys(r["method"] for r in rows))
    spec = list(C.PASS_NULLS) + list(C.S2_RATES)
    md = ["# Phase 6 C2 development tuning (development seeds; P*/G1 0-29, N1-N6/S2 0-49)", "",
          "AND detections. Screen = the preregistered rule applied to development windows "
          "(limits floor(0.07 N) on N1-N6 and S2; power loss <= floor(0.03 N) = 0 at N = 30 vs BASELINE-K; "
          "the power screen is therefore stricter than on test).", "",
          "| method | P3 Henon | P3 logistic | pooled P3 | " + " | ".join(C.POWER_CHECK) + " | "
          + " | ".join(spec) + " | max null+S2 AND | screen |",
          "|---|---|---|---|" + "---|" * (len(C.POWER_CHECK) + len(spec)) + "---|---|"]
    for m in methods:
        g = lambda c: cell.get((m, c))  # noqa: E731
        if any(g(c) is None for c in list(C.TREND_POSITIVES) + list(C.POWER_CHECK) + spec):
            continue
        p3 = [g(c)["and_k"] for c in C.TREND_POSITIVES]
        pw = []
        ok = True
        for c in C.POWER_CHECK:
            b = cell[("baseline_k", c)]["and_k"]
            k = g(c)["and_k"]
            lim = math.floor(0.03 * g(c)["N"] + 1e-9)
            ok &= b - k <= lim
            pw.append(f"{k}" + ("" if b - k <= lim else f" (base {b})"))
        sp = []
        for c in spec:
            k, n = g(c)["and_k"], g(c)["N"]
            ok &= k <= math.floor(0.07 * n + 1e-9)
            sp.append(f"{k}/{n}")
        mx = max(g(c)["and_k"] for c in spec)
        md.append(f"| {m} | {p3[0]}/30 | {p3[1]}/30 | {sum(p3)}/60 | " + " | ".join(pw) + " | " + " | ".join(sp)
                  + f" | {mx} | {'pass' if ok else 'fail'} |")
    out = A.HERE / "results" / "dev" / "tables"
    out.mkdir(parents=True, exist_ok=True)
    (out / "tuning.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
