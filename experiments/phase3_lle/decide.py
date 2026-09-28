"""
Phase 3, B4 -- apply the PREREGISTERED primary decision rule (PREREGISTRATION.md
Section 3) to the test-seed results.  Written before any test result was seen.

    python -m experiments.phase3_lle.decide [--phase test]

Writes results/<phase>/tables/primary.md.
"""
from __future__ import annotations

import argparse
import math

from experiments.phase3_lle import analysis as A
from experiments.phase3_lle import config as C

CANDIDATES = ["c1_rosenstein_m2_iaaft", "c2_kantz_m3_sat_iaaft", "c3_eps_m2_iaaft", "c4_zero_one_iaaft"]
ALPHA = 0.05
N = C.PRIMARY_WINDOW


def cell(rows, cond):
    r = [x for x in rows if x["window"] == N and x["condition"] == cond]
    return r[0] if r else None


def evaluate(phase):
    out = []
    for i, m in enumerate(["baseline"] + CANDIDATES):
        rows = A.summarize(A.load(phase, m))
        if not rows:
            out.append({"method": m, "missing": True})
            continue
        wn, ar = cell(rows, "white_noise"), cell(rows, "ar1")
        lo_wn = A.wilson(wn["detected_k"], wn["N"])[0]
        lo_ar = A.wilson(ar["detected_k"], ar["N"])[0]
        passes = lo_wn <= ALPHA and lo_ar <= ALPHA
        l20, h20 = cell(rows, "logistic@20dB"), cell(rows, "henon@20dB")
        pooled_k, pooled_n = l20["detected_k"] + h20["detected_k"], l20["N"] + h20["N"]
        lg, hn = cell(rows, "logistic"), cell(rows, "henon")
        tb = (abs(lg["bias"]) + abs(hn["bias"])) if (lg["lle_n"] and hn["lle_n"]) else math.inf
        out.append({"method": m, "index": i, "eligible": m != "baseline",
                    "wn": f"{wn['detected_k']}/{wn['N']}", "wn_lo": lo_wn,
                    "ar": f"{ar['detected_k']}/{ar['N']}", "ar_lo": lo_ar, "passes": passes,
                    "pooled20": pooled_k, "pooled20_n": pooled_n, "tiebreak_abs_bias": tb,
                    "null_fp": wn["detected_k"] + ar["detected_k"]})
    elig = [r for r in out if r.get("eligible") and r.get("passes")]
    elig.sort(key=lambda r: (-r["pooled20"], r["tiebreak_abs_bias"], r["null_fp"], r["index"]))
    winner = elig[0]["method"] if elig else None
    return out, elig, winner


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", default="test", choices=("dev", "test"))
    a = ap.parse_args(argv)
    out, ranked, winner = evaluate(a.phase)
    md = [f"# Primary decision rule ({a.phase} seeds, {N} samples)", "",
          "PASS: lower 95 % Wilson bound of the false-positive rate <= 0.05 on white_noise AND ar1. "
          "WINNER: highest pooled detection on logistic@20dB + henon@20dB; tie-break: lowest "
          "|bias logistic| + |bias henon| (clean), then fewer null false positives, then method index. "
          "The baseline is reported but not eligible.", "",
          "| method | WN FP (Wilson lower) | AR1 FP (Wilson lower) | PASS | 20 dB pooled | clean |bias| sum |",
          "|---|---|---|---|---|---|"]
    for r in out:
        if r.get("missing"):
            md.append(f"| {r['method']} | missing | | | | |")
            continue
        tb = "NA (no LLE)" if not math.isfinite(r["tiebreak_abs_bias"]) else f"{r['tiebreak_abs_bias']:.4f}"
        md.append(f"| {r['method']} | {r['wn']} ({r['wn_lo']:.3f}) | {r['ar']} ({r['ar_lo']:.3f}) | "
                  f"{'yes' if r['passes'] else 'no'}{'' if r['eligible'] else ' (not eligible)'} | "
                  f"{r['pooled20']}/{r['pooled20_n']} | {tb} |")
    md += ["", "Ranking of passing candidates: " + (", ".join(r["method"] for r in ranked) or "none"), "",
           f"**WINNER: {winner}**" if winner else "**No candidate passes; there is no winner.**"]
    tdir = A.RESULTS / a.phase / "tables"
    tdir.mkdir(parents=True, exist_ok=True)
    (tdir / "primary.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
