"""
Phase 4 -- the PREREGISTERED primary decision rule (PREREGISTRATION.md
Section 3), written before any test result was seen.

    python -m experiments.phase4_upo.decide --phase test --file test

At 256 samples:
  PASS    detections on EACH of white_noise, ar1, sinusoid, two_tone
          <= floor(0.07 N) windows (7/100, 3/50, ...), stated as counts
  WINNER  among passing candidates: highest pooled detection count on
          henon@30dB + henon@20dB
  tie     median localization error on clean henon (nearest detected peak to
          the analytic on-attractor fixed point; no detection -> worst), then
          lower pooled false-positive count on the four PASS systems, then
          lower candidate index.
The baseline is reported but not eligible.
"""
from __future__ import annotations

import argparse
import math

from experiments.phase4_upo import analysis as A
from experiments.phase4_upo import config as C
from experiments.phase4_upo import methods as MM

N = C.PRIMARY_WINDOW
FP_FRACTION = 0.07


def evaluate(phase, file):
    rows = A.summarize(A.load(phase, file), list(MM.METHODS))
    cell = {(r["method"], r["condition"]): r for r in rows if r["window"] == N}
    out = []
    for m, spec in MM.METHODS.items():
        fps = {s: cell[(m, s)] for s in C.SPECIFICITY_SYSTEMS}
        limits = {s: math.floor(FP_FRACTION * fps[s]["N"] + 1e-9) for s in fps}
        passes = all(fps[s]["det_k"] <= limits[s] for s in fps)
        h30, h20 = cell[(m, "henon@30dB")], cell[(m, "henon@20dB")]
        loc = cell[(m, "henon")]["loc_median"]
        out.append({"method": m, "eligible": spec["eligible"], "index": spec["index"],
                    "fp": {s: (fps[s]["det_k"], fps[s]["N"], limits[s]) for s in fps}, "passes": passes,
                    "pooled": h30["det_k"] + h20["det_k"], "pooled_n": h30["N"] + h20["N"],
                    "loc": math.inf if loc is None else loc,
                    "fp_total": sum(fps[s]["det_k"] for s in fps)})
    ranked = sorted([r for r in out if r["eligible"] and r["passes"]],
                    key=lambda r: (-r["pooled"], r["loc"], r["fp_total"], r["index"]))
    return out, ranked, (ranked[0]["method"] if ranked else None)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", default="test", choices=("dev", "test"))
    ap.add_argument("--file", default="test")
    a = ap.parse_args(argv)
    out, ranked, winner = evaluate(a.phase, a.file)
    md = [f"# Phase 4 primary decision ({a.phase} seeds, {N} samples, file {a.file})", "",
          "PASS: detections <= floor(0.07 N) on EACH of white_noise, ar1, sinusoid, two_tone. "
          "WINNER: most pooled detections on henon@30dB + henon@20dB; tie-break: median localization "
          "error on clean henon, then fewer pooled false positives, then candidate index. "
          "The baseline is not eligible.", "",
          "| method | WN | AR(1) | sinusoid | two-tone | PASS | Hénon 30+20 dB | clean Hénon loc. error |",
          "|---|---|---|---|---|---|---|---|"]
    for r in out:
        f = r["fp"]
        cells = [f"{f[s][0]}/{f[s][1]} (≤{f[s][2]})" for s in C.SPECIFICITY_SYSTEMS]
        loc = "NA" if not math.isfinite(r["loc"]) else f"{r['loc']:.4f}"
        md.append(f"| {r['method']} | " + " | ".join(cells) +
                  f" | {'yes' if r['passes'] else 'no'}{'' if r['eligible'] else ' (not eligible)'} | "
                  f"{r['pooled']}/{r['pooled_n']} | {loc} |")
    md += ["", "Ranking of passing candidates: " + (", ".join(r["method"] for r in ranked) or "none"), "",
           f"**WINNER: {winner}**" if winner else "**No candidate passes; there is no winner.**"]
    tdir = A.HERE / "results" / a.phase / "tables"
    tdir.mkdir(parents=True, exist_ok=True)
    (tdir / f"{a.file}_primary.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
