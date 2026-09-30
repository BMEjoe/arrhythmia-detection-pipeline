"""
Descriptive breakdown of the stored test records for the report (no new
detector runs, nothing changed): per condition at m = 2, the fraction of
windows with >= 1 Level-B UPO peak, the fraction passing the instability gate
(= UPO detection), LLE-test p-value quartiles and the fraction with p <= 0.05,
and the mean production-LLE (analyze_segment) TDMI tau / Cao m.

    python -m experiments.phase5_rr.report_extras
"""
from __future__ import annotations

import pathlib

import numpy as np

from experiments.phase5_rr import analysis as A
from experiments.phase5_rr import config as C

HERE = pathlib.Path(__file__).resolve().parent


def main():
    recs = A.load("test", "test")
    md = ["# Phase 5 test records: component breakdown at m = 2 (descriptive)", "",
          "| condition | N | >= 1 Level-B peak | UPO gate passed | LLE p quartiles | LLE p <= 0.05 | "
          "AND | production TDMI tau / Cao m (median) |", "|---|---|---|---|---|---|---|---|"]
    for cond in C.CONDITIONS:
        ds = [r["m2"] for r in recs if r["task"]["condition"] == cond]
        ok = [d for d in ds if d.get("error") is None]
        lb = sum((d.get("upo_n_level_b") or 0) > 0 for d in ds)
        gate = sum(d["upo_n_gated"] > 0 for d in ds)
        ps = np.array([d["lle_p"] for d in ds if d["lle_p"] is not None])
        q = np.percentile(ps, [25, 50, 75]) if len(ps) else [np.nan] * 3
        lle = sum(p <= 0.05 for p in ps)
        both = sum(d["upo_n_gated"] > 0 and d["lle_p"] is not None and d["lle_p"] <= 0.05 for d in ds)
        taus = [d["tdmi_tau"] for d in ok if "tdmi_tau" in d]
        ms = [d["cao_m"] for d in ok if "cao_m" in d]
        md.append(f"| {cond} | {len(ds)} | {lb} | {gate} | {q[0]:.2f} / {q[1]:.2f} / {q[2]:.2f} | {lle} | {both} | "
                  f"{np.median(taus):.0f} / {np.median(ms):.0f} |")
    md += ["", "Counts include the 3 windows where analyze_segment raised (component decisions from the "
           "direct calls); the preregistered AND counts those windows as not detected."]
    (HERE / "results" / "test" / "tables" / "component_breakdown.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
