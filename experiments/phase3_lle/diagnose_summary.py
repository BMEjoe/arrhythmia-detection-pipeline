"""
Phase 3, B2 -- summary tables and plots of results/dev/diagnosis.jsonl.

Writes results/dev/tables/B2_*.md and plots/B2_*.png.  Development seeds only.
"""
from __future__ import annotations

import json
import pathlib

import numpy as np

from experiments.phase3_lle import config as C

HERE = pathlib.Path(__file__).resolve().parent
SRC = HERE / "results" / "dev" / "diagnosis.jsonl"
TABLES = HERE / "results" / "dev" / "tables"
PLOTS = HERE / "plots"


def load():
    return [json.loads(line) for line in open(SRC)]


def med(v):
    v = [x for x in v if x is not None and np.isfinite(x)]
    return float(np.median(v)) if v else float("nan")


def iqr(v):
    v = [x for x in v if x is not None and np.isfinite(x)]
    if not v:
        return "NA"
    q1, q2, q3 = np.percentile(v, [25, 50, 75])
    return f"{q2:.3f} [{q1:.3f}, {q3:.3f}]"


def groups(R):
    out = {}
    for r in R:
        out.setdefault((C.condition_label(r["system"], r["snr_db"]), r["window_length"]), []).append(r)
    order = [C.condition_label(s, snr) for s, snr in C.conditions()]
    return sorted(out.items(), key=lambda kv: (kv[0][1], order.index(kv[0][0])))


def table(rows, header):
    lines = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    lines += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return "\n".join(lines)


def ref(label):
    return C.LLE_REFERENCE.get(label.split("@")[0])


def main():
    R = load()
    TABLES.mkdir(parents=True, exist_ok=True)
    PLOTS.mkdir(parents=True, exist_ok=True)
    G = groups(R)
    md = ["# B2 diagnosis (development seeds; first 30 seeds per condition)", "",
          "LLE = slope of the frozen fp.rosenstein_lle mean log-divergence curve, per sample. "
          "Cells: median [Q1, Q3] over windows. Reference: logistic ln 2 = 0.693, Henon 0.419.", ""]

    # D1 delay and D2 dimension: production fit on each embedding
    rows = []
    for (lab, n), g in G:
        r0 = ref(lab)
        cells = [lab, n, len(g)]
        for v in ("V0", "V1", "V2", "V3"):
            vals = [x["variants"][v]["prod"] for x in g]
            cells.append(iqr(vals))
        if r0:
            cells.append(" / ".join(f"{med([x['variants'][v]['prod'] for x in g]) - r0:+.3f}"
                                    for v in ("V0", "V1", "V2", "V3")))
        else:
            cells.append("")
        cells.append(f"{med([x['tau_prod'] for x in g]):.0f}")
        cells.append(f"{med([x['m_prod'] for x in g]):.0f}")
        cells.append(f"{med([x['m_cao_lag1'] for x in g]):.0f}")
        rows.append(cells)
    md += ["## D1/D2 -- delay and dimension (production early fit on each embedding)", "",
           "V0 = production (TDMI tau*, Cao m* at tau*); V1 = tau 1, m*; V2 = tau 1, Cao m at lag 1; "
           "V3 = tau 1, true m (2 for non-map systems).", "",
           table(rows, ["condition", "N", "windows", "V0", "V1", "V2", "V3", "median bias V0/V1/V2/V3",
                        "tau*", "m*", "m Cao lag 1"]), ""]

    # D3 fit region: saturation vs the fixed 6-point early fit (k = 0..5)
    rows = []
    for (lab, n), g in G:
        for v in ("V0", "V3"):
            k95 = [x["variants"][v].get("k95") for x in g]
            kh = [x["variants"][v].get("k_half") for x in g]
            within = sum(1 for k in k95 if k is not None and k <= 5)
            rows.append([lab, n, v, f"{med(kh):.0f}", f"{med(k95):.0f}", f"{within}/{len(g)}",
                         iqr([x["variants"][v]["prod"] for x in g]),
                         iqr([x["variants"][v].get("k1_6") for x in g]),
                         iqr([x["variants"][v].get("sat") for x in g])])
    md += ["## D3 -- fit region", "",
           "k_half = first step at which the curve has risen half-way from y[1] to its plateau "
           "(mean of y[25:40]); k95 = first step within 95 % of the rise from y[0]. The production "
           "fit uses steps 0..5. 'k95 <= 5' = windows whose curve saturates inside the production fit.",
           "", table(rows, ["condition", "N", "emb", "k_half median", "k95 median", "k95 <= 5",
                            "slope prod (0..5)", "slope 1..6", "slope 1..k_half"]), ""]

    # D4 noise floor: k = 0 -> 1 jump and early slope vs SNR (V3)
    rows = []
    for (lab, n), g in G:
        if lab.split("@")[0] not in C.CHAOTIC_SYSTEMS and lab not in C.NULL_SYSTEMS:
            continue
        r0 = ref(lab)
        prod = [x["variants"]["V3"]["prod"] for x in g]
        s16 = [x["variants"]["V3"].get("k1_6") for x in g]
        rows.append([lab, n, iqr([x["variants"]["V3"].get("jump01") for x in g]),
                     iqr([x["variants"]["V3"].get("jump12") for x in g]),
                     iqr([x["variants"]["V3"].get("y0") for x in g]),
                     iqr(prod), f"{med(prod) - r0:+.3f}" if r0 else "",
                     iqr(s16), f"{med(s16) - r0:+.3f}" if r0 else ""])
    md += ["## D4 -- noise floor (tau 1, true m)", "",
           "jump01 = y[1] - y[0], jump12 = y[2] - y[1] (log-distance increments over the first two "
           "steps); y0 = mean log initial nearest-neighbour distance.", "",
           table(rows, ["condition", "N", "jump01", "jump12", "y0", "slope prod (0..5)", "bias",
                        "slope 1..6", "bias"]), ""]
    (TABLES / "B2_diagnosis.md").write_text("\n".join(md))
    plots(R)
    print("\n".join(md))


def plots(R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    show = [("logistic", None), ("logistic", 30.0), ("logistic", 20.0), ("logistic", 10.0),
            ("henon", None), ("henon", 30.0), ("henon", 20.0), ("henon", 10.0),
            ("white_noise", None), ("ar1", None), ("sinusoid", None), ("logistic_p4", None)]
    for n in C.WINDOW_LENGTHS:
        fig, axes = plt.subplots(3, 4, figsize=(13, 8.5), sharex=True)
        for ax, (sysn, snr) in zip(axes.ravel(), show):
            recs = [r for r in R if r["system"] == sysn and r["snr_db"] == snr and r["window_length"] == n]
            for v, col in (("V0", "#b2182b"), ("V1", "#ef8a62"), ("V2", "#67a9cf"), ("V3", "#2166ac")):
                for i, r in enumerate(recs[:5]):
                    y = r["variants"][v]["curve"]
                    if y is None:
                        continue
                    y = np.array([np.nan if t is None else t for t in y])
                    ax.plot(np.arange(len(y)), y - y[0], color=col, lw=0.8, alpha=0.7,
                            label=f"{v} tau={r['variants'][v]['tau']} m={r['variants'][v]['m']}" if i == 0 else None)
            ax.axvspan(0, 5, color="0.85", zorder=0, label="production fit k=0..5")
            r0 = C.LLE_REFERENCE.get(sysn)
            if r0:
                k = np.arange(0, 12)
                ax.plot(k, r0 * k, "k--", lw=1, label=f"reference slope {r0:.3f}")
            ax.set_title(C.condition_label(sysn, snr), fontsize=9)
            ax.set_xlim(0, 40)
            ax.legend(fontsize=6, loc="lower right")
        for ax in axes[-1]:
            ax.set_xlabel("step k")
        for ax in axes[:, 0]:
            ax.set_ylabel("mean log divergence - y[0]")
        fig.suptitle(f"B2: Rosenstein divergence curves, N = {n}, dev seeds 0-4 (grey = production fit region)")
        fig.tight_layout()
        fig.savefig(PLOTS / f"B2_divergence_curves_N{n}.png", dpi=100)
        plt.close(fig)


if __name__ == "__main__":
    main()
