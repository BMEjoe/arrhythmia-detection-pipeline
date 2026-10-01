"""
Phase 8 Part A realism table and example plots for every selected regime (no detector).
Realism seeds 900000-900009 (never used elsewhere).

    python -m experiments.phase8_cardiac.realism
"""
from __future__ import annotations

import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from experiments.phase5_rr.realism import measures  # noqa: E402
from experiments.phase8_cardiac import series as SR  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "results" / "realism"
PLOTS = HERE / "plots"
SEEDS = range(900000, 900010)
INK, MUTED, SERIES = "#1f1f1e", "#6b6a63", "#2a78d6"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    PLOTS.mkdir(exist_ok=True)
    regimes = SR.load_regimes()
    rows = []
    for r in regimes:
        for variant in ("i_clean", "iv_q", "vi_E3_10"):
            ms = [measures(SR.window(r["family"], r["regime_idx"], r["params"], 1024, variant, s)) for s in SEEDS]
            m = pd.DataFrame(ms).median()
            rows.append({"family": r["family"], "split": r["split"], "regime": r["name"], "label": r["label"],
                         "lam": r["lam"], "variant": variant, **m.to_dict()})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "realism.csv", index=False)
    cols = ["family", "split", "regime", "label", "variant", "mean_rr_s", "min_rr_s", "max_rr_s", "sdnn_ms",
            "rmssd_ms", "lf_hf", "n_distinct"]
    with open(OUT / "realism.md", "w") as fh:
        fh.write("# Phase 8 realism (median over 10 realism seeds, 1024 beats; no detector)\n\n")
        fh.write("Targets (Phase 5): mean RR 0.6-1.0 s; short-term SDNN 50 +/- 16 ms, RMSSD 42 +/- 15 ms "
                 "(Nunan et al. 2010); LF/HF 1.5-2.0 (Task Force 1996).\n\n")
        fh.write(df[cols].round(3).to_markdown(index=False))
        fh.write("\n")
    # example plots: one figure per family, one panel per regime (clean, first 256 beats)
    for fam, g in pd.DataFrame(regimes).groupby("family"):
        n = len(g)
        fig, axes = plt.subplots(n, 2, figsize=(11, 1.8 * n), squeeze=False, sharex="col")
        for row, (_, r) in enumerate(g.iterrows()):
            for col, variant in enumerate(("i_clean", "vi_E3_10")):
                x = SR.window(fam, r["regime_idx"], r["params"], 256, variant, 900000)
                ax = axes[row, col]
                ax.plot(np.arange(len(x)), x, color=SERIES, lw=1.0)
                ax.set_title(f"{r['name']}  [{r['label']}, LE {r['lam']:.3g}]  {variant}", fontsize=8, color=INK, loc="left")
                ax.grid(True, color="#e6e5df", lw=0.6)
                ax.tick_params(labelsize=7, colors=MUTED)
                for s in ("top", "right"):
                    ax.spines[s].set_visible(False)
                ax.set_ylabel("RR (s)", fontsize=7, color=MUTED)
        for ax in axes[-1]:
            ax.set_xlabel("beat", fontsize=7, color=MUTED)
        fig.tight_layout()
        fig.savefig(PLOTS / f"A_examples_{fam}.png", dpi=110)
        plt.close(fig)
    print(df[cols].round(3).to_string())


if __name__ == "__main__":
    main()
