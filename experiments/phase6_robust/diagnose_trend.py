"""
Phase 6 Part C1: mechanism of the trend power loss (DEVELOPMENT seeds 0-29 only).

For P1 / P3 (Henon, logistic; the same seed = the same map realization and RR
scale, P3 adds the N3-type trend) the C2 UPO analysis (combined_chaos_config)
is run on
  clean   P1 window,
  trend   P3 window,
  trend + linear detrend (diagnostic only: shows that removing the trend
          restores the peak).
Recorded per window: Level-B peak count, gate decision, So-transform
histogram peak height (max bin / points), and the location of the
So-transform mode computed separately from the transformed points whose
source beat lies in the first, middle and last third of the window.  The
mode's drift across thirds is compared with the trend's change over the same
span.  Writes results/dev/diagnose_trend.{json,md} and plots/C1_*.png.

    python -m experiments.phase6_robust.diagnose_trend
"""
from __future__ import annotations

import json
import multiprocessing as mp
import pathlib
import sys
import warnings
from dataclasses import replace

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import final_pipeline as fp  # noqa: E402
from experiments.phase5_rr import systems as S5  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
CFG = fp.combined_chaos_config()
LIN = replace(CFG, rr_detrend="linear")
SEEDS = range(30)
PAIRS = (("henon", "P1_henon_rr", "P3_henon_rr_trend"), ("logistic", "P2_logistic_rr", "P3_logistic_rr_trend"))


def thirds_modes(p1, n_beats, bins=60):
    sc, src = np.asarray(p1["scalar"]), np.asarray(p1["source_indices"])
    lo, hi = np.percentile(sc, [1, 99])
    edges = np.linspace(lo, hi, bins + 1)
    modes = []
    for a, b in ((0, n_beats / 3), (n_beats / 3, 2 * n_beats / 3), (2 * n_beats / 3, n_beats)):
        sel = (src >= a) & (src < b)
        h, _ = np.histogram(sc[sel], bins=edges)
        modes.append(float(0.5 * (edges[np.argmax(h)] + edges[np.argmax(h) + 1])) if sel.any() else None)
    return modes


def analyse(x):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out = fp.run_upo_analysis(x, config=CFG, rng=np.random.default_rng(CFG.random_seed))
    p1 = out["periods"][1]
    h = np.asarray(p1["histogram"], dtype=float)
    return {"level_b": len(out.get("significant_uop_candidates") or []),
            "gate": bool(out.get("instability_gate_detected")),
            "peak_height": float(h.max() / max(h.sum(), 1)),
            "modes_thirds": thirds_modes(p1, len(x)),
            "status": out["status"]}


def one(task):
    name, clean_c, trend_c, seed = task
    xc = S5.generate(clean_c, seed)
    pt = S5.parts(trend_c, seed)
    xt = pt["rr"]
    xd = fp.detrend_rr(xt, LIN)
    frac = pt["params"]["trend_fraction"]
    # trend change between the centres of the first and last thirds (2/3 of the window)
    trend_span = frac * pt["mean_rr"] * (2.0 / 3.0)
    return {"system": name, "seed": seed, "trend_fraction": frac, "trend_change_thirds_s": trend_span,
            "clean": analyse(xc), "trend": analyse(xt), "trend_linear_detrended": analyse(xd)}


def plots(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    pdir = HERE / "plots"
    pdir.mkdir(exist_ok=True)
    seed = 0
    fig, axes = plt.subplots(2, 4, figsize=(18, 8))
    for r, (name, clean_c, trend_c) in enumerate(PAIRS):
        xc, pt = S5.generate(clean_c, seed), S5.parts(trend_c, seed)
        xt = pt["rr"]
        ax = axes[r, 0]
        ax.plot(xc, lw=0.7, label=clean_c)
        ax.plot(xt, lw=0.7, label=trend_c, alpha=0.8)
        ax.set_title(f"{name}: RR series (dev seed {seed})", fontsize=9)
        ax.legend(fontsize=7)
        ax = axes[r, 1]
        k = np.arange(len(xt) - 1)
        sc = ax.scatter(xt[:-1], xt[1:], c=k, s=4, cmap="viridis")
        lim = [min(xt.min(), xc.min()), max(xt.max(), xc.max())]
        ax.plot(lim, lim, "k--", lw=0.6)
        ax.set_title("trend: return map coloured by beat (fixed point = diagonal crossing)", fontsize=8)
        fig.colorbar(sc, ax=ax, fraction=0.04)
        for c, (x, lab) in enumerate(((xc, "clean"), (xt, "trend"))):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                p1 = fp.run_upo_analysis(x, config=CFG, rng=np.random.default_rng(CFG.random_seed))["periods"][1]
            ax = axes[r, 2 + c]
            s, src = np.asarray(p1["scalar"]), np.asarray(p1["source_indices"])
            lo, hi = np.percentile(s, [1, 99])
            for a, b, col in ((0, 85, "#2b6cb0"), (85, 171, "#38a169"), (171, 256, "#c05621")):
                sel = (src >= a) & (src < b)
                ax.hist(s[sel], bins=np.linspace(lo, hi, 61), histtype="step", color=col, label=f"beats {a}-{b}")
            ax.set_title(f"{name} {lab}: So-transform scalar by window third", fontsize=8)
            ax.legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(pdir / "C1_trend_mechanism.png", dpi=85)
    plt.close(fig)


def main():
    tasks = [(n, c, t, s) for n, c, t in PAIRS for s in SEEDS]
    with mp.get_context("fork").Pool(4) as pool:
        rows = pool.map(one, tasks)
    out = HERE / "results" / "dev"
    out.mkdir(parents=True, exist_ok=True)
    (out / "diagnose_trend.json").write_text(json.dumps(rows, indent=1))
    md = ["# Part C1: trend mechanism (development seeds 0-29, C2 UPO analysis)", "",
          "Mode drift = |So-transform mode in the last third - mode in the first third| (s); trend change = "
          "|trend difference between the centres of the first and last thirds| (s).", "",
          "| system | variant | Level-B windows | gate (UPO) windows | median peak height | median mode drift (s) | "
          "median trend change (s) | median drift / trend change |", "|---|---|---|---|---|---|---|---|"]
    for name, _, _ in PAIRS:
        rs = [r for r in rows if r["system"] == name]
        tc = np.array([abs(r["trend_change_thirds_s"]) for r in rs])
        for v in ("clean", "trend", "trend_linear_detrended"):
            drift = np.array([abs(r[v]["modes_thirds"][2] - r[v]["modes_thirds"][0]) for r in rs])
            md.append(f"| {name} | {v} | {sum(r[v]['level_b'] > 0 for r in rs)}/30 | {sum(r[v]['gate'] for r in rs)}/30 | "
                      f"{np.median([r[v]['peak_height'] for r in rs]):.3f} | {np.median(drift):.4f} | "
                      f"{np.median(tc):.4f} | {np.median(drift / tc):.2f} |")
    (out / "diagnose_trend.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))
    plots(rows)


if __name__ == "__main__":
    main()
