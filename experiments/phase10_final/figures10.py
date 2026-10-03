"""
Phase 10 figures (static PNGs for docs/PHASE10_FINAL_RESULTS.md) from results/<phase>/analysis/analysis.json
and results/verification/titration_verification.json.

    python -m experiments.phase10_final.figures10 --phase conf

Style: reference categorical palette (fixed slot order), 2 px lines, >= 8 px markers, recessive grid,
one y-axis per panel, legends for >= 2 series, text in ink colours.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
FAMILY_TITLE = {"henon": "Hénon map", "logistic": "Logistic map", "coupled_vdp": "Coupled van der Pol (Phase 8)",
                "phase_reset": "Phase-resetting map (Phase 8)"}


def _style(ax):
    ax.set_facecolor(SURF)
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)
    ax.title.set_color(INK)


def fig_verification(out):
    p = HERE / "results" / "verification" / "titration_verification.json"
    if not p.exists():
        return
    v = json.load(open(p))["verification"]["V4"]["scan"]
    r = np.array([a for a, _, _ in v])
    le = np.array([b for _, b, _ in v])
    nl = np.array([c for _, _, c in v])
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.5, 5.2), sharex=True, facecolor=SURF)
    a1.plot(r, 100 * nl, "-o", color=SERIES[0], lw=2, ms=4)
    a1.set_ylabel("noise limit NL (% of SD)")
    a1.set_title("Titration of the logistic map (1,000 points; PNAS 2001 Fig. 1 design)", fontsize=11, loc="left")
    a2.plot(r, le, "-", color=SERIES[1], lw=2)
    a2.axhline(0, color=INK2, lw=1)
    a2.set_ylabel("Lyapunov exponent")
    a2.set_xlabel("r")
    for a in (a1, a2):
        _style(a)
    fig.tight_layout()
    fig.savefig(out / "fig_titration_verification.png", dpi=150)
    plt.close(fig)


def fig_curves(A, out, det):
    D = A.get("Q1", {}).get("detectors", {}).get(det)
    if not D:
        return
    f = np.round(np.arange(0.05, 0.9001, 0.01), 2)
    fams = ["henon", "logistic", "coupled_vdp", "phase_reset"]
    fig, axes = plt.subplots(2, 4, figsize=(15, 7), facecolor=SURF, sharex=True)
    for col, fam in enumerate(fams):
        keys = [k for k in D["curves"] if k.startswith(fam + "|")]
        keys.sort(key=lambda k: D["curves"][k]["lam"])
        ax = axes[0, col]
        ax2 = axes[1, col]
        for i, k in enumerate(keys):
            c = D["curves"][k]
            lab = f"{c['param']} (λ {c['lam']:.3g})"
            ax.fill_between(f, c["p_ci95"][0], c["p_ci95"][1], color=SERIES[i], alpha=0.12, lw=0)
            ax.plot(f, c["p_hat"], color=SERIES[i], lw=2, label=lab)
            ef = [float(x) for x in c["empirical"]]
            ep = [c["empirical"][x]["p"] for x in c["empirical"]]
            ax.plot(ef, ep, "o", color=SERIES[i], ms=5, mec=SURF, mew=1)
            pi = D["exclusion_L0"][k]["p_hat"]["pi_upper"]
            ax2.plot(f, pi, color=SERIES[i], lw=2, label=lab)
        ax.set_ylim(-0.02, 1.02)
        ax.set_title(FAMILY_TITLE[fam], fontsize=10, loc="left")
        ax2.set_ylim(0, 1.02)
        for th in (0.05, 0.20):
            ax2.axhline(th, color=INK2, lw=1, ls="--")
        ax2.set_xlabel("chaos fraction f")
        if col == 0:
            ax.set_ylabel(f"{det} detection probability")
            ax2.set_ylabel("π upper (95 %)")
        ax.legend(fontsize=7.5, frameon=False, loc="upper left")
        _style(ax)
        _style(ax2)
    U = D["real"]["U"]
    fig.suptitle(f"{det}: detection probability of additive chaos in real windows (top) and exclusion bound "
                 f"(bottom; real-rate bound U = {U:.4f})", fontsize=11, color=INK, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(out / f"fig_q1_{det}.png", dpi=150)
    plt.close(fig)


def fig_p2(A, out):
    Q = A.get("Q2", {})
    if not Q:
        return
    arms = [("P2_raw_by_group", "raw"), ("masked_by_group", "masked"), ("edited_by_group", "edited"),
            ("wu_by_group", "Wu NN")]
    fig, ax = plt.subplots(figsize=(7.5, 4), facecolor=SURF)
    for gi, g in enumerate(("NSR", "CHF")):
        xs, ys, lo, hi = [], [], [], []
        for ai, (key, name) in enumerate(arms):
            r = Q.get(key, {}).get(g)
            if not r:
                continue
            xs.append(ai + (gi - 0.5) * 0.25)
            ys.append(100 * r["mean_subject_rate"])
            lo.append(100 * (r["mean_subject_rate"] - r["ci95"][0]))
            hi.append(100 * (r["ci95"][1] - r["mean_subject_rate"]))
        ax.errorbar(xs, ys, yerr=[lo, hi], fmt="o", ms=8, color=SERIES[gi], lw=2, capsize=0, label=g)
    ax.set_xticks(range(len(arms)), [n for _, n in arms])
    ax.set_ylabel("titration-positive 12-min segments (%)\nmean of subject rates, 95 % CI")
    ax.set_ylim(0, 100)
    ax.legend(frameon=False)
    ax.set_title("Noise titration detection rate by group and preprocessing", fontsize=11, loc="left")
    _style(ax)
    fig.tight_layout()
    fig.savefig(out / "fig_p2_titration.png", dpi=150)
    plt.close(fig)


def fig_p3(A, out):
    Q = A.get("Q3", {})
    if not Q:
        return
    bt = Q["burden_table"]
    labels = [b["burden"] for b in bt]
    fig, ax = plt.subplots(figsize=(7.5, 4.2), facecolor=SURF)
    x = np.arange(len(labels))
    for i, (c, name) in enumerate((("TIT", "titration"), ("LLE", "LLE alone"), ("UPO", "UPO alone"),
                                   ("K1", "K1"), ("K3", "K3"))):
        y = [100 * (b[c]["rate"] or 0) for b in bt]
        ax.plot(x, y, "-o", color=SERIES[i], lw=2, ms=7, label=name)
    ax.set_xticks(x, [f"{l}\n(n={b['n']})" for l, b in zip(labels, bt)])
    ax.set_xlabel("ectopy burden (non-normal beats in the 513-beat window)")
    ax.set_ylabel("positive windows (%)")
    ax.legend(frameon=False, fontsize=8.5, ncol=2)
    ax.set_title("Positive rate by ectopy burden (512-interval windows)", fontsize=11, loc="left")
    _style(ax)
    fig.tight_layout()
    fig.savefig(out / "fig_p3_burden.png", dpi=150)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", choices=["dev", "conf"], required=True)
    a = ap.parse_args(argv)
    out = HERE / "results" / a.phase / "analysis" / "plots"
    out.mkdir(parents=True, exist_ok=True)
    A = json.load(open(HERE / "results" / a.phase / "analysis" / "analysis.json"))
    fig_verification(out)
    for det in ("K3", "K1", "LLE", "UPO"):
        fig_curves(A, out, det)
    fig_p2(A, out)
    fig_p3(A, out)
    print("figures in", out)


if __name__ == "__main__":
    main()
