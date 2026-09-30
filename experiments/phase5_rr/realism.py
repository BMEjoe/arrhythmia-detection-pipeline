"""
Phase 5 Step 1: realism check of the RR generators (DEVELOPMENT seeds only; no
detector is run here).

    python -m experiments.phase5_rr.realism

Per window: mean RR, SDNN, RMSSD (ms); LF (0.04-0.15 Hz) and HF (0.15-0.40 Hz)
power (ms^2), LF/HF and LF n.u. = 100 LF / (LF + HF) (Task Force 1996 band
definitions), from the RR series interpolated at 4 Hz (cubic spline over the
beat times), linearly detrended, Welch periodogram (64 s Hann segments, 50 %
overlap); sample skewness; number of distinct values.  Targets and sources:
README.md.  Writes results/realism/realism_dev.{csv,md} and example plots in
plots/.
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
from scipy import interpolate, signal, stats

ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.phase5_rr import config as C  # noqa: E402
from experiments.phase5_rr import systems as S  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
FS = 4.0
LF_BAND, HF_BAND = (0.04, 0.15), (0.15, 0.40)

# Targets (README.md): mean RR (user), SDNN / RMSSD short-term (Nunan et al. 2010,
# mean +- 1 SD), LF/HF (Task Force 1996, 5-min supine 1.5-2.0; Nunan 2.8 +- 2.6).
TARGETS = {"mean_rr_s": (0.6, 1.0), "sdnn_ms": (34.0, 66.0), "rmssd_ms": (27.0, 57.0),
           "lf_hf_taskforce": (1.5, 2.0), "lf_hf_nunan": (0.2, 5.4)}


def psd(rr):
    t = np.cumsum(rr)
    t = t - t[0]
    grid = np.arange(0.0, t[-1], 1.0 / FS)
    y = interpolate.CubicSpline(t, rr)(grid) * 1000.0
    f, p = signal.welch(y, fs=FS, window="hann", nperseg=int(64 * FS), noverlap=int(32 * FS),
                        detrend="linear")
    return f, p


def band(f, p, lo, hi):
    sel = (f >= lo) & (f < hi)
    return float(np.trapezoid(p[sel], f[sel]))


def measures(rr):
    rr = np.asarray(rr, dtype=float)
    f, p = psd(rr)
    lf, hf = band(f, p, *LF_BAND), band(f, p, *HF_BAND)
    return {"mean_rr_s": float(rr.mean()), "min_rr_s": float(rr.min()), "max_rr_s": float(rr.max()),
            "sdnn_ms": float(rr.std(ddof=1) * 1000), "rmssd_ms": float(np.sqrt(np.mean(np.diff(rr) ** 2)) * 1000),
            "lf_ms2": lf, "hf_ms2": hf, "lf_hf": lf / hf if hf > 0 else np.nan,
            "lf_nu": 100 * lf / (lf + hf) if lf + hf > 0 else np.nan,
            "skewness": float(stats.skew(rr)), "n_distinct": int(len(np.unique(rr))),
            "duration_s": float(rr.sum())}


def main():
    rows = []
    for cond in C.CONDITIONS:
        for seed in C.seeds("dev", cond):
            m = measures(S.generate(cond, seed))
            m.update(condition=cond, seed=seed)
            rows.append(m)
    outdir = HERE / "results" / "realism"
    outdir.mkdir(parents=True, exist_ok=True)
    keys = ["condition", "seed", "mean_rr_s", "min_rr_s", "max_rr_s", "sdnn_ms", "rmssd_ms", "lf_ms2", "hf_ms2",
            "lf_hf", "lf_nu", "skewness", "n_distinct", "duration_s"]
    with open(outdir / "realism_dev.csv", "w") as fh:
        fh.write(",".join(keys) + "\n")
        for r in rows:
            fh.write(",".join(str(r[k]) for k in keys) + "\n")

    def q(v):
        return f"{np.median(v):.3g} [{np.percentile(v, 10):.3g}, {np.percentile(v, 90):.3g}]"

    def ok(v, key):
        lo, hi = TARGETS[key]
        return "yes" if lo <= np.median(v) <= hi else "no"

    md = ["# Phase 5 realism check (development seeds; no detector run)", "",
          "Median [10th, 90th percentile] over development windows (nulls and secondary 100, others 30). "
          "Targets (README.md): mean RR 0.6-1.0 s for every window; SDNN 34-66 ms and RMSSD 27-57 ms "
          "(Nunan et al. 2010, short-term mean +- 1 SD); LF/HF 1.5-2.0 (Task Force 1996, 5-min supine) or "
          "0.2-5.4 (Nunan et al. 2010, mean +- 1 SD).  'in' columns test the median.", "",
          "| condition | N | mean RR (s) | all windows 0.6-1.0 s | SDNN (ms) | in | RMSSD (ms) | in | LF/HF | "
          "in TF / Nunan | LF n.u. | skewness | min RR (s) | distinct values |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for cond in C.CONDITIONS:
        rs = [r for r in rows if r["condition"] == cond]
        g = {k: np.array([r[k] for r in rs]) for k in keys[2:]}
        allin = all(0.6 <= r["mean_rr_s"] <= 1.0 for r in rs)
        md.append(f"| {cond} | {len(rs)} | {q(g['mean_rr_s'])} | {'yes' if allin else 'no'} | {q(g['sdnn_ms'])} | "
                  f"{ok(g['sdnn_ms'], 'sdnn_ms')} | {q(g['rmssd_ms'])} | {ok(g['rmssd_ms'], 'rmssd_ms')} | "
                  f"{q(g['lf_hf'])} | {ok(g['lf_hf'], 'lf_hf_taskforce')} / {ok(g['lf_hf'], 'lf_hf_nunan')} | "
                  f"{q(g['lf_nu'])} | {q(g['skewness'])} | {q(g['min_rr_s'])} | {q(g['n_distinct'])} |")
    (outdir / "realism_dev.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))
    plots()


def plots():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    pdir = HERE / "plots"
    pdir.mkdir(exist_ok=True)
    conds = list(C.CONDITIONS)
    fig, axes = plt.subplots(8, 3, figsize=(15, 20), sharex=True)
    for ax, cond in zip(axes.T.ravel(), conds):
        rr = S.generate(cond, 0)
        ax.plot(rr, lw=0.8, color="#2b6cb0")
        ax.set_title(cond, fontsize=9)
        ax.set_ylabel("RR (s)", fontsize=8)
        ax.tick_params(labelsize=7)
    for ax in axes[-1]:
        ax.set_xlabel("beat")
    fig.suptitle("Phase 5 generators: development seed 0 of each condition (256 RR intervals)")
    fig.tight_layout()
    fig.savefig(pdir / "examples_timeseries.png", dpi=90)
    plt.close(fig)
    fig, axes = plt.subplots(8, 3, figsize=(15, 20), sharex=True)
    for ax, cond in zip(axes.T.ravel(), conds):
        ps = []
        for seed in range(10):
            f, p = psd(S.generate(cond, seed))
            ps.append(p)
        ax.semilogy(f, np.median(ps, axis=0), color="#c05621", lw=1)
        for lo, hi, col in ((0.04, 0.15, "#bee3f8"), (0.15, 0.40, "#c6f6d5")):
            ax.axvspan(lo, hi, color=col, alpha=0.5, lw=0)
        ax.set_xlim(0, 0.5)
        ax.set_title(cond, fontsize=9)
        ax.tick_params(labelsize=7)
    for ax in axes[-1]:
        ax.set_xlabel("Hz")
    fig.suptitle("Median Welch PSD (ms$^2$/Hz) over development seeds 0-9; LF (blue) and HF (green) bands")
    fig.tight_layout()
    fig.savefig(pdir / "examples_psd.png", dpi=90)
    plt.close(fig)
    # return maps (x_k, x_k+1) show the deterministic structure (or its absence)
    fig, axes = plt.subplots(4, 6, figsize=(18, 12))
    for ax, cond in zip(axes.ravel(), conds):
        rr = S.generate(cond, 0)
        ax.plot(rr[:-1], rr[1:], ".", ms=2, color="#2d3748")
        ax.set_title(cond, fontsize=8)
        ax.tick_params(labelsize=6)
    fig.suptitle("Return maps RR_k+1 vs RR_k (development seed 0)")
    fig.tight_layout()
    fig.savefig(pdir / "examples_return_maps.png", dpi=80)
    plt.close(fig)


if __name__ == "__main__":
    main()
