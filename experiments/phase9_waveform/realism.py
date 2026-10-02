"""
Phase 9 Part A realism table and example plots for every family and variant (no detector, no
chaos measure).  Realism seeds 980000-980004 (never used elsewhere), 512-beat windows.

    python -m experiments.phase9_waveform.realism --workers 4

Variants per source:
  clean            no added ectopy, no noise, true beat times
  ect              models: + E3_10 ectopy (KTz regimes carry their own input), clean, true beats
  real             ect + nstdb mix 12 dB + detected beats + jitter (the most realistic variant)
Noise grid (N1 null, phase_reset tau=1.2, ktz P=92 constant pacing): mix24/12/6, bw12, ma12, em12,
detected beats.
Columns: RR statistics of the window RR (Phase 5 realism.measures), beat detection (missed,
extra), ECG amplitudes (median R and T deflection of normal beats on the clean record, relative
to the PR-segment level), model QT (Q event - 2 sigma to T event + 2 sigma) and nstdb gain.
"""
from __future__ import annotations

import argparse
import os
import pathlib
from multiprocessing import Pool

os.environ.setdefault("OMP_NUM_THREADS", "1")
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from experiments.phase5_rr.realism import measures  # noqa: E402
from experiments.phase9_waveform import library as L  # noqa: E402
from experiments.phase9_waveform import regimes9 as R9  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "results" / "realism"
PLOTS = HERE / "plots"
SEEDS = range(980000, 980005)
N = 512
FS = 360
GRID_SOURCES = (("null", "N1_linear_rr", "none"), ("model", "phase_reset:tau=1.2", "none"), ("ktz", "ktz:P=92", "none"))
GRID_NOISE = ("mix24", "mix12", "mix6", "bw12", "ma12", "em12")


def variants(s):
    ect_fixed = s["ectopy_fixed"]
    if s["source"] == "null":
        return [("clean", "none", "none", "true", False), ("real", "none", "mix12", "detected", True)]
    if s["source"] == "ktz":
        return [("clean", ect_fixed, "none", "true", False), ("real", ect_fixed, "mix12", "detected", True)]
    return [("clean", "none", "none", "true", False), ("ect", "E3_10", "none", "true", False),
            ("real", "E3_10", "mix12", "detected", True)]


def amplitudes(w):
    e = w["ecg_clean"]
    idx, ty = w["true_idx"], w["types_true"]
    R, T = [], []
    for k in range(L.LEAD, min(len(idx) - 1, L.LEAD + 200)):
        if ty[k] != "N" or ty[k + 1] == "V":
            continue
        i, rr = idx[k], (idx[k + 1] - idx[k]) / FS
        base = np.median(e[max(0, i - int(0.12 * FS)):i - int(0.08 * FS)])
        R.append(e[i - 5:i + 6].max() - base)
        seg = e[i + int(0.1 * FS):i + int(0.7 * rr * FS)]
        if len(seg):
            T.append(seg.max() - base)
    return float(np.median(R)), float(np.median(T))


def one(task):
    s, vname, ect, noise, beats, jit, seed = task
    w = L.build(s["source"], s["name"], ect, N, seed, noise, "DEV", beats, jit, keep_clean=True)
    row = {"source": s["source"], "family": s["family"], "regime": R9.key(s), "label": s["label"], "split": s["split"],
           "variant": vname, "noise": noise, "seed": seed, "error": w.get("error")}
    if w.get("error"):
        return row
    row.update(measures(w["rr"]))
    row["R_mV"], row["T_mV"] = amplitudes(w)
    row["n_V"] = int(np.sum(w["beat_types"] == "V"))
    row["n_A"] = int(np.sum(w["beat_types"] == "A"))
    if beats == "detected":
        row["missed"], row["extra"] = w["n_missed"], w["n_extra"]
        row["rr_abs_err_ms_med"] = float(np.median(np.abs(w["rr"] - w["rr_true"])) * 1e3)
    if "noise" in w:
        row["noise_gain"] = w["noise"]["gain"]
    if "apd_true_s" in w:
        row["apd_sd_ms"] = float(np.std(w["apd_true_s"]) * 1e3)
    return row


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    PLOTS.mkdir(exist_ok=True)
    src = R9.sources()
    tasks = [(s, *v, seed) for s in src for v in variants(s) for seed in SEEDS]
    by = {(s["source"], s["name"], s.get("ectopy_fixed") or "none"): s for s in src}
    for g in GRID_SOURCES:
        s = by[g]
        for nz in GRID_NOISE:
            tasks += [(s, f"grid_{nz}", g[2], nz, "detected", False, seed) for seed in SEEDS]
    tasks.sort(key=lambda t: t[0]["family"] != "coupled_vdp")          # slow ones first
    with Pool(a.workers) as pool:
        rows = list(pool.imap_unordered(one, tasks, chunksize=2))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "realism_raw.csv", index=False)
    num = [c for c in df.columns if df[c].dtype.kind in "fi" and c != "seed"]
    med = df.groupby(["family", "split", "regime", "label", "variant"], sort=False)[num].median().reset_index()
    med.to_csv(OUT / "realism.csv", index=False)
    cols = ["family", "split", "regime", "label", "variant", "mean_rr_s", "min_rr_s", "max_rr_s", "sdnn_ms", "rmssd_ms",
            "lf_hf", "R_mV", "T_mV", "n_V", "n_A", "missed", "extra", "rr_abs_err_ms_med", "apd_sd_ms"]
    cols = [c for c in cols if c in med.columns]
    with open(OUT / "realism.md", "w") as fh:
        fh.write("# Phase 9 realism (median over 5 realism seeds 980000-980004, 512 beats; no detector or measure)\n\n")
        fh.write("Variants: clean = true beats, no noise; ect = + E3 10 % couplets (models); real = ect + nstdb mix "
                 "12 dB + Pan-Tompkins detected beats + Phase 7 V jitter; grid_* = noise grid, detected beats.\n"
                 "Targets (Phase 5): mean RR 0.6-1.0 s; short-term SDNN 50 +/- 16 ms, RMSSD 42 +/- 15 ms "
                 "(Nunan et al. 2010); LF/HF 1.5-2.0 (Task Force 1996). ECGSYN normal beat at 60 bpm: R 1.2 mV "
                 "(peak), T 0.41 mV (verify_ecgsyn.py). Errors: " + str(int(df["error"].notna().sum())) + " windows.\n\n")
        fh.write(med[cols].round(3).to_markdown(index=False))
        fh.write("\n")
    plots(src)
    print(med[cols].round(3).to_string())


def _strip(ax, w, secs, title):
    s = int(w["true_idx"][L.LEAD] - 0.3 * FS)
    e = s + int(secs * FS)
    t = np.arange(e - s) / FS
    ax.plot(t, w["ecg"][s:e], color="#2a78d6", lw=0.7)
    for i, ty in zip(w["true_idx"], w["types_true"]):
        if s <= i < e and ty != "N":
            ax.text((i - s) / FS, ax.get_ylim()[1] if False else w["ecg"][i] + 0.25, ty, fontsize=6, color="#c23b22", ha="center")
    ax.set_title(title, fontsize=7, loc="left")
    ax.tick_params(labelsize=6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def plots(src):
    fams = {}
    for s in src:
        fams.setdefault(s["family"], []).append(s)
    for fam, ss in fams.items():
        pick = ss if len(ss) <= 6 else [x for i, x in enumerate(ss) if i % max(1, len(ss) // 6) == 0][:6]
        fig, axes = plt.subplots(len(pick), 3, figsize=(14, 1.6 * len(pick)), squeeze=False)
        for r, s in enumerate(pick):
            vs = variants(s)
            for c, v in enumerate(vs[:3] if len(vs) == 3 else [vs[0], vs[0], vs[1]]):
                vname, ect, noise, beats, jit = v
                w = L.build(s["source"], s["name"], ect, N, SEEDS[0], noise, "DEV", beats, jit)
                _strip(axes[r, c], w, 8.0, f"{R9.key(s)} [{s['label']}] {vname}")
        fig.tight_layout()
        fig.savefig(PLOTS / f"A_ecg_{fam}.png", dpi=100)
        plt.close(fig)
    # noise grid and RR tachograms
    by = {(s["source"], s["name"], s.get("ectopy_fixed") or "none"): s for s in src}
    fig, axes = plt.subplots(len(GRID_NOISE) + 1, 1, figsize=(12, 1.5 * (len(GRID_NOISE) + 1)))
    s = by[GRID_SOURCES[1]]
    w = L.build("model", s["name"], "E3_10", N, SEEDS[0], "none", "DEV", "true")
    _strip(axes[0], w, 10.0, "phase_reset tau=1.2, E3_10 ectopy, clean")
    for ax, nz in zip(axes[1:], GRID_NOISE):
        w = L.build("model", s["name"], "E3_10", N, SEEDS[0], nz, "DEV", "detected", True)
        _strip(ax, w, 10.0, f"... + nstdb {nz} (DEV signal 0)")
    fig.tight_layout()
    fig.savefig(PLOTS / "A_noise_grid.png", dpi=100)
    plt.close(fig)
    fig, axes = plt.subplots(3, 1, figsize=(12, 6))
    for ax, (g, nz) in zip(axes, ((GRID_SOURCES[0], "mix12"), (GRID_SOURCES[1], "mix6"), (GRID_SOURCES[2], "mix12"))):
        s = by[g]
        wt = L.build(s["source"], s["name"], g[2], N, SEEDS[0], "none", "DEV", "true")
        wd = L.build(s["source"], s["name"], g[2], N, SEEDS[0], nz, "DEV", "detected", True)
        ax.plot(wt["rr"], color="0.5", lw=1.5, label="true beats")
        ax.plot(wd["rr"], color="#c23b22", lw=0.7, label=f"detected, {nz}")
        ax.set_title(R9.key(s), fontsize=8, loc="left")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(PLOTS / "A_rr_true_vs_detected.png", dpi=100)
    plt.close(fig)


if __name__ == "__main__":
    main()
