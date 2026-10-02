"""
Phase 9 Part C1: accuracy of the fixed-window repolarization features (beatfeat.py) against the
generator's ground truth, on DEVELOPMENT material only (seeds 9000-9004).

    python -m experiments.phase9_waveform.verify_beatfeat

Predeclared criteria (written before this script was run):
 F1 T apex time (R_used + RTp) vs the analytic T event time t_T (library.true_fiducials), normal beats,
    nstdb mix 12 dB, detected beats (Pan-Tompkins), sources N1_linear_rr, E3_couplets_10pct,
    E5_atrial_10pct, mackey_glass tau=23 + E1: |mean| <= 13.9 ms and SD <= 27.8 ms (the published T-peak
    SD of Martinez et al. 2004 on QTDB, and 2 x that SD) and >= 95 % of beats within 150 ms.
 F2 the same on clean windows with true beats: |mean| <= 5 ms, SD <= 5 ms.
 F3 Ta at mix 12 dB vs Ta on the clean ECG of the same window (true beats, normal beats): median over
    windows of (median |difference| / median |Ta_clean|) <=
    0.1.
 F4 APD tracking: N1 RR with an imposed APD sequence APD_k = APD_ref (1 + 0.05 sin(2 pi k / 7.3))
    (APD_ref = 0.25 s; morphology.beat_params A2 rule), clean and mix 12: regression slope of
    measured RTp deviation on the imposed APD deviation within [0.8, 1.2], correlation >= 0.9 (clean)
    and >= 0.5 (mix 12).
Ventricular beats are reported without a criterion.
"""
from __future__ import annotations

import json
import pathlib

import numpy as np

from experiments.phase9_waveform import beatfeat as BF
from experiments.phase9_waveform import ecgsyn as E
from experiments.phase9_waveform import families as F
from experiments.phase9_waveform import library as L
from experiments.phase9_waveform import morphology as MO
from experiments.phase9_waveform import noise as NZ

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "results" / "verification"
SRC = [("null", "N1_linear_rr", "none"), ("null", "E3_couplets_10pct", "none"),
       ("null", "E5_atrial_10pct", "none"), ("model", "mackey_glass:tau=23.0", "E1")]
SEEDS = range(9000, 9005)


def apex_errors(w, X):
    """T apex absolute time error (ms) per window beat matched to a true beat, and its type."""
    tf = w["truth"]
    true_idx = w["true_idx"]
    bi = np.asarray(w["beat_idx"])
    errs, types = [], []
    mp = L.match(bi, true_idx, int(round(L.MATCH_TOL_S * L.FS)))
    for k in range(1, len(bi)):
        j = mp[k]
        if j < 0 or not np.isfinite(X[k - 1, 1]) or not np.isfinite(tf["t_peak"][j]):
            continue
        apex = bi[k] / L.FS + X[k - 1, 1]
        errs.append(1000 * (apex - tf["t_peak"][j]))
        types.append(w["types_true"][j])
    return np.array(errs), np.array(types)


def stats(e):
    return {"n": int(len(e)), "mean_ms": float(np.mean(e)), "sd_ms": float(np.std(e, ddof=1)),
            "frac_within_150ms": float(np.mean(np.abs(e) <= 150))}


def apd_check(seed, noise):
    """F4: manual window as library.build, with an imposed sinusoidal APD sequence."""
    n = 512
    rr, types = F.null_window("N1_linear_rr", seed, n + L.MARGIN)
    rng = np.random.default_rng(seed + 31)
    m = float(np.mean(rr[:n]))
    rr_full = np.concatenate([np.full(L.LEAD, m), rr])
    types_full = np.array(["N"] * len(rr_full) + ["N"])
    K = len(rr_full) + 1
    apd_ref = 0.25
    apd = apd_ref * (1 + 0.05 * np.sin(2 * np.pi * np.arange(K) / 7.3))
    bt = np.concatenate([[0.0], np.cumsum(rr_full)])
    th, a, b = MO.beat_params(rr_full, types_full[:K], apd, apd_ref)
    bt_ext = np.concatenate([bt, [bt[-1] + m]])
    th = np.vstack([th, th[-1:]]); a = np.vstack([a, a[-1:]]); b = np.vstack([b, b[-1:]])
    z = E.synthesize(bt_ext, th, a, b, fs=L.FS, q=2, t0=0.0, t1=bt[-1] + 0.5 * m, resp_phase=0.3)
    ecg = E.to_mV(z)
    idx = np.round(bt * L.FS).astype(int)
    idx = idx[idx < len(ecg)]
    if noise != "none":
        ecg, _ = NZ.add_noise(ecg, idx[L.LEAD:], "mix", 12, rng, 0)
    ecg = NZ.quantize_mitbih(ecg)
    sel = idx[L.LEAD:L.LEAD + n + 1]
    X = BF.beat_features(ecg, sel, L.FS)
    dev_meas = X[:, 1] - np.nanmean(X[:, 1])
    dev_true = apd[L.LEAD + 1:L.LEAD + n + 1] - apd_ref
    ok = np.isfinite(dev_meas)
    slope = float(np.polyfit(dev_true[ok], dev_meas[ok], 1)[0])
    return slope, float(np.corrcoef(dev_true[ok], dev_meas[ok])[0, 1])


def main():
    res = {"F1": {}, "F2": {}, "F3": {}, "F4": {}}
    allN1, allV1, allN2 = [], [], []
    ta_ratio = []
    for src, name, ect in SRC:
        for s in SEEDS:
            wc = L.build(src, name, ect, 512, s, "none", "DEV", "true", False, keep_truth=True)
            Xc = BF.beat_features(wc["ecg"], wc["beat_idx"], L.FS)
            e2, t2 = apex_errors(wc, Xc)
            allN2 += list(e2[t2 == "N"])
            wn = L.build(src, name, ect, 512, s, "mix12", "DEV", "detected", False, keep_truth=True)
            Xn = BF.beat_features(wn["ecg"], wn["beat_idx"], L.FS)
            e1, t1 = apex_errors(wn, Xn)
            allN1 += list(e1[t1 == "N"])
            allV1 += list(e1[t1 == "V"])
            # F3: same window, true beats, noisy vs clean ECG
            wt = L.build(src, name, ect, 512, s, "mix12", "DEV", "true", False, keep_clean=True)
            Xa = BF.beat_features(wt["ecg"], wt["beat_idx"], L.FS)
            Xb = BF.beat_features(NZ.quantize_mitbih(wt["ecg_clean"]), wt["beat_idx"], L.FS)
            nb = wt["beat_types"][1:] == "N"
            d = np.abs(Xa[nb, 2] - Xb[nb, 2])
            ta_ratio.append(float(np.nanmedian(d) / np.nanmedian(np.abs(Xb[nb, 2]))))
    res["F1"] = {"N_mix12_detected": stats(np.array(allN1)), "V_mix12_detected": stats(np.array(allV1))}
    res["F2"] = {"N_clean_true": stats(np.array(allN2))}
    res["F3"] = {"median_abs_diff_over_median_abs_Ta_per_window_max": float(np.max(ta_ratio)),
                 "median_over_windows": float(np.median(ta_ratio))}
    f4 = {nz: [apd_check(s, nz) for s in SEEDS] for nz in ("none", "mix12")}
    res["F4"] = {nz: {"slopes": [v[0] for v in f4[nz]], "corrs": [v[1] for v in f4[nz]]} for nz in f4}
    f1 = res["F1"]["N_mix12_detected"]
    f2 = res["F2"]["N_clean_true"]
    ok = (abs(f1["mean_ms"]) <= 13.9 and f1["sd_ms"] <= 27.8 and f1["frac_within_150ms"] >= 0.95
          and abs(f2["mean_ms"]) <= 5 and f2["sd_ms"] <= 5
          and np.median(ta_ratio) <= 0.1
          and all(0.8 <= v[0] <= 1.2 and v[1] >= 0.9 for v in f4["none"])
          and all(0.8 <= v[0] <= 1.2 and v[1] >= 0.5 for v in f4["mix12"]))
    res["VERIFIED"] = bool(ok)
    OUT.mkdir(parents=True, exist_ok=True)
    json.dump(res, open(OUT / "beatfeat.json", "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
