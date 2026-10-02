"""
Phase 9 Part C1 (revised): per-beat repolarization features WITHOUT a delineator.

Why: the only published delineator available as verified open-source code, NeuroKit2 0.2.13
ecg_delineate(method='dwt') (Martinez et al. 2004), FAILED its predeclared QTDB verification
(results/verification/delineation_qtdb.json: best-lead QRS onset -15.5 +/- 37.2 ms vs published
4.6 +/- 7.7; T end -20.9 +/- 34.9 ms, sensitivity 90 %), so it is not used (METHODS.md).  ecgpuwave
(Laguna et al. 1994) needs a Fortran compiler, which this container does not have.  Features that
need QRS onset / offset or T end (QT, QRSd) are therefore DROPPED.

Fixed-window measurements (predeclared; no fitted parameter, the classical R-to-T-apex and T-amplitude
measurements of the ECG):
  baseline_k  median ECG over [R_k - 0.5 RR_prev, R_k - 0.3 RR_prev] (the TP segment of ECGSYN
              beats: previous T end ~0.39 RR after its R, P onset ~0.27 RR before R)
  T apex      sample of max |ECG - baseline_k| in [R_k + 80 ms, R_k + 0.6 RR_next]
  RTp_k       T apex - R_k (s)
  Ta_k        ECG(T apex) - baseline_k (mV, signed; inverted for the discordant T of PVCs)
  QRSa_k      max - min of the ECG within +/- 60 ms of R_k (mV)
with R_k the beat positions used for the RR series (true or detected) and RR_prev / RR_next the
intervals before and after R_k.  Beat k of the feature series is the beat that ENDS RR interval k
(as multivariate.features), so the feature series are aligned with the RR window.
Accuracy is verified against generator ground truth (verify_beatfeat.py) before use.
"""
from __future__ import annotations

import numpy as np

FEATS = ("RR", "RTp", "Ta", "QRSa")


def beat_features(ecg, beat_idx, fs=360):
    """(n x 4) array (RR, RTp, Ta, QRSa) for the n = len(beat_idx) - 1 intervals; RTp/Ta of the
    last beat use the mean RR as RR_next."""
    ecg = np.asarray(ecg, float)
    bi = np.asarray(beat_idx, int)
    K = len(bi)
    rr = np.diff(bi)
    w60 = int(round(0.06 * fs))
    w80 = int(round(0.08 * fs))
    mean_rr = float(np.mean(rr))
    rtp = np.full(K, np.nan)
    ta = np.full(K, np.nan)
    qa = np.full(K, np.nan)
    for k in range(1, K):
        i = bi[k]
        rp = rr[k - 1]
        rn = rr[k] if k < K - 1 else mean_rr
        b0, b1 = int(i - 0.5 * rp), int(i - 0.3 * rp)
        if b0 < 0 or b1 <= b0:
            continue
        base = float(np.median(ecg[b0:b1]))
        s, e = i + w80, int(i + 0.6 * rn)
        if e <= s or e > len(ecg):
            continue
        seg = ecg[s:e] - base
        j = int(np.argmax(np.abs(seg)))
        rtp[k] = (s + j - i) / fs
        ta[k] = seg[j]
        qa[k] = float(np.ptp(ecg[max(0, i - w60):i + w60 + 1]))
    X = np.stack([rr / fs, rtp[1:], ta[1:], qa[1:]], axis=1)
    return X
