"""
Phase 9 Part C: per-beat features from the delineated ECG, multivariate embedding and multivariate
IAAFT surrogates.

C1 delineation (predeclared): NeuroKit2 0.2.13 ecg_delineate(method='dwt') (wavelet delineator of
Martinez et al. 2004; verify_delineation.py checks it against the published QT-database accuracy),
run on the recorded ECG (360 Hz) with the window's beat sample indices as R peaks.  Per beat k:
  RR_k     interval from beat k-1 to beat k (s)
  QT_k     QRS onset (ECG_R_Onsets) to T end (ECG_T_Offsets) (s)
  RT_k     R peak to T peak (s) (secondary; robust when onsets fail)
  QRSd_k   QRS onset to QRS offset (ECG_R_Offsets) (s)
  QRSa_k   QRS amplitude: max - min of the ECG within +/- 60 ms of the R peak (mV)
  Ta_k     T amplitude: ECG at T peak minus the median ECG in the 40 ms before QRS onset (mV)
A feature value that is NaN (delineation failure) is replaced by linear interpolation over beats
(at most 10 % of beats per feature, otherwise the feature is flagged unusable for the window).

Ground truth on the synthetic ECG (for the accuracy report only): from the generator's event
parameters (ecgsyn.py, morphology.py): event time t_i = t_R + theta_i / omega, Gaussian width in time
sigma_i = b_i / omega (omega of the half-cycle containing the event); QRS onset = t_Q - 2 sigma_Q,
QRS offset = t_S + 2 sigma_S, T end = t_T + 2 sigma_T (the tangent-method end of a Gaussian wave),
T peak = t_T.

C2 multivariate IAAFT (Schreiber T, Schmitz A, Physica D 142:346 (2000), arXiv:chao-dyn/9909037,
Sect. 4.6, Eqs. 18-20): per channel the amplitude step is rank ordering (exact distributions); in
the spectral step every channel gets the data's Fourier amplitudes |S_k,m| and phases
rho_k,m + alpha_k with tan alpha_k = sum_m sin(psi_k,m - rho_k,m) / sum_m cos(psi_k,m - rho_k,m)
(psi = current surrogate phases), which keeps the data's relative phases (cross-spectrum).
Iterated until the rank orders no longer change (max 500 iterations); the last step is the rank
step.  Channels are standardized before and restored after.
NULL HYPOTHESIS: the multivariate series is a stationary multivariate linear Gaussian process
(given auto- and cross-spectra) observed through an invertible static (monotone) transformation of
each channel.  Rejection rules out only that null: it does not establish chaos (nonlinear stochastic
processes, ectopic beat patterns, static nonlinear couplings between features, nonstationarity all
reject it).
"""
from __future__ import annotations

import numpy as np


def mv_iaaft(X, rng, max_iter=500):
    """Multivariate IAAFT surrogate of X (N x M)."""
    X = np.asarray(X, float)
    N, Mc = X.shape
    mu, sd = X.mean(0), X.std(0)
    sd[sd == 0] = 1.0
    Z = (X - mu) / sd
    srt = np.sort(Z, axis=0)
    F = np.fft.rfft(Z, axis=0)
    amp = np.abs(F)
    rho = np.angle(F)
    # start: independent random permutation of each channel
    S = np.stack([rng.permutation(Z[:, m]) for m in range(Mc)], axis=1)
    prev = None
    for _ in range(max_iter):
        G = np.fft.rfft(S, axis=0)
        psi = np.angle(G)
        alpha = np.arctan2(np.sum(np.sin(psi - rho), axis=1), np.sum(np.cos(psi - rho), axis=1))
        Snew = np.fft.irfft(amp * np.exp(1j * (rho + alpha[:, None])), n=N, axis=0)
        ranks = np.argsort(np.argsort(Snew, axis=0), axis=0)
        S = np.take_along_axis(srt, ranks, axis=0)
        if prev is not None and np.array_equal(ranks, prev):
            break
        prev = ranks
    return S * sd + mu


def embed_mv(X, m, tau=1):
    """Multivariate delay embedding: each channel with dimension m (same for all), delay tau;
    row t holds (X[t], X[t+tau], ..., X[t+(m-1)tau]) for every channel."""
    X = np.asarray(X, float)
    n = len(X) - (m - 1) * tau
    return np.concatenate([X[k * tau:k * tau + n] for k in range(m)], axis=1)


def mv_nlp_error(X, m=2, tau=1, k=5, theiler=5, target=0):
    """Leave-one-out local-average prediction of channel `target` one beat ahead from the
    multivariate delay vector (channels standardized), robust normalized error (median |e| / MAD), as
    the Phase 8 univariate statistic."""
    X = np.asarray(X, float)
    Z = (X - X.mean(0)) / np.maximum(X.std(0), 1e-12)
    E = embed_mv(Z, m, tau)
    n = len(E) - 1
    A = E[:n]
    y = Z[(m - 1) * tau + 1:(m - 1) * tau + 1 + n, target]
    d = ((A[:, None, :] - A[None, :, :]) ** 2).sum(-1)
    idx = np.arange(n)
    d[np.abs(idx[:, None] - idx[None, :]) <= theiler] = np.inf
    nn = np.argpartition(d, k, axis=1)[:, :k]
    e = y - y[nn].mean(1)
    zt = Z[:, target]
    return float(np.median(np.abs(e)) / max(np.median(np.abs(zt - np.median(zt))), 1e-12))


# ------------------------------------------------------------------------------- C1 features
FEATURES = ("RR", "QT", "RT", "QRSd", "QRSa", "Ta")


def delineate(ecg, beat_idx, fs=360):
    """NeuroKit2 DWT delineation at the given R positions; returns per-beat sample positions (float,
    NaN where missing) of QRS onset/offset, T peak, T end."""
    import warnings
    import neurokit2 as nk
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _, w = nk.ecg_delineate(np.asarray(ecg, float), np.asarray(beat_idx, int), sampling_rate=fs, method="dwt")
    key = {"qrs_on": "ECG_R_Onsets", "qrs_off": "ECG_R_Offsets", "t_peak": "ECG_T_Peaks", "t_end": "ECG_T_Offsets"}
    return {k: np.asarray(w[v], dtype=float) for k, v in key.items()}


def _interp_nan(v, max_frac=0.10):
    v = np.asarray(v, float).copy()
    bad = ~np.isfinite(v)
    if bad.mean() > max_frac or (~bad).sum() < 2:
        return v, False
    if bad.any():
        idx = np.arange(len(v))
        v[bad] = np.interp(idx[bad], idx[~bad], v[~bad])
    return v, True


def features(ecg, beat_idx, fs=360, return_raw=False):
    """Per-beat feature matrix (n x 6) for beats 1..n of beat_idx (n = len(beat_idx) - 1): RR (s), QT (s),
    RT (s), QRSd (s), QRSa (mV), Ta (mV); dict of per-feature usability flags; optionally the
    delineation."""
    ecg = np.asarray(ecg, float)
    bi = np.asarray(beat_idx, int)
    d = delineate(ecg, bi, fs)
    n = len(bi) - 1
    w60 = int(round(0.06 * fs))
    w40 = int(round(0.04 * fs))
    qrsa = np.array([np.ptp(ecg[max(0, i - w60):i + w60 + 1]) for i in bi])
    base = np.full(len(bi), np.nan)
    ta = np.full(len(bi), np.nan)
    for k, i in enumerate(bi):
        on = d["qrs_on"][k]
        o = int(on) if np.isfinite(on) else i - int(0.08 * fs)
        base[k] = np.median(ecg[max(0, o - w40):max(1, o)]) if o - w40 >= 0 else np.nan
        tp = d["t_peak"][k]
        if np.isfinite(tp) and np.isfinite(base[k]) and int(tp) < len(ecg):
            ta[k] = ecg[int(tp)] - base[k]
    raw = {"RR": np.diff(bi) / fs,
           "QT": (d["t_end"] - d["qrs_on"])[1:] / fs,
           "RT": (d["t_peak"] - bi)[1:] / fs,
           "QRSd": (d["qrs_off"] - d["qrs_on"])[1:] / fs,
           "QRSa": qrsa[1:],
           "Ta": ta[1:]}
    cols, ok = [], {}
    for f in FEATURES:
        v, good = _interp_nan(raw[f])
        cols.append(v)
        ok[f] = good
    X = np.stack(cols, axis=1)
    if return_raw:
        return X, ok, d
    return X, ok
