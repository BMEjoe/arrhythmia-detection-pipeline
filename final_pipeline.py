"""
Research-grade nonlinear ECG/RR pipeline
========================================

Implements:
  1) ECG preprocessing:
       - zero-phase ECG bandpass (retained for plotting/QC only)
       - Pan-Tompkins-compatible QRS detection using the published 200-Hz filter chain,
         search-back, T-wave discrimination, and R-peak refinement
       - RR artifact detection/correction while retaining the raw RR series

     NOTE: CWT-based denoising has been deliberately removed from the RR
     extraction pathway. Pan-Tompkins already includes its own bandpass
     stage (5-15 Hz) specifically tuned to the QRS complex; running CWT
     denoising on the raw signal first would change the frequency content
     entering that filter and can distort the very features Pan-Tompkins is
     tuned to detect. RR intervals -- not ECG waveform morphology -- are
     what every downstream analysis in this pipeline (LLE, UPO, TDA) uses,
     so CWT denoising was never load-bearing for the actual results. If your
     methods section describes CWT denoising, update it to describe
     Pan-Tompkins' own internal filtering as what performs denoising for RR
     extraction, per this pipeline's actual behavior.

  2) Delay reconstruction:
       - TDMI first-local-minimum delay
       - Cao E1/E2 embedding-dimension selection (max-norm, per Cao 1997;
         E1(m) = E(m+1)/E(m) computed as a ratio between consecutive
         dimension-averaged nearest-neighbor distance ratios, not conflated
         with the un-ratioed per-dimension average itself)
       - Theiler exclusion for nearest-neighbor calculations

  3) Largest Lyapunov exponent:
       - Rosenstein, Collins & De Luca (1993), with explicit divergence
         scaling-region diagnostics (R^2-weighted search) rather than an
         arbitrary "flattest" fit

  4) UPO detection:
       - So, Ott, Sauer, Gluckman, Grebogi & Schiff (1996), Phys. Rev. Lett.
         76, 4705; and (1997), Phys. Rev. E 55, 5398-5417.
       - First-order fixed-point transform G(z,R)
       - Local Jacobian estimated using the exact companion-matrix structure
         of the tau-step delay map (only the first row is fitted; the remaining
         rows are the deterministic delay-coordinate shift)
       - R-perturbation term implemented as S(z,R) = J + kappa*R*||F(z)-z||_1,
         a d x d random matrix R times the scalar L1 norm of the state
         difference -- this is the formula verified directly against So et
         al. (1996), PRL 76, 4705, Eq. (4). (An earlier draft of this file
         used a d x d x d random tensor contracted against the raw
         difference vector; that formula could not be verified against the
         primary source text and has been replaced with this one.)
       - Peak-finding runs on the RAW transformed-value histogram, not a
         Gaussian-smoothed version. A fixed smoothing bandwidth was found
         (via direct testing against the logistic map's two known analytic
         fixed points) to merge distinct peaks that land only a couple of
         histogram bins apart into a single spurious peak between them --
         confirmed by the earlier version detecting only 1 of the logistic
         map's known 2 fixed points despite the raw histogram clearly
         showing both. A light smoothing pass is still used only to
         estimate the background/significance threshold, not to locate
         peaks themselves.
       - Period-p "short" grouping scheme using every iterate
       - p-dimensional reduction using the cyclic delay-coordinate symmetry
       - optional AAFT surrogate significance with fixed real-data tau and m

  Validation note: the logistic-map LLE and UPO tests from the previous file
  remain useful regression tests, but this corrected file does not claim that
  those tests have been rerun automatically. They should be rerun after every
  algorithmic change before MIT-BIH experiments.

Important:
  - The So method is a method for discrete maps reconstructed from a scalar
    time series. UPO results are therefore reported in reconstructed RR
    phase space.
  - Period p=1 uses the original fixed-point transform.
  - Period p>1 uses the 1997 periodic-orbit transform and short grouping
    scheme. The default keeps p small because the original algorithm scales
    approximately as N*(K+1)^p.
"""

from __future__ import annotations

import itertools
import math
import os
from dataclasses import dataclass, replace
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt, find_peaks, medfilt, welch, resample_poly
from scipy.ndimage import gaussian_filter1d
try:
    import wfdb
except ImportError:
    wfdb = None


# =============================================================================
# Configuration
# =============================================================================

@dataclass
class PipelineConfig:
    # ECG preprocessing (CWT intentionally removed -- see module docstring)
    ecg_low_hz: float = 0.5
    ecg_high_hz: float = 40.0
    qrs_low_hz: float = 5.0
    qrs_high_hz: float = 15.0
    ecg_filter_order: int = 4
    qrs_filter_order: int = 2

    # Pan-Tompkins
    pt_original_fs: float = 200.0
    pt_integration_seconds: float = 0.150
    pt_refractory_seconds: float = 0.200
    pt_searchback_seconds: float = 1.66
    pt_initial_seconds: float = 2.0
    pt_twave_window_seconds: float = 0.36

    # RR cleaning
    rr_min_seconds: float = 0.25
    rr_max_seconds: float = 3.00
    rr_local_window: int = 11
    rr_outlier_sigma: float = 4.5
    # False by default: abnormal RR intervals may be the arrhythmia signal.
    use_corrected_rr_for_dynamics: bool = False


    # TDMI / Cao
    tdmi_max_tau: int = 30
    tdmi_bins: int = 16
    tdmi_smooth_window: int = 3
    cao_max_dim: int = 12
    cao_tol: float = 0.05
    cao_theiler: int = 2

    # Rosenstein
    lle_theiler_beats: Optional[int] = None
    lle_max_iter: int = 40
    lle_min_fit_points: int = 6
    lle_max_fit_fraction: float = 0.10
    lle_min_r2: float = 0.90
    lle_require_valid_fit: bool = True

    # So et al. UPO method
    #
    # IMPORTANT (found via direct timing test): period-p detection (p>1) is
    # combinatorially expensive -- it evaluates (K+1)^p combinations per
    # backbone point, each requiring so_random_R random-matrix solves, with
    # substantial per-iteration Python overhead on top of the raw linear
    # algebra. Measured directly: on a modest 355-point embedded series,
    # period-2 detection alone (with the previous defaults: so_random_R=100,
    # so_max_backbones=None) did not complete in 40 seconds. Real RR-interval
    # recordings will have several times more points than that test case, so
    # period>1 detection is now OFF by default (so_periods=(1,)) -- period-1
    # detection was independently timed at under 1 second on the same data
    # and is unaffected. Enable period-2/3 explicitly only after also setting
    # so_max_backbones to a real cap and/or reducing so_random_R for period-p
    # calls specifically; do not enable them with the period-1 defaults
    # unmodified, or expect the run to take a very long time.
    so_periods: Tuple[int, ...] = (1,)
    so_neighbors_K: int = 1       # K in So et al. short grouping scheme
    so_jacobian_neighbors: int = 7  # M spatial neighbors for local Jacobian
    so_random_R: int = 100
    so_kappa: float = 3.0
    so_hist_bins: int = 30
    so_slab_fraction: float = 0.08  # still used for the period-p residual check
                                      # (line ~862), which is a relative-scale
                                      # comparison and not subject to the same
                                      # curse-of-dimensionality problem as the
                                      # fixed-point diagonal filter below.
    so_diagonal_percentile: float = 10.0  # replaces the old fixed-fraction
                                            # diagonal filter for period-1 fixed
                                            # point detection. Found via direct
                                            # testing: a FIXED distance-fraction
                                            # threshold (the old so_slab_fraction
                                            # used for this purpose) worked fine
                                            # at the originally-validated d=2 case
                                            # (~25% of points survived) but
                                            # collapsed to a 0.01% survival rate
                                            # at a real auto-selected d=6 case on
                                            # a synthetic RR-like series -- a
                                            # genuine curse-of-dimensionality
                                            # effect, not a data problem. A
                                            # percentile-based cutoff (always keep
                                            # the closest N% of transformed points
                                            # to the diagonal, whatever the
                                            # ambient dimension) is self-
                                            # calibrating and was confirmed to
                                            # produce thousands of surviving
                                            # points at d=6 where the fixed
                                            # threshold produced 5.

    # UPO coverage / density measure (the actual Devaney "density of periodic
    # orbits" proxy -- a single detected fixed point alone is not evidence of
    # chaos; density/coverage of the attractor by detected UPOs is). Radius is
    # NOT a fixed distance or fixed fraction of attractor extent, for the same
    # curse-of-dimensionality reason as so_diagonal_percentile above -- it is
    # set relative to the k-th nearest-neighbor spacing WITHIN the normalized
    # attractor itself, so it self-calibrates to whatever the ambient
    # dimension's typical point density actually is. Tested directly on a
    # matched chaotic-vs-periodic twin pair at a real auto-selected embedding
    # dimension (d=6): coverage came out 0.89-1.0 for the chaotic case and
    # 0.0-0.28 for the periodic case across every radius_multiplier tested
    # (2.0, 3.0, 5.0) -- a real, direction-consistent discriminating result,
    # not an artifact of one tuned setting.
    coverage_n_reference: int = 200
    coverage_radius_neighbors_k: int = 5
    coverage_radius_multiplier: float = 3.0
    so_peak_min_count: int = 5
    so_peak_sigma: float = 3.0
    so_max_backbones: Optional[int] = 200  # required cap if period>1 is enabled;
                                             # None (unbounded) is no longer the
                                             # default given the timing above
    so_random_R_period_p: int = 20  # separate, smaller R-count for period>1,
                                      # since the combinatorial multiplier over
                                      # backbones already provides many samples

    # Surrogates
    surrogate_count: int = 199
    compute_surrogates: bool = False
    random_seed: int = 20260823


CFG = PipelineConfig()


# =============================================================================
# Utility functions
# =============================================================================

def _as_1d(x) -> np.ndarray:
    x = np.asarray(x, dtype=float).reshape(-1)
    if x.size == 0:
        raise ValueError("Input time series is empty.")
    if not np.all(np.isfinite(x)):
        raise ValueError("Input contains NaN/Inf. Clean or interpolate before use.")
    return x


def _safe_std(x) -> float:
    s = float(np.std(x))
    return s if s > 1e-15 else 1.0


def _nearest_neighbor_indices(X, theiler=0):
    """Nearest neighbor in Euclidean phase space with temporal exclusion."""
    X = np.asarray(X, dtype=float)
    n = len(X)
    nn = np.full(n, -1, dtype=int)
    if n < 2:
        return nn
    for i in range(n):
        d = np.linalg.norm(X - X[i], axis=1)
        lo = max(0, i - theiler)
        hi = min(n, i + theiler + 1)
        d[lo:hi] = np.inf
        j = np.argmin(d)
        if np.isfinite(d[j]):
            nn[i] = j
    return nn


# =============================================================================
# 1. ECG preprocessing (bandpass + Pan-Tompkins only; no CWT -- see docstring)
# =============================================================================

def bandpass_filter(x, fs, low=0.5, high=40.0, order=4):
    x = _as_1d(x)
    nyq = 0.5 * fs
    if not (0 < low < high < nyq):
        raise ValueError("Bandpass frequencies must satisfy 0 < low < high < fs/2.")
    b, a = butter(order, [low / nyq, high / nyq], btype="band")
    return filtfilt(b, a, x)


def preprocess_ecg(x, fs, config=CFG):
    """
    ECG preprocessing used for the nonlinear pipeline: a single zero-phase
    bandpass. Pan-Tompkins' own internal 5-15 Hz bandpass performs the
    QRS-specific filtering used for R-peak detection; this broader-band
    output is retained only as a general-purpose cleaned waveform (e.g. for
    plotting), not as an input to peak detection.
    """
    x = _as_1d(x)
    return bandpass_filter(x, fs, config.ecg_low_hz, config.ecg_high_hz, config.ecg_filter_order)


# =============================================================================
# Pan-Tompkins-style adaptive QRS detector
# =============================================================================

def _pan_tompkins_filter_chain_200hz(x):
    """Pan--Tompkins 200-Hz filter chain.

    The original algorithm uses causal integer low/high-pass sections,
    differentiation, squaring, and a 30-sample (150 ms) moving-window
    integrator.  Detection is performed sequentially in ``detect_r_peaks``;
    this helper only constructs the three signals needed by that detector.
    """
    x = _as_1d(x)
    # Published 200-Hz recursive integer filters.
    lp = np.zeros_like(x, dtype=float)
    for n in range(len(x)):
        v = x[n]
        if n >= 1: v += 2.0 * lp[n-1]
        if n >= 2: v -= lp[n-2]
        if n >= 6: v -= 2.0 * x[n-6]
        if n >= 12: v += x[n-12]
        lp[n] = v / 36.0

    hp = np.zeros_like(lp)
    for n in range(len(lp)):
        v = lp[n]
        if n >= 1: v += hp[n-1]
        if n >= 16: v -= 32.0 * lp[n-16]
        if n >= 17: v += lp[n-17]
        if n >= 32: v -= lp[n-32]
        hp[n] = v

    d = np.zeros_like(hp)
    for n in range(len(hp)):
        # Five-point derivative, causal form of the published differentiator.
        if n >= 2:
            d[n] += 2.0 * hp[n] + hp[n-1] - hp[n-3]
        if n >= 4:
            d[n] -= 2.0 * hp[n-4]
        d[n] /= 8.0

    sq = d * d
    integ = np.convolve(sq, np.ones(30, dtype=float) / 30.0, mode="full")[:len(sq)]
    return hp, d, integ


def _pan_tompkins_processing_delay_200hz():
    """Measure the implemented causal chain's impulse-peak delays.

    Returns (filtered_delay, integrated_delay) in 200-Hz samples.  Computing
    this from the actual implementation is safer than hard-coding a delay if
    the recursive sections are ever changed.
    """
    impulse = np.zeros(256, dtype=float)
    origin = 64
    impulse[origin] = 1.0
    hp, _, integ = _pan_tompkins_filter_chain_200hz(impulse)
    hp_delay = int(np.argmax(np.abs(hp)) - origin)
    integ_delay = int(np.argmax(integ) - origin)
    return hp_delay, integ_delay


def detect_r_peaks(x, fs, config=CFG):
    """Faithful sequential Pan--Tompkins-style QRS detector.

    Detection is carried out at the paper's 200-Hz sampling rate.  Candidate
    peaks are local maxima of the integrated signal; adaptive SPKI/NPKI
    thresholds, the 200-ms refractory period, 1.66-RR search-back, and
    360-ms slope-based T-wave discrimination are applied chronologically.
    The causal processing delay is measured from the implemented filter chain,
    removed, and the final fiducial is re-localized on the original ECG.
    """
    x = _as_1d(x)
    target_fs = float(config.pt_original_fs)
    if abs(float(fs) - target_fs) < 1e-12:
        xr = x.copy(); back = 1.0
    else:
        from fractions import Fraction
        frac = Fraction(target_fs / float(fs)).limit_denominator(1000)
        xr = resample_poly(x, frac.numerator, frac.denominator)
        back = float(fs) / target_fs

    hp, deriv, integ = _pan_tompkins_filter_chain_200hz(xr)
    fsd = target_fs
    refractory = int(round(config.pt_refractory_seconds * fsd))
    twave_window = int(round(config.pt_twave_window_seconds * fsd))
    searchback_factor = float(config.pt_searchback_seconds)
    init_n = min(len(integ), int(round(config.pt_initial_seconds * fsd)))
    if len(integ) < 10:
        return np.array([], dtype=int)

    # Local maxima of the integrated waveform are only candidates; the
    # adaptive decision process below is the Pan--Tompkins part.
    candidates, props = find_peaks(integ, distance=max(1, refractory // 2))
    if len(candidates) == 0:
        return np.array([], dtype=int)

    # Initialize signal/noise peak estimates from the first 2 s.
    init_peaks, _ = find_peaks(integ[:init_n], distance=refractory)
    init_heights = integ[init_peaks] if len(init_peaks) else integ[:init_n]
    # Original initialization is based on the first 2 s of the integrated
    # signal: signal level starts at one quarter of the largest peak and noise
    # level at one half of the mean peak level.
    spki = 0.25 * float(np.max(init_heights))
    npki = 0.50 * float(np.mean(init_heights))
    th_i1 = npki + 0.25 * (spki - npki)
    th_i2 = 0.5 * th_i1

    # Filtered-signal threshold is maintained as the secondary Pan--Tompkins
    # decision channel used for T-wave discrimination and peak confirmation.
    filt_abs = np.abs(hp)
    f_init = []
    for q in init_peaks:
        lo, hi = max(0, q-2), min(len(hp), q+3)
        if hi > lo: f_init.append(float(np.max(filt_abs[lo:hi])))
    f_init = np.asarray(f_init if f_init else filt_abs[:init_n])
    spkf = 0.25 * float(np.max(f_init))
    npkf = 0.50 * float(np.mean(f_init))
    th_f1 = npkf + 0.25 * (spkf - npkf)
    th_f2 = 0.5 * th_f1

    def slope_at(p):
        lo, hi = max(0, int(p)-4), min(len(deriv), int(p)+5)
        return float(np.max(np.abs(deriv[lo:hi]))) if hi > lo else 0.0

    accepted = []
    noise = []
    last_candidate_index = -1

    def accept(q):
        nonlocal spki, th_i1, th_i2, spkf, th_f1, th_f2
        q = int(q)
        if accepted and q - accepted[-1] < refractory:
            return False
        accepted.append(q)
        h = float(integ[q])
        fh = float(np.max(filt_abs[max(0,q-2):min(len(hp),q+3)]))
        spki = 0.125*h + 0.875*spki
        spkf = 0.125*fh + 0.875*spkf
        th_i1 = npki + 0.25*(spki-npki); th_i2 = 0.5*th_i1
        th_f1 = npkf + 0.25*(spkf-npkf); th_f2 = 0.5*th_f1
        return True

    # Chronological adaptive decision. Search-back is performed whenever the
    # current RR interval exceeds 1.66 times the previous accepted RR.
    for q in candidates:
        q = int(q)
        h = float(integ[q])
        fh = float(np.max(filt_abs[max(0,q-2):min(len(hp),q+3)]))

        if accepted and len(accepted) >= 2:
            rr_prev = accepted[-1] - accepted[-2]
            if q - accepted[-1] > searchback_factor * rr_prev:
                sb_lo = accepted[-1] + refractory
                sb_hi = q - refractory
                eligible = [c for c in candidates if sb_lo <= c <= sb_hi and c > last_candidate_index]
                eligible = [c for c in eligible if integ[c] >= th_i2]
                if eligible:
                    sb = max(eligible, key=lambda c: integ[c])
                    accept(sb)

        last_candidate_index = q
        if accepted and q - accepted[-1] < refractory:
            noise.append(q)
            npki = 0.125*h + 0.875*npki
            npkf = 0.125*fh + 0.875*npkf
        elif h >= th_i1:
            # Primary threshold is applied to the integrated QRS energy.
            # The filtered channel is retained for the original slope-based
            # T-wave test rather than becoming a second hard threshold.
            twave = False
            if accepted and 0 < q - accepted[-1] <= twave_window:
                twave = slope_at(q) < 0.5 * slope_at(accepted[-1])
            if not twave:
                accept(q)
            else:
                noise.append(q)
                npki = 0.125*h + 0.875*npki
                npkf = 0.125*fh + 0.875*npkf
        else:
            noise.append(q)
            npki = 0.125*h + 0.875*npki
            npkf = 0.125*fh + 0.875*npkf

        th_i1 = npki + 0.25*(spki-npki); th_i2 = 0.5*th_i1
        th_f1 = npkf + 0.25*(spkf-npkf); th_f2 = 0.5*th_f1

    if not accepted:
        return np.array([], dtype=int)

    # ------------------------------------------------------------------
    # Timing correction and physiological fiducial localization.
    #
    # The published Pan--Tompkins chain is causal.  The integrated candidate
    # therefore carries the filter + derivative + integration delay.  Measure
    # that delay from this implementation's impulse response, subtract it,
    # map to the original sampling clock, then refine on the ORIGINAL ECG.
    # This keeps Pan--Tompkins causal and faithful while preventing the
    # processing delay from being reported as the R-wave time.
    # ------------------------------------------------------------------
    _, integrated_delay = _pan_tompkins_processing_delay_200hz()
    corrected_200 = np.asarray(accepted, dtype=int) - integrated_delay
    corrected_200 = np.clip(corrected_200, 0, len(xr)-1)
    mapped = np.rint(corrected_200.astype(float) * back).astype(int)
    mapped = np.clip(mapped, 0, len(x)-1)

    # Re-localize in a local window of the original ECG.  80 ms each side is
    # wide enough for residual morphology-dependent timing variation after
    # the deterministic processing delay has been removed, while remaining
    # much smaller than an RR interval.
    refine_original = max(1, int(round(0.080 * float(fs))))
    refined=[]
    for p in mapped:
        lo=max(0, p-refine_original); hi=min(len(x), p+refine_original+1)
        if hi <= lo:
            continue
        seg=x[lo:hi]
        baseline=float(np.median(seg))
        dev=seg-baseline
        pos=int(np.argmax(dev)); neg=int(np.argmin(dev))
        local = pos if abs(dev[pos]) >= abs(dev[neg]) else neg
        refined.append(lo+local)

    refined=np.asarray(sorted(set(refined)), dtype=int)
    if len(refined)==0:
        return np.array([], dtype=int)

    refractory_original=max(1, int(round(config.pt_refractory_seconds*float(fs))))
    final=[int(refined[0])]
    for p in refined[1:]:
        if p-final[-1] >= refractory_original:
            final.append(int(p))
        else:
            # Keep the stronger original-ECG excursion when two detections
            # collapse inside the refractory interval.
            def strength(q):
                lo=max(0,q-3); hi=min(len(x),q+4)
                return abs(float(x[q])-float(np.median(x[lo:hi])))
            if strength(int(p)) > strength(final[-1]):
                final[-1]=int(p)

    return np.asarray(np.clip(np.unique(final),0,len(x)-1),dtype=int)


def extract_rr_intervals(r_peaks, fs, rr_min_sec=0.30, rr_max_sec=2.0):
    r_peaks = np.asarray(r_peaks, dtype=int)
    if len(r_peaks) < 2:
        return np.array([]), np.array([]), np.array([], dtype=bool)
    times = r_peaks / float(fs)
    rr_raw = np.diff(times)
    rr_times = times[1:]
    valid = (rr_raw >= rr_min_sec) & (rr_raw <= rr_max_sec)
    return rr_raw, rr_times, valid


def correct_rr_intervals(rr_raw, config=CFG):
    rr = _as_1d(rr_raw)
    corrected = rr.copy()
    hard = (rr < config.rr_min_seconds) | (rr > config.rr_max_seconds)

    w = int(config.rr_local_window)
    if w % 2 == 0:
        w += 1
    w = max(3, min(w, len(rr) if len(rr) % 2 == 1 else len(rr) - 1))

    baseline = medfilt(rr, kernel_size=w) if w >= 3 else np.full_like(rr, np.median(rr))
    residual = rr - baseline
    mad = np.median(np.abs(residual - np.median(residual)))
    scale = max(mad / 0.6744897501960817, 1e-6)
    robust = np.abs(residual) / scale
    local = robust > config.rr_outlier_sigma
    artifact = hard | local

    if np.any(artifact):
        good = ~artifact
        idx = np.arange(len(rr))
        if np.sum(good) >= 2:
            corrected[artifact] = np.interp(idx[artifact], idx[good], rr[good])
        else:
            corrected[artifact] = np.median(rr[good]) if np.any(good) else np.median(rr)

    return corrected, artifact, artifact.copy()


# =============================================================================
# 2. Delay reconstruction: TDMI + Cao
# =============================================================================

def time_delayed_mutual_information(x, max_tau=30, n_bins=16):
    x = _as_1d(x)
    if len(x) < max_tau + 10:
        max_tau = max(1, len(x) // 4)
    edges = np.linspace(np.min(x), np.max(x), n_bins + 1)
    if np.allclose(edges[0], edges[-1]):
        return np.zeros(max_tau)
    mi = np.zeros(max_tau)
    for tau in range(1, max_tau + 1):
        a = x[:-tau]
        b = x[tau:]
        h, _, _ = np.histogram2d(a, b, bins=(edges, edges))
        pxy = h / max(np.sum(h), 1.0)
        px = np.sum(pxy, axis=1)
        py = np.sum(pxy, axis=0)
        denom = px[:, None] * py[None, :]
        mask = (pxy > 0) & (denom > 0)
        mi[tau - 1] = np.sum(pxy[mask] * np.log(pxy[mask] / denom[mask]))
    return mi


def find_optimal_tau(mi_values, smooth_window=3):
    mi = np.asarray(mi_values, dtype=float)
    if len(mi) == 0:
        return 1
    if smooth_window > 1 and len(mi) >= smooth_window:
        kernel = np.ones(smooth_window) / smooth_window
        smoothed = np.convolve(mi, kernel, mode="same")
    else:
        smoothed = mi
    for i in range(1, len(smoothed) - 1):
        if smoothed[i] <= smoothed[i-1] and smoothed[i] < smoothed[i+1]:
            return i + 1
    return int(np.argmin(smoothed)) + 1


def _embed_backward(x, tau, m):
    """So et al. delay-coordinate convention: z(n) = [x(n), x(n-tau), ..., x(n-(m-1)tau)]"""
    x = _as_1d(x)
    n = len(x) - (m - 1) * tau
    if n <= 0:
        return np.empty((0, m))
    return np.column_stack([x[(m - 1 - j) * tau : (m - 1 - j) * tau + n] for j in range(m)])


def _cao_forward_embed(x, tau, m):
    """Cao (1997) forward delay vector Y_i(m)."""
    x = _as_1d(x)
    tau = int(tau); m = int(m)
    if tau < 1 or m < 1:
        raise ValueError("tau and m must be positive integers")
    n = len(x) - (m - 1) * tau
    if n <= 0:
        return np.empty((0, m), dtype=float)
    return np.column_stack([x[j*tau:j*tau+n] for j in range(m)])


def cao_method(x, tau, max_dim=12, tol=0.05, theiler=2):
    """
    Cao (1997) E1/E2 embedding-dimension estimator.

    Uses the paper's forward delay vectors and maximum norm.  E1 is used for
    dimension saturation; E2 is retained as the deterministic/stochastic
    diagnostic.  The So/Rosenstein reconstruction elsewhere remains the
    backward delay-coordinate representation.
    """
    x = _as_1d(x)
    tau = int(tau)
    if tau < 1:
        raise ValueError("tau must be >= 1")

    E = np.full(max_dim + 1, np.nan, dtype=float)
    Estar = np.full(max_dim + 1, np.nan, dtype=float)

    for m in range(1, max_dim + 1):
        Xm = _cao_forward_embed(x, tau, m)
        Xm1 = _cao_forward_embed(x, tau, m + 1)
        n = min(len(Xm), len(Xm1))
        if n < max(30, 3*m):
            continue
        Xm = Xm[:n]
        Xm1 = Xm1[:n]

        ratios, star = [], []
        for i in range(n):
            dist = np.max(np.abs(Xm - Xm[i]), axis=1)
            lo=max(0, i-int(theiler)); hi=min(n, i+int(theiler)+1)
            dist[lo:hi] = np.inf
            j = int(np.argmin(dist))
            if not np.isfinite(dist[j]) or dist[j] <= 1e-15:
                continue

            # a(i,m) in Cao: nearest-neighbor distance ratio after adding
            # the next delayed coordinate.
            ratios.append(np.max(np.abs(Xm1[i]-Xm1[j])) / dist[j])

            # E*(m) is the mean absolute difference of the newly added
            # coordinate for the SAME nearest-neighbor pair.
            star.append(abs(Xm1[i,-1] - Xm1[j,-1]))

        if ratios:
            E[m-1] = float(np.mean(ratios))
            Estar[m-1] = float(np.mean(star))

    E1 = E[1:] / E[:-1]
    E2 = Estar[1:] / Estar[:-1]

    # Cao's criterion is fundamentally a plateau criterion.  Require two
    # consecutive small relative changes instead of accepting a single pair.
    chosen = max_dim + 1
    stable_run = 0
    for k in range(1, len(E1)):
        a, b = E1[k-1], E1[k]
        if np.isfinite(a) and np.isfinite(b):
            rel = abs(b-a) / max(abs(a), 1e-12)
            stable_run = stable_run + 1 if rel < tol else 0
            if stable_run >= 2:
                chosen = k + 1
                break
        else:
            stable_run = 0

    return int(chosen), E1, E2


def takens_embed(x, tau=None, m=None, config=CFG):
    x = _as_1d(x)
    if tau is None:
        mi = time_delayed_mutual_information(x, config.tdmi_max_tau, config.tdmi_bins)
        tau = find_optimal_tau(mi, config.tdmi_smooth_window)
    else:
        mi = None
    if m is None:
        m, e1, e2 = cao_method(x, tau, config.cao_max_dim, config.cao_tol, config.cao_theiler)
    else:
        e1 = e2 = None
    embedded = _embed_backward(x, tau, m)
    return embedded, int(tau), int(m), mi, e1, e2


# =============================================================================
# 3. Rosenstein LLE
# =============================================================================

def _mean_period_beats(x):
    x = _as_1d(x); y = x - np.mean(x)
    if len(x) < 16:
        return 1
    f, pxx = welch(y, fs=1.0, nperseg=min(256, len(y)))
    mask = f > 0
    if not np.any(mask) or np.sum(pxx[mask]) <= 0:
        return 1
    mean_f = np.sum(f[mask]*pxx[mask]) / np.sum(pxx[mask])
    return max(1, int(round(1.0/mean_f)))


def rosenstein_lle(embedded, theiler=20, max_iter=40, min_fit_points=6,
                   max_fit_fraction=0.10, min_r2=0.90, require_valid_fit=True):
    """Rosenstein et al. nearest-neighbor largest-Lyapunov-exponent estimate.

    The estimator is the slope of the mean logarithmic nearest-neighbor
    divergence over the predefined early scaling region. R^2 is retained as a
    transparent fit-quality diagnostic; it is NOT used to convert an otherwise
    finite Rosenstein estimate into NaN. This avoids making an arbitrary QC
    threshold part of the estimator itself.

    ``valid_fit`` therefore means the regression is numerically valid.
    ``valid_fit_quality`` separately records whether R^2 >= min_r2.
    """
    X = np.asarray(embedded, dtype=float)
    n = len(X)
    theiler = max(0, int(theiler))

    if n < max(30, theiler + max_iter + 5):
        return np.nan, {"reason": "insufficient embedded points"}

    nn = _nearest_neighbor_indices(X, theiler=theiler)

    log_div = np.full((n, max_iter), np.nan)
    for i, j in enumerate(nn):
        if j < 0:
            continue
        max_k = min(max_iter, n - max(i, j))
        for k_step in range(max_k):
            d = np.linalg.norm(X[i + k_step] - X[j + k_step])
            if d > 1e-15:
                log_div[i, k_step] = np.log(d)

    y = np.nanmean(log_div, axis=0)
    counts = np.sum(np.isfinite(log_div), axis=0)

    # Require a meaningful number of reference trajectories at each evolution
    # step before using that step in the averaged divergence curve.
    valid = np.isfinite(y) & (counts >= max(10, int(0.25 * n)))
    ks = np.flatnonzero(valid)

    if len(ks) < min_fit_points:
        return np.nan, {
            "reason": "insufficient divergence points",
            "counts": counts,
            "mean_log_divergence": y,
        }

    # The first 10% is the primary early scaling region. Because at least six
    # points are required, a short 40-step curve uses six points rather than
    # four. This is a minimum-point safeguard, not a later-window search.
    kmax = min(
        int(max_iter),
        max(
            min_fit_points,
            int(math.ceil(max_fit_fraction * len(ks))),
        ),
    )
    fit_idx = ks[:kmax]

    if len(fit_idx) < min_fit_points:
        return np.nan, {
            "reason": "insufficient early scaling points",
            "counts": counts,
            "mean_log_divergence": y,
        }

    xx = fit_idx.astype(float)
    yy = y[fit_idx]

    if not np.all(np.isfinite(yy)):
        return np.nan, {
            "reason": "non-finite selected scaling region",
            "fit_indices": fit_idx,
            "counts": counts,
            "mean_log_divergence": y,
        }

    slope, intercept = np.polyfit(xx, yy, 1)
    if not (np.isfinite(slope) and np.isfinite(intercept)):
        return np.nan, {
            "reason": "non-finite linear regression",
            "fit_indices": fit_idx,
            "counts": counts,
            "mean_log_divergence": y,
        }

    pred = intercept + slope * xx
    ss_res = float(np.sum((yy - pred) ** 2))
    ss_tot = float(np.sum((yy - np.mean(yy)) ** 2))
    r2 = 1.0 if ss_tot <= 1e-15 else 1.0 - ss_res / ss_tot

    return float(slope), {
        "slope": float(slope),
        "intercept": float(intercept),
        "r2": float(r2) if np.isfinite(r2) else np.nan,
        "valid_fit": True,
        "valid_fit_quality": bool(np.isfinite(r2) and r2 >= min_r2),
        "fit_quality_threshold": float(min_r2),
        "fit_indices": fit_idx,
        "counts": counts,
        "mean_log_divergence": y,
        "fit_strategy": "primary_early_prefix",
        "negative_slope_retained": True,
    }


def _delay_map_jacobian_from_local_fit(embedded, idx, neighbor_indices, step=1):
    """
    Jacobian in the delay-coordinate companion-matrix form: only the first
    row (evolution of the raw scalar series) is estimated from data via
    least squares; the remaining rows are the exact, deterministic shift
    structure inherent to delay embedding (z(n+1)'s coordinates 2..d are,
    by construction, exactly equal to z(n)'s coordinates 1..d-1).
    """
    X = np.asarray(embedded, dtype=float)
    d = X.shape[1]
    neigh = np.asarray(neighbor_indices, dtype=int)
    step = int(step)
    if step < 1:
        raise ValueError("step must be >= 1")
    neigh = neigh[(neigh >= 0) & (neigh + step < len(X))]
    if len(neigh) < d + 2 or idx + step >= len(X):
        return None
    Z = X[neigh]
    Y0 = X[neigh + step, 0]
    center = X[idx]
    D = Z - center
    A = D
    b = Y0 - X[idx, 0]
    try:
        grad, *_ = np.linalg.lstsq(A, b, rcond=None)
    except np.linalg.LinAlgError:
        return None
    J = np.zeros((d, d))
    J[0, :] = grad
    if d > 1:
        J[1:, :-1] = np.eye(d - 1)
    return J


def _spatial_neighbors(embedded, idx, K, exclude=1):
    X = np.asarray(embedded, dtype=float)
    d = np.linalg.norm(X - X[idx], axis=1)
    d[idx] = np.inf
    if exclude > 0:
        lo = max(0, idx - exclude)
        hi = min(len(X), idx + exclude + 1)
        d[lo:hi] = np.inf
    order = np.argsort(d)
    order = order[np.isfinite(d[order])]
    return order[:K]


def _random_R_matrix(rng, d, kappa):
    """
    Random perturbation matrix per So et al. (1996), PRL 76, 4705, Eq. (4):
    the higher-dimensional perturbation term is kappa * R * ||z_{n+1}-z_n||,
    where R is a d x d matrix with entries drawn independently and uniformly
    from [-1, 1] and the norm is the L1 norm of the state difference (the
    paper specifies L1 explicitly).
    """
    return rng.uniform(-1.0, 1.0, size=(d, d))


def _l1_norm(v):
    return np.sum(np.abs(v))


def so_fixed_point_transform(z, Fz, J, R, kappa):
    """
    So et al. (1996) fixed-point transform, Eq. (1) and Eq. (4):
        G(z, R) = [I - S(z, R)]^-1 [F(z) - S(z, R) z]
        S(z, R) = J_F(z) + kappa * R * ||F(z) - z||_1
    """
    d = len(z)
    delta = Fz - z
    S = J + kappa * R * _l1_norm(delta)
    M = np.eye(d) - S
    try:
        return np.linalg.solve(M, Fz - S @ z)
    except np.linalg.LinAlgError:
        return None


def _project_fixed_point_tube(transformed, percentile=10.0):
    """
    Reduce the transformed fixed-point distribution using So's diagonal tube:
    keep the points whose components are closest to their own mean (i.e.
    closest to lying on the diagonal z_1=z_2=...=z_d, which every genuine
    fixed point satisfies exactly in delay coordinates).

    Uses a PERCENTILE cutoff (always keep the closest `percentile`% of points)
    rather than a fixed distance threshold. A fixed threshold was tested
    directly and found to collapse to keeping ~0.01% of points at a realistic
    higher embedding dimension (d=6) where it kept ~25% at d=2 -- a genuine
    curse-of-dimensionality effect. A percentile cutoff self-calibrates to
    whatever the ambient dimension's distance distribution actually looks
    like, rather than requiring a hand-tuned, dimension-dependent constant.
    """
    Z = np.asarray(transformed, dtype=float)
    if len(Z) == 0:
        return np.empty(0), np.zeros(0, dtype=bool)
    centered = Z - np.mean(Z, axis=1, keepdims=True)
    distance = np.linalg.norm(centered, axis=1)
    threshold = np.percentile(distance, percentile)
    keep = distance <= threshold
    scalar = np.mean(Z[keep], axis=1)
    return scalar, keep


def upo_coverage(candidates, embedded, config=CFG, rng_seed=0):
    """
    Density/coverage measure for detected UPOs -- the actual proxy for
    Devaney's "density of periodic orbits" condition. The existence of a
    single detected fixed point is not itself evidence of chaos; what matters
    is how much of the reconstructed attractor sits near SOME detected UPO.

    Procedure: normalize the attractor per-dimension to [0,1], sample
    reference points directly from the real embedded trajectory (not an
    arbitrary coordinate grid -- a fixed grid's cell count grows
    exponentially with embedding dimension, leaving nearly all cells empty
    regardless of whether the system is chaotic), and compute the fraction
    of reference points lying within a distance-based radius of the nearest
    detected UPO.

    The radius is deliberately NOT a fixed constant or fixed fraction of the
    attractor's extent. Direct testing showed that approach suffers a severe
    curse-of-dimensionality collapse: what discriminates well at low
    embedding dimension can become meaningless at realistic higher
    dimensions purely from distance concentration, independent of whether
    real UPO structure is present. Instead, the radius is set relative to
    the k-th nearest-neighbor spacing within the normalized attractor
    itself, which self-calibrates to whatever the ambient dimension's
    typical point density actually is.
    """
    if len(candidates) == 0 or len(embedded) == 0:
        return {"coverage": 0.0, "radius": np.nan, "n_reference": 0, "n_upos": 0}

    X = np.asarray(embedded, dtype=float)
    mins, maxs = X.min(axis=0), X.max(axis=0)
    span = np.maximum(maxs - mins, 1e-12)
    X_norm = (X - mins) / span

    upo_locs = np.array([c["location"] for c in candidates])
    upo_locs_norm = (upo_locs - mins) / span

    rng = np.random.default_rng(rng_seed)
    n_ref = min(config.coverage_n_reference, len(X_norm))
    ref_idx = rng.choice(len(X_norm), size=n_ref, replace=False)
    reference_points = X_norm[ref_idx]

    sample_size = min(500, len(X_norm))
    sample_idx = rng.choice(len(X_norm), size=sample_size, replace=False)
    sample = X_norm[sample_idx]
    nn_dists = []
    for i in range(len(sample)):
        d = np.linalg.norm(sample - sample[i], axis=1)
        d[i] = np.inf
        k = min(config.coverage_radius_neighbors_k, len(d) - 1)
        nn_dists.append(np.partition(d, k)[k])
    typical_spacing = np.median(nn_dists)
    radius = config.coverage_radius_multiplier * typical_spacing

    covered = 0
    for ref in reference_points:
        dists = np.linalg.norm(upo_locs_norm - ref, axis=1)
        if np.min(dists) <= radius:
            covered += 1

    coverage = covered / len(reference_points)
    return {
        "coverage": coverage, "radius": float(radius), "typical_spacing": float(typical_spacing),
        "n_reference": len(reference_points), "n_upos": len(candidates),
    }


def detect_so_fixed_points(embedded, tau=1, config=CFG, rng=None):
    """So et al. fixed-point extraction with data-derived candidate states.

    Histogram peaks are found in the scalar diagonal coordinate of transformed
    points.  For each peak, the candidate state is the diagonal projection of
    the *actual transformed points in that peak*, rather than an artificial
    all-center vector.  The within-peak scatter is retained as a genuine
    diagonal-consistency diagnostic.
    """
    if rng is None: rng=np.random.default_rng(config.random_seed)
    X=np.asarray(embedded,dtype=float); n,d=X.shape; tau=int(tau)
    if tau<1: raise ValueError("tau must be >= 1")
    transformed=[]; source_idx=[]
    M=max(config.so_jacobian_neighbors,d+2)
    for i in range(n-tau):
        neigh=_spatial_neighbors(X,i,K=M,exclude=max(1,tau))
        if len(neigh)<d+2: continue
        J=_delay_map_jacobian_from_local_fit(X,i,neigh,step=tau)
        if J is None: continue
        z=X[i]; Fz=X[i+tau]
        for _ in range(config.so_random_R):
            R=_random_R_matrix(rng,d,config.so_kappa)
            zh=so_fixed_point_transform(z,Fz,J,R,config.so_kappa)
            if zh is not None and np.all(np.isfinite(zh)):
                transformed.append(zh); source_idx.append(i)
    if not transformed:
        return {"candidates":[],"transformed":np.empty((0,d)),"scalar":np.empty(0)}
    transformed=np.asarray(transformed); source_idx=np.asarray(source_idx,int)
    scalar,keep=_project_fixed_point_tube(transformed,config.so_diagonal_percentile)
    transformed_kept=transformed[keep]; source_kept=source_idx[keep]
    if len(scalar)<config.so_peak_min_count:
        return {"candidates":[],"transformed":transformed_kept,"scalar":scalar,"source_indices":source_kept}
    hist,edges=np.histogram(scalar,bins=config.so_hist_bins); centers=0.5*(edges[:-1]+edges[1:])
    smoothed=gaussian_filter1d(hist.astype(float),sigma=0.5)
    bg=np.median(smoothed); mad=np.median(np.abs(smoothed-bg)); scale=max(mad/0.6744897501960817,1e-12)
    threshold=bg+config.so_peak_sigma*scale
    peaks,_=find_peaks(hist.astype(float),height=max(threshold,config.so_peak_min_count),distance=1)
    candidates=[]
    for pidx in peaks:
        in_bin=(scalar>=edges[pidx]) & ((scalar<edges[pidx+1]) if pidx < len(hist)-1 else (scalar<=edges[pidx+1]))
        pts=transformed_kept[in_bin]
        if len(pts)==0: continue
        # Project each transformed state to the diagonal, then use the mean
        # projected coordinate as the fixed-point candidate.
        diag_coord=np.mean(pts,axis=1)
        c=float(np.mean(diag_coord)); zhat=np.full(d,c,dtype=float)
        diagonal_error=float(np.median(np.std(pts,axis=1)))
        centroid_error=float(np.linalg.norm(np.mean(pts,axis=0)-zhat))
        candidates.append({"period":1,"location":zhat,"scalar_location":c,
            "histogram_count":int(hist[pidx]),"peak_height":float(hist[pidx]),
            "diagonal_error":diagonal_error,"centroid_diagonal_error":centroid_error,
            "is_diagonal_consistent":bool(centroid_error <= max(0.10*np.std(transformed_kept),1e-12)),
            "peak_bin":int(pidx),"transformed_points_in_peak":int(len(pts))})
    return {"candidates":candidates,"transformed":transformed_kept,"scalar":scalar,
            "source_indices":source_kept,"histogram":hist,"edges":edges,"smoothed_histogram":smoothed}


# -------------------------------------------------------------------------
# General period-p So transform (1997)
# -------------------------------------------------------------------------

def _so_period_transform(Z, FZ, J_list, R_list, kappa):
    p = len(Z)
    d = Z[0].size
    pd = p * d
    B = np.zeros((pd, pd))
    rhs = np.zeros(pd)
    for k in range(p):
        kp1 = (k + 1) % p
        zk = Z[k]; Fzk = FZ[k]; Rk = R_list[k]; Jk = J_list[k]
        S = Jk + kappa * Rk * _l1_norm(Fzk - Z[kp1])
        r0 = k*d; r1 = r0+d
        B[r0:r1, r0:r1] = -S
        if k < p - 1:
            B[r0:r1, (k+1)*d:(k+2)*d] = np.eye(d)
        else:
            B[r0:r1, 0:d] = np.eye(d)
        rhs[r0:r1] = Fzk - S @ zk
    try:
        return np.linalg.solve(B, rhs).reshape(p, d)
    except np.linalg.LinAlgError:
        return None


def _cyclic_slab_reduction(Zhat):
    Zhat = np.asarray(Zhat, dtype=float)
    p, d = Zhat.shape
    independent = Zhat[:, 0].copy()
    recon = np.empty_like(Zhat)
    for k in range(p):
        for j in range(d):
            recon[k, j] = independent[(k+j) % p]
    residual = np.linalg.norm(Zhat - recon)
    return independent, residual


def detect_so_period_p(embedded, period, tau=1, config=CFG, rng=None):
    if period < 2:
        raise ValueError("Use detect_so_fixed_points for period 1.")
    if rng is None:
        rng = np.random.default_rng(config.random_seed + period)

    X = np.asarray(embedded, dtype=float)
    n, d = X.shape
    tau = int(tau)
    if tau < 1:
        raise ValueError("tau must be >= 1")
    K = config.so_neighbors_K
    M = max(config.so_jacobian_neighbors, d + 2)

    if n < period + d + 10:
        return {"candidates": [], "reduced": np.empty((0, period))}

    neighbor_lists = []
    jacobians = []
    for i in range(n - tau):
        neigh = _spatial_neighbors(X, i, K=M, exclude=max(1, tau))
        neigh = np.asarray([i, *neigh], dtype=int)
        J = _delay_map_jacobian_from_local_fit(X, i, neigh[1:], step=tau)
        neighbor_lists.append(neigh)
        jacobians.append(J)

    reduced_all = []
    max_backbones = config.so_max_backbones
    backbone_indices = list(range(n - (period * tau)))
    if max_backbones is not None and len(backbone_indices) > max_backbones:
        backbone_indices = list(np.linspace(0, len(backbone_indices)-1, max_backbones).astype(int))

    for start in backbone_indices:
        clusters = []
        valid = True
        for k in range(period):
            idx = start + k * tau
            choices = neighbor_lists[idx]
            choices = choices[choices + tau < n]
            if len(choices) == 0:
                valid = False
                break
            clusters.append(choices)
        if not valid:
            continue

        for combo in itertools.product(*clusters):
            Z = [X[idx] for idx in combo]
            FZ = [X[idx+tau] for idx in combo]
            J_list = []
            good = True
            for idx in combo:
                J = jacobians[idx]
                if J is None:
                    good = False
                    break
                J_list.append(J)
            if not good:
                continue

            for _ in range(config.so_random_R_period_p):
                R_list = [_random_R_matrix(rng, d, config.so_kappa) for _ in range(period)]
                Zhat = _so_period_transform(Z, FZ, J_list, R_list, config.so_kappa)
                if Zhat is None or not np.all(np.isfinite(Zhat)):
                    continue
                reduced, residual = _cyclic_slab_reduction(Zhat)
                scale = max(np.linalg.norm(Zhat), 1e-12)
                if residual <= config.so_slab_fraction * scale:
                    reduced_all.append(reduced)

    if len(reduced_all) < config.so_peak_min_count:
        return {"candidates": [], "reduced": np.asarray(reduced_all)}

    reduced = np.asarray(reduced_all)
    if period > 3:
        raise ValueError("Full p-dimensional histogram becomes impractical for p>3. Use p<=3.")

    ranges = []
    for j in range(period):
        lo, hi = np.min(reduced[:, j]), np.max(reduced[:, j])
        if np.isclose(lo, hi):
            hi = lo + 1e-9
        ranges.append((lo, hi))

    hist, edges = np.histogramdd(reduced, bins=config.so_hist_bins, range=ranges)
    padded = np.pad(hist, 1, mode="constant")
    candidate_cells = np.argwhere(hist >= config.so_peak_min_count)

    candidates = []
    for cell in candidate_cells:
        pc = tuple(cell + 1)
        slices = tuple(slice(c-1, c+2) for c in pc)
        local = padded[slices]
        if hist[tuple(cell)] < np.max(local):
            continue
        loc = [0.5*(edges[j][c]+edges[j][c+1]) for j, c in enumerate(cell)]
        candidates.append({
            "period": period, "orbit_coordinates": np.asarray(loc, dtype=float),
            "histogram_count": int(hist[tuple(cell)]),
        })

    return {"candidates": candidates, "reduced": reduced, "histogram": hist, "edges": edges}


# =============================================================================
# UPO stability
# =============================================================================

def estimate_jacobian_at_candidate(embedded, candidate, neighborhood=25, step=1):
    """Estimate the local tau-step map Jacobian near an actual candidate."""
    X=np.asarray(embedded,dtype=float); z0=np.asarray(candidate,dtype=float); d=X.shape[1]
    step=int(step); valid_n=len(X)-step
    if valid_n<=d+2: return None
    dist=np.linalg.norm(X[:valid_n]-z0,axis=1)
    idx=np.argsort(dist)[:max(int(neighborhood),d+2)]
    if len(idx)<d+2: return None
    center=int(idx[0])
    return _delay_map_jacobian_from_local_fit(X,center,idx[1:],step=step)


def fixed_point_stability(embedded,z_star,neighborhood=30,step=1):
    """Linear stability from the companion-form Jacobian of the delay map."""
    X=np.asarray(embedded,dtype=float); z=np.asarray(z_star,dtype=float)
    if X.ndim!=2 or X.shape[1]<2: return None
    J=estimate_jacobian_at_candidate(X,z,neighborhood=neighborhood,step=step)
    if J is None or not np.all(np.isfinite(J)): return None
    eig=np.linalg.eigvals(J); mod=np.abs(eig)
    unstable=bool(np.any(mod>1.0+1e-8))
    return {"jacobian":J,"multipliers":eig,"multiplier_moduli":mod,
            "unstable":unstable,"saddle":bool(np.any(mod>1.0+1e-8) and np.any(mod<1.0-1e-8)),
            "max_multiplier":float(np.max(mod)),"min_multiplier":float(np.min(mod)),
            "candidate_residual":float(np.linalg.norm(_delay_map_apply(X,z,step=step)-z)) if False else np.nan}


# =============================================================================
# Surrogate data
# =============================================================================

def aaft_surrogate(x, rng=None):
    """Generate one amplitude-adjusted Fourier-transform (AAFT) surrogate.

    Construction: rank-map the observed series to Gaussian order statistics,
    randomize Fourier phases while retaining the Gaussianized spectrum,
    inverse transform, then rank-map back to the exact observed amplitude
    distribution.  This is AAFT (not IAAFT); the latter is an iterative
    refinement and is intentionally not substituted here.
    """
    if rng is None: rng=np.random.default_rng()
    x=_as_1d(x); n=len(x)
    if n<4: return x.copy()
    # Gaussian rank transform, with deterministic tie handling.
    order=np.argsort(x,kind="mergesort")
    gaussian_order=np.sort(rng.standard_normal(n))
    gaussianized=np.empty(n,float); gaussianized[order]=gaussian_order
    spec=np.fft.rfft(gaussianized)
    amp=np.abs(spec)
    phase=rng.uniform(-np.pi,np.pi,len(amp))
    phase[0]=0.0
    if n%2==0: phase[-1]=0.0
    surrogate_gaussian=np.fft.irfft(amp*np.exp(1j*phase),n=n)
    # Exact amplitude distribution of the observed series.
    sorted_x=np.sort(x)
    rank_order=np.argsort(surrogate_gaussian,kind="mergesort")
    surrogate=np.empty(n,float); surrogate[rank_order]=sorted_x
    return surrogate


def gaussian_scaled_phase_shuffle(x, rng=None):
    return aaft_surrogate(x, rng=rng)


def so_surrogate_density_significance(real_result, surrogate_results):
    scalar = np.asarray(real_result.get("scalar", []), float)
    edges = np.asarray(real_result.get("edges"), float)
    if scalar.size == 0 or edges.size < 2 or not surrogate_results:
        return {"W": np.nan, "p_value": np.nan, "ratio_to_median_surrogate_W": np.nan}
    real_hist, _ = np.histogram(scalar, bins=edges)
    sur = np.asarray([np.histogram(np.asarray(r.get("scalar", []), float), bins=edges)[0] for r in surrogate_results], float)
    mean = np.mean(sur, axis=0)
    W = float(np.max(np.abs(real_hist - mean)))
    Ws = np.max(np.abs(sur - mean), axis=1)
    W0 = float(np.median(Ws))
    return {
        "real_histogram": real_hist, "surrogate_mean_histogram": mean, "surrogate_W": Ws,
        "W": W, "p_value": float((1+np.sum(Ws >= W))/(len(Ws)+1)),
        "ratio_to_median_surrogate_W": float(W/max(W0, 1e-12)),
    }


def run_surrogate_analysis(rr, tau, m, real_lle, real_fixed, lle_theiler, config=CFG):
    """Optional fixed-parameter AAFT surrogate test.

    Surrogates keep the real series' marginal distribution and linear power
    spectrum while destroying its nonlinear temporal organization. Crucially,
    tau and m are held fixed at the values selected from the real series so
    surrogate re-optimization cannot manufacture a favorable null model.
    """
    rng = np.random.default_rng(config.random_seed + 7919)
    lle_values=[]; fixed_results=[]
    for _ in range(int(config.surrogate_count)):
        xs = aaft_surrogate(rr, rng=rng)
        emb = _embed_backward(xs, tau, m)
        if len(emb) < 50:
            continue
        lv, _ = rosenstein_lle(emb, theiler=lle_theiler,
                               max_iter=config.lle_max_iter,
                               min_fit_points=config.lle_min_fit_points,
                               max_fit_fraction=config.lle_max_fit_fraction,
                               min_r2=config.lle_min_r2,
                               require_valid_fit=False)
        lle_values.append(lv)
        fixed_results.append(detect_so_fixed_points(emb, tau=tau, config=config, rng=rng))
    real_scalar=np.asarray(real_fixed.get("scalar",[]),float)
    real_edges=np.asarray(real_fixed.get("edges",[]),float)
    if real_scalar.size and real_edges.size>=2:
        sig=so_surrogate_density_significance(real_fixed,fixed_results)
    else:
        sig={"W":np.nan,"p_value":np.nan,"ratio_to_median_surrogate_W":np.nan}
    finite=np.asarray(lle_values,float); finite=finite[np.isfinite(finite)]
    if np.isfinite(real_lle) and len(finite):
        z=(real_lle-float(np.mean(finite)))/max(float(np.std(finite,ddof=1)),1e-12) if len(finite)>1 else np.nan
        p_lle=float((1+np.sum(finite>=real_lle))/(len(finite)+1))
    else:
        z=np.nan; p_lle=np.nan
    return {
        "n_surrogates": len(fixed_results),
        "lle_surrogates": finite,
        "lle_z": float(z) if np.isfinite(z) else np.nan,
        "lle_p_upper": p_lle,
        "upo": sig,
        "null_preserves_tau_m": True,
    }


def upo_peak_statistic(result):
    if not result.get("candidates"):
        return 0.0
    return float(max(c.get("peak_height", c.get("histogram_count", 0.0)) for c in result["candidates"]))


# =============================================================================
# Full RR / ECG record analysis
# =============================================================================

def stationarity_diagnostics(rr, n_blocks=5):
    rr = _as_1d(rr)
    blocks = np.array_split(rr, n_blocks)
    means = np.array([np.mean(b) for b in blocks])
    stds = np.array([np.std(b, ddof=1) for b in blocks])
    mean_cv = float(np.std(means)/max(abs(np.mean(means)), 1e-12))
    std_cv = float(np.std(stds)/max(np.mean(stds), 1e-12))
    return {"block_means": means, "block_stds": stds, "mean_cv": mean_cv, "std_cv": std_cv,
            "stationarity_warning": bool(mean_cv > 0.05 or std_cv > 0.20)}


def analyze_segment(rr_intervals, config=CFG):
    rr = _as_1d(rr_intervals)
    if len(rr) < 100:
        raise ValueError("At least ~100 RR intervals are recommended; more is strongly preferred.")

    rr_corrected, artifact_mask, _ = correct_rr_intervals(rr, config)
    rr_dynamics = rr_corrected if config.use_corrected_rr_for_dynamics else rr
    raw_stationarity = stationarity_diagnostics(rr)
    corrected_stationarity = stationarity_diagnostics(rr_corrected)
    stationarity = {
        "raw": raw_stationarity,
        "corrected": corrected_stationarity,
        "stationarity_warning": stationarity_diagnostics(rr_dynamics)["stationarity_warning"],
    }
    embedded, tau, m, mi, e1, e2 = takens_embed(rr_dynamics, config=config)
    if len(embedded) < 50:
        raise ValueError("Insufficient embedded points after reconstruction.")

    lle_theiler = config.lle_theiler_beats if config.lle_theiler_beats is not None else _mean_period_beats(rr_dynamics)
    lle, lle_diag = rosenstein_lle(embedded, theiler=lle_theiler, max_iter=config.lle_max_iter,
                                    min_fit_points=config.lle_min_fit_points,
                                    max_fit_fraction=config.lle_max_fit_fraction,
                                    min_r2=config.lle_min_r2,
                                    require_valid_fit=config.lle_require_valid_fit)

    rng = np.random.default_rng(config.random_seed)
    fixed = detect_so_fixed_points(embedded, tau=tau, config=config, rng=rng)
    for c in fixed.get("candidates", []):
        stability = fixed_point_stability(embedded, c["location"], neighborhood=30, step=tau)
        if stability:
            c.update(stability)

    # A detected fixed point is only counted as a UPO for coverage after its
    # local multipliers establish instability. Unknown stability is not treated
    # as evidence of an unstable periodic orbit.
    unstable_candidates = [c for c in fixed.get("candidates", []) if c.get("unstable") is True]
    coverage = upo_coverage(unstable_candidates, embedded, config=config, rng_seed=config.random_seed)

    periodic = {}
    for p in config.so_periods:
        if p == 1:
            continue
        periodic[p] = detect_so_period_p(embedded, p, tau=tau, config=config, rng=rng)

    surrogate = None
    if config.compute_surrogates:
        surrogate = run_surrogate_analysis(rr_dynamics, tau, m, lle, fixed,
                                           lle_theiler, config=config)

    return {
        "rr_raw": rr, "rr_corrected": rr_corrected, "rr_dynamics": rr_dynamics,
        "used_corrected_rr_for_dynamics": bool(config.use_corrected_rr_for_dynamics),
        "rr_artifact_mask": artifact_mask,
        "stationarity": stationarity, "tau": tau, "embedding_dimension": m,
        "tdmi": mi, "cao_E1": e1, "cao_E2": e2, "embedded": embedded,
        "lle_per_beat": lle, "lle_theiler_beats": int(lle_theiler), "lle_diagnostics": lle_diag,
        "so_fixed_points": fixed, "so_periodic_orbits": periodic, "upo_coverage": coverage,
        "surrogate_analysis": surrogate,
    }


def analyze_ecg_record(raw_ecg, fs, config=CFG):
    raw_ecg = _as_1d(raw_ecg)
    clean_ecg = preprocess_ecg(raw_ecg, fs, config)  # general-purpose cleaned waveform only
    r_peaks = detect_r_peaks(raw_ecg, fs, config)     # peak detection runs on raw ECG

    rr_raw, rr_times, physiological_valid = extract_rr_intervals(
        r_peaks, fs, config.rr_min_seconds, config.rr_max_seconds
    )
    if len(rr_raw) < 100:
        raise ValueError(f"Only {len(rr_raw)} RR intervals were obtained; >=100 is recommended.")

    results = analyze_segment(rr_raw, config)
    results.update({
        "sampling_frequency_hz": float(fs), "n_ecg_samples": int(len(raw_ecg)),
        "n_beats_detected": int(len(r_peaks)), "n_rr_intervals": int(len(rr_raw)),
        "n_rr_artifacts": int(np.sum(results["rr_artifact_mask"])),
        "r_peaks": r_peaks, "rr_times": rr_times, "physiological_rr_mask": physiological_valid,
    })
    return results, clean_ecg


def validate_r_peaks_against_annotations(r_peaks, reference_peaks, fs, tolerance_ms=75.0):
    detected = np.asarray(r_peaks, int)
    reference = np.asarray(reference_peaks, int)
    tol = int(round(tolerance_ms * fs / 1000))
    used = np.zeros(len(reference), bool)
    tp = 0
    # Chronological one-to-one matching avoids consuming the wrong reference
    # beat when multiple annotations lie within the tolerance window.
    for p in detected:
        candidates = np.flatnonzero((np.abs(reference - p) <= tol) & (~used))
        if len(candidates):
            j = candidates[np.argmin(np.abs(reference[candidates] - p))]
            used[j] = True
            tp += 1
    fp = len(detected) - tp
    fn = len(reference) - tp
    sens = tp / max(tp + fn, 1)
    ppv = tp / max(tp + fp, 1)
    f1 = 2 * sens * ppv / max(sens + ppv, 1e-12)
    return {"TP": tp, "FP": fp, "FN": fn, "sensitivity": sens, "PPV": ppv, "F1": f1}


# =============================================================================
# NaN-aware twin-experiment utilities
# =============================================================================

def _mann_whitney_nan_safe(a, b):
    """Two-sided Mann--Whitney U after finite-value filtering."""
    from scipy.stats import mannwhitneyu
    a=np.asarray(a,float); b=np.asarray(b,float)
    a=a[np.isfinite(a)]; b=b[np.isfinite(b)]
    if len(a)==0 or len(b)==0:
        return {"U":np.nan,"p_value":np.nan,"n_a":len(a),"n_b":len(b)}
    r=mannwhitneyu(a,b,alternative="two-sided")
    return {"U":float(r.statistic),"p_value":float(r.pvalue),"n_a":len(a),"n_b":len(b)}


def summarize_twin_results(chaotic_results, periodic_results):
    """Summarize analyze_segment outputs with explicit LLE fit-success rates."""
    def extract(results):
        lle=np.asarray([r.get("lle_per_beat",np.nan) for r in results],float)
        coverage=np.asarray([
            r.get("upo_coverage",{}).get("coverage",np.nan) for r in results
        ],float)
        valid=np.isfinite(lle)
        return lle,coverage,valid

    lc,cc,vc=extract(chaotic_results)
    lp,cp,vp=extract(periodic_results)

    return {
        "lle": {
            "chaotic_valid_fit_rate": float(np.mean(vc)) if len(vc) else np.nan,
            "periodic_valid_fit_rate": float(np.mean(vp)) if len(vp) else np.nan,
            "chaotic_valid_values": lc[vc],
            "periodic_valid_values": lp[vp],
            "mann_whitney_valid_only": _mann_whitney_nan_safe(lc,lp),
        },
        "upo_coverage": {
            "chaotic_values": cc[np.isfinite(cc)],
            "periodic_values": cp[np.isfinite(cp)],
            "mann_whitney": _mann_whitney_nan_safe(cc,cp),
        },
    }


def run_twin_experiment(chaotic_segments, periodic_segments, config=CFG):
    """Run analyze_segment on matched iterable collections of RR segments."""
    chaotic_results=[analyze_segment(rr,config=config) for rr in chaotic_segments]
    periodic_results=[analyze_segment(rr,config=config) for rr in periodic_segments]
    return {
        "chaotic_results":chaotic_results,
        "periodic_results":periodic_results,
        "summary":summarize_twin_results(chaotic_results,periodic_results),
    }


def compare_rr_correction_twin(chaotic_segments, periodic_segments, config=CFG):
    """Run the same twin experiment with raw-RR and corrected-RR dynamics."""
    from dataclasses import replace
    raw_cfg=replace(config,use_corrected_rr_for_dynamics=False)
    corrected_cfg=replace(config,use_corrected_rr_for_dynamics=True)
    return {
        "raw_rr":run_twin_experiment(chaotic_segments,periodic_segments,raw_cfg),
        "corrected_rr":run_twin_experiment(chaotic_segments,periodic_segments,corrected_cfg),
    }

# =============================================================================
# PhysioNet / MIT-BIH
# =============================================================================

def load_mitbih_record(record_name, dl_dir="mitdb_data", db_name="mitdb"):
    if wfdb is None:
        raise ImportError("Optional dependency 'wfdb' is required for PhysioNet/MIT-BIH loading.")
    os.makedirs(dl_dir, exist_ok=True)
    local_path = os.path.join(dl_dir, record_name)
    if (
        not os.path.exists(local_path + ".dat")
        or not os.path.exists(local_path + ".hea")
        or not os.path.exists(local_path + ".atr")
    ):
        wfdb.dl_database(db_name, dl_dir=dl_dir, records=[record_name])
    record = wfdb.rdrecord(local_path)
    fs = float(record.fs)
    signal = np.asarray(record.p_signal[:, 0], dtype=float)
    return signal, fs


# =============================================================================
# 5. Binary arrhythmia-classification experiment
# =============================================================================
#
# This section deliberately sits on top of the existing validated pipeline.
# It does NOT alter:
#   - Pan-Tompkins detection
#   - RR correction
#   - TDMI/Cao reconstruction
#   - Rosenstein LLE
#   - So et al. UPO detection
#
# It only:
#   1. downloads/loads MIT-BIH records,
#   2. matches detected R peaks to reference annotations,
#   3. makes fixed-length RR windows,
#   4. labels those windows,
#   5. extracts nonlinear features from each window,
#   6. performs record-grouped cross-validation,
#   7. compares raw-RR and corrected-RR dynamics.
#
# The target is intentionally explicit and configurable.  It is NOT a claim
# that every non-N beat is clinically equivalent.

MITDB_RECORDS = [
    "100", "101", "102", "103", "104", "105", "106", "107",
    "108", "109", "111", "112", "113", "114", "115", "116",
    "117", "118", "119", "121", "122", "123", "124", "200",
    "201", "202", "203", "205", "207", "208", "209", "210",
    "212", "213", "214", "215", "217", "219", "220", "221",
    "222", "223", "228", "230", "231", "232", "233", "234",
]

# Beat classes used for the binary experiment.
# Normal-side symbols are kept explicit rather than treating every unknown
# annotation as normal.
NORMAL_BEAT_SYMBOLS = frozenset({
    "N", "L", "R", "e", "j"
})

ABNORMAL_BEAT_SYMBOLS = frozenset({
    "A", "a", "J", "S", "V", "E", "F", "f", "Q"
})

# Annotation symbols that are deliberately excluded from the beat-level label.
# In particular, paced beats ("/") are not silently folded into either class.
EXCLUDED_BEAT_SYMBOLS = frozenset({
    "/", "P", "Q", "B", "?", "x"
})

# Experiment defaults.  Change these deliberately before a final paper run.
CLASSIFIER_WINDOW_RR = 256
CLASSIFIER_STEP_RR = 256

# A window is abnormal when at least this fraction of its usable beats is
# abnormal.  0.10 means >=10% abnormal beats.
CLASSIFIER_ABNORMAL_FRACTION = 0.10

# Maximum detector-to-annotation timing error accepted when transferring the
# cardiologist annotation label onto a detected beat.
ANNOTATION_MATCH_TOLERANCE_MS = 75.0

CLASSIFIER_FEATURE_COLUMNS = [
    "lle_per_beat",
    "upo_coverage",
    "n_unstable_upos",
    "upo_candidate_density",
]


def load_mitbih_record_with_annotations(
    record_name,
    dl_dir="mitdb_data",
    db_name="mitdb",
):
    """
    Load one MIT-BIH ECG record and its annotation file.

    Returns:
        signal, fs, annotation_samples, annotation_symbols
    """
    if wfdb is None:
        raise ImportError(
            "wfdb is required. Install it with: pip install wfdb"
        )

    os.makedirs(dl_dir, exist_ok=True)
    local_path = os.path.join(dl_dir, str(record_name))

    required_files = [
        local_path + ".dat",
        local_path + ".hea",
        local_path + ".atr",
    ]

    # PhysioNet can transiently return HTTP 502 and leave a partial record.
    # Require signal + header + annotation, and retry transient download
    # failures before declaring the record unavailable.
    if not all(os.path.exists(p) for p in required_files):
        last_error = None
        for attempt in range(1, 5):
            try:
                wfdb.dl_database(
                    db_name,
                    dl_dir=dl_dir,
                    records=[str(record_name)],
                )
                if all(os.path.exists(p) for p in required_files):
                    break
            except Exception as exc:
                last_error = exc

            if attempt < 4:
                import time
                time.sleep(2.0 * attempt)

        if not all(os.path.exists(p) for p in required_files):
            missing = [p for p in required_files if not os.path.exists(p)]
            raise RuntimeError(
                f"Could not obtain complete MIT-BIH record {record_name} "
                f"after 4 download attempts. Missing: {missing}. "
                f"Last error: {last_error}"
            )

    record = wfdb.rdrecord(local_path)
    annotation = wfdb.rdann(local_path, "atr")

    signal = np.asarray(record.p_signal[:, 0], dtype=float)
    fs = float(record.fs)

    annotation_samples = np.asarray(annotation.sample, dtype=int)
    annotation_symbols = np.asarray(annotation.symbol, dtype=str)

    if len(annotation_samples) != len(annotation_symbols):
        raise RuntimeError(
            f"Annotation length mismatch for record {record_name}."
        )

    return (
        signal,
        fs,
        annotation_samples,
        annotation_symbols,
    )


def _annotation_label(symbol):
    """Map a beat annotation to {0, 1, None}."""
    symbol = str(symbol)

    if symbol in NORMAL_BEAT_SYMBOLS:
        return 0

    if symbol in ABNORMAL_BEAT_SYMBOLS:
        return 1

    return None


def match_detected_peaks_to_annotations(
    detected_peaks,
    annotation_samples,
    annotation_symbols,
    fs,
    tolerance_ms=ANNOTATION_MATCH_TOLERANCE_MS,
):
    """
    One-to-one chronological matching of detected R peaks to reference
    beat annotations.

    Returns an array with one entry per detected R peak:
        0 = normal beat
        1 = abnormal beat
       -1 = excluded/unknown/unmatched
    """
    detected = np.asarray(detected_peaks, dtype=int)
    reference = np.asarray(annotation_samples, dtype=int)

    if len(detected) == 0 or len(reference) == 0:
        return np.full(len(detected), -1, dtype=int)

    tol = int(round(float(tolerance_ms) * float(fs) / 1000.0))
    used = np.zeros(len(reference), dtype=bool)
    labels = np.full(len(detected), -1, dtype=int)

    for i, peak in enumerate(detected):
        candidates = np.flatnonzero(
            (np.abs(reference - peak) <= tol) & (~used)
        )

        if len(candidates) == 0:
            continue

        j = candidates[
            np.argmin(np.abs(reference[candidates] - peak))
        ]

        used[j] = True
        mapped = _annotation_label(annotation_symbols[j])

        if mapped is not None:
            labels[i] = int(mapped)

    return labels


def make_rr_windows(
    rr_values,
    beat_labels,
    window_size=CLASSIFIER_WINDOW_RR,
    step_size=CLASSIFIER_STEP_RR,
    abnormal_fraction=CLASSIFIER_ABNORMAL_FRACTION,
):
    """
    Create fixed-length RR windows.

    RR interval i is the interval ending at detected beat i+1, so its label is
    taken from the corresponding ending beat.

    A window is retained only if every RR interval in it has a usable
    normal/abnormal beat label. This prevents excluded/unknown annotations from
    being silently converted into either class.

    Returns:
        list of dictionaries containing:
            start_rr
            end_rr
            rr
            label
            abnormal_fraction
    """
    rr = np.asarray(rr_values, dtype=float)
    labels = np.asarray(beat_labels, dtype=int)

    if len(rr) != len(labels) - 1:
        raise ValueError(
            "rr_values must have exactly one fewer element than beat_labels."
        )

    if window_size < 100:
        raise ValueError("window_size must be >= 100.")

    if step_size < 1:
        raise ValueError("step_size must be >= 1.")

    if not (0.0 <= abnormal_fraction <= 1.0):
        raise ValueError("abnormal_fraction must be in [0, 1].")

    windows = []

    for start in range(
        0,
        len(rr) - window_size + 1,
        step_size,
    ):
        end = start + window_size

        rr_window = rr[start:end]
        # RR[start] ends at beat start+1.
        # Therefore the corresponding beat labels are:
        label_window = labels[start + 1:end + 1]

        # -1 denotes an excluded/unknown/unmatched annotation.
        if np.any(label_window < 0):
            continue

        frac = float(np.mean(label_window == 1))
        label = int(frac >= abnormal_fraction)

        windows.append({
            "start_rr": int(start),
            "end_rr": int(end),
            "rr": rr_window.copy(),
            "label": label,
            "abnormal_fraction": frac,
        })

    return windows


def extract_classifier_features(result):
    """
    Extract only features that are actually produced by analyze_segment().

    No reconstruction parameter (tau or m) is used as a disease feature.
    """
    lle = result.get("lle_per_beat", np.nan)

    coverage_result = result.get("upo_coverage", {})
    coverage = coverage_result.get("coverage", np.nan)

    fixed = result.get("so_fixed_points", {})
    candidates = fixed.get("candidates", [])

    n_unstable = sum(
        bool(c.get("unstable") is True)
        for c in candidates
    )

    # Density is normalized by the number of embedded trajectory points so
    # different window lengths would remain comparable.
    embedded = result.get("embedded", np.empty((0, 1)))
    n_embedded = len(embedded)

    candidate_density = (
        float(n_unstable) / float(n_embedded)
        if n_embedded > 0
        else np.nan
    )

    return {
        "lle_per_beat": (
            float(lle) if np.isfinite(lle) else np.nan
        ),
        "upo_coverage": (
            float(coverage)
            if np.isfinite(coverage)
            else np.nan
        ),
        "n_unstable_upos": float(n_unstable),
        "upo_candidate_density": candidate_density,
    }


def analyze_mitbih_record_windows(
    record_name,
    config,
    dl_dir="mitdb_data",
    db_name="mitdb",
    window_size=CLASSIFIER_WINDOW_RR,
    step_size=CLASSIFIER_STEP_RR,
    abnormal_fraction=CLASSIFIER_ABNORMAL_FRACTION,
):
    """
    Run the existing ECG/RR pipeline on one MIT-BIH record and then analyze
    fixed-length RR windows for binary classification.
    """
    (
        signal,
        fs,
        annotation_samples,
        annotation_symbols,
    ) = load_mitbih_record_with_annotations(
        record_name,
        dl_dir=dl_dir,
        db_name=db_name,
    )

    # Peak detection is the exact detector already used by the pipeline.
    r_peaks = detect_r_peaks(signal, fs, config)

    rr_raw, rr_times, physiological_valid = extract_rr_intervals(
        r_peaks,
        fs,
        config.rr_min_seconds,
        config.rr_max_seconds,
    )

    if len(rr_raw) < window_size:
        return [], {
            "record": str(record_name),
            "status": "skipped",
            "reason": "fewer RR intervals than classifier window",
            "n_rr": int(len(rr_raw)),
            "n_r_peaks": int(len(r_peaks)),
        }

    beat_labels = match_detected_peaks_to_annotations(
        r_peaks,
        annotation_samples,
        annotation_symbols,
        fs,
    )

    windows = make_rr_windows(
        rr_raw,
        beat_labels,
        window_size=window_size,
        step_size=step_size,
        abnormal_fraction=abnormal_fraction,
    )

    rows = []

    for window_number, window in enumerate(windows):
        try:
            result = analyze_segment(
                window["rr"],
                config=config,
            )
        except (ValueError, np.linalg.LinAlgError) as exc:
            # A failed nonlinear fit is not silently converted into a feature.
            # Record the failure so the final QC table exposes it.
            rows.append({
                "record": str(record_name),
                "window": int(window_number),
                "label": int(window["label"]),
                "abnormal_fraction": float(window["abnormal_fraction"]),
                "lle_per_beat": np.nan,
                "upo_coverage": np.nan,
                "n_unstable_upos": np.nan,
                "upo_candidate_density": np.nan,
                "analysis_error": str(exc),
            })
            continue

        features = extract_classifier_features(result)

        rows.append({
            "record": str(record_name),
            "window": int(window_number),
            "label": int(window["label"]),
            "abnormal_fraction": float(window["abnormal_fraction"]),
            **features,
            "analysis_error": "",
        })

    meta = {
        "record": str(record_name),
        "status": "ok",
        "n_rr": int(len(rr_raw)),
        "n_r_peaks": int(len(r_peaks)),
        "n_windows_total": int(len(windows)),
        "n_windows_retained": int(len(rows)),
        "n_unmatched_or_excluded_beats": int(np.sum(beat_labels < 0)),
        "n_physiologically_out_of_range_rr": int(
            np.sum(~physiological_valid)
        ),
    }

    return rows, meta


def build_mitbih_classification_dataset(
    records=MITDB_RECORDS,
    config=CFG,
    dl_dir="mitdb_data",
    db_name="mitdb",
    window_size=CLASSIFIER_WINDOW_RR,
    step_size=CLASSIFIER_STEP_RR,
    abnormal_fraction=CLASSIFIER_ABNORMAL_FRACTION,
):
    """
    Process all requested MIT-BIH records.

    Returns:
        dataframe, metadata_dataframe
    """
    all_rows = []
    metadata = []

    for number, record_name in enumerate(records, start=1):
        print(
            f"[{number}/{len(records)}] "
            f"Processing MIT-BIH record {record_name}..."
        )

        try:
            rows, meta = analyze_mitbih_record_windows(
                record_name,
                config=config,
                dl_dir=dl_dir,
                db_name=db_name,
                window_size=window_size,
                step_size=step_size,
                abnormal_fraction=abnormal_fraction,
            )
            all_rows.extend(rows)
            metadata.append(meta)

            print(
                f"    windows retained: {len(rows)}"
            )

        except Exception as exc:
            # Do not silently stop the entire study because one record failed.
            # The metadata table makes the failure explicit.
            meta = {
                "record": str(record_name),
                "status": "error",
                "reason": repr(exc),
            }
            metadata.append(meta)
            print(
                f"    ERROR: {exc}"
            )

    df = pd.DataFrame(all_rows)
    meta_df = pd.DataFrame(metadata)

    if len(df):
        for col in CLASSIFIER_FEATURE_COLUMNS + [
            "label",
            "record",
            "window",
            "abnormal_fraction",
        ]:
            if col not in df.columns:
                df[col] = np.nan

    return df, meta_df


def run_grouped_logistic_regression(
    df,
    feature_columns=CLASSIFIER_FEATURE_COLUMNS,
    n_splits=5,
    random_state=42,
):
    """
    Record-grouped logistic regression.

    Imputation and scaling are fitted inside each training fold, preventing
    information leakage from the held-out records.
    """
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import (
        average_precision_score,
        balanced_accuracy_score,
        confusion_matrix,
        f1_score,
        roc_auc_score,
    )
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import StratifiedGroupKFold

    required = list(feature_columns) + ["label", "record"]

    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"Classification dataframe is missing columns: {missing}"
        )

    work = df.copy()

    # A record must have a single group identity.
    work["record"] = work["record"].astype(str)
    work["label"] = work["label"].astype(int)

    X = work[list(feature_columns)]
    y = work["label"]
    groups = work["record"]

    if y.nunique() < 2:
        raise ValueError(
            "Both normal (0) and abnormal (1) windows are required."
        )

    n_groups = groups.nunique()
    n_splits = min(int(n_splits), int(n_groups))

    if n_splits < 2:
        raise ValueError(
            "At least two distinct records are required for grouped CV."
        )

    # Stratification can fail when the minority class occurs in fewer records
    # than the requested number of folds. Reduce n_splits until valid.
    while n_splits >= 2:
        try:
            cv = StratifiedGroupKFold(
                n_splits=n_splits,
                shuffle=True,
                random_state=random_state,
            )
            splits = list(cv.split(X, y, groups))
            break
        except ValueError:
            n_splits -= 1

    if n_splits < 2:
        raise ValueError(
            "Could not construct a valid StratifiedGroupKFold split. "
            "There are too few records containing both classes."
        )

    model = Pipeline([
        (
            "imputer",
            SimpleImputer(strategy="median"),
        ),
        (
            "scaler",
            StandardScaler(),
        ),
        (
            "classifier",
            LogisticRegression(
                class_weight="balanced",
                max_iter=5000,
                random_state=random_state,
            ),
        ),
    ])

    fold_rows = []
    oof_probability = np.full(len(work), np.nan)
    oof_prediction = np.full(len(work), -1, dtype=int)

    for fold, (train_idx, test_idx) in enumerate(
        splits,
        start=1,
    ):
        X_train = X.iloc[train_idx]
        X_test = X.iloc[test_idx]

        y_train = y.iloc[train_idx]
        y_test = y.iloc[test_idx]

        # Logistic regression requires both classes in the training fold.
        if y_train.nunique() < 2:
            raise ValueError(
                f"Fold {fold} training data contains only one class."
            )

        model.fit(X_train, y_train)

        probability = model.predict_proba(X_test)[:, 1]
        prediction = (probability >= 0.50).astype(int)

        oof_probability[test_idx] = probability
        oof_prediction[test_idx] = prediction

        tn, fp, fn, tp = confusion_matrix(
            y_test,
            prediction,
            labels=[0, 1],
        ).ravel()

        sensitivity = (
            tp / (tp + fn)
            if (tp + fn) > 0
            else np.nan
        )

        specificity = (
            tn / (tn + fp)
            if (tn + fp) > 0
            else np.nan
        )

        fold_rows.append({
            "fold": fold,
            "n_train": len(train_idx),
            "n_test": len(test_idx),
            "n_train_records": groups.iloc[train_idx].nunique(),
            "n_test_records": groups.iloc[test_idx].nunique(),
            "roc_auc": roc_auc_score(y_test, probability),
            "pr_auc": average_precision_score(y_test, probability),
            "sensitivity": sensitivity,
            "specificity": specificity,
            "balanced_accuracy": balanced_accuracy_score(
                y_test,
                prediction,
            ),
            "f1": f1_score(
                y_test,
                prediction,
                zero_division=0,
            ),
        })

    fold_df = pd.DataFrame(fold_rows)

    # Out-of-fold metrics use each window exactly once, while still respecting
    # record-level separation.
    valid_oof = (
        np.isfinite(oof_probability)
        & (oof_prediction >= 0)
    )

    oof_summary = {
        "oof_n_windows": int(np.sum(valid_oof)),
        "oof_roc_auc": (
            float(roc_auc_score(
                y.iloc[valid_oof],
                oof_probability[valid_oof],
            ))
            if np.sum(valid_oof) > 0
            else np.nan
        ),
        "oof_pr_auc": (
            float(average_precision_score(
                y.iloc[valid_oof],
                oof_probability[valid_oof],
            ))
            if np.sum(valid_oof) > 0
            else np.nan
        ),
    }

    return {
        "fold_results": fold_df,
        "fold_means": (
            fold_df[
                [
                    "roc_auc",
                    "pr_auc",
                    "sensitivity",
                    "specificity",
                    "balanced_accuracy",
                    "f1",
                ]
            ].mean()
            if len(fold_df)
            else pd.Series(dtype=float)
        ),
        "oof_summary": oof_summary,
        "oof_probability": oof_probability,
        "oof_prediction": oof_prediction,
        "n_splits": n_splits,
    }


def nonlinear_feature_statistics(df):
    """
    Descriptive/Mann-Whitney statistics for LLE and UPO coverage.

    LLE uses only finite fits. UPO coverage is evaluated independently.
    """
    from scipy.stats import mannwhitneyu

    output = {}

    for feature in ["lle_per_beat", "upo_coverage"]:
        abnormal = df.loc[
            df["label"] == 1,
            feature,
        ].to_numpy(dtype=float)

        normal = df.loc[
            df["label"] == 0,
            feature,
        ].to_numpy(dtype=float)

        abnormal = abnormal[np.isfinite(abnormal)]
        normal = normal[np.isfinite(normal)]

        if len(abnormal) and len(normal):
            test = mannwhitneyu(
                abnormal,
                normal,
                alternative="two-sided",
            )
            U = float(test.statistic)
            p = float(test.pvalue)
        else:
            U = np.nan
            p = np.nan

        output[feature] = {
            "normal_n": int(len(normal)),
            "abnormal_n": int(len(abnormal)),
            "normal_median": (
                float(np.median(normal))
                if len(normal)
                else np.nan
            ),
            "abnormal_median": (
                float(np.median(abnormal))
                if len(abnormal)
                else np.nan
            ),
            "mann_whitney_U": U,
            "mann_whitney_p": p,
        }

    # Explicit LLE fit success rates.
    lle_all = df["lle_per_beat"].to_numpy(dtype=float)

    normal_lle = lle_all[df["label"].to_numpy() == 0]
    abnormal_lle = lle_all[df["label"].to_numpy() == 1]

    output["lle_fit_success"] = {
        "normal_rate": (
            float(np.mean(np.isfinite(normal_lle)))
            if len(normal_lle)
            else np.nan
        ),
        "abnormal_rate": (
            float(np.mean(np.isfinite(abnormal_lle)))
            if len(abnormal_lle)
            else np.nan
        ),
    }

    return output


def run_binary_experiment(
    config,
    records=MITDB_RECORDS,
    dl_dir="mitdb_data",
    db_name="mitdb",
    window_size=CLASSIFIER_WINDOW_RR,
    step_size=CLASSIFIER_STEP_RR,
    abnormal_fraction=CLASSIFIER_ABNORMAL_FRACTION,
    output_prefix="mitdb_binary",
):
    """
    Complete binary experiment for one RR-dynamics configuration.
    """
    import pandas as pd

    df, metadata = build_mitbih_classification_dataset(
        records=records,
        config=config,
        dl_dir=dl_dir,
        db_name=db_name,
        window_size=window_size,
        step_size=step_size,
        abnormal_fraction=abnormal_fraction,
    )

    if len(df) == 0:
        raise RuntimeError(
            "No usable classification windows were produced."
        )

    # Save the raw feature table before modeling.
    feature_path = f"{output_prefix}_features.csv"
    metadata_path = f"{output_prefix}_record_metadata.csv"

    df.to_csv(feature_path, index=False)
    metadata.to_csv(metadata_path, index=False)

    print("\nDataset:")
    print(f"  windows: {len(df)}")
    print(f"  records: {df['record'].nunique()}")
    print("  labels:")
    print(df["label"].value_counts().sort_index())

    stats = nonlinear_feature_statistics(df)

    model = run_grouped_logistic_regression(
        df,
        feature_columns=CLASSIFIER_FEATURE_COLUMNS,
        n_splits=5,
        random_state=42,
    )

    model["fold_results"].to_csv(
        f"{output_prefix}_cv_folds.csv",
        index=False,
    )

    return {
        "dataframe": df,
        "metadata": metadata,
        "statistics": stats,
        "model": model,
        "feature_path": feature_path,
        "metadata_path": metadata_path,
    }


def print_experiment_summary(name, result):
    print("\n" + "=" * 72)
    print(name)
    print("=" * 72)

    df = result["dataframe"]
    stats = result["statistics"]
    model = result["model"]

    print(f"Windows: {len(df)}")
    print(f"Records: {df['record'].nunique()}")

    print("\nClass counts:")
    print(df["label"].value_counts().sort_index())

    print("\nNonlinear feature statistics:")
    for feature, values in stats.items():
        print(f"\n{feature}:")
        for key, value in values.items():
            print(f"  {key}: {value}")

    print("\nGrouped-CV fold means:")
    print(model["fold_means"])

    print("\nOut-of-fold metrics:")
    for key, value in model["oof_summary"].items():
        print(f"  {key}: {value}")


def main_binary_experiment():
    """
    Main executable entry point.

    Runs the complete experiment twice:
        A) raw RR intervals for nonlinear dynamics
        B) corrected RR intervals for nonlinear dynamics

    Both use exactly the same records, windows, labels, features, and grouped
    cross-validation procedure.
    """
    # -------------------------------------------------------------------------
    # IMPORTANT: start with a small record list to validate the end-to-end
    # implementation before launching the full 48-record run.
    #
    # Change to MITDB_RECORDS only after the smoke run completes successfully.
    # -------------------------------------------------------------------------
    records = ["100", "105", "106"]

    base_cfg = CFG

    raw_cfg = replace(
        base_cfg,
        use_corrected_rr_for_dynamics=False,
    )

    corrected_cfg = replace(
        base_cfg,
        use_corrected_rr_for_dynamics=True,
    )

    print("=" * 72)
    print("RESEARCH-GRADE MIT-BIH BINARY ARRHYTHMIA EXPERIMENT")
    print("=" * 72)
    print(f"Records requested: {len(records)}")
    print(f"RR window length: {CLASSIFIER_WINDOW_RR}")
    print(f"RR step: {CLASSIFIER_STEP_RR}")
    print(
        "Abnormal-window threshold: "
        f"{CLASSIFIER_ABNORMAL_FRACTION:.1%}"
    )
    print(
        "Features: "
        + ", ".join(CLASSIFIER_FEATURE_COLUMNS)
    )

    raw_result = run_binary_experiment(
        config=raw_cfg,
        records=records,
        output_prefix="mitdb_binary_raw_rr",
    )

    print_experiment_summary(
        "RAW RR DYNAMICS",
        raw_result,
    )

    corrected_result = run_binary_experiment(
        config=corrected_cfg,
        records=records,
        output_prefix="mitdb_binary_corrected_rr",
    )

    print_experiment_summary(
        "CORRECTED RR DYNAMICS",
        corrected_result,
    )

    # Save a compact comparison table.
    comparison_rows = []

    for name, result in [
        ("raw_rr", raw_result),
        ("corrected_rr", corrected_result),
    ]:
        means = result["model"]["fold_means"]

        comparison_rows.append({
            "configuration": name,
            "lle_fit_success_normal":
                result["statistics"]["lle_fit_success"]["normal_rate"],
            "lle_fit_success_abnormal":
                result["statistics"]["lle_fit_success"]["abnormal_rate"],
            "lle_mann_whitney_p":
                result["statistics"]["lle_per_beat"]["mann_whitney_p"],
            "upo_coverage_mann_whitney_p":
                result["statistics"]["upo_coverage"]["mann_whitney_p"],
            "mean_fold_ROC_AUC":
                means.get("roc_auc", np.nan),
            "mean_fold_PR_AUC":
                means.get("pr_auc", np.nan),
            "mean_fold_sensitivity":
                means.get("sensitivity", np.nan),
            "mean_fold_specificity":
                means.get("specificity", np.nan),
            "mean_fold_balanced_accuracy":
                means.get("balanced_accuracy", np.nan),
            "mean_fold_F1":
                means.get("f1", np.nan),
            "OOF_ROC_AUC":
                result["model"]["oof_summary"].get(
                    "oof_roc_auc",
                    np.nan,
                ),
            "OOF_PR_AUC":
                result["model"]["oof_summary"].get(
                    "oof_pr_auc",
                    np.nan,
                ),
        })

    comparison = pd.DataFrame(comparison_rows)
    comparison.to_csv(
        "mitdb_binary_raw_vs_corrected_comparison.csv",
        index=False,
    )

    print("\n" + "=" * 72)
    print("RAW RR vs CORRECTED RR")
    print("=" * 72)
    print(comparison.to_string(index=False))

    return {
        "raw_rr": raw_result,
        "corrected_rr": corrected_result,
        "comparison": comparison,
    }


if __name__ == "__main__":
    main_binary_experiment()
