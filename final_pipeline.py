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

  4) UPO detection (Phase 2A-2D; see section 4 for full provenance):
       - So, Ott, Schiff, Kaplan, Sauer & Grebogi, Phys. Rev. Lett. 76, 4705
         (1996) [PRL]; So, Ott, Sauer, Gluckman, Grebogi & Schiff, Phys. Rev.
         E 55, 5398 (1997) [PRE].
       - Local delay-map Jacobian by least squares over M >= d spatial
         neighbours (PRE Eq. 20), companion structure of the delay map
       - Fixed-point transform G(z,R) = [I-S]^-1 [F(z) - S z] with either the
         PRL norm randomization S = J + kappa R ||F(z)-z||_1 (R d x d, U[-1,1])
         or the PRE tensor randomization S = J + R.[F(z)-z] (PRE Eq. 4).
         In BOTH forms J is the PRE Eq. 20 spatial-neighbour fit; "prl_norm"
         names only the PRL randomization term, not PRL's temporal S_n.
       - Diagonal tube -> histogram -> SOURCE PEAKS (Level A).  Peaks are
         located on the RAW histogram; light smoothing is used only for the
         background threshold (smoothing merged the logistic map's peaks).
       - Period-p block-cyclic transform (PRE Eq. 30), "short" grouping scheme,
         cyclic slab with recon[k,j] = s[(k-j) % p], p-dimensional histogram
       - So et al. surrogate significance (Gaussian-scaled phase shuffle,
         signed W, W0, J(W), r_J) -> Level B significant_uop_candidates
       - Period-1 stability: HYBRID PRL/PRE -- PRL averaging of kappa=0
         Jacobians estimated by the PRE Eq. 20 spatial-neighbour fit (not the
         literal PRL temporal S_n construction)
       - PROJECT EXTENSIONS, labelled as such: tau-step map, verification
         gates (Level C), candidate-centred/monodromy stability, coverage,
         per-point peak rates, classifier features

  Validation: tests/test_upo_*.py (synthetic logistic / skewed Henon systems,
  including the PRE Table I skewed-Henon orbits).

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
from scipy.signal import butter, filtfilt, find_peaks, lfilter, medfilt, welch, resample_poly
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
    pt_rr_average_beats: int = 8       # RR_AVERAGE1/2 window (paper Eqs. 24-25)
    pt_rr_low_factor: float = 0.92     # RR LOW LIMIT  (Eq. 26)
    pt_rr_high_factor: float = 1.16    # RR HIGH LIMIT (Eq. 27)
    # pt_searchback_seconds is the RR MISSED LIMIT factor (1.66, Eq. 28)

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
                                    # (PRE Eq. 20 requires M >= d; PRE Ikeda: M=7)
    so_random_R: int = 100
    so_kappa: float = 3.0
    so_hist_bins: int = 30

    # UPO map mode.  "one_sample" is the source-faithful So et al. one-sample
    # delay map (lag-1 embedding, Cao dimension at lag 1, tau must be 1).
    # "tau_step" is a PROJECT EXTENSION reusing the TDMI tau and embedding.
    upo_map_mode: str = "one_sample"
    # "prl_norm": S = J + kappa R ||F(z)-z||_1 (PRL randomization term);
    # "pre_tensor": S = J + R.[F(z)-z] (PRE Eq. 4).  In both, J is the PRE
    # Eq. 20 spatial-neighbour fit -- "prl_norm" is NOT PRL's temporal S_n.
    so_randomization: str = "prl_norm"
    # PROJECT choice: histogram axes span the attractor range padded by this
    # fraction of its span (see _so_histogram_range).
    so_hist_range_pad: float = 0.10
    # PROJECT choice: period-p slab keeps this percentile of the transformed
    # points closest to the cyclic hyperplane (same rationale as the tube).
    so_slab_percentile: float = 10.0
    so_min_embedded_points: int = 30

    # So et al. surrogate significance for UPO peaks (independent of the LLE
    # surrogate analysis controlled by compute_surrogates).
    so_assess_significance: bool = False
    so_surrogate_count: int = 50     # PRE Sec. IV used 50 surrogates
    so_significance_alpha: float = 0.05  # PROJECT convention on J(deviation)

    # Level C verification gates -- PROJECT EXTENSION, never source detection.
    so_verify_peaks: bool = True
    verify_neighbors: int = 30
    verify_max_residual: float = 0.05     # ||F_hat(z*) - z*|| / std(x)
    verify_min_r2: float = 0.90
    verify_min_support_ratio: float = 2.0  # peak count / median occupied-cell count

    # Classifier feature contract: "source", "significant_source", "extended"
    upo_feature_mode: str = "source"

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

    # Peak coverage -- a PROJECT-derived metric (fraction of 200 sampled
    # attractor points lying near detected peaks); it is NOT a density of
    # periodic orbits.  Radius is
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

# Pan & Tompkins (1985) filters, realised as exact FIR equivalents of the
# published recursive integer filters (all gains normalised to 1 at the
# filter output, never inside the recursion).
#
#   Low-pass  Eq.1/3 : H(z) = (1 - z^-6)^2 / (1 - z^-1)^2 / 36
#                      = [ (1+z^-1+...+z^-5) ]^2 / 36        (11-tap triangle)
#   High-pass Eq.4-6 : intended H(z) = z^-16 - (1/32)(1 - z^-32)/(1 - z^-1)
#                      = z^-16 - (1/32) * sum_{k=0}^{31} z^-k  (32-tap FIR)
#       (The printed Eq.4/6 place the pole of the first-order low-pass at
#       z = -1, i.e. y[n] = 32x[n-16] - [y[n-1] + x[n] - x[n-32]], which has
#       a (1+z^-1) denominator and does not cancel the (1 - z^-32) zeros at
#       z = 1; that is inconsistent with the stated 5 Hz cutoff / gain 32 /
#       delay 16.  The intended filter has a (1 - z^-1) denominator.)
#   Derivative Eq.7/9: (1/8)(-z^-2 - 2z^-1 + 2z + z^2), realised causally
#                      with the paper's two-sample delay:
#                      d[n] = (1/8)(x[n] + 2x[n-1] - 2x[n-3] - x[n-4]).
#                      (T=1 sample; the paper's 1/T scale factor only rescales
#                      the slope and every threshold is data-adaptive.)
_PT_LP_TAPS = np.convolve(np.ones(6), np.ones(6)) / 36.0
_PT_HP_TAPS = np.full(32, -1.0 / 32.0)
_PT_HP_TAPS[16] += 1.0
_PT_DERIV_TAPS = np.array([1.0, 2.0, 0.0, -2.0, -1.0]) / 8.0
# Startup handling (implementation choice; NOT specified by Pan & Tompkins):
# the record is extended on the left by 128 samples of edge replication of
# x[0], filtered from zero state, and the padded outputs are discarded.  This
# does not change any filter transfer function.  It changes the assumed
# pre-record signal (x[0] held constant, instead of the zero-state
# initialisation's implicit x = 0), which affects only the first ~40 (band-pass),
# ~44 (derivative) and ~73 (integrated) output samples; beyond that the output
# equals plain zero-state filtering.  Longest filter memory is 10 (LP) + 31 (HP)
# + 4 (D) + 29 (MWI) = 74 samples, so 128 samples of padding lets the padding's
# own zero-state transient decay before the record begins.
_PT_PAD_SAMPLES = 128


def _pan_tompkins_filter_chain_200hz(x, integration_samples=30):
    """Pan--Tompkins 200-Hz filter chain (causal, deterministic startup).

    Returns (hp, d, integ): band-passed signal, derivative, and the
    moving-window integral of the squared derivative, all aligned with ``x``
    (i.e. each retains its own causal group delay).

    Startup: the input is extended on the left by edge replication of x[0]
    and filtered from zero state; the padded outputs (which contain the
    zero-state transient) are discarded.  This is an implementation choice, not
    something Pan & Tompkins specify: it leaves the transfer functions
    unchanged but replaces the zero-state assumption about the pre-record
    signal with "x[0] held constant".  No sample is ever read at a negative
    index, and no ECG sample is altered.
    """
    x = _as_1d(x)
    n_int = max(1, int(integration_samples))
    xp = np.concatenate([np.full(_PT_PAD_SAMPLES, x[0]), x])
    lp = lfilter(_PT_LP_TAPS, [1.0], xp)
    hp = lfilter(_PT_HP_TAPS, [1.0], lp)
    d = lfilter(_PT_DERIV_TAPS, [1.0], hp)
    integ = lfilter(np.full(n_int, 1.0 / n_int), [1.0], d * d)
    p = _PT_PAD_SAMPLES
    return hp[p:], d[p:], integ[p:]


def _pan_tompkins_processing_delay_200hz(integration_samples=30):
    """Nominal (approximate) processing-delay compensation, in 200-Hz samples.

    Returns (filtered_delay, integrated_delay).  ``filtered_delay`` is the
    band-pass impulse-response peak delay (LP 5 + HP ~16 samples).
    ``integrated_delay`` is derived from the linear filter chain plus the
    moving-window integration: the centroid of the squared-derivative impulse
    response (~23 samples: LP + HP + derivative) plus half the integration
    window ((N-1)/2 = 14.5), rounded (37.486 -> 37 at N = 30).

    This is NOT an exact fixed detector delay.  The integrator acts on the
    squared signal, so the position of the integrated peak relative to the R
    wave depends on QRS morphology; on synthetic beats it was measured at
    roughly 30-45 samples.  The value is only a nominal centre for the
    subsequent +-80 ms re-localization on the original ECG, which is what
    determines the reported fiducial; very broad QRS complexes can fall
    outside that window.
    """
    n_int = max(1, int(integration_samples))
    impulse = np.zeros(256, dtype=float)
    origin = 64
    impulse[origin] = 1.0
    hp, d, _ = _pan_tompkins_filter_chain_200hz(impulse, n_int)
    hp_delay = int(np.argmax(np.abs(hp)) - origin)
    e = d * d
    centroid = float(np.sum(np.arange(len(e)) * e) / np.sum(e))
    integ_delay = int(round(centroid - origin + 0.5 * (n_int - 1)))
    return hp_delay, integ_delay


class _RRTracker:
    """RR_AVERAGE1 / RR_AVERAGE2 bookkeeping (Pan & Tompkins Eqs. 24-29)."""

    def __init__(self, n_avg, low_factor, high_factor, missed_factor):
        self.n = int(n_avg)
        self.low_f, self.high_f, self.missed_f = low_factor, high_factor, missed_factor
        self.recent = []      # RR_AVERAGE1 list: last n RRs regardless of value
        self.selected = []    # RR_AVERAGE2 list: last n RRs within limits

    @property
    def avg1(self):
        return float(np.mean(self.recent)) if self.recent else None

    @property
    def avg2(self):
        return float(np.mean(self.selected)) if self.selected else None

    def limits(self):
        a2 = self.avg2
        if a2 is None:
            return None
        return self.low_f * a2, self.high_f * a2, self.missed_f * a2

    def add(self, rr):
        rr = float(rr)
        lim = self.limits()
        self.recent = (self.recent + [rr])[-self.n:]
        if lim is None or lim[0] <= rr <= lim[1]:
            self.selected = (self.selected + [rr])[-self.n:]

    def irregular(self):
        """True when the last n RRs are not all inside [LOW, HIGH] (Eq. 29).

        With fewer than n RRs the rhythm is not classified as irregular.
        """
        lim = self.limits()
        if lim is None or len(self.recent) < self.n:
            return False
        return not all(lim[0] <= r <= lim[1] for r in self.recent)


def _pan_tompkins_decide(integ, hp, deriv, config, fsd, integration_samples,
                         trace=None):
    """Chronological Pan--Tompkins decision process on the 200-Hz signals.

    Returns accepted QRS indices (in the integrated-signal clock, i.e. still
    carrying the filter delay).  If ``trace`` is a list, one dict per accepted
    beat (index, searchback flag, SPKI before/after, peak height) is appended
    for testing/diagnostics.
    """
    refractory = int(round(config.pt_refractory_seconds * fsd))
    twave_window = int(round(config.pt_twave_window_seconds * fsd))
    init_n = min(len(integ), int(round(config.pt_initial_seconds * fsd)))
    n_int = int(integration_samples)

    # Candidates are separated by the full refractory period.
    candidates, _ = find_peaks(integ, distance=max(1, refractory))
    if len(candidates) == 0:
        return []

    filt_abs = np.abs(hp)
    abs_d = np.abs(deriv)

    def peak_f(q):    # filtered-signal peak inside the integration window
        return float(np.max(filt_abs[max(0, q - n_int + 1):q + 1]))

    def slope_at(q):  # max |slope| inside the integration window ending at q
        return float(np.max(abs_d[max(0, q - n_int + 1):q + 1]))

    init_peaks, _ = find_peaks(integ[:init_n], distance=refractory)
    init_heights = integ[init_peaks] if len(init_peaks) else integ[:init_n]
    spki = 0.25 * float(np.max(init_heights))
    npki = 0.50 * float(np.mean(init_heights))
    f_init = [peak_f(int(q)) for q in init_peaks] or [float(np.max(filt_abs[:init_n]))]
    spkf = 0.25 * float(np.max(f_init))
    npkf = 0.50 * float(np.mean(f_init))

    rr = _RRTracker(config.pt_rr_average_beats, config.pt_rr_low_factor,
                    config.pt_rr_high_factor, float(config.pt_searchback_seconds))

    accepted = []
    accepted_slope = []
    rejected = []        # (index, integ height) rejected as noise since last QRS
                         # (refractory- and T-wave-rejected peaks are not eligible)

    def th_i1():
        return npki + 0.25 * (spki - npki)

    def accept(q, h, fh, searchback):
        nonlocal spki, spkf
        if accepted:
            rr.add(q - accepted[-1])
        accepted.append(q)
        accepted_slope.append(slope_at(q))
        w = 0.25 if searchback else 0.125          # Eq. 16/21 vs Eq. 12/17
        spki_before = spki
        spki = w * h + (1.0 - w) * spki
        if trace is not None:
            trace.append(dict(index=q, searchback=searchback, height=h,
                              spki_before=spki_before, spki_after=spki))
        spkf = w * fh + (1.0 - w) * spkf
        rejected.clear()

    def searchback(q_now):
        """Recover the maximal reserved peak once RR_MISSED has elapsed."""
        while accepted:
            lim = rr.limits()
            if lim is None:
                return
            t_missed = accepted[-1] + int(math.ceil(lim[2]))
            if q_now <= t_missed:
                return
            th2 = 0.5 * th_i1()
            pool = [(c, h) for c, h in rejected
                    if accepted[-1] + refractory <= c <= t_missed and h >= th2]
            if not pool:
                return
            c, h = max(pool, key=lambda t: t[1])
            accept(c, h, peak_f(c), searchback=True)

    for q in candidates:
        q = int(q)
        h = float(integ[q])
        fh = peak_f(q)

        searchback(q)

        thr = th_i1()
        if rr.irregular():
            thr *= 0.5                              # Eq. 22: only THRESHOLD I1 used
        if accepted and q - accepted[-1] < refractory:
            is_signal, keep = False, False
        elif h >= thr:
            twave = (accepted and q - accepted[-1] <= twave_window
                     and slope_at(q) < 0.5 * accepted_slope[-1])
            is_signal, keep = (not twave), False
        else:
            is_signal, keep = False, True

        if is_signal:
            accept(q, h, fh, searchback=False)
        else:
            npki = 0.125 * h + 0.875 * npki
            npkf = 0.125 * fh + 0.875 * npkf
            if keep:
                rejected.append((q, h))

    # Missed-beat check for the tail of the record.
    searchback(len(integ) - 1)
    return accepted


def detect_r_peaks(x, fs, config=CFG):
    """Sequential Pan--Tompkins-style QRS detector (not a literal reproduction).

    Detection is carried out at the paper's 200-Hz sampling rate.  Candidate
    peaks are local maxima of the integrated signal; adaptive SPKI/NPKI
    thresholds, the 200-ms refractory period, 1.66-RR search-back, and
    360-ms slope-based T-wave discrimination are applied chronologically.
    The nominal processing delay is subtracted (see the timing note below) and
    the final fiducial is re-localized on the original ECG.

    Known source-level discrepancy (to be evaluated separately): the actual
    accept/reject decision uses the INTEGRATED channel only.  The 1985 paper
    describes a second set of thresholds (THRESHOLD F1/F2) on the band-pass
    filtered channel and requires a peak to be recognised in BOTH the
    integration and band-pass waveforms.  The filtered-channel running
    estimates (SPKF/NPKF) are maintained here but are not used for any
    decision, and Eq. 23 (halving THRESHOLD F1 in irregular rhythm) is not
    implemented.  Other places where the paper is ambiguous (irregular-rhythm
    threshold halving, search-back window, the 1985 high-pass equation) are
    documented at the corresponding code.  This implementation must not be
    described as reproducing every detail of the original algorithm.
    """
    x = _as_1d(x)
    target_fs = float(config.pt_original_fs)
    if abs(float(fs) - target_fs) < 1e-12:
        xr = x.copy(); back = 1.0
    else:
        from fractions import Fraction
        frac = Fraction(target_fs / float(fs)).limit_denominator(1000)
        # 'line' padding: the default zero padding would create a step
        # (and a large startup transient) whenever the ECG has a DC offset.
        xr = resample_poly(x, frac.numerator, frac.denominator, padtype="line")
        back = float(fs) / target_fs

    fsd = target_fs
    integration_samples = max(1, int(round(config.pt_integration_seconds * fsd)))
    hp, deriv, integ = _pan_tompkins_filter_chain_200hz(xr, integration_samples)
    if len(integ) < 10:
        return np.array([], dtype=int)
    accepted = _pan_tompkins_decide(integ, hp, deriv, config, fsd, integration_samples)

    if not accepted:
        return np.array([], dtype=int)

    # ------------------------------------------------------------------
    # Timing correction and physiological fiducial localization.
    #
    # The published Pan--Tompkins chain is causal, so an integrated-signal
    # candidate lags the R wave.  Subtract a NOMINAL delay (linear filter-chain
    # delay plus the moving-window delay; approximate, because the integrated
    # peak shifts with QRS morphology -- see
    # _pan_tompkins_processing_delay_200hz), map to the original sampling
    # clock, then refine on the ORIGINAL ECG within +-80 ms.  The refined
    # ECG location, not the nominal delay, is the reported fiducial.
    # ------------------------------------------------------------------
    _, integrated_delay = _pan_tompkins_processing_delay_200hz(integration_samples)
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


# =============================================================================
# 4. UPO detection -- So et al. (PRL 1996; PRE 1997)
# =============================================================================
#
# Primary sources (both read in full for this implementation):
#   [PRL] P. So, E. Ott, S. J. Schiff, D. T. Kaplan, T. Sauer, C. Grebogi,
#         "Detecting unstable periodic orbits in chaotic experimental data",
#         Phys. Rev. Lett. 76, 4705 (1996).
#   [PRE] P. So, E. Ott, T. Sauer, B. J. Gluckman, C. Grebogi, S. J. Schiff,
#         "Extracting unstable periodic orbits from chaotic time series data",
#         Phys. Rev. E 55, 5398 (1997).
#
# Provenance is tracked explicitly.  SOURCE components follow the papers:
#   - local delay-map Jacobian by least squares over M spatial neighbours
#     [PRE Eq. (20)], companion structure of the delay map [PRE Eq. (1)]
#   - fixed-point transform G(z,R) = [I-S]^-1 [F(z) - S z]  [PRE Eq. (3)]
#       PRL norm form  : S = J + kappa R ||F(z)-z||_1, R d x d, U[-1,1]
#       PRE tensor form: S = J + R.[F(z)-z], R a d x d x d tensor [PRE Eq. (4)]
#     (J is always the PRE Eq. (20) fit; the PRL's own temporal-difference
#     S_n is not implemented, so "prl_norm" is a PRL randomization on a PRE
#     Jacobian.)
#   - diagonal tube reduction and histogram peaks [PRE Sec. II D]
#   - period-p block transform, "short" grouping scheme, cyclic slab and
#     p-dimensional first-component histogram [PRE Eqs. (27)-(30), Sec. III]
#   - Gaussian-scaled phase-shuffle surrogates, signed W, W0 (J(W0)=0.5),
#     J(W) and r_J = W/W0 [PRL; PRE Sec. IV]
#
# HYBRID PRL/PRE (neither purely source nor a project extension):
#   - period-1 stability: the PRL procedure of averaging kappa=0 S_n over the
#     data points whose transformed values form the peak, but with each S_n
#     estimated by the PRE Eq. (20) spatial-neighbour fit instead of the PRL
#     temporal-difference construction.  The literal PRL variant is not
#     implemented and gives materially different values on the Henon map.
#
# PROJECT EXTENSIONS (never described as the So et al. method):
#   - tau-step delay map (map_mode="tau_step")
#   - candidate residual, R^2 and support-ratio verification gates
#   - candidate-centred / monodromy stability (incl. p > 1)
#   - peak coverage, per-point peak rates, classifier feature construction
#   - histogram range fixed to the padded attractor range, percentile tube
#     and slab widths, and the median+MAD peak threshold (the papers state
#     only "a small cross-section tube", "a thin slab" and "look for peaks")
#
# Result API.  Detection results never contain a generic "candidates" field:
#   Level A  source_peak_candidates         histogram peaks (SOURCE)
#   Level B  significant_uop_candidates     None = significance NOT assessed;
#                                           []   = assessed, none significant
#   Level C  verification_gated_candidates  PROJECT EXTENSION gates
#            verification_failed_peaks
# A Level-A source peak remains a source peak whatever Level C decides.

SO_PRL_CITATION = ("So, Ott, Schiff, Kaplan, Sauer & Grebogi, "
                   "Phys. Rev. Lett. 76, 4705 (1996)")
SO_PRE_CITATION = ("So, Ott, Sauer, Gluckman, Grebogi & Schiff, "
                   "Phys. Rev. E 55, 5398 (1997)")

UPO_MAP_MODES = ("one_sample", "tau_step")
SO_RANDOMIZATIONS = ("prl_norm", "pre_tensor")
UPO_FEATURE_MODES = ("source", "significant_source", "extended")

UPO_STATUS_OK = "ok"
# "No peak" outcomes: the analysis ran and found no histogram peak, so peak
# counts / rates / coverages are 0.  KNOWN LIMITATION (open for Phase 2E):
# tube_empty and too_few_in_tube are treated as "no peaks" here, but whether
# they are genuinely peak-free or a data-sufficiency failure has not been
# resolved scientifically.  No threshold has been invented to decide it.
UPO_NO_PEAK_STATUSES = frozenset({"no_peaks", "tube_empty", "too_few_in_tube"})
# Genuine pipeline failures: UPO features are NaN, never 0.
UPO_FAILURE_STATUSES = frozenset({
    "constant_data", "nonfinite_input", "too_few_points",
    "embedding_not_saturated", "no_valid_transforms", "analysis_error",
})

LEVEL_A_FIELD = "source_peak_candidates"
LEVEL_B_FIELD = "significant_uop_candidates"
LEVEL_C_PASS_FIELD = "verification_gated_candidates"
LEVEL_C_FAIL_FIELD = "verification_failed_peaks"

UPO_LEVELS = {
    "A": {"field": LEVEL_A_FIELD, "provenance": "SOURCE",
          "meaning": "histogram peaks of the So-transformed data; not by "
                     "themselves significant, unstable, verified or clinical"},
    "B": {"field": LEVEL_B_FIELD, "provenance": "SOURCE",
          "meaning": "source peaks whose deviation above the surrogate mean "
                     "is significant under the So et al. J(W) distribution"},
    "C": {"field": [LEVEL_C_PASS_FIELD, LEVEL_C_FAIL_FIELD],
          "provenance": "PROJECT EXTENSION",
          "meaning": "source peaks passing / failing project residual, R^2 "
                     "and support-ratio gates; never redefines Level A"},
}

UPO_FEATURE_CONTRACT = {
    "source": (
        "lle_per_beat",
        "source_peak_count",
        "source_peaks_per_point",
        "source_peak_coverage",
    ),
    "significant_source": (
        "lle_per_beat",
        "significant_peak_count",
        "significant_peaks_per_point",
        "significant_peak_coverage",
        "source_rJ",
    ),
    "extended": (
        "lle_per_beat",
        "verified_unstable_count",
        "verified_unstable_per_point",
        "verified_unstable_coverage",
    ),
}

UPO_FEATURE_DEFINITIONS = {
    "lle_per_beat": "Rosenstein largest Lyapunov exponent (independent of UPO analysis)",
    "source_peak_count": "number of So histogram peaks (Level A), all analysed periods",
    "source_peaks_per_point": "source_peak_count / number of UPO-embedded points (a rate, not a density)",
    "source_peak_coverage": "PROJECT metric: fraction of sampled attractor points near a source peak",
    "significant_peak_count": "number of Level B (surrogate-significant) source peaks",
    "significant_peaks_per_point": "significant_peak_count / number of UPO-embedded points",
    "significant_peak_coverage": "PROJECT metric: coverage by Level B peaks",
    "source_rJ": "So et al. r_J = W / W0 for the period-1 histogram",
    "verified_unstable_count": "PROJECT EXTENSION: Level C gated peaks with unstable monodromy",
    "verified_unstable_per_point": "verified_unstable_count / number of UPO-embedded points",
    "verified_unstable_coverage": "PROJECT metric: coverage by verified unstable peaks",
}


def validate_upo_config(config=CFG):
    """Reject ambiguous or contradictory UPO settings instead of guessing."""
    if config.upo_map_mode not in UPO_MAP_MODES:
        raise ValueError(f"upo_map_mode must be one of {UPO_MAP_MODES}, "
                         f"got {config.upo_map_mode!r}")
    if config.so_randomization not in SO_RANDOMIZATIONS:
        raise ValueError(f"so_randomization must be one of {SO_RANDOMIZATIONS}, "
                         f"got {config.so_randomization!r}")
    if config.upo_feature_mode not in UPO_FEATURE_MODES:
        raise ValueError(f"upo_feature_mode must be one of {UPO_FEATURE_MODES}, "
                         f"got {config.upo_feature_mode!r}")
    if config.upo_feature_mode == "significant_source" and not config.so_assess_significance:
        raise ValueError("upo_feature_mode='significant_source' requires "
                         "so_assess_significance=True")
    if any(int(p) < 1 for p in config.so_periods):
        raise ValueError("so_periods must contain positive integers")
    return True


def resolve_upo_map(map_mode, tau):
    """Map-mode metadata.  one_sample is the So et al. map and requires tau=1."""
    if map_mode not in UPO_MAP_MODES:
        raise ValueError(f"map_mode must be one of {UPO_MAP_MODES}, got {map_mode!r}")
    tau = int(tau)
    if tau < 1:
        raise ValueError("tau must be >= 1")
    if map_mode == "one_sample":
        if tau != 1:
            raise ValueError(
                "map_mode='one_sample' is the source-faithful So et al. one-sample "
                f"delay map and requires lag-1 embedded vectors (tau=1); got tau={tau}. "
                "Use map_mode='tau_step' (PROJECT EXTENSION) for tau > 1.")
        return {"map_mode": "one_sample", "tau": 1, "step": 1,
                "source_faithful_map": True,
                "map_label": "SOURCE: one-sample delay map z(n) -> z(n+1) (So et al.)"}
    return {"map_mode": "tau_step", "tau": tau, "step": tau,
            "source_faithful_map": False,
            "map_label": (f"PROJECT EXTENSION: tau-step delay map z(n) -> z(n+{tau}) "
                          "on a lag-tau embedding; NOT the So et al. source map")}


def _detection_label(period, map_info, randomization):
    form = "PRL norm" if randomization == "prl_norm" else "PRE tensor"
    kind = "fixed-point transform" if period == 1 else f"period-{period} block transform"
    return (f"So et al. histogram source peaks, period {period} "
            f"({kind}, {form} randomization, PRE Eq. 20 Jacobian) [{map_info['map_label']}]")


# -------------------------------------------------------------------------
# Local delay-map Jacobian [PRE Eq. (20)]
# -------------------------------------------------------------------------

def _companion(grad):
    grad = np.asarray(grad, dtype=float)
    d = grad.size
    J = np.zeros((d, d))
    J[0, :] = grad
    if d > 1:
        J[1:, :-1] = np.eye(d - 1)
    return J


def _delay_map_jacobian_from_local_fit(embedded, idx, neighbor_indices, step=1):
    """
    Local delay-map Jacobian [PRE Eq. (20)]:

        w_1^k(n+s) - z_1(n+s) = grad f(z(n)) . [w^k - z(n)],   k = 1..M, M >= d

    i.e. the response X[nbr+s, 0] - X[idx+s, 0] is regressed (least squares)
    on X[nbr] - X[idx].  Only the first row of the companion matrix is
    estimated; rows 2..d are the exact delay shift [PRE Eq. (1)].
    """
    X = np.asarray(embedded, dtype=float)
    d = X.shape[1]
    step = int(step)
    if step < 1:
        raise ValueError("step must be >= 1")
    neigh = np.asarray(neighbor_indices, dtype=int)
    neigh = neigh[(neigh >= 0) & (neigh + step < len(X))]
    if len(neigh) < d or idx + step >= len(X):
        return None
    A = X[neigh] - X[idx]
    b = X[neigh + step, 0] - X[idx + step, 0]
    try:
        grad, _, rank, _ = np.linalg.lstsq(A, b, rcond=None)
    except np.linalg.LinAlgError:
        return None
    if rank < d or not np.all(np.isfinite(grad)):
        return None
    return _companion(grad)


def _spatial_neighbors(embedded, idx, K, exclude=1, valid_limit=None):
    """K nearest phase-space neighbours of X[idx], excluding |j-idx| <= exclude."""
    X = np.asarray(embedded, dtype=float)
    lim = len(X) if valid_limit is None else int(valid_limit)
    d = np.linalg.norm(X[:lim] - X[idx], axis=1)
    lo = max(0, idx - max(int(exclude), 0))
    hi = min(lim, idx + max(int(exclude), 0) + 1)
    d[lo:hi] = np.inf
    if idx < lim:
        d[idx] = np.inf
    order = np.argsort(d, kind="stable")
    order = order[np.isfinite(d[order])]
    return order[:int(K)]


def _local_jacobians(X, step, M, exclude):
    """Jacobian for every point that has an image z(n+step)."""
    n, d = X.shape
    M = max(int(M), d)
    limit = n - step
    Js = [None] * n
    for i in range(limit):
        neigh = _spatial_neighbors(X, i, K=M, exclude=exclude, valid_limit=limit)
        Js[i] = _delay_map_jacobian_from_local_fit(X, i, neigh, step=step)
    return Js


# -------------------------------------------------------------------------
# So transforms
# -------------------------------------------------------------------------

def _random_R_matrix(rng, d, kappa):
    """PRL randomization: d x d matrix, entries i.i.d. uniform on [-1, 1]."""
    return rng.uniform(-1.0, 1.0, size=(d, d))


def _random_R_tensor(rng, d, kappa):
    """PRE randomization: d x d x d tensor, entries i.i.d. uniform on [-1, 1]."""
    return rng.uniform(-1.0, 1.0, size=(d, d, d))


def _draw_R(rng, count, d, randomization):
    if randomization == "prl_norm":
        return rng.uniform(-1.0, 1.0, size=(count, d, d))
    if randomization == "pre_tensor":
        return rng.uniform(-1.0, 1.0, size=(count, d, d, d))
    raise ValueError(f"unknown randomization {randomization!r}")


def _l1_norm(v):
    return np.sum(np.abs(v))


def so_perturbation(R, delta, kappa, randomization="prl_norm"):
    """
    Randomizing term added to the Jacobian.

      prl_norm  : kappa * R * ||delta||_1         (R: [..., d, d])   [PRL]
      pre_tensor: kappa * (R . delta)_ij = kappa * sum_k R_ijk delta_k
                                                  (R: [..., d, d, d]) [PRE Eq. (4)]

    delta is F(z) - z for fixed points and F(z_k) - z_{k+1} for period p.
    "prl_norm" refers to this randomization term only: the J it is added to
    is the PRE Eq. 20 fit, never the PRL temporal-difference S_n.
    PRE writes R.[F(z)-z] with R free; the tensor here is kappa times a
    U[-1,1] tensor, matching the PRE 1-D example R = k*eta, eta ~ U[-1,1].
    """
    R = np.asarray(R, dtype=float)
    delta = np.asarray(delta, dtype=float)
    if randomization == "prl_norm":
        return kappa * R * _l1_norm(delta)
    if randomization == "pre_tensor":
        return kappa * np.einsum("...ijk,k->...ij", R, delta)
    raise ValueError(f"unknown randomization {randomization!r}")


def so_fixed_point_transform(z, Fz, J, R, kappa, randomization="prl_norm"):
    """
    So fixed-point transform [PRE Eq. (3)]:
        z_hat = G(z, R) = [I - S]^-1 [F(z) - S z]
    with S = J + so_perturbation(R, F(z) - z).
    """
    z = np.asarray(z, dtype=float)
    Fz = np.asarray(Fz, dtype=float)
    S = np.asarray(J, dtype=float) + so_perturbation(R, Fz - z, kappa, randomization)
    try:
        return np.linalg.solve(np.eye(len(z)) - S, Fz - S @ z)
    except np.linalg.LinAlgError:
        return None


def so_fixed_point_transform_tensor(z, Fz, J, R_tensor, kappa):
    """PRE tensor-form fixed-point transform (see so_fixed_point_transform)."""
    return so_fixed_point_transform(z, Fz, J, R_tensor, kappa, randomization="pre_tensor")


def _batched_solve(A, b):
    """Solve A x = b for stacked systems; singular systems give NaN rows."""
    try:
        return np.linalg.solve(A, b[..., None])[..., 0]
    except np.linalg.LinAlgError:
        out = np.full(b.shape, np.nan)
        for i in range(len(A)):
            try:
                out[i] = np.linalg.solve(A[i], b[i])
            except np.linalg.LinAlgError:
                pass
        return out


def _batched_fixed_point_transform(z, Fz, J, Rs, kappa, randomization):
    d = len(z)
    S = J[None, :, :] + so_perturbation(Rs, Fz - z, kappa, randomization)
    A = np.eye(d)[None, :, :] - S
    rhs = Fz[None, :] - np.einsum("rij,j->ri", S, z)
    return _batched_solve(A, rhs)


def _so_period_block_system(Z, FZ, S_list):
    """
    Block-cyclic system of PRE Eq. (30).  Block row k (0-based) reads
        -S_k zhat_k + zhat_{k+1 mod p} = F(z_k) - S_k z_k,
    the linearization F(z_k) = z*(k+1) + S_k [z_k - z*(k)] of PRE Eq. (29).
    """
    p = len(Z)
    d = np.asarray(Z[0]).size
    B = np.zeros((p * d, p * d))
    rhs = np.zeros(p * d)
    for k in range(p):
        kp1 = (k + 1) % p
        S = S_list[k]
        B[k*d:(k+1)*d, k*d:(k+1)*d] = -S
        B[k*d:(k+1)*d, kp1*d:(kp1+1)*d] += np.eye(d)
        rhs[k*d:(k+1)*d] = FZ[k] - S @ Z[k]
    return B, rhs


def _so_period_transform(Z, FZ, J_list, R_list, kappa, randomization="prl_norm"):
    """
    Period-p So transform [PRE Eq. (30)] for one combination of test points,
    with S(z_k, z_{k+1}, R_k) = J(z_k) + R_k.[F(z_k) - z_{k+1}].
    Returns the (p, d) array (zhat_1, ..., zhat_p) or None if singular.
    """
    p = len(Z)
    S_list = [np.asarray(J_list[k], float)
              + so_perturbation(R_list[k], np.asarray(FZ[k]) - np.asarray(Z[(k+1) % p]),
                                kappa, randomization)
              for k in range(p)]
    B, rhs = _so_period_block_system(Z, FZ, S_list)
    try:
        return np.linalg.solve(B, rhs).reshape(p, -1)
    except np.linalg.LinAlgError:
        return None


def _batched_period_transform(Z, FZ, J_list, Rs, kappa, randomization):
    """Rs: list of p stacks of random matrices/tensors, each of length r."""
    p = len(Z)
    d = Z[0].size
    r = len(Rs[0])
    B = np.zeros((r, p*d, p*d))
    rhs = np.zeros((r, p*d))
    eye = np.eye(d)
    for k in range(p):
        kp1 = (k + 1) % p
        S = J_list[k][None] + so_perturbation(Rs[k], FZ[k] - Z[kp1], kappa, randomization)
        B[:, k*d:(k+1)*d, k*d:(k+1)*d] = -S
        B[:, k*d:(k+1)*d, kp1*d:(kp1+1)*d] += eye
        rhs[:, k*d:(k+1)*d] = FZ[k][None] - np.einsum("rij,j->ri", S, Z[k])
    return _batched_solve(B, rhs).reshape(r, p, d)


# -------------------------------------------------------------------------
# Delay-coordinate cyclic structure [PRE Eqs. (27)-(28)]
# -------------------------------------------------------------------------

def cyclic_orbit_matrix(s, d):
    """
    Delay-coordinate period-p orbit from its p independent values.

    With z(n) = (x(n), x(n-1), ..., x(n-d+1)) and z*(k+1) = F(z*(k)), the
    j-th component (0-based) of the k-th orbit point is x at time k-j, so
        recon[k, j] = s[(k - j) % p].
    (The earlier (k + j) % p indexing was wrong for p >= 3.)
    """
    s = np.asarray(s, dtype=float)
    p = len(s)
    k = np.arange(p)[:, None]
    j = np.arange(int(d))[None, :]
    return s[(k - j) % p]


def _cyclic_slab_reduction(Zhat):
    """
    First-component reduction of a transformed period-p point [PRE Sec. III B]:
    s = (zhat_1(1), ..., zhat_1(p)), reconstructed with cyclic_orbit_matrix;
    the residual is the distance from the cyclic hyperplane (0 on it).
    """
    Zhat = np.asarray(Zhat, dtype=float)
    p, d = Zhat.shape
    s = Zhat[:, 0].copy()
    return s, float(np.linalg.norm(Zhat - cyclic_orbit_matrix(s, d)))


def _batched_cyclic_residual(Zhat):
    r, p, d = Zhat.shape
    k = np.arange(p)[:, None]
    j = np.arange(d)[None, :]
    recon = Zhat[:, :, 0][:, (k - j) % p]
    return np.linalg.norm((Zhat - recon).reshape(r, -1), axis=1)


def _project_fixed_point_tube(transformed, percentile=10.0):
    """
    Diagonal tube reduction [PRE Sec. II D]: every delay-coordinate fixed
    point lies on z_1 = ... = z_d.  The perpendicular distance to the
    diagonal is ||z - mean(z) 1||; the PROJECT choice here keeps the closest
    `percentile` % (a self-calibrating tube width -- a fixed width collapses
    at higher d).  Returns (scalar diagonal coordinate, keep mask).
    """
    Z = np.asarray(transformed, dtype=float)
    if len(Z) == 0:
        return np.empty(0), np.zeros(0, dtype=bool)
    distance = np.linalg.norm(Z - np.mean(Z, axis=1, keepdims=True), axis=1)
    keep = distance <= np.percentile(distance, percentile)
    return np.mean(Z[keep], axis=1), keep


def canonical_cyclic_rotation(s):
    """Rotation of a cyclic sequence that starts at its largest element."""
    s = np.asarray(s, dtype=float)
    return np.roll(s, -int(np.argmax(s)))


def minimal_cyclic_period(s, tol):
    """Smallest q dividing p with s invariant under rotation by q."""
    s = np.asarray(s, dtype=float)
    p = len(s)
    for q in range(1, p + 1):
        if p % q == 0 and np.max(np.abs(s - np.roll(s, -q))) <= tol:
            return q
    return p


# -------------------------------------------------------------------------
# Histogram and peaks
# -------------------------------------------------------------------------

def _so_histogram_range(embedded, pad_fraction):
    """
    PROJECT choice: histogram axes span the observed attractor coordinate
    range, padded by pad_fraction of its span on each side.  Transformed
    points far outside (near-singular I - S) are excluded; without this the
    min/max of the transformed data can stretch the grid so that distinct
    peaks merge.  Periodic orbits slightly off the attractor remain inside
    the padding [PRE: the Henon fixed point -1.903 lies off the attractor].
    """
    x = np.asarray(embedded, dtype=float)[:, 0]
    lo, hi = float(np.min(x)), float(np.max(x))
    span = hi - lo
    if not np.isfinite(span) or span <= 1e-12:
        return None
    pad = float(pad_fraction) * span
    return (lo - pad, hi + pad)


def _histogram_peak_threshold(hist, config):
    """
    PROJECT peak rule (the papers say only "look for peaks"): threshold
    = max(median + so_peak_sigma * MAD-scale of a lightly smoothed
    histogram, so_peak_min_count).  Smoothing is used for the background
    estimate only; peaks are located on the raw histogram.
    """
    from scipy.ndimage import gaussian_filter
    h = np.asarray(hist, dtype=float)
    sm = gaussian_filter(h, sigma=0.5, mode="constant")
    bg = float(np.median(sm))
    scale = max(float(np.median(np.abs(sm - bg))) / 0.6744897501960817, 1e-12)
    return max(bg + config.so_peak_sigma * scale, float(config.so_peak_min_count))


def _histogram_peaks(hist, threshold):
    """Local maxima (3^p neighbourhood) with count >= threshold; one per plateau."""
    from scipy.ndimage import maximum_filter, label
    h = np.asarray(hist, dtype=float)
    local = (h == maximum_filter(h, size=3, mode="constant", cval=-np.inf)) & (h >= threshold)
    lab, nlab = label(local, structure=np.ones((3,) * h.ndim))
    cells = []
    for lbl in range(1, nlab + 1):
        members = np.argwhere(lab == lbl)
        cells.append(tuple(int(v) for v in members[0]))
    cells.sort(key=lambda c: -h[c])
    return cells


def _in_cell(points, edges, cell):
    """Boolean mask of points (n, p) inside histogram cell `cell`."""
    pts = np.asarray(points, dtype=float).reshape(len(points), -1)
    mask = np.ones(len(pts), dtype=bool)
    for j, c in enumerate(cell):
        e = edges[j]
        upper = pts[:, j] <= e[c + 1] if c == len(e) - 2 else pts[:, j] < e[c + 1]
        mask &= (pts[:, j] >= e[c]) & upper
    return mask


def _histogram_on_edges(points, edges):
    pts = np.asarray(points, dtype=float).reshape(len(points), -1) if len(points) else \
        np.empty((0, len(edges)))
    hist, _ = np.histogramdd(pts, bins=[np.asarray(e) for e in edges])
    return hist


# -------------------------------------------------------------------------
# Detection results
# -------------------------------------------------------------------------

def _new_detection_result(period, map_info, config, status, n_points, d):
    return {
        "detection_label": _detection_label(period, map_info, config.so_randomization),
        "period": int(period),
        "map_mode": map_info["map_mode"],
        "source_faithful_map": bool(map_info["source_faithful_map"]),
        "map_label": map_info["map_label"],
        "tau": int(map_info["tau"]),
        "step": int(map_info["step"]),
        "randomization": config.so_randomization,
        "status": status,
        "n_points": int(n_points),
        "embedding_dimension": int(d),
        "levels": UPO_LEVELS,
        LEVEL_A_FIELD: [],
        "significance_assessed": False,
        "significance": None,
        LEVEL_B_FIELD: None,
        "verification_assessed": False,
        LEVEL_C_PASS_FIELD: None,
        LEVEL_C_FAIL_FIELD: None,
        "citations": (SO_PRL_CITATION, SO_PRE_CITATION),
    }


def _input_status(X, config, step, period):
    n, d = X.shape
    if not np.all(np.isfinite(X)):
        return "nonfinite_input"
    if n and float(np.ptp(X[:, 0])) <= 1e-12:
        return "constant_data"
    need = max(int(config.so_min_embedded_points), period * step + max(config.so_jacobian_neighbors, d) + 2)
    if n < need:
        return "too_few_points"
    return None


SOURCE_STABILITY_PROVENANCE = (
    "HYBRID PRL/PRE: PRL averaging of kappa=0 Jacobians estimated using the "
    "PRE Eq. 20 spatial-neighbor least-squares fit; hybrid PRL/PRE "
    "implementation, not the literal PRL temporal S_n construction")


def source_period1_stability(jacobians, member_indices):
    """
    HYBRID PRL/PRE period-1 stability.

    Procedure from the PRL: average S_n with kappa = 0 over the distinct data
    points whose transformed values lie in the peak, and take eigenvalues
    ("Lyapunov numbers" per map step).  Here each S_n (kappa = 0) is the
    local Jacobian estimated by the PRE Eq. 20 spatial-neighbour least-squares
    fit -- NOT the literal PRL construction of S_n from temporal differences
    of consecutive delay vectors.  The result therefore does not reproduce
    the literal PRL stability procedure.  (The result key remains
    "source_stability" for API stability; its "provenance" field is
    authoritative.)
    """
    Js = [jacobians[i] for i in np.unique(np.asarray(member_indices, dtype=int))
          if jacobians[i] is not None]
    if not Js:
        return None
    S_bar = np.mean(Js, axis=0)
    eig = np.linalg.eigvals(S_bar)
    mod = np.sort(np.abs(eig))[::-1]
    return {
        "provenance": SOURCE_STABILITY_PROVENANCE,
        "jacobian_method": "PRE Eq. 20 spatial-neighbor least squares",
        "literal_prl_temporal_S_n": False,
        "mean_S": S_bar,
        "eigenvalues": eig,
        "lyapunov_numbers": mod,
        "n_points_averaged": int(len(Js)),
        "unstable": bool(mod[0] > 1.0 + 1e-8),
        "saddle": bool(mod[0] > 1.0 + 1e-8 and mod[-1] < 1.0 - 1e-8),
    }


def detect_so_fixed_points(embedded, tau=1, config=CFG, rng=None, map_mode=None,
                           hist_range=None):
    """
    So et al. period-1 detection:
        So transformation -> transformed points -> diagonal tube
        -> histogram -> peaks (Level A source peaks).

    Stability is NOT part of the Level A detection criterion: after a peak
    has been detected, a separate diagnostic ("source_stability", HYBRID
    PRL/PRE, see source_period1_stability) is attached to it.  It never adds,
    removes or relabels a peak and may be None.

    `tau` is the embedding lag of `embedded`.  map_mode defaults to
    config.upo_map_mode; one_sample (SOURCE) requires tau=1.
    """
    if rng is None:
        rng = np.random.default_rng(config.random_seed)
    map_info = resolve_upo_map(config.upo_map_mode if map_mode is None else map_mode, tau)
    step = map_info["step"]
    X = np.asarray(embedded, dtype=float)
    if X.ndim != 2:
        raise ValueError("embedded must be a 2-D array (points x dimension)")
    n, d = X.shape
    status = _input_status(X, config, step, 1)
    result = _new_detection_result(1, map_info, config, status or UPO_STATUS_OK, n, d)
    result.update({"transformed": np.empty((0, d)), "scalar": np.empty(0),
                   "source_indices": np.empty(0, dtype=int), "histogram": None, "edges": None})
    if status is not None:
        return result
    randomization = config.so_randomization
    Js = _local_jacobians(X, step, config.so_jacobian_neighbors, max(1, step))
    out, src = [], []
    for i in range(n - step):
        if Js[i] is None:
            continue
        zh = _batched_fixed_point_transform(
            X[i], X[i + step], Js[i], _draw_R(rng, config.so_random_R, d, randomization),
            config.so_kappa, randomization)
        ok = np.all(np.isfinite(zh), axis=1)
        out.append(zh[ok])
        src.append(np.full(int(ok.sum()), i, dtype=int))
    result["jacobians"] = Js
    transformed = np.concatenate(out) if out else np.empty((0, d))
    source_idx = np.concatenate(src) if src else np.empty(0, dtype=int)
    if len(transformed) == 0:
        result["status"] = "no_valid_transforms"
        return result
    scalar, keep = _project_fixed_point_tube(transformed, config.so_diagonal_percentile)
    rng_range = hist_range if hist_range is not None else _so_histogram_range(X, config.so_hist_range_pad)
    edges = [np.linspace(rng_range[0], rng_range[1], config.so_hist_bins + 1)]
    inside = (scalar >= edges[0][0]) & (scalar <= edges[0][-1])
    transformed_kept = transformed[keep][inside]
    scalar = scalar[inside]
    source_kept = source_idx[keep][inside]
    hist = _histogram_on_edges(scalar[:, None], edges)
    result.update({"transformed": transformed_kept, "scalar": scalar,
                   "source_indices": source_kept, "histogram": hist, "edges": edges,
                   "reduced": scalar[:, None], "histogram_range": tuple(rng_range)})
    if len(scalar) == 0:
        result["status"] = "tube_empty"
        return result
    if len(scalar) < config.so_peak_min_count:
        result["status"] = "too_few_in_tube"
        return result
    threshold = _histogram_peak_threshold(hist, config)
    result["peak_threshold"] = threshold
    peaks = []
    for cell in _histogram_peaks(hist, threshold):
        m = _in_cell(scalar[:, None], edges, cell)
        pts = transformed_kept[m]
        if len(pts) == 0:
            continue
        c = float(np.mean(np.mean(pts, axis=1)))
        members = np.unique(source_kept[m])
        peaks.append({
            "period": 1,
            "provenance": "SOURCE",
            "level": "A",
            "location": np.full(d, c),
            "scalar_location": c,
            "orbit_coordinates": np.array([c]),
            "orbit_points": np.full((1, d), c),
            "minimal_period": 1,
            "histogram_count": int(hist[cell]),
            "peak_cell": cell,
            "transformed_points_in_peak": int(len(pts)),
            "diagonal_error": float(np.median(np.std(pts, axis=1))),
            "source_member_indices": members,
            # Diagnostic attached to an already-detected peak; not a detection criterion.
            "source_stability": source_period1_stability(Js, members),
        })
    result[LEVEL_A_FIELD] = peaks
    result["status"] = UPO_STATUS_OK if peaks else "no_peaks"
    return result


def detect_so_period_p(embedded, period, tau=1, config=CFG, rng=None, map_mode=None,
                       hist_range=None):
    """
    So et al. period-p detection (p >= 2) [PRE Sec. III]:
      - backbone of p temporally consecutive points, K nearest spatial
        neighbours of each -> p clusters of K+1 test points ("short" scheme)
      - every combination through the block-cyclic transform, PRE Eq. (30),
        with many random R per combination
      - thin slab about the cyclic p-plane, first-component reduction
        s = (zhat_1(1), ..., zhat_1(p)), p-dimensional histogram, peaks.
    Each histogram peak is a Level A source peak.  Cyclic rotations of one
    orbit appear as separate peaks; `canonical_orbit` groups them.
    """
    period = int(period)
    if period < 2:
        raise ValueError("Use detect_so_fixed_points for period 1.")
    if period > 3:
        raise ValueError("Full p-dimensional histogram becomes impractical for p>3. Use p<=3.")
    if rng is None:
        rng = np.random.default_rng(config.random_seed + period)
    map_info = resolve_upo_map(config.upo_map_mode if map_mode is None else map_mode, tau)
    step = map_info["step"]
    X = np.asarray(embedded, dtype=float)
    if X.ndim != 2:
        raise ValueError("embedded must be a 2-D array (points x dimension)")
    n, d = X.shape
    status = _input_status(X, config, step, period)
    result = _new_detection_result(period, map_info, config, status or UPO_STATUS_OK, n, d)
    result.update({"reduced": np.empty((0, period)), "histogram": None, "edges": None})
    if status is not None:
        return result
    randomization = config.so_randomization
    limit = n - step
    Js = _local_jacobians(X, step, config.so_jacobian_neighbors, max(1, step))
    clusters_of = {}

    def cluster(i):
        if i not in clusters_of:
            nb = _spatial_neighbors(X, i, K=config.so_neighbors_K, exclude=max(1, step),
                                    valid_limit=limit)
            members = np.asarray([i, *nb], dtype=int)
            clusters_of[i] = members[[Js[m] is not None for m in members]]
        return clusters_of[i]

    starts = list(range(n - period * step))
    if config.so_max_backbones is not None and len(starts) > config.so_max_backbones:
        starts = list(np.unique(np.linspace(0, len(starts) - 1,
                                            config.so_max_backbones).astype(int)))
    reduced_all, resid_all = [], []
    r = int(config.so_random_R_period_p)
    for start in starts:
        cl = [cluster(start + k * step) for k in range(period)]
        if any(len(c) == 0 for c in cl):
            continue
        for combo in itertools.product(*cl):
            Z = [X[c] for c in combo]
            FZ = [X[c + step] for c in combo]
            J_list = [Js[c] for c in combo]
            Rs = [_draw_R(rng, r, d, randomization) for _ in range(period)]
            Zhat = _batched_period_transform(Z, FZ, J_list, Rs, config.so_kappa, randomization)
            ok = np.all(np.isfinite(Zhat.reshape(r, -1)), axis=1)
            if np.any(ok):
                reduced_all.append(Zhat[ok][:, :, 0])
                resid_all.append(_batched_cyclic_residual(Zhat[ok]))
    if not reduced_all:
        result["status"] = "no_valid_transforms"
        return result
    reduced = np.concatenate(reduced_all)
    resid = np.concatenate(resid_all)
    keep = resid <= np.percentile(resid, config.so_slab_percentile)
    reduced = reduced[keep]
    rng_range = hist_range if hist_range is not None else _so_histogram_range(X, config.so_hist_range_pad)
    edges = [np.linspace(rng_range[0], rng_range[1], config.so_hist_bins + 1)] * period
    inside = np.all((reduced >= rng_range[0]) & (reduced <= rng_range[1]), axis=1)
    reduced = reduced[inside]
    hist = _histogram_on_edges(reduced, edges)
    result.update({"reduced": reduced, "histogram": hist, "edges": edges,
                   "histogram_range": tuple(rng_range), "jacobians": Js})
    if len(reduced) == 0:
        result["status"] = "tube_empty"
        return result
    if len(reduced) < config.so_peak_min_count:
        result["status"] = "too_few_in_tube"
        return result
    threshold = _histogram_peak_threshold(hist, config)
    result["peak_threshold"] = threshold
    bin_width = (rng_range[1] - rng_range[0]) / config.so_hist_bins
    peaks = []
    for cell in _histogram_peaks(hist, threshold):
        m = _in_cell(reduced, edges, cell)
        if not np.any(m):
            continue
        s = np.mean(reduced[m], axis=0)
        peaks.append({
            "period": period,
            "provenance": "SOURCE",
            "level": "A",
            "orbit_coordinates": s,
            "canonical_orbit": canonical_cyclic_rotation(s),
            "minimal_period": minimal_cyclic_period(s, tol=bin_width),
            "orbit_points": cyclic_orbit_matrix(s, d),
            "location": cyclic_orbit_matrix(s, d)[0],
            "histogram_count": int(hist[cell]),
            "peak_cell": cell,
            "transformed_points_in_peak": int(np.sum(m)),
            # The papers define no separate stability procedure for p > 1.
            "source_stability": None,
        })
    result[LEVEL_A_FIELD] = peaks
    result["status"] = UPO_STATUS_OK if peaks else "no_peaks"
    return result


# -------------------------------------------------------------------------
# Level C -- PROJECT EXTENSION verification and stability
# -------------------------------------------------------------------------

def _local_affine_model(X, q, neighborhood, step):
    """Least-squares affine model of the first delay coordinate near q."""
    n, d = X.shape
    valid_n = n - step
    k = max(int(neighborhood), d + 2)
    if valid_n < k:
        return None
    dist = np.linalg.norm(X[:valid_n] - q, axis=1)
    idx = np.argsort(dist, kind="stable")[:k]
    A = np.column_stack([np.ones(k), X[idx] - q])
    y = X[idx + step, 0]
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ coef
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 if ss_tot <= 1e-15 else 1.0 - ss_res / ss_tot
    Fq = np.empty(d)
    Fq[0] = coef[0]
    Fq[1:] = q[:-1]
    return {"F_q": Fq, "jacobian": _companion(coef[1:]), "r2": float(r2)}


def estimate_jacobian_at_candidate(embedded, candidate, neighborhood=25, step=1):
    """PROJECT EXTENSION: companion Jacobian of a local fit centred on a candidate."""
    X = np.asarray(embedded, dtype=float)
    model = _local_affine_model(X, np.asarray(candidate, dtype=float), neighborhood, int(step))
    return None if model is None else model["jacobian"]


def monodromy_stability(embedded, orbit_points, neighborhood=30, step=1):
    """
    PROJECT EXTENSION (not the So et al. stability method): product of
    candidate-centred local Jacobians around the orbit, J(z_p)...J(z_1).
    For p = 1 this is candidate-centred fixed-point stability.
    """
    X = np.asarray(embedded, dtype=float)
    pts = np.atleast_2d(np.asarray(orbit_points, dtype=float))
    Mono = np.eye(X.shape[1])
    for q in pts:
        J = estimate_jacobian_at_candidate(X, q, neighborhood=neighborhood, step=step)
        if J is None or not np.all(np.isfinite(J)):
            return None
        Mono = J @ Mono
    eig = np.linalg.eigvals(Mono)
    mod = np.sort(np.abs(eig))[::-1]
    return {
        "provenance": "PROJECT EXTENSION (candidate-centred monodromy)",
        "monodromy": Mono,
        "multipliers": eig,
        "multiplier_moduli": mod,
        "unstable": bool(mod[0] > 1.0 + 1e-8),
        "saddle": bool(mod[0] > 1.0 + 1e-8 and mod[-1] < 1.0 - 1e-8),
    }


def fixed_point_stability(embedded, z_star, neighborhood=30, step=1):
    """PROJECT EXTENSION: candidate-centred period-1 stability."""
    X = np.asarray(embedded, dtype=float)
    if X.ndim != 2 or X.shape[1] < 2:
        return None
    return monodromy_stability(X, np.asarray(z_star, dtype=float)[None, :],
                               neighborhood=neighborhood, step=step)


def candidate_verification(embedded, peak, histogram, config=CFG, step=1):
    """
    PROJECT EXTENSION verification of one source peak:
      residual      max_k ||F_hat(z*(k)) - z*(k+1)|| / std(x), local affine F_hat
      r2            min_k R^2 of those local fits
      support_ratio peak count / median count of occupied histogram cells
    """
    X = np.asarray(embedded, dtype=float)
    pts = np.atleast_2d(peak["orbit_points"])
    p = len(pts)
    scale = max(float(np.std(X[:, 0])), 1e-12)
    residuals, r2s = [], []
    for k in range(p):
        model = _local_affine_model(X, pts[k], config.verify_neighbors, step)
        if model is None:
            return {"residual": np.nan, "r2": np.nan, "support_ratio": np.nan,
                    "passes": False, "failed_gates": ["model_unavailable"]}
        residuals.append(float(np.linalg.norm(model["F_q"] - pts[(k + 1) % p])) / scale)
        r2s.append(model["r2"])
    h = np.asarray(histogram, dtype=float)
    occupied = h[h > 0]
    support = float(peak["histogram_count"]) / max(float(np.median(occupied)), 1e-12) \
        if occupied.size else np.nan
    residual, r2 = max(residuals), min(r2s)
    failed = []
    if not residual <= config.verify_max_residual:
        failed.append("residual")
    if not r2 >= config.verify_min_r2:
        failed.append("r2")
    if not support >= config.verify_min_support_ratio:
        failed.append("support_ratio")
    return {"residual": residual, "r2": r2, "support_ratio": support,
            "passes": not failed, "failed_gates": failed,
            "provenance": "PROJECT EXTENSION"}


def apply_verification_gates(result, embedded, config=CFG):
    """
    Fill Level C.  Source peaks (Level A) are never removed or modified in
    meaning: failing peaks are listed in verification_failed_peaks and
    remain in source_peak_candidates.
    """
    passed, failed = [], []
    hist = result.get("histogram")
    for peak in result.get(LEVEL_A_FIELD, []):
        ver = candidate_verification(embedded, peak, hist, config=config, step=result["step"])
        stab = monodromy_stability(embedded, peak["orbit_points"],
                                   neighborhood=config.verify_neighbors, step=result["step"])
        entry = dict(peak)
        entry.update({"level": "C", "verification": ver, "extension_stability": stab,
                      "verified_unstable": bool(ver["passes"] and stab is not None
                                                and stab["unstable"])})
        (passed if ver["passes"] else failed).append(entry)
    result[LEVEL_C_PASS_FIELD] = passed
    result[LEVEL_C_FAIL_FIELD] = failed
    result["verification_assessed"] = True
    return result


# -------------------------------------------------------------------------
# Coverage -- PROJECT metric (not UPO density)
# -------------------------------------------------------------------------

def peak_coverage(locations, embedded, config=CFG, rng_seed=0):
    """
    PROJECT-derived coverage (NOT a density of periodic orbits):

      1. normalise the embedded attractor to [0, 1] per coordinate;
      2. sample n_ref = min(coverage_n_reference (=200), N) attractor points
         without replacement (rng_seed);
      3. radius = coverage_radius_multiplier x median distance to the
         coverage_radius_neighbors_k-th nearest neighbour within a sample
         of min(500, N) attractor points;
      4. coverage = fraction of the n_ref points within `radius` of at least
         one peak location (d-dimensional orbit points).

    No peaks -> 0.0.
    """
    X = np.asarray(embedded, dtype=float)
    locs = np.asarray(locations, dtype=float).reshape(-1, X.shape[1]) if len(X) else np.empty((0, 1))
    if len(locs) == 0 or len(X) == 0:
        return {"coverage": 0.0, "radius": np.nan, "n_reference": 0, "n_locations": 0}
    mins, maxs = X.min(axis=0), X.max(axis=0)
    span = np.maximum(maxs - mins, 1e-12)
    X_norm = (X - mins) / span
    locs_norm = (locs - mins) / span
    rng = np.random.default_rng(rng_seed)
    n_ref = min(config.coverage_n_reference, len(X_norm))
    reference = X_norm[rng.choice(len(X_norm), size=n_ref, replace=False)]
    sample = X_norm[rng.choice(len(X_norm), size=min(500, len(X_norm)), replace=False)]
    nn = []
    for i in range(len(sample)):
        dd = np.linalg.norm(sample - sample[i], axis=1)
        dd[i] = np.inf
        k = min(config.coverage_radius_neighbors_k, len(dd) - 1)
        nn.append(np.partition(dd, k)[k])
    spacing = float(np.median(nn))
    radius = config.coverage_radius_multiplier * spacing
    dmin = np.min(np.linalg.norm(reference[:, None, :] - locs_norm[None, :, :], axis=2), axis=1)
    return {"coverage": float(np.mean(dmin <= radius)), "radius": float(radius),
            "typical_spacing": spacing, "n_reference": int(n_ref),
            "n_locations": int(len(locs)),
            "definition": "PROJECT: fraction of sampled attractor points near a peak"}


def _population_points(peaks, d):
    if not peaks:
        return np.empty((0, d))
    return np.vstack([np.atleast_2d(p["orbit_points"]) for p in peaks])


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

def so_J(surrogate_W, w):
    """So et al. J(W'): fraction of surrogates whose maximum deviation W exceeds W'."""
    Ws = np.asarray(surrogate_W, dtype=float)
    return float(np.mean(Ws > w)) if Ws.size else np.nan


def so_surrogate_significance(observed_result, surrogate_results, alpha=0.05):
    """
    So et al. surrogate test [PRL; PRE Sec. IV].

    Every histogram is taken on the OBSERVED edges.  With rho_bar_sur the
    mean surrogate histogram:
        w(z)  = rho(z) - rho_bar_sur(z)          (signed; never abs())
        W     = max_z w(z)                       observed maximum deviation
        W_i   = max_z [rho_sur,i(z) - rho_bar_sur(z)]
        J(W') = fraction of W_i exceeding W'
        W0    = median(W_i), so J(W0) = 0.5
        r_J   = W / W0
    The finite-sample p-value (1 + #{W_i >= W}) / (n + 1) is reported
    separately.  A source peak is significant when J(deviation at its cell)
    < alpha (alpha is a PROJECT convention; the papers report J and r_J).
    """
    edges = observed_result.get("edges")
    obs = observed_result.get("histogram")
    if edges is None or obs is None:
        raise ValueError("observed result has no histogram; significance cannot be assessed")
    obs = np.asarray(obs, dtype=float)
    sur = np.asarray([_histogram_on_edges(np.asarray(r.get("reduced", np.empty((0, len(edges))))), edges)
                      for r in surrogate_results], dtype=float)
    if len(sur) < 2:
        raise ValueError("at least two surrogate histograms are required")
    mean = np.mean(sur, axis=0)
    flat = (sur - mean).reshape(len(sur), -1)
    Ws = np.max(flat, axis=1)
    deviation = obs - mean
    W = float(np.max(deviation))
    W0 = float(np.median(Ws))
    per_peak = []
    for peak in observed_result.get(LEVEL_A_FIELD, []):
        dev = float(deviation[peak["peak_cell"]])
        J_peak = so_J(Ws, dev)
        per_peak.append({"deviation": dev, "J": J_peak,
                         "p_value_finite": float((1 + np.sum(Ws >= dev)) / (len(Ws) + 1)),
                         "rJ": float(dev / W0) if W0 > 0 else np.nan,
                         "significant": bool(J_peak < alpha)})
    return {
        "provenance": "SOURCE (So et al. surrogate test)",
        "n_surrogates": int(len(sur)),
        "observed_histogram": obs,
        "surrogate_mean_histogram": mean,
        "deviation": deviation,
        "surrogate_W": Ws,
        "W": W,
        "W0": W0,
        "J_W": so_J(Ws, W),
        "J_W0": so_J(Ws, W0),
        "rJ": float(W / W0) if W0 > 0 else np.nan,
        "p_value_finite": float((1 + np.sum(Ws >= W)) / (len(Ws) + 1)),
        "alpha": float(alpha),
        "per_peak": per_peak,
    }


def assess_so_significance(result, x, config=CFG, rng=None, embedding_dimension=None):
    """
    Fill Level B for a detection result.  Gaussian-scaled phase-shuffle
    surrogates of the scalar series x are embedded with the observed lag
    and dimension and run through the identical detector on the observed
    histogram range.  On success significant_uop_candidates is a list
    (possibly empty); if the test cannot be carried out it stays None and
    significance_assessed stays False.
    """
    if result["status"] in UPO_FAILURE_STATUSES or result.get("edges") is None:
        result["significance"] = {"not_assessed_reason": f"detection status {result['status']}"}
        return result
    if rng is None:
        rng = np.random.default_rng(config.random_seed + 104729)
    x = _as_1d(x)
    tau = result["tau"]
    m = int(embedding_dimension or result["embedding_dimension"])
    period = result["period"]
    sur_results = []
    for _ in range(int(config.so_surrogate_count)):
        emb = _embed_backward(gaussian_scaled_phase_shuffle(x, rng=rng), tau, m)
        kw = dict(tau=tau, config=config, rng=rng, map_mode=result["map_mode"],
                  hist_range=result["histogram_range"])
        r = detect_so_fixed_points(emb, **kw) if period == 1 else detect_so_period_p(emb, period, **kw)
        if r["status"] not in UPO_FAILURE_STATUSES:
            sur_results.append(r)
    try:
        sig = so_surrogate_significance(result, sur_results, alpha=config.so_significance_alpha)
    except ValueError as exc:
        result["significance"] = {"not_assessed_reason": str(exc)}
        return result
    significant = []
    for peak, info in zip(result[LEVEL_A_FIELD], sig["per_peak"]):
        if info["significant"]:
            entry = dict(peak)
            entry.update({"level": "B", "significance": info})
            significant.append(entry)
    result["significance"] = sig
    result["significance_assessed"] = True
    result[LEVEL_B_FIELD] = significant
    return result


def run_surrogate_analysis(rr, tau, m, real_lle, lle_theiler, config=CFG):
    """Optional fixed-parameter AAFT surrogate test for the LLE.

    Surrogates keep the real series' marginal distribution and linear power
    spectrum while destroying its nonlinear temporal organization. Crucially,
    tau and m are held fixed at the values selected from the real series so
    surrogate re-optimization cannot manufacture a favorable null model.

    UPO surrogate significance is separate (so_assess_significance /
    assess_so_significance) and does not depend on this LLE analysis.
    """
    rng = np.random.default_rng(config.random_seed + 7919)
    lle_values=[]
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
    finite=np.asarray(lle_values,float); finite=finite[np.isfinite(finite)]
    if np.isfinite(real_lle) and len(finite):
        z=(real_lle-float(np.mean(finite)))/max(float(np.std(finite,ddof=1)),1e-12) if len(finite)>1 else np.nan
        p_lle=float((1+np.sum(finite>=real_lle))/(len(finite)+1))
    else:
        z=np.nan; p_lle=np.nan
    return {
        "n_surrogates": len(lle_values),
        "lle_surrogates": finite,
        "lle_z": float(z) if np.isfinite(z) else np.nan,
        "lle_p_upper": p_lle,
        "null_preserves_tau_m": True,
    }


def _aggregate_upo_status(statuses):
    for s in statuses:
        if s in UPO_FAILURE_STATUSES:
            return s
    if UPO_STATUS_OK in statuses:
        return UPO_STATUS_OK
    return statuses[0] if statuses else "no_peaks"


# Significance status of a UPO analysis.  Keeps "not requested", "requested
# but could not be assessed", "assessed" (Level B is then a list, possibly
# empty) and "detector failed" distinct.
SIGNIFICANCE_NOT_REQUESTED = "not_requested"
SIGNIFICANCE_ASSESSED = "assessed"
SIGNIFICANCE_NOT_ASSESSED = "significance_not_assessed"
SIGNIFICANCE_DETECTOR_FAILURE = "not_assessed_detector_failure"
SIGNIFICANCE_STATUSES = (SIGNIFICANCE_NOT_REQUESTED, SIGNIFICANCE_ASSESSED,
                         SIGNIFICANCE_NOT_ASSESSED, SIGNIFICANCE_DETECTOR_FAILURE)


def _initial_significance_status(config, status):
    if not config.so_assess_significance:
        return SIGNIFICANCE_NOT_REQUESTED
    if status in UPO_FAILURE_STATUSES:
        return SIGNIFICANCE_DETECTOR_FAILURE
    return SIGNIFICANCE_NOT_ASSESSED


def _empty_upo_analysis(map_info, config, status, tau, m, n_points=0, embedded=None):
    return {
        "map_mode": map_info["map_mode"] if map_info else config.upo_map_mode,
        "source_faithful_map": bool(map_info["source_faithful_map"]) if map_info else
                               config.upo_map_mode == "one_sample",
        "map_label": map_info["map_label"] if map_info else None,
        "tau": tau, "embedding_dimension": m, "n_embedded_points": int(n_points),
        "status": status,
        "significance_requested": bool(config.so_assess_significance),
        "significance_assessed": False,
        "significance_status": _initial_significance_status(config, status),
        "significance_not_assessed_reasons": {},
        "verification_assessed": False,
        "randomization": config.so_randomization,
        "periods": {},
        "levels": UPO_LEVELS,
        LEVEL_A_FIELD: [], LEVEL_B_FIELD: None,
        LEVEL_C_PASS_FIELD: None, LEVEL_C_FAIL_FIELD: None,
        "verified_unstable": None,
        "coverage": {"source": None, "significant": None, "verified_unstable": None},
        "source_rJ": np.nan,
        "embedded": embedded,
    }


def run_upo_analysis(x, config=CFG, lle_tau=None, lle_m=None, rng=None):
    """
    UPO analysis of a scalar series.

    one_sample (SOURCE, default): lag-1 embedding with the Cao dimension
    selected at lag 1; the detector receives lag-1 vectors.
    tau_step (PROJECT EXTENSION): reuses the TDMI tau and the LLE's Cao
    dimension / embedding.

    Returns status-labelled metadata plus the Level A/B/C populations of
    all analysed periods.  Failures are reported as statuses, not raised.
    """
    validate_upo_config(config)
    xs = np.asarray(x, dtype=float).reshape(-1)
    if not np.all(np.isfinite(xs)):
        return _empty_upo_analysis(None, config, "nonfinite_input", None, None)
    if xs.size == 0 or float(np.ptp(xs)) <= 1e-12:
        return _empty_upo_analysis(None, config, "constant_data", None, None)
    if config.upo_map_mode == "one_sample":
        tau_u = 1
        m_u, _, _ = cao_method(xs, 1, config.cao_max_dim, config.cao_tol, config.cao_theiler)
    else:
        if lle_tau is None:
            mi = time_delayed_mutual_information(xs, config.tdmi_max_tau, config.tdmi_bins)
            lle_tau = find_optimal_tau(mi, config.tdmi_smooth_window)
        tau_u = int(lle_tau)
        m_u = int(lle_m) if lle_m is not None else \
            cao_method(xs, tau_u, config.cao_max_dim, config.cao_tol, config.cao_theiler)[0]
    map_info = resolve_upo_map(config.upo_map_mode, tau_u)
    if m_u > config.cao_max_dim:
        return _empty_upo_analysis(map_info, config, "embedding_not_saturated", tau_u, int(m_u))
    emb = _embed_backward(xs, tau_u, m_u)
    if rng is None:
        rng = np.random.default_rng(config.random_seed)
    sig_rng = np.random.default_rng(config.random_seed + 104729)
    periods = {}
    for p in sorted(set(int(v) for v in config.so_periods)):
        if p == 1:
            r = detect_so_fixed_points(emb, tau=tau_u, config=config, rng=rng)
        else:
            r = detect_so_period_p(emb, p, tau=tau_u, config=config, rng=rng)
        if config.so_assess_significance:
            assess_so_significance(r, xs, config=config, rng=sig_rng, embedding_dimension=m_u)
        if config.so_verify_peaks and r["status"] not in UPO_FAILURE_STATUSES:
            apply_verification_gates(r, emb, config=config)
        periods[p] = r
    out = _empty_upo_analysis(map_info, config,
                              _aggregate_upo_status([r["status"] for r in periods.values()]),
                              tau_u, int(m_u), n_points=len(emb), embedded=emb)
    out["periods"] = periods
    d = emb.shape[1]
    source = [c for r in periods.values() for c in r[LEVEL_A_FIELD]]
    out[LEVEL_A_FIELD] = source
    if out["status"] in UPO_FAILURE_STATUSES:
        return out
    out["significance_assessed"] = bool(periods) and all(r["significance_assessed"] for r in periods.values())
    if config.so_assess_significance:
        # Requested but not computable (e.g. fewer than two usable surrogates)
        # is reported explicitly; Level B then stays None, never [].
        out["significance_status"] = (SIGNIFICANCE_ASSESSED if out["significance_assessed"]
                                      else SIGNIFICANCE_NOT_ASSESSED)
        out["significance_not_assessed_reasons"] = {
            p: (r["significance"] or {}).get("not_assessed_reason")
            for p, r in periods.items() if not r["significance_assessed"]}
    if out["significance_assessed"]:
        out[LEVEL_B_FIELD] = [c for r in periods.values() for c in r[LEVEL_B_FIELD]]
        p1 = periods.get(1)
        out["source_rJ"] = p1["significance"]["rJ"] if p1 is not None else np.nan
    out["verification_assessed"] = bool(periods) and all(r["verification_assessed"] for r in periods.values())
    if out["verification_assessed"]:
        out[LEVEL_C_PASS_FIELD] = [c for r in periods.values() for c in r[LEVEL_C_PASS_FIELD]]
        out[LEVEL_C_FAIL_FIELD] = [c for r in periods.values() for c in r[LEVEL_C_FAIL_FIELD]]
        out["verified_unstable"] = [c for c in out[LEVEL_C_PASS_FIELD] if c["verified_unstable"]]
    seed = config.random_seed
    out["coverage"] = {
        "source": peak_coverage(_population_points(source, d), emb, config, rng_seed=seed),
        "significant": (peak_coverage(_population_points(out[LEVEL_B_FIELD], d), emb, config, rng_seed=seed)
                        if out[LEVEL_B_FIELD] is not None else None),
        "verified_unstable": (peak_coverage(_population_points(out["verified_unstable"], d), emb, config,
                                            rng_seed=seed)
                              if out["verified_unstable"] is not None else None),
    }
    return out


def upo_summary_features(upo):
    """Per-population peak counts, per-point rates and coverages."""
    n = int(upo.get("n_embedded_points") or 0)
    status = upo.get("status")
    failed = status in UPO_FAILURE_STATUSES

    def population(peaks, coverage):
        if failed or peaks is None:
            return np.nan, np.nan, np.nan
        count = float(len(peaks))
        rate = count / n if n > 0 else np.nan
        cov = 0.0 if not peaks else (coverage or {}).get("coverage", np.nan)
        return count, rate, float(cov)

    cov = upo.get("coverage") or {}
    s = population(upo.get(LEVEL_A_FIELD), cov.get("source"))
    g = population(upo.get(LEVEL_B_FIELD), cov.get("significant"))
    v = population(upo.get("verified_unstable"), cov.get("verified_unstable"))
    rJ = upo.get("source_rJ", np.nan)
    return {
        "source_peak_count": s[0], "source_peaks_per_point": s[1], "source_peak_coverage": s[2],
        "significant_peak_count": g[0], "significant_peaks_per_point": g[1],
        "significant_peak_coverage": g[2],
        "source_rJ": float(rJ) if (not failed and upo.get("significance_assessed")
                                   and rJ is not None and np.isfinite(rJ)) else np.nan,
        "verified_unstable_count": v[0], "verified_unstable_per_point": v[1],
        "verified_unstable_coverage": v[2],
    }


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

    # UPO analysis: separate from the LLE.  In one_sample mode it uses its own
    # lag-1 embedding (Cao at lag 1); tau_step reuses tau, m and the embedding.
    try:
        upo = run_upo_analysis(rr_dynamics, config=config, lle_tau=tau, lle_m=m,
                               rng=np.random.default_rng(config.random_seed))
    except np.linalg.LinAlgError as exc:
        upo = _empty_upo_analysis(None, config, "analysis_error", None, None)
        upo["error"] = repr(exc)
    periodic = {p: r for p, r in upo["periods"].items() if p != 1}

    surrogate = None
    if config.compute_surrogates:
        surrogate = run_surrogate_analysis(rr_dynamics, tau, m, lle,
                                           lle_theiler, config=config)

    return {
        "rr_raw": rr, "rr_corrected": rr_corrected, "rr_dynamics": rr_dynamics,
        "used_corrected_rr_for_dynamics": bool(config.use_corrected_rr_for_dynamics),
        "rr_artifact_mask": artifact_mask,
        "stationarity": stationarity, "tau": tau, "embedding_dimension": m,
        "tdmi": mi, "cao_E1": e1, "cao_E2": e2, "embedded": embedded,
        "lle_per_beat": lle, "lle_theiler_beats": int(lle_theiler), "lle_diagnostics": lle_diag,
        "upo": upo,
        "so_fixed_points": upo["periods"].get(1), "so_periodic_orbits": periodic,
        "upo_map_mode": upo["map_mode"],
        "upo_source_faithful_map": upo["source_faithful_map"],
        "upo_tau": upo["tau"],
        "upo_embedding_dimension": upo["embedding_dimension"],
        "upo_n_embedded_points": upo["n_embedded_points"],
        "upo_status": upo["status"],
        "upo_significance_assessed": upo["significance_assessed"],
        "upo_significance_status": upo["significance_status"],
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
            upo_summary_features(r["upo"])["source_peak_coverage"] if "upo" in r else np.nan
            for r in results
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
        "source_peak_coverage": {
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

UPO_ROW_METADATA_COLUMNS = (
    "upo_feature_mode",
    "upo_map_mode",
    "upo_status",
    "upo_significance_assessed",
    "upo_significance_status",
)


def classifier_feature_columns(config=CFG):
    """
    Explicit classifier feature list for config.upo_feature_mode, taken
    from UPO_FEATURE_CONTRACT.  There is no implicit default feature list.
    """
    validate_upo_config(config)
    return list(UPO_FEATURE_CONTRACT[config.upo_feature_mode])


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


def extract_classifier_features(result, feature_mode):
    """
    Features of one analyze_segment() result for an explicit feature mode
    (UPO_FEATURE_CONTRACT).  Uses only the analysed RR window: no labels,
    abnormal fraction, record identity, other windows or model outputs.
    No reconstruction parameter (tau or m) is used as a disease feature.

    No peaks -> counts, rates and coverages are 0.
    UPO pipeline failure (UPO_FAILURE_STATUSES) -> UPO features are NaN.
    Significance requested but not assessable (significance_status ==
    "significance_not_assessed") -> significance-dependent features are NaN
    (never 0); source-peak features are unaffected.

    Using significant_source on a result for which significance was never
    requested is a configuration error and raises ValueError.
    """
    if feature_mode not in UPO_FEATURE_CONTRACT:
        raise ValueError(f"feature_mode must be one of {tuple(UPO_FEATURE_CONTRACT)}, "
                         f"got {feature_mode!r}")
    upo = result["upo"]
    if feature_mode == "significant_source" and not upo.get("significance_requested", False) \
            and upo["status"] not in UPO_FAILURE_STATUSES:
        raise ValueError("feature_mode='significant_source' requires UPO surrogate "
                         "significance to have been requested (so_assess_significance=True)")
    lle = result.get("lle_per_beat", np.nan)
    values = upo_summary_features(upo)
    values["lle_per_beat"] = float(lle) if lle is not None and np.isfinite(lle) else np.nan
    return {col: values[col] for col in UPO_FEATURE_CONTRACT[feature_mode]}


def upo_row_metadata(result, feature_mode):
    """Per-row UPO provenance recorded next to the features."""
    upo = result["upo"]
    return {
        "upo_feature_mode": feature_mode,
        "upo_map_mode": upo["map_mode"],
        "upo_status": upo["status"],
        "upo_significance_assessed": bool(upo["significance_assessed"]),
        "upo_significance_status": upo["significance_status"],
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
    feature_mode = config.upo_feature_mode
    feature_columns = classifier_feature_columns(config)

    for window_number, window in enumerate(windows):
        # Features are computed from the RR window alone; the label and
        # abnormal fraction are attached only after feature extraction.
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
                **{col: np.nan for col in feature_columns},
                "upo_feature_mode": feature_mode,
                "upo_map_mode": config.upo_map_mode,
                "upo_status": "analysis_error",
                "upo_significance_assessed": False,
                "upo_significance_status": (SIGNIFICANCE_DETECTOR_FAILURE
                                            if config.so_assess_significance
                                            else SIGNIFICANCE_NOT_REQUESTED),
                "analysis_error": str(exc),
            })
            continue

        features = extract_classifier_features(result, feature_mode)

        rows.append({
            "record": str(record_name),
            "window": int(window_number),
            "label": int(window["label"]),
            "abnormal_fraction": float(window["abnormal_fraction"]),
            **features,
            **upo_row_metadata(result, feature_mode),
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
        for col in classifier_feature_columns(config) + list(UPO_ROW_METADATA_COLUMNS) + [
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
    feature_columns,
    n_splits=5,
    random_state=42,
):
    """
    Record-grouped logistic regression.

    feature_columns is required and explicit (see classifier_feature_columns);
    there is no hard-coded default feature list.

    Imputation and scaling are fitted inside each training fold, preventing
    information leakage from the held-out records.
    """
    feature_columns = list(feature_columns)
    if not feature_columns:
        raise ValueError("feature_columns must be a non-empty explicit list")
    leaking = {"label", "record", "abnormal_fraction", "window"} & set(feature_columns)
    if leaking:
        raise ValueError(f"feature_columns must not include {sorted(leaking)}")
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


def nonlinear_feature_statistics(df, feature_columns):
    """
    Descriptive/Mann-Whitney statistics for the explicit feature columns.

    Each feature uses only its own finite values.
    """
    from scipy.stats import mannwhitneyu

    output = {}

    for feature in feature_columns:
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

    feature_columns = classifier_feature_columns(config)

    stats = nonlinear_feature_statistics(df, feature_columns)

    model = run_grouped_logistic_regression(
        df,
        feature_columns=feature_columns,
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
        "feature_mode": config.upo_feature_mode,
        "feature_columns": feature_columns,
        "upo_map_mode": config.upo_map_mode,
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
    print(f"UPO feature mode: {base_cfg.upo_feature_mode}")
    print(f"UPO map mode: {base_cfg.upo_map_mode}")
    print(
        "Features: "
        + ", ".join(classifier_feature_columns(base_cfg))
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

        coverage_col = [c for c in result["feature_columns"] if c.endswith("_coverage")][0]
        comparison_rows.append({
            "configuration": name,
            "upo_feature_mode": result["feature_mode"],
            "upo_map_mode": result["upo_map_mode"],
            "lle_fit_success_normal":
                result["statistics"]["lle_fit_success"]["normal_rate"],
            "lle_fit_success_abnormal":
                result["statistics"]["lle_fit_success"]["abnormal_rate"],
            "lle_mann_whitney_p":
                result["statistics"]["lle_per_beat"]["mann_whitney_p"],
            f"{coverage_col}_mann_whitney_p":
                result["statistics"][coverage_col]["mann_whitney_p"],
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
