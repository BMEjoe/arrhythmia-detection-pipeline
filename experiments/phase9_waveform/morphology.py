"""
Phase 9 Part A1/A2: per-beat ECGSYN extrema parameters (predeclared; fixed before any detector
or measure was run on a Phase 9 series).

Normal beats (type 'N'): paper Table I / ecgsyn.c defaults, adjusted for the window's mean
heart rate exactly as ecgsyn.c (ecgsyn.hr_adjusted).

Ventricular ectopic beats (type 'V'; used for Phase 5 S2, Phase 6 E1-E4, E6 and the Phase 8
ectopy variants S2_5, E1, E3_10, all of which are ventricular-type with a full compensatory
pause).  Basis: LITFL "Premature Ventricular Complex" (litfl.com/premature-ventricular-
complex-pvc-ecg-library, accessed 2026-10-02): "Broad QRS complex (>= 120 ms) with abnormal
morphology", "Discordant ST segment and T wave changes", "Usually followed by a full
compensatory pause"; no conducted sinus P wave precedes the complex.  Implementation:
  - P event removed (a_P = 0);
  - QRS widened by W_QRS = 1.5: b_Q, b_R, b_S and the angles theta_Q, theta_S multiplied by 1.5
    (the model QRS, Q event - 2 sigma to S event + 2 sigma, is ~131 ms at 75 bpm for a normal
    beat and ~197 ms for this PVC);
  - QRS amplitudes a_Q, a_R, a_S multiplied by A_QRS = 0.68, which gives an R deflection about
    1.5 x the normal one (calibrated numerically at 75 bpm: 1.53 x);
  - T discordant (inverted): a_T -> -K_T a_T with K_T = 0.50 and b_T x W_T = 1.5, which gives an
    inverted T of the same magnitude as the normal upright T (calibrated at 75 bpm: 1.00 x).
Atrial ectopic beats (type 'A'; Phase 6 E5).  Basis: LITFL "Premature Atrial Complex"
(litfl.com/premature-atrial-complex-pac): "Abnormal (non-sinus) P wave usually followed by a
normal QRS complex (< 120 ms)"; PACs from low atrial areas produce "an inverted P wave".
Implementation: a_P -> -a_P (inverted P); QRS and T normal.

Repolarization (APD) morphology (A2, KTz family): beat k's T event is moved by the APD
deviation and its amplitude scaled by the APD ratio:
  theta_T,k = theta_T + 2 pi (APD_k - APD_ref) / RR_{k+1}  (the T event lies in the interval
             after R_k, whose angular speed is 2 pi / RR_{k+1}), so the T peak (and hence QT)
             is shifted in time by exactly APD_k - APD_ref;
  a_T,k     = a_T * APD_k / APD_ref  (T amplitude proportional to APD: T-wave alternans is the
             ECG manifestation of APD alternans, Qu et al. Phys Rep 2014, PMC4175480, Sect. 1).
APD_ref is the regime's long-run mean APD (ground_truth_ktz.py).
"""
from __future__ import annotations

import numpy as np

from experiments.phase9_waveform import ecgsyn as E

W_QRS = 1.5
A_QRS = 0.68
K_T = 0.50
W_T = 1.5
P, Q, R, S, T = range(5)


def beat_params(rr, types, apd_s=None, apd_ref_s=None):
    """Per-beat (theta, a, b), each (K, 5), for K = len(rr) + 1 beats; beat k ends interval
    k - 1 (beat 0 starts the window).  types: length K array of 'N', 'V', 'A'.
    apd_s / apd_ref_s (seconds): optional APD per beat (length K) for the T-wave mapping."""
    rr = np.asarray(rr, float)
    K = len(rr) + 1
    hr = 60.0 / float(np.mean(rr))
    th0, a0, b0 = E.hr_adjusted(hr)
    th = np.tile(th0, (K, 1))
    a = np.tile(a0, (K, 1))
    b = np.tile(b0, (K, 1))
    types = np.asarray(types)
    v = types == "V"
    th[v, Q] *= W_QRS
    th[v, S] *= W_QRS
    b[v, Q:S + 1] *= W_QRS
    a[v, Q:S + 1] *= A_QRS
    a[v, P] = 0.0
    a[v, T] = -K_T * a0[T]
    b[v, T] *= W_T
    at = types == "A"
    a[at, P] = -a0[P]
    if apd_s is not None:
        apd_s = np.asarray(apd_s, float)
        rr_next = np.concatenate([rr, [rr[-1]]])           # interval after R_k
        th[:, T] = th[:, T] + 2.0 * np.pi * (apd_s - apd_ref_s) / rr_next
        a[:, T] = a[:, T] * (apd_s / apd_ref_s)
    return th, a, b
