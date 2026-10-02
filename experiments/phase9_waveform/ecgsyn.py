"""
Phase 9 Part A1: the McSharry et al. (2003) dynamical ECG model (ECGSYN), driven by an
EXTERNAL beat-time series.

Source (primary, open): McSharry PE, Clifford GD, Tarassenko L, Smith LA, "A dynamical model
for generating synthetic electrocardiogram signals", IEEE Trans Biomed Eng 50(3):289-294
(2003); the authors' HTML version is on PhysioNet (physionet.org/content/ecgsyn/1.0.0/paper/),
Eqs. (1)-(4) and Table I.  Reference implementation: the authors' ecgsyn.c (PhysioNet
ECGSYN 1.0.0), used for the heart-rate adjustment of the extrema parameters and for the
numerical verification (verify_ecgsyn.py).

Model (paper Eq. 1, 2, 4):
    xdot = alpha x - omega y,  ydot = alpha y + omega x,  alpha = 1 - sqrt(x^2 + y^2)
    zdot = - sum_{i in P,Q,R,S,T} a_i dtheta_i exp(-dtheta_i^2 / (2 b_i^2)) - (z - z0(t))
    dtheta_i = theta - theta_i,  theta = atan2(y, x),  z0(t) = A sin(2 pi f2 t)
    omega(t) = 2 pi / T(t)  (T(t) = the current RR interval, piecewise constant, as ecgsyn.c's rrpc)
On the limit cycle (start (x, y) = (1, 0)) theta advances exactly as the integral of omega,
so theta is evaluated analytically: with input beat times t_k, theta(t) = 2 pi (t - t_k) /
(t_{k+1} - t_k) on [t_k, t_{k+1}), wrapped to (-pi, pi].  The R event (theta_R = 0) therefore
falls exactly on each input beat time, and only z is integrated (fourth-order Runge-Kutta,
fixed step 1/fs_int, as the paper and ecgsyn.c).

Per-beat parameters (the only extension, needed for ectopic and repolarization
morphology): beat k owns the half-cycle before t_k (its P, Q and the rise of R) and the
half-cycle after it (fall of R, S, T): on [t_k, t_{k+1}) the parameters of beat k are used for
theta >= 0 and those of beat k+1 for theta < 0.  With identical parameters for every beat this
is exactly Eq. (1).

Units: z is in model units; to_mV() applies ONE fixed affine map (FIXED_SCALE, FIXED_OFFSET),
determined once from the reference normal beat at 60 bpm so that its range is
[-0.4, 1.2] mV as in ecgsyn.c.  (ecgsyn.c rescales each record by its own min/max; a fixed map
keeps amplitudes comparable across windows and morphologies.)
"""
from __future__ import annotations

import math

import numpy as np
from numba import njit

# Table I (paper) == ecgsyn.c defaults:          P       Q      R     S      T
THETA_DEG = np.array([-60.0, -15.0, 0.0, 15.0, 90.0])
A_I = np.array([1.2, -5.0, 30.0, -7.5, 0.75])
B_I = np.array([0.25, 0.1, 0.1, 0.1, 0.4])
A_RESP = 0.005          # ecgsyn.c zbase amplitude (model units) == 0.15 mV after scaling (paper Eq. 2)
F_RESP = 0.25           # Hz (paper f2, ecgsyn.c fhi)
Z_INIT = 0.04           # ecgsyn.c zinitial


def hr_adjusted(hrmean):
    """Extrema parameters adjusted for mean heart rate exactly as ecgsyn.c (dorun()):
    b_i *= sqrt(hr/60); theta_P *= (hr/60)^(1/4); theta_Q, theta_S *= sqrt(hr/60);
    theta_R, theta_T unchanged.  Returns (theta [rad], a, b), each shape (5,)."""
    hrfact = math.sqrt(hrmean / 60.0)
    hrfact2 = math.sqrt(hrfact)
    th = THETA_DEG * math.pi / 180.0
    th = th * np.array([hrfact2, hrfact, 1.0, hrfact, 1.0])
    return th, A_I.copy(), B_I * hrfact


@njit(cache=True)
def _zdot(t, z, beat_t, th, a, b, A, f, ph0):
    # locate interval k with beat_t[k] <= t < beat_t[k+1] (binary search)
    lo, hi = 0, beat_t.shape[0] - 1
    if t <= beat_t[0]:
        k = 0
    elif t >= beat_t[hi]:
        k = hi - 1
    else:
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if beat_t[mid] <= t:
                lo = mid
            else:
                hi = mid
        k = lo
    frac = (t - beat_t[k]) / (beat_t[k + 1] - beat_t[k])
    if frac < 0.5:
        theta = 2.0 * math.pi * frac
        j = k
    else:
        theta = 2.0 * math.pi * (frac - 1.0)
        j = k + 1
    s = 0.0
    for i in range(5):
        dth = theta - th[j, i]
        bi = b[j, i]
        s += a[j, i] * dth * math.exp(-0.5 * dth * dth / (bi * bi))
    z0 = A * math.sin(2.0 * math.pi * f * t + ph0)
    return -s - (z - z0)


@njit(cache=True)
def _integrate(beat_t, th, a, b, t0, n_out, q, h, A, f, ph0, zinit):
    out = np.empty(n_out)
    z = zinit
    t = t0
    for n in range(n_out):
        out[n] = z
        for _ in range(q):
            k1 = _zdot(t, z, beat_t, th, a, b, A, f, ph0)
            k2 = _zdot(t + 0.5 * h, z + 0.5 * h * k1, beat_t, th, a, b, A, f, ph0)
            k3 = _zdot(t + 0.5 * h, z + 0.5 * h * k2, beat_t, th, a, b, A, f, ph0)
            k4 = _zdot(t + h, z + h * k3, beat_t, th, a, b, A, f, ph0)
            z = z + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
            t = t + h
    return out


def synthesize(beat_t, theta, a, b, fs=360, q=2, t0=None, t1=None, A=A_RESP, f=F_RESP,
               resp_phase=0.0, zinit=Z_INIT):
    """z(t) sampled at fs on [t0, t1) (default: first to last beat time).
    beat_t: (K,) increasing beat (R-event) times [s]; theta, a, b: (K, 5) per-beat extrema
    parameters (rows for beats). fs_int = q * fs (ecgsyn.c integrates at sf and keeps every
    q-th sample)."""
    beat_t = np.ascontiguousarray(beat_t, dtype=float)
    th = np.ascontiguousarray(np.broadcast_to(theta, (len(beat_t), 5)), dtype=float)
    aa = np.ascontiguousarray(np.broadcast_to(a, (len(beat_t), 5)), dtype=float)
    bb = np.ascontiguousarray(np.broadcast_to(b, (len(beat_t), 5)), dtype=float)
    t0 = float(beat_t[0]) if t0 is None else float(t0)
    t1 = float(beat_t[-1]) if t1 is None else float(t1)
    n_out = int(math.floor((t1 - t0) * fs + 1e-9))
    return _integrate(beat_t, th, aa, bb, t0, n_out, int(q), 1.0 / (q * fs), float(A), float(f),
                      float(resp_phase), float(zinit))


# ------------------------------------------------------------------ fixed amplitude scaling
def _reference_scale():
    th, a, b = hr_adjusted(60.0)
    bt = np.arange(0.0, 41.0)
    z = synthesize(bt, th, a, b, fs=360, q=2, A=0.0)
    seg = z[int(20 * 360):int(40 * 360)]          # after the initial transient
    zmin, zmax = float(seg.min()), float(seg.max())
    scale = 1.6 / (zmax - zmin)
    return scale, -0.4 - zmin * scale


FIXED_SCALE, FIXED_OFFSET = None, None


def to_mV(z):
    global FIXED_SCALE, FIXED_OFFSET
    if FIXED_SCALE is None:
        FIXED_SCALE, FIXED_OFFSET = _reference_scale()
    return np.asarray(z) * FIXED_SCALE + FIXED_OFFSET
