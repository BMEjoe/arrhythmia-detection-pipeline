"""Deterministic synthetic ECG generator for Pan-Tompkins regression tests."""
import numpy as np


def _g(t, mu, sig, amp):
    return amp * np.exp(-0.5 * ((t - mu) / sig) ** 2)


def make_ecg(fs, beat_times, *, qrs_amp=1.0, qrs_sigma=0.010, t_amp=0.25,
             t_delay=0.25, t_sigma=0.045, amps=None, t_amps=None,
             sigmas=None, duration=None, noise=0.0, wander=0.0,
             offset=0.0, seed=0):
    """Sum of Gaussian R (+ small Q/S) and T waves at ``beat_times`` (s).

    Returns (ecg, true_r_sample_indices).
    """
    beat_times = np.asarray(beat_times, dtype=float)
    dur = duration if duration is not None else beat_times[-1] + 1.0
    n = int(round(dur * fs))
    t = np.arange(n) / fs
    x = np.zeros(n)
    for i, b in enumerate(beat_times):
        a = qrs_amp * (amps[i] if amps is not None else 1.0)
        ta = t_amp * (t_amps[i] if t_amps is not None else 1.0)
        s = qrs_sigma * (sigmas[i] if sigmas is not None else 1.0)
        x += _g(t, b, s, a)
        x += _g(t, b - 2.5 * s, s, -0.10 * a) + _g(t, b + 2.5 * s, s, -0.10 * a)
        x += _g(t, b + t_delay, t_sigma, ta)
    rng = np.random.default_rng(seed)
    if noise:
        x += noise * rng.standard_normal(n)
    if wander:
        x += wander * np.sin(2 * np.pi * 0.25 * t)
    x += offset
    return x, np.rint(beat_times * fs).astype(int)


def score(detected, truth, fs, tol_s=0.075):
    """Greedy one-to-one matching. Returns dict of TP/FP/FN/Se/PPV/F1/timing."""
    detected = np.asarray(detected, dtype=int)
    truth = np.asarray(truth, dtype=int)
    tol = tol_s * fs
    used = set()
    errs = []
    tp = 0
    for tr in truth:
        best, bd = None, None
        for j, d in enumerate(detected):
            if j in used:
                continue
            e = abs(d - tr)
            if e <= tol and (bd is None or e < bd):
                best, bd = j, e
        if best is not None:
            used.add(best); tp += 1; errs.append((detected[best] - tr) / fs)
    fp = len(detected) - tp
    fn = len(truth) - tp
    se = tp / max(1, tp + fn)
    ppv = tp / max(1, tp + fp)
    f1 = 2 * se * ppv / max(1e-12, se + ppv)
    return dict(TP=tp, FP=fp, FN=fn, Se=se, PPV=ppv, F1=f1,
                mean_err_ms=1000 * float(np.mean(errs)) if errs else float("nan"),
                max_abs_err_ms=1000 * float(np.max(np.abs(errs))) if errs else float("nan"))
