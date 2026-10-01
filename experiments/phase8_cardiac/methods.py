"""
Phase 8 detector methods (Part B).  Each method maps an RR window to a JSON-able record
with a boolean "detected" (and diagnostics).  Errors count as not detected.

  frozen      the Phase 6 recommended detector (Phase 7 DETECTOR), unchanged: AND of
              lle_chaos_test and the instability-gated UPO test; also its components.
Candidates are added below as they are designed (development only).
"""
from __future__ import annotations

import traceback

import numpy as np


def _frozen(rr):
    from experiments.phase7_mitbih import detector as DT
    r = DT.evaluate(rr)
    keep = ("and_detected", "lle_detected", "upo_detected", "lle_p", "lle_z", "lle_stat", "upo_score",
            "upo_status", "error", "detrend_applied")
    out = {k: r.get(k) for k in keep}
    out["detected"] = bool(r["and_detected"])
    return out


METHODS = {"frozen": _frozen}


def run(name, rr):
    try:
        return METHODS[name](np.asarray(rr, dtype=float))
    except Exception as exc:                       # noqa: BLE001
        return {"detected": False, "error": f"{type(exc).__name__}: {exc}", "traceback": traceback.format_exc(limit=2)}


# ----------------------------------------------------------------------------- development statistics
def embed(x, m, tau=1):
    n = len(x) - (m - 1) * tau
    return np.stack([x[i * tau: i * tau + n] for i in range(m)][::-1], axis=1)   # row t: (x_t, x_{t-1}, ...)


def nl_pred_error(x, m=3, k=5, theiler=5, horizon=1, tau=1, robust=True):
    """Leave-one-out local-average prediction of x[t+h] from the k nearest neighbours of the
    delay vector at t (Theiler window excluded).  Returns the normalized error:
    robust: median |e| / median |x - median(x)|; else RMSE / SD."""
    x = np.asarray(x, dtype=float)
    E = embed(x, m, tau)                         # E[i] corresponds to time t_i = i + (m-1)tau
    t0 = (m - 1) * tau
    nvec = len(E) - horizon
    A = E[:nvec]
    target = x[t0 + horizon: t0 + horizon + nvec]
    d = ((A[:, None, :] - A[None, :, :]) ** 2).sum(-1)
    idx = np.arange(nvec)
    d[np.abs(idx[:, None] - idx[None, :]) <= theiler] = np.inf
    nn = np.argpartition(d, k, axis=1)[:, :k]
    pred = target[nn].mean(1)
    e = target - pred
    if robust:
        return float(np.median(np.abs(e)) / max(np.median(np.abs(x - np.median(x))), 1e-12))
    return float(np.sqrt(np.mean(e ** 2)) / max(np.std(x), 1e-12))


def surrogate_test(stat_fn, x, n_sur, rng, smaller_is_signal=True):
    """One-sided rank test against IAAFT surrogates (fp.iaaft_surrogate)."""
    import final_pipeline as fp
    obs = stat_fn(x)
    sur = np.array([stat_fn(fp.iaaft_surrogate(x, rng)) for _ in range(n_sur)])
    if smaller_is_signal:
        exceed = np.sum(sur <= obs)
    else:
        exceed = np.sum(sur >= obs)
    return obs, sur, float((1 + exceed) / (n_sur + 1))


def _rng_for(x, salt):
    import zlib
    h = zlib.crc32(np.ascontiguousarray(np.asarray(x, dtype=float)).tobytes())
    return np.random.default_rng(np.random.SeedSequence([20261004, int(h), int(salt)]))


def _nlpred_scan(rr):
    """Diagnostic (development only): nonlinear-prediction p-values for m = 2..5 with 39
    IAAFT surrogates, robust and RMSE statistics."""
    out = {}
    for m in (2, 3, 4, 5):
        rng = _rng_for(rr, m)
        obs, sur, p = surrogate_test(lambda z: nl_pred_error(z, m=m), rr, 39, rng)
        out[f"np_m{m}_p"] = p
        out[f"np_m{m}_z"] = float((np.mean(sur) - obs) / max(np.std(sur, ddof=1), 1e-12))
    out["detected"] = bool(min(out[f"np_m{m}_p"] for m in (2, 3, 4, 5)) <= 0.0125)
    return out


METHODS["nlpred_scan"] = _nlpred_scan


# ----------------------------------------------------------------------------- noise titration
def _vwk_design(y, kappa, d):
    """Volterra-Wiener-Korenberg (polynomial autoregressive) design matrix: constant and all
    monomials of degree 1..d of (y[n-1], ..., y[n-kappa]); targets y[n], n = kappa..N-1."""
    from itertools import combinations_with_replacement
    N = len(y)
    lags = np.stack([y[kappa - j - 1: N - j - 1] for j in range(kappa)], axis=1)
    cols = [np.ones(N - kappa)]
    for deg in range(1, d + 1):
        for comb in combinations_with_replacement(range(kappa), deg):
            cols.append(np.prod(lags[:, comb], axis=1))
    return np.stack(cols, axis=1), y[kappa:]


def vwk_nonlinearity(y, kmax=6, dmax=3, alpha=0.01):
    """Barahona-Poon style test (as described in Poon & Barahona PNAS 2001 and Wysocki et al.,
    arXiv:nlin/0606032): best linear model (d = 1, kappa = 1..kmax) and best nonlinear model
    (d = 2..dmax, kappa = 1..kmax) by C(r) = log(eps(r)) + r/N, eps = normalized residual
    variance, r = number of terms; nonlinearity iff C_nl < C_lin and the nested F-test of the
    best nonlinear model against the linear model with the same kappa (its linear part) has
    p < alpha.  Returns (detected, info)."""
    from scipy import stats
    y = (np.asarray(y, float) - np.mean(y)) / max(np.std(y), 1e-12)
    best = {}
    for kappa in range(1, kmax + 1):
        for d in range(1, dmax + 1):
            X, t = _vwk_design(y, kappa, d)
            n_eff = len(t)
            if X.shape[1] >= n_eff - 5:
                continue
            coef, *_ = np.linalg.lstsq(X, t, rcond=None)
            rss = float(np.sum((t - X @ coef) ** 2))
            eps = rss / max(float(np.sum((t - t.mean()) ** 2)), 1e-300)
            C = np.log(max(eps, 1e-300)) + X.shape[1] / n_eff
            best[(kappa, d)] = (C, rss, X.shape[1], n_eff)
    lin = min((v[0], k) for k, v in best.items() if k[1] == 1)
    nl = min((v[0], k) for k, v in best.items() if k[1] > 1)
    kn = nl[1]
    C_nl, rss_nl, r_nl, n_eff = best[kn]
    C_l, rss_l, r_l, _ = best[(kn[0], 1)]
    F = ((rss_l - rss_nl) / max(r_nl - r_l, 1)) / max(rss_nl / max(n_eff - r_nl, 1), 1e-300)
    p = float(stats.f.sf(F, max(r_nl - r_l, 1), max(n_eff - r_nl, 1)))
    det = bool(nl[0] < lin[0] and p < alpha)
    return det, {"C_lin": float(lin[0]), "C_nl": float(nl[0]), "k_lin": lin[1][0], "k_nl": kn[0], "d_nl": kn[1],
                 "F_p": p}


def noise_titration(x, rng, max_level=2.0, n_bisect=6):
    """Noise limit NL: the largest added-white-noise SD (fraction of the series SD) at which
    vwk_nonlinearity still detects nonlinearity (bisection on [0, max_level]); NL = 0 when
    the raw series is not detected.  One noise realization per level (scaled draws of a
    single Gaussian sequence, so detection is evaluated on a nested noise family)."""
    x = np.asarray(x, float)
    det0, info = vwk_nonlinearity(x)
    if not det0:
        return 0.0, info
    z = rng.standard_normal(len(x)) * np.std(x)
    lo, hi = 0.0, max_level
    if vwk_nonlinearity(x + hi * z)[0]:
        return hi, info
    for _ in range(n_bisect):
        mid = 0.5 * (lo + hi)
        if vwk_nonlinearity(x + mid * z)[0]:
            lo = mid
        else:
            hi = mid
    return lo, info


def _titration(rr):
    nl, info = noise_titration(rr, _rng_for(rr, 77))
    return {"detected": bool(nl > 0), "NL": float(nl), **info}


METHODS["titration"] = _titration


# ----------------------------------------------------------------------------- ectopy-preserving surrogates
def outlier_mask(x):
    """Intervals flagged by the pipeline's own robust local outlier rule
    (fp.correct_rr_intervals with CFG: hard limits and > 4.5 robust SD from an 11-beat
    running median).  No annotations are used."""
    import final_pipeline as fp
    _, mask, _ = fp.correct_rr_intervals(np.asarray(x, float), fp.CFG)
    return np.asarray(mask, bool)


def ep_iaaft(x, rng, mask=None):
    """Ectopy-preserving IAAFT surrogate: IAAFT of the series with flagged intervals
    replaced by linear interpolation (fp.correct_rr_intervals), then the flagged original
    intervals re-inserted at their original positions, shifted by the local level change
    (surrogate minus cleaned original at that index).  Null hypothesis: a monotone transform
    of a linear Gaussian process plus the observed ectopic/outlier events at the observed
    times."""
    import final_pipeline as fp
    x = np.asarray(x, float)
    if mask is None:
        mask = outlier_mask(x)
    clean, _, _ = fp.correct_rr_intervals(x, fp.CFG)
    s = fp.iaaft_surrogate(clean, rng)
    if mask.any():
        s = s.copy()
        s[mask] = x[mask] + (s[mask] - clean[mask])
    return s


def surrogate_test_gen(stat_fn, x, n_sur, rng, gen, smaller_is_signal=True):
    obs = stat_fn(x)
    sur = np.array([stat_fn(gen(x, rng)) for _ in range(n_sur)])
    exceed = np.sum(sur <= obs) if smaller_is_signal else np.sum(sur >= obs)
    return obs, sur, float((1 + exceed) / (n_sur + 1))


def _nlpred_ep_scan(rr):
    """Diagnostic: robust nonlinear-prediction error vs ectopy-preserving IAAFT, m = 2..5."""
    mask = outlier_mask(rr)
    out = {"n_flagged": int(mask.sum())}
    for m in (2, 3, 4, 5):
        rng = _rng_for(rr, 100 + m)
        obs, sur, p = surrogate_test_gen(lambda z: nl_pred_error(z, m=m), rr, 39, rng,
                                         lambda z, r: ep_iaaft(z, r, mask))
        out[f"np_m{m}_p"] = p
        out[f"np_m{m}_z"] = float((np.mean(sur) - obs) / max(np.std(sur, ddof=1), 1e-12))
    out["detected"] = bool(max(out[f"np_m{m}_z"] for m in (2, 3, 4, 5)) > 3.0)
    return out


METHODS["nlpred_ep_scan"] = _nlpred_ep_scan
