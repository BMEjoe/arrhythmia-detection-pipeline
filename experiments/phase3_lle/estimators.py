"""
Phase 3 LLE estimators and tests.

Every method maps one scalar window x to a result dict with at least
    lle        LLE estimate, natural log per sample (NaN if undefined)
    p          one-sided surrogate p-value (small = more divergence than the null)
    detected   bool, the method's chaos decision
    status     "ok" or a failure reason (failures count as "not detected")
Pipeline functions are imported from final_pipeline, never copied.
"""
from __future__ import annotations

import numpy as np

import final_pipeline as fp

CFG = fp.CFG
ALPHA = 0.05


def _f(v):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    return v if np.isfinite(v) else None


# -----------------------------------------------------------------------------
# B1 baseline: existing Rosenstein path + existing AAFT test, unchanged
# -----------------------------------------------------------------------------

def baseline(x, config=CFG, with_curve=False):
    """
    Exactly the LLE part of final_pipeline.analyze_segment with
    compute_surrogates=True (default CFG: TDMI tau, Cao m at tau, Theiler =
    mean period, Rosenstein early fit, 199 AAFT surrogates at fixed tau/m with
    rng = random_seed + 7919).  Decision: lle_p_upper <= 0.05, whose nominal size
    under exchangeability is 10/200 = 0.05.  The UPO part of analyze_segment is
    not run (it does not affect the LLE).
    """
    x = np.asarray(x, dtype=float)
    emb, tau, m, *_ = fp.takens_embed(x, config=config)
    out = {"tau": int(tau), "m": int(m), "n_embedded": int(len(emb))}
    if len(emb) < 50:          # analyze_segment raises here
        out.update(lle=None, p=None, detected=False, status="insufficient embedded points (analyze_segment raises)")
        return out
    th = config.lle_theiler_beats if config.lle_theiler_beats is not None else fp._mean_period_beats(x)
    lle, diag = fp.rosenstein_lle(emb, theiler=th, max_iter=config.lle_max_iter,
                                  min_fit_points=config.lle_min_fit_points,
                                  max_fit_fraction=config.lle_max_fit_fraction,
                                  min_r2=config.lle_min_r2, require_valid_fit=config.lle_require_valid_fit)
    sur = fp.run_surrogate_analysis(x, tau, m, lle, th, config=config)
    p = _f(sur["lle_p_upper"])
    fi = diag.get("fit_indices")
    out.update({
        "lle": _f(lle), "r2": _f(diag.get("r2")), "theiler": int(th),
        "lle_reason": diag.get("reason"),
        "fit_k": None if fi is None else [int(fi[0]), int(fi[-1])],
        "p": p, "z": _f(sur["lle_z"]), "n_surrogates": int(sur["n_surrogates"]),
        "surrogate_median": _f(np.median(sur["lle_surrogates"])) if len(sur["lle_surrogates"]) else None,
        "detected": bool(p is not None and p <= ALPHA),
        "status": "ok" if np.isfinite(lle) else f"lle undefined: {diag.get('reason')}",
    })
    if with_curve and diag.get("mean_log_divergence") is not None:
        out["curve"] = [_f(v) for v in diag["mean_log_divergence"]]
    return out


METHODS = {"baseline": baseline}


# =============================================================================
# Building blocks for Phase 3 candidates (B3).  New estimator code; the delay
# embedding, AAFT surrogates and Cao are the pipeline's own functions.
# =============================================================================

import zlib  # noqa: E402

SURROGATE_ENTROPY = 20260927


def window_rng(x, stream=0):
    """Deterministic RNG unique to this window (hash of the data), so windows do
    not share surrogate streams (cf. Phase 2E finding 8)."""
    h = zlib.crc32(np.ascontiguousarray(np.asarray(x, dtype=float)).tobytes())
    return np.random.default_rng(np.random.SeedSequence([SURROGATE_ENTROPY, int(h), int(stream)]))


def divergence_curve(X, theiler, max_iter=40, n_neighbors=1, min_dist=0.0):
    """
    Mean log nearest-neighbour divergence y[k], k = 0..max_iter-1.

    n_neighbors = 1, min_dist = 0 reproduces fp.rosenstein_lle's curve (same
    Theiler exclusion, same argmin neighbour, same 1e-15 floor, same count rule;
    checked in tests/test_phase3_estimators.py).  n_neighbors > 1 averages the
    distances of the n nearest admissible neighbours before taking the log
    (Kantz 1994 style); min_dist > 0 excludes neighbours closer than min_dist
    (noise-scale exclusion).
    """
    X = np.asarray(X, dtype=float)
    n = len(X)
    D = np.sqrt(np.maximum(((X[:, None, :] - X[None, :, :]) ** 2).sum(-1), 0.0))
    idx = np.arange(n)
    D[np.abs(idx[:, None] - idx[None, :]) <= int(theiler)] = np.inf
    if min_dist > 0:
        D[D <= min_dist] = np.inf
    nb = np.argsort(D, axis=1, kind="stable")[:, :n_neighbors]
    ok = np.isfinite(np.take_along_axis(D, nb, axis=1))
    logd = np.full((n, max_iter), np.nan)
    for k in range(max_iter):
        i = idx[:, None].repeat(n_neighbors, 1)
        ii, jj = i + k, nb + k
        valid = ok & (ii < n) & (jj < n)
        d = np.full((n, n_neighbors), np.nan)
        vi = np.nonzero(valid)
        d[vi] = np.linalg.norm(X[ii[vi]] - X[jj[vi]], axis=1)
        if n_neighbors == 1:
            dd = d[:, 0]
        else:
            with np.errstate(invalid="ignore"):
                dd = np.where(valid.any(1), np.nanmean(np.where(valid, d, np.nan), axis=1), np.nan)
        good = np.isfinite(dd) & (dd > 1e-15)
        logd[good, k] = np.log(dd[good])
    with np.errstate(invalid="ignore"), warnings_ignored():
        y = np.nanmean(logd, axis=0)
    counts = np.sum(np.isfinite(logd), axis=0)
    y[~(np.isfinite(y) & (counts >= max(10, int(0.25 * n))))] = np.nan
    return y, counts


class warnings_ignored:
    def __enter__(self):
        import warnings
        self._c = warnings.catch_warnings()
        self._c.__enter__()
        warnings.simplefilter("ignore")

    def __exit__(self, *a):
        return self._c.__exit__(*a)


def saturation_fit(y, k_start=1, frac=0.5, min_points=3, max_end=12, plateau_from=25):
    """
    Slope of y over k_start..k_end, where k_end is the first step at which y has
    risen `frac` of the way from y[k_start] to the plateau (mean of y[plateau_from:]),
    clipped to [k_start + min_points - 1, max_end].  Returns (slope, k_end) or
    (nan, None) if undefined.
    """
    y = np.asarray(y, float)
    if len(y) <= max_end or not np.all(np.isfinite(y[k_start:k_start + min_points])):
        return np.nan, None
    tail = y[plateau_from:]
    plateau = float(np.nanmean(tail)) if np.any(np.isfinite(tail)) else np.nan
    k_lo = k_start + min_points - 1
    k_end = k_lo
    if np.isfinite(plateau) and plateau > y[k_start]:
        target = y[k_start] + frac * (plateau - y[k_start])
        above = np.flatnonzero(np.isfinite(y[k_start:max_end + 1]) & (y[k_start:max_end + 1] >= target))
        k_end = int(above[0] + k_start) if len(above) else max_end
        k_end = min(max(k_end, k_lo), max_end)
    ks = np.arange(k_start, k_end + 1)
    yy = y[ks]
    good = np.isfinite(yy)
    if good.sum() < min_points:
        return np.nan, None
    return float(np.polyfit(ks[good], yy[good], 1)[0]), int(k_end)


def iaaft_surrogate(x, rng, max_iter=200):
    """Iterative AAFT (Schreiber & Schmitz 1996): exact amplitude distribution,
    power spectrum matched iteratively; stops when the rank order is stable."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    amp = np.abs(np.fft.rfft(x))
    sx = np.sort(x)
    s = rng.permutation(x)
    prev = None
    for _ in range(max_iter):
        ph = np.angle(np.fft.rfft(s))
        s_spec = np.fft.irfft(amp * np.exp(1j * ph), n=n)
        rank = np.argsort(np.argsort(s_spec, kind="stable"), kind="stable")
        s = sx[rank]
        if prev is not None and np.array_equal(rank, prev):
            break
        prev = rank
    return s


def surrogates(x, kind, n_sur, rng):
    if kind == "aaft":
        return [fp.aaft_surrogate(x, rng=rng) for _ in range(n_sur)]
    if kind == "iaaft":
        return [iaaft_surrogate(x, rng) for _ in range(n_sur)]
    raise ValueError(kind)


def zero_one_K(x, rng, n_c=100, frac_n=10):
    """0-1 test for chaos, modified correlation method (Gottwald & Melbourne
    2009): median over c ~ U(pi/5, 4pi/5) of corr(n, D_c(n)), n <= N/frac_n."""
    x = np.asarray(x, dtype=float)
    N = len(x)
    ncut = max(3, N // frac_n)
    j = np.arange(1, N + 1)
    nn = np.arange(1, ncut + 1)
    Ks = []
    ex2 = np.mean(x) ** 2
    for c in rng.uniform(np.pi / 5, 4 * np.pi / 5, n_c):
        p = np.cumsum(x * np.cos(j * c))
        q = np.cumsum(x * np.sin(j * c))
        M = np.array([np.mean((p[k:] - p[:-k]) ** 2 + (q[k:] - q[:-k]) ** 2) for k in nn])
        Dc = M - ex2 * (1 - np.cos(nn * c)) / (1 - np.cos(c))
        if np.std(Dc) > 0:
            Ks.append(np.corrcoef(nn, Dc)[0, 1])
    return float(np.median(Ks)) if Ks else np.nan
