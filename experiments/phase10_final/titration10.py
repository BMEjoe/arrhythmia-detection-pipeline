"""
Phase 10 Question 2: noise titration of chaos (Poon & Barahona 2001) with the Barahona-Poon (1996)
Volterra-Wiener-Korenberg (VWK) nonlinearity indicator, implemented from the accessible sources
(METHODS.md Section 2 lists every source, choice and deviation):

  [PB01]  Poon & Barahona, PNAS 98:7107 (2001), PMC34630 (text read via the PMC page): algorithm,
          Fig. 1 legend (bisection on alpha, y + alpha*xi with unit-variance white xi, 1 % level,
          NL = alpha / sigma_y, 10 noise realizations averaged, 1,000-point series), Fig. 2 values.
  [Wu09]  Wu et al., PLoS ONE 4:e4323 (2009), Methods, Eq. 1 (VWK model) and Eq. 2
          (C(r) = log eps(r) + r / N, r = number of LEADING terms of Eq. 1), F-test at 1 %.
  [Wy06]  Wysocki et al. 2006, arXiv:nlin/0606032, Methods 2.4: best linear model = d = 1 with kappa
          minimising the Akaike criterion; best nonlinear model with d > 1; F-test / Mann-Whitney at
          1 %; default routine settings kappa = 6, d = 3.
  [BP96]  Barahona & Poon, Nature 381:215 (1996): NOT accessible here (paywall); its definitions
          are taken from the restatements above.

Model (Eq. 1, [Wu09]):  y_n = a_0 + a_1 y_{n-1} + ... + a_k y_{n-k} + a_{k+1} y_{n-1}^2
                         + a_{k+2} y_{n-1} y_{n-2} + ... + a_{M-1} y_{n-k}^d,
terms ordered by degree, and within a degree lexicographically in the lag indices
(itertools.combinations_with_replacement order), M = (k + d)! / (k! d!).
Coefficients: least squares on the leading r terms for every r, computed by modified Gram-Schmidt
orthogonalisation of the ordered terms (the recursive, orthogonal estimation of Korenberg that
[Wu09] calls "recursively estimated"; identical to ordinary least squares on each leading set).
eps(r) = sqrt(RSS(r) / sum (y_n - mean y)^2)  (normalised RMS one-step error; with this reading
C(r) = AIC / (2 N) + const, the "Akaike criterion" named by [PB01], [Wy06], [Wu09]).
"""
from __future__ import annotations

from itertools import combinations_with_replacement
from math import comb

import numpy as np
from scipy import stats

KAPPA_MAX = 6      # [Wy06] default of the Poon-lab routine (nonlinear models)
# Memory limit of the LINEAR family: linear models may use as many terms as the largest nonlinear
# model (M = 84 terms -> kappa_lin <= 83), so a limit cycle with up to 41 harmonics has an exact
# linear description ([PB01]: periodic signals are "described both by nonlinear models and by linear
# models of enough memory").  kappa_lin <= 6 fails [PB01] Fig. 2's periodic example (METHODS.md 2.3).
D_MAX = 3          # [Wy06] default of the Poon-lab routine
KAPPA_LIN_MAX = comb(KAPPA_MAX + D_MAX, D_MAX) - 1
ALPHA = 0.01       # [PB01], [Wu09], [Wy06]
N_REAL = 10        # [PB01] Fig. 1 legend (10 noise realisations); [Wu09] "5-10 times"
HI0 = 2.0          # initial bisection bracket for alpha / sigma_y (doubled while still detected)
HI_MAX = 16.0
N_BISECT = 10      # resolution HI0 / 2**10 ~ 0.2 % of sigma_y
MIN_ROWS = 300     # minimum number of usable regression rows (masked variant)


def n_terms(kappa, d):
    return comb(kappa + d, d)


def design(y, kappa, d, rows):
    """Ordered VWK terms (Eq. 1) for target indices `rows` (each row n uses y[n-1] .. y[n-kappa])."""
    lags = np.stack([y[rows - j] for j in range(1, kappa + 1)], axis=1)
    cols = [np.ones(len(rows))]
    for deg in range(1, d + 1):
        for c in combinations_with_replacement(range(kappa), deg):
            cols.append(np.prod(lags[:, c], axis=1))
    return np.stack(cols, axis=1)


def leading_rss(X, t, tol=1e-10):
    """RSS of the least-squares fit on the leading r columns, r = 1..M (modified Gram-Schmidt with
    re-orthogonalisation).  A column that is numerically in the span of the preceding ones adds
    nothing (Korenberg's procedure skips such terms)."""
    n, M = X.shape
    Q = np.zeros((n, M))
    rss = np.empty(M)
    res = t.astype(float).copy()
    cur = float(res @ res)
    k = 0
    for j in range(M):
        v = X[:, j].astype(float).copy()
        nv0 = np.linalg.norm(v)
        if k:
            for _ in range(2):
                v -= Q[:, :k] @ (Q[:, :k].T @ v)
        nv = np.linalg.norm(v)
        if nv0 > 0 and nv > tol * nv0:
            q = v / nv
            Q[:, k] = q
            k += 1
            b = float(q @ res)
            res -= b * q
            cur = float(res @ res)
        rss[j] = cur
    return rss, res


EPS_POWER = 0.5    # C(r) = EPS_POWER * log(RSS / SST) + r / N; 0.5 <=> eps = normalised RMS (primary)


def _fits(y, rows, kappa_max, d_max, kappa_lin_max):
    """All candidate models: {(kappa, d_family, r): (C, rss, residuals?)} reduced to the best linear
    and best nonlinear model (with residual vectors)."""
    t = y[rows]
    sst = float(np.sum((t - t.mean()) ** 2))
    N = len(rows)
    # linear family: leading terms of the kappa_lin_max, d = 1 model (r = 1 .. kappa_lin_max + 1)
    Xl = design(y, kappa_lin_max, 1, rows)
    rss_l, _ = leading_rss(Xl, t)
    r = np.arange(1, Xl.shape[1] + 1)
    C_l = EPS_POWER * np.log(np.maximum(rss_l / sst, 1e-300)) + r / N
    il = int(np.argmin(C_l))
    best_lin = {"C": float(C_l[il]), "r": int(il + 1), "kappa": int(il), "rss": float(rss_l[il])}
    best_nl = None
    for kappa in range(1, kappa_max + 1):
        X = design(y, kappa, d_max, rows)
        rss, _ = leading_rss(X, t)
        rr = np.arange(1, X.shape[1] + 1)
        C = EPS_POWER * np.log(np.maximum(rss / sst, 1e-300)) + rr / N
        lo = kappa + 1                       # first nonlinear term is column index kappa + 1
        if lo >= X.shape[1]:
            continue
        j = lo + int(np.argmin(C[lo:]))
        cand = {"C": float(C[j]), "r": int(j + 1), "kappa": kappa,
                "d": int(next(dd for dd in range(2, d_max + 1) if j + 1 <= n_terms(kappa, dd))),
                "rss": float(rss[j])}
        if best_nl is None or cand["C"] < best_nl["C"]:
            best_nl = cand
    return best_lin, best_nl, N, sst


def _residuals(y, rows, r, kappa, d):
    X = design(y, kappa, d, rows)[:, :r]
    coef, *_ = np.linalg.lstsq(X, y[rows], rcond=None)
    return y[rows] - X @ coef


def vwk_test(y, valid_rows=None, kappa_max=KAPPA_MAX, d_max=D_MAX, alpha=ALPHA, test="F_var",
             kappa_lin_max=None):
    """Barahona-Poon nonlinearity detection on series y.
    valid_rows: optional boolean array over targets n (True = row usable); default all n >= kappa_max.
    test: 'F_var'   one-sided F-test of the residual variances of the best linear vs best nonlinear
                    model, df (N - r_lin, N - r_nl)  [PRIMARY; see METHODS.md 2.3]
          'F_nested' partial F-test (RSS_lin - RSS_nl)/(r_nl - r_lin) / (RSS_nl/(N - r_nl))
          'MWW'     one-sided Mann-Whitney U on the squared residuals ([PB01] "Whitney-Mann").
    Nonlinear iff C_nl < C_lin and p < alpha."""
    kl = KAPPA_LIN_MAX if kappa_lin_max is None else kappa_lin_max
    y = np.asarray(y, float)
    s = np.std(y)
    y = (y - y.mean()) / (s if s > 0 else 1.0)
    n = len(y)
    rows = np.arange(max(kappa_max, kl), n)
    if valid_rows is not None:
        rows = rows[np.asarray(valid_rows, bool)[rows]]
    if len(rows) < max(MIN_ROWS if valid_rows is not None else 0, 3 * n_terms(kappa_max, d_max)):
        return {"nonlinear": False, "analysable": False, "n_rows": int(len(rows))}
    lin, nl, N, sst = _fits(y, rows, kappa_max, d_max, kl)
    out = {"analysable": True, "n_rows": int(N), "C_lin": lin["C"], "C_nl": nl["C"], "r_lin": lin["r"],
           "r_nl": nl["r"], "kappa_lin": lin["kappa"], "kappa_nl": nl["kappa"], "d_nl": nl["d"],
           "eps_lin": float(np.sqrt(lin["rss"] / sst)), "eps_nl": float(np.sqrt(nl["rss"] / sst))}
    better = nl["C"] < lin["C"]
    if test == "F_var":
        df1, df2 = N - lin["r"], N - nl["r"]
        F = (lin["rss"] / df1) / max(nl["rss"] / df2, 1e-300)
        p = float(stats.f.sf(F, df1, df2))
    elif test == "F_nested":
        dr = nl["r"] - lin["r"]
        if dr <= 0:
            p = 0.0 if nl["rss"] < lin["rss"] else 1.0
        else:
            F = ((lin["rss"] - nl["rss"]) / dr) / max(nl["rss"] / (N - nl["r"]), 1e-300)
            p = float(stats.f.sf(F, dr, N - nl["r"]))
    elif test == "MWW":
        rl = _residuals(y, rows, lin["r"], max(lin["kappa"], 1), 1) if lin["r"] > 1 else y[rows] - y[rows].mean()
        rn = _residuals(y, rows, nl["r"], nl["kappa"], d_max)
        p = float(stats.mannwhitneyu(rl ** 2, rn ** 2, alternative="greater").pvalue)
    else:
        raise ValueError(test)
    out["p"] = p
    out["nonlinear"] = bool(better and p < alpha)
    return out


def titrate(y, rng, valid_rows=None, n_real=N_REAL, test="F_var", kappa_max=KAPPA_MAX, d_max=D_MAX,
            scheme="bisection", step=0.01, kappa_lin_max=None):
    """Noise titration [PB01]: for each of n_real unit-variance white Gaussian xi, find the largest alpha
    such that y + alpha * xi is still nonlinear (bisection, as [PB01] Fig. 1 legend); NL_i = alpha /
    sigma_y; NL = mean NL_i.  NL = 0 if y itself is not nonlinear.
    scheme='increment': instead add noise in increments of `step` * sigma_y until the first
    non-detection ([Wu09], [Wy06] wording); NL_i = last detected level (verification only)."""
    y = np.asarray(y, float)
    sy = float(np.std(y))
    kw = dict(valid_rows=valid_rows, test=test, kappa_max=kappa_max, d_max=d_max, kappa_lin_max=kappa_lin_max)
    base = vwk_test(y, **kw)
    out = {"base": base, "analysable": base.get("analysable", False)}
    if not base.get("nonlinear"):
        out.update(NL=0.0, NL_i=[0.0] * n_real, positive=False, NL_power=0.0)
        return out
    nls = []
    for _ in range(n_real):
        xi = rng.standard_normal(len(y))
        xi = (xi - xi.mean()) / xi.std()
        det = lambda a: vwk_test(y + a * sy * xi, **kw)["nonlinear"]  # noqa: E731
        if scheme == "bisection":
            lo, hi = 0.0, HI0
            while det(hi) and hi < HI_MAX:
                lo, hi = hi, 2 * hi
            if lo >= HI_MAX:
                nls.append(HI_MAX)
                continue
            for _ in range(N_BISECT):
                mid = 0.5 * (lo + hi)
                if det(mid):
                    lo = mid
                else:
                    hi = mid
            nls.append(lo)
        else:
            a = 0.0
            while a + step <= HI_MAX and det(a + step):
                a += step
            nls.append(a)
    NL = float(np.mean(nls))
    out.update(NL=NL, NL_i=[float(v) for v in nls], positive=bool(NL > 0),
               NL_power=float(np.mean(np.square(nls))))
    return out


# ------------------------------------------------------------------------------- masks / editing
def rows_from_mask(mask_intervals, kappa_max=None):
    """K3-style masking for the VWK regression: row n (target interval n, lags n-1..n-kappa_max, with
    kappa_max the longest lag of any candidate model) is usable iff none of the intervals
    n - kappa_max .. n is ectopy-related (the K3 rule: a delay vector is used only if its whole span
    is free of masked intervals)."""
    if kappa_max is None:
        kappa_max = max(KAPPA_MAX, KAPPA_LIN_MAX)
    m = np.asarray(mask_intervals, bool)
    n = len(m)
    c = np.concatenate([[0], np.cumsum(m.astype(int))])
    ok = np.zeros(n, bool)
    idx = np.arange(kappa_max, n)
    ok[idx] = (c[idx + 1] - c[idx - kappa_max]) == 0
    return ok


def ectopy_mask(labels):
    """K3 rule: interval k (beat k -> k+1) is ectopy-related iff beat k or beat k+1 is not 'N'."""
    bad = np.asarray(labels).astype(str) != "N"
    return bad[:-1] | bad[1:]
