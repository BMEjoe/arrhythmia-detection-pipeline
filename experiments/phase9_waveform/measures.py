"""
Phase 9 Part D: noise-robust, scale-resolved chaos measures (implementations).

Sources (METHODS.md Section 5 has the verification):
 D1 SDLE   Gao J, Hu J, Tung WW, Cao Y et al., "Multiscale analysis of biological data by
           scale-dependent Lyapunov exponent", Front Physiol 2:110 (2011), PMC3264951, Eqs. 1-11
           and the lambda_k_curves procedure (pairs in shells eps_k <= ||Vi - Vj|| <= eps_k + d eps_k,
           |i - j| >= t_uncorrelated, Lambda(t) = <ln ||V_{i+t} - V_{j+t}||> per shell, SDLE = local
           slope (Lambda(k+1) - Lambda(k-1)) / 2 at eps_t = exp(<ln eps_t>)); first defined in Gao et al.
           Phys Rev E 74:066204 (2006).  Euclidean norm (Eq. 6).
 D2 FSLE   Aurell E, Boffetta G, Crisanti A, Paladin G, Vulpiani A, "Predictability in the large: an
           extension of the concept of Lyapunov exponent", J Phys A 30:1 (1997), arXiv:chao-dyn/9606014;
           Boffetta et al., Phys Rep 356:367 (2002), arXiv:nlin/0101029; Cencini et al., PRE 62:427
           (2000), arXiv:nlin/0002018, Eq. 16: lambda(delta) = ln r / <T_r(delta)>, T_r = time for a
           separation to grow from delta to r delta.  Data version: nearest neighbour (max norm,
           Theiler window) of every reference delay vector, separation of the two delay vectors
           followed in time, first-passage times of delta_n = delta_0 r^n.
 D3 eps-entropy  (eps, tau)-entropy per unit time from Grassberger-Procaccia correlation sums
           (Cencini et al. 2000, Eqs. 10-13): h_m(eps) = (1/tau) ln(C_m(eps) / C_{m+1}(eps)), max norm,
           Theiler window; Gaspard & Wang, Phys Rep 235:291 (1993).
 D4 PE and complexity-entropy plane: Bandt C, Pompe B, PRL 88:174102 (2002); Rosso OA et al., PRL
           99:154102 (2007); as implemented and documented in ordpy (Pessa & Ribeiro, Chaos 31:063110
           (2021), arXiv:2102.06786): normalized permutation entropy H and statistical complexity
           C = Q_J[P, P_e] H with the normalized Jensen-Shannon divergence.
 D5 RQA DET  Marwan N, "How to avoid potential pitfalls in recurrence plot based data analysis",
           Int J Bifurcat Chaos 21:1003 (2011), arXiv:1007.2215, Eq. (1): DET = sum_{l>=lmin} l P(l) /
           sum_ij R_ij; recurrence threshold set adaptively to a fixed recurrence rate (Sect. 3.2).
 D6 Local projective noise reduction (GHKSS): Grassberger P, Hegger R, Kantz H, Schaffrath C,
           Schreiber T, Chaos 3:127 (1993) (blocked); algorithm as described in Hegger, Kantz &
           Schreiber, Chaos 9:413 (1999), arXiv:chao-dyn/9810005, Sect. V B (metric tensor Eq. 14 with
           the first and last coordinates fixed, projection Eqs. 15-17 onto the Q-dimensional local
           subspace, curvature correction Theta - <Theta>), and the authors' TISEAN 3.0.1 program
           ghkss.c (used for verification).
"""
from __future__ import annotations

import math

import numpy as np
from numba import njit
from scipy.spatial import cKDTree


def embed(x, m, L=1):
    """Forward delay vectors V_i = (x_i, x_{i+L}, ..., x_{i+(m-1)L}), i = 0 .. N-(m-1)L-1."""
    x = np.asarray(x, float)
    n = len(x) - (m - 1) * L
    return np.stack([x[k * L:k * L + n] for k in range(m)], axis=1)


def standardize(x):
    x = np.asarray(x, float)
    return (x - x.mean()) / max(x.std(), 1e-300)


# ============================================================================= D1 SDLE
@njit(cache=True)
def _lambda_curves(V, pi, pj, T):
    """Sum and count of ln ||V_{i+t} - V_{j+t}|| over pairs, t = 0..T (Euclidean)."""
    n, m = V.shape
    s = np.zeros(T + 1)
    c = np.zeros(T + 1)
    for p in range(pi.shape[0]):
        i, j = pi[p], pj[p]
        for t in range(T + 1):
            if i + t >= n or j + t >= n:
                break
            d = 0.0
            for k in range(m):
                q = V[i + t, k] - V[j + t, k]
                d += q * q
            if d > 0.0:
                s[t] += 0.5 * math.log(d)
                c[t] += 1.0
    return s, c


def sdle(x, m=4, L=1, shells=None, t_unc=None, T=30, max_pairs=20000, rng=None, n_shells=4, eps_top=None,
         standardize_input=True):
    """Scale-dependent Lyapunov exponent (D1).

    shells: list of (eps_lo, eps_hi); default (predeclared): n_shells shells, the first
    [eps_top / sqrt 2, eps_top] with eps_top = sqrt(m / 10) (i.e. SD / sqrt(10) per coordinate of the
    standardized series, Gao et al. Front Bioeng Biotechnol 2020's recommendation, in m
    dimensions), each next shell 1/sqrt(2) smaller.  t_unc default (m - 1) L + 1.
    Returns dict: per-shell (eps_t, lambda) curves, the pooled curve sorted by eps, and pair counts."""
    rng = np.random.default_rng(0) if rng is None else rng
    z = standardize(x) if standardize_input else np.asarray(x, float)
    V = embed(z, m, L)
    n = len(V)
    t_unc = (m - 1) * L + 1 if t_unc is None else t_unc
    if shells is None:
        top = math.sqrt(m / 10.0) if eps_top is None else eps_top
        shells = []
        for k in range(n_shells):
            hi = top / (math.sqrt(2.0) ** k)
            shells.append((hi / math.sqrt(2.0), hi))
    tree = cKDTree(V)
    rmax = max(h for _, h in shells)
    pairs = tree.query_pairs(rmax, output_type="ndarray")
    if len(pairs):
        pairs = pairs[np.abs(pairs[:, 1] - pairs[:, 0]) >= t_unc]
        pairs = pairs[(pairs[:, 0] < n - T) & (pairs[:, 1] < n - T)]
    d0 = np.linalg.norm(V[pairs[:, 0]] - V[pairs[:, 1]], axis=1) if len(pairs) else np.array([])
    out = {"shells": [], "eps": [], "lam": [], "w": []}
    for lo, hi in shells:
        sel = np.nonzero((d0 >= lo) & (d0 <= hi))[0]
        npair = len(sel)
        if npair > max_pairs:
            sel = rng.choice(sel, max_pairs, replace=False)
        if npair < 10:
            out["shells"].append({"lo": lo, "hi": hi, "n_pairs": int(npair)})
            continue
        s, c = _lambda_curves(V, pairs[sel, 0].astype(np.int64), pairs[sel, 1].astype(np.int64), T)
        Lam = s / np.maximum(c, 1)
        eps_t = np.exp(Lam)
        lam = (Lam[2:] - Lam[:-2]) / 2.0
        out["shells"].append({"lo": lo, "hi": hi, "n_pairs": int(npair), "Lambda": Lam, "eps_t": eps_t[1:-1],
                              "lam": lam})
        out["eps"].append(eps_t[1:-1])
        out["lam"].append(lam)
        out["w"].append(np.full(len(lam), float(len(sel))))
    if out["eps"]:
        e = np.concatenate(out["eps"])
        l_ = np.concatenate(out["lam"])
        w = np.concatenate(out["w"])
        o = np.argsort(e)
        out["eps_all"], out["lam_all"], out["w_all"] = e[o], l_[o], w[o]
    return out


def sdle_on_grid(res, grid):
    """Pool the shell curves onto a log-spaced eps grid: pair-count-weighted mean of lambda in
    each grid cell (Gao et al. 2020's recommended combination, simplified to binning)."""
    if "eps_all" not in res:
        return np.full(len(grid) - 1, np.nan)
    e, l_, w = res["eps_all"], res["lam_all"], res["w_all"]
    idx = np.digitize(e, grid) - 1
    out = np.full(len(grid) - 1, np.nan)
    for k in range(len(grid) - 1):
        s = idx == k
        if s.any():
            out[k] = np.sum(l_[s] * w[s]) / np.sum(w[s])
    return out


# ============================================================================= D2 FSLE
@njit(cache=True)
def _fsle_core(V, nn, deltas, maxt):
    """For each pair, the level n is the highest delta_n <= current separation; when the separation
    first exceeds delta_{n+1} (possibly several levels at once in one sampling step), level n gets
    time (t - t_enter) and log growth ln(d(t) / d_enter) (TISEAN fsle's estimator, which equals
    ln r per passage when growth is continuous)."""
    n, m = V.shape
    nl = deltas.shape[0]
    tsum = np.zeros(nl)
    gsum = np.zeros(nl)
    tcnt = np.zeros(nl)
    for i in range(n):
        j = nn[i]
        if j < 0:
            continue
        lev = -1
        t_enter = 0
        d_enter = 0.0
        for t in range(maxt):
            if i + t >= n or j + t >= n:
                break
            d = 0.0
            for k in range(m):
                q = abs(V[i + t, k] - V[j + t, k])
                if q > d:
                    d = q
            if t == 0:
                if d <= 0.0 or d >= deltas[nl - 1]:
                    break
                lev = 0
                while lev + 1 < nl and deltas[lev + 1] <= d:
                    lev += 1
                if deltas[lev] > d:
                    break
                t_enter = 0
                d_enter = d
                continue
            if lev + 1 < nl and d >= deltas[lev + 1]:
                tsum[lev] += t - t_enter
                gsum[lev] += math.log(d / d_enter)
                tcnt[lev] += 1.0
                while lev + 1 < nl and deltas[lev + 1] <= d:
                    lev += 1
                t_enter = t
                d_enter = d
                if lev + 1 >= nl:
                    break
    return tsum, gsum, tcnt


def fsle(x, m=1, L=1, delta0=None, r=math.sqrt(2.0), n_levels=20, theiler=None, maxt=500,
         standardize_input=True):
    """Finite-size Lyapunov exponent (D2): lambda(delta_n) = ln r / <T(delta_n -> delta_{n+1})>
    (Boffetta et al. Eq. 3.37) estimated from sampled data with TISEAN's overshoot correction:
    lambda(delta_n) = sum ln(d_exit / d_enter) / sum T over the passages out of level n.
    Reference points: every delay vector and its nearest neighbour (max norm, |i - j| > theiler,
    searched among the 2 theiler + 2 nearest); the pair is followed from the first level delta_n
    above its initial separation (as TISEAN fsle starts each pair at the level of its initial
    distance).  delta_n = delta0 r^n, default delta0 = 1e-3 SD.  Returns (delta_n, lambda_n,
    counts, n_pairs)."""
    z = standardize(x) if standardize_input else np.asarray(x, float)
    V = embed(z, m, L)
    theiler = (m - 1) * L if theiler is None else theiler
    delta0 = 1e-3 if delta0 is None else delta0
    deltas = delta0 * r ** np.arange(n_levels + 1)
    tree = cKDTree(V)
    k = min(len(V), 2 * theiler + 3)
    dist, idx = tree.query(V, k=k, p=np.inf)
    nn = -np.ones(len(V), dtype=np.int64)
    for i in range(len(V)):
        for d, j in zip(dist[i, 1:], idx[i, 1:]):
            if abs(j - i) > theiler and 0 < d < deltas[-1]:
                nn[i] = j
                break
    tsum, gsum, tcnt = _fsle_core(V, nn, deltas, maxt)
    with np.errstate(invalid="ignore", divide="ignore"):
        lam = np.where(tcnt > 0, gsum / np.maximum(tsum, 1e-300), np.nan)
    return deltas[:-1], lam[:-1], tcnt[:-1], int((nn >= 0).sum())


# ============================================================================= D3 eps-entropy
@njit(cache=True)
def _corr_counts(x, m_max, L, eps, theiler):
    """Pair counts C_m(eps) (max norm) for m = 1..m_max for every eps (sorted ascending): count of
    pairs i < j, j - i > theiler, whose m-dimensional max-norm distance < eps, using a common
    set of pairs (vectors defined for m_max)."""
    n = x.shape[0] - (m_max - 1) * L
    ne = eps.shape[0]
    cnt = np.zeros((m_max, ne))
    npairs = 0.0
    for i in range(n):
        for j in range(i + theiler + 1, n):
            npairs += 1.0
            d = 0.0
            for k in range(m_max):
                q = abs(x[i + k * L] - x[j + k * L])
                if q > d:
                    d = q
                # pair contributes to C_{k+1}(eps) for eps > d
                if d >= eps[ne - 1]:
                    break
                lo = 0
                hi = ne
                while lo < hi:
                    mid = (lo + hi) // 2
                    if eps[mid] > d:
                        hi = mid
                    else:
                        lo = mid + 1
                cnt[k, lo] += 1.0
    for k in range(m_max):
        for e in range(1, ne):
            cnt[k, e] += cnt[k, e - 1]
    return cnt, npairs


def correlation_sums(x, m_max=6, L=1, eps=None, theiler=None, standardize_input=True):
    """C_m(eps), m = 1..m_max, max norm, Theiler window (default (m_max-1) L)."""
    z = standardize(x) if standardize_input else np.asarray(x, float)
    eps = np.geomspace(1e-3, 4.0, 60) if eps is None else np.asarray(eps, float)
    theiler = (m_max - 1) * L if theiler is None else theiler
    cnt, npairs = _corr_counts(z, m_max, L, np.sort(eps), theiler)
    return np.sort(eps), cnt / max(npairs, 1.0)


def eps_entropy(C, tau=1.0):
    """h_m(eps) = ln(C_m / C_{m+1}) / tau for m = 1 .. m_max - 1 (rows)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log(C[:-1] / C[1:]) / tau


# ============================================================================= D4 PE / CECP
def ordinal_distribution(x, d=6, tau=1, return_missing=True):
    """Bandt-Pompe ordinal-pattern distribution (patterns = argsort of each window, as ordpy;
    ties broken by order of appearance).  Returns probabilities of all d! patterns (zeros
    included) in lexicographic order of the permutations."""
    from itertools import permutations
    x = np.asarray(x, float)
    n = len(x) - (d - 1) * tau
    W = np.stack([x[k * tau:k * tau + n] for k in range(d)], axis=1)
    pats = np.argsort(W, axis=1, kind="stable")
    allp = {p: i for i, p in enumerate(permutations(range(d)))}
    code = np.array([allp[tuple(p)] for p in pats])
    counts = np.bincount(code, minlength=len(allp))
    return counts / counts.sum()


def _shannon(p):
    p = p[p > 0]
    return float(-np.sum(p * np.log(p)))


def pe_cecp(x, d=6, tau=1):
    """Normalized permutation entropy H and statistical complexity C (Rosso et al. 2007)."""
    p = ordinal_distribution(x, d, tau)
    n = len(p)
    S = _shannon(p)
    H = S / math.log(n)
    pe = np.full(n, 1.0 / n)
    js = _shannon(0.5 * (p + pe)) - 0.5 * S - 0.5 * math.log(n)
    q0 = -2.0 / (((n + 1.0) / n) * math.log(n + 1.0) - 2.0 * math.log(2.0 * n) + math.log(n))
    C = q0 * js * H
    return H, C


# ============================================================================= D5 RQA DET
@njit(cache=True)
def _det_from_matrix(dist, eps, lmin, theiler):
    n = dist.shape[0]
    rec = 0.0
    inl = 0.0
    for k in range(theiler + 1, n):             # upper triangle diagonals, offset k
        run = 0
        for i in range(n - k):
            if dist[i, i + k] <= eps:
                run += 1
                rec += 1.0
            else:
                if run >= lmin:
                    inl += run
                run = 0
        if run >= lmin:
            inl += run
    return inl / rec if rec > 0 else np.nan, rec


def rqa_det(x, m=3, tau=1, rr=0.05, lmin=2, theiler=0, standardize_input=True):
    """DET with the recurrence threshold set so the recurrence rate (off the Theiler band) equals
    rr (Marwan 2011, Sect. 3.2 'adaptive recurrence threshold').  Euclidean norm (Marwan's
    examples).  Returns (DET, eps)."""
    z = standardize(x) if standardize_input else np.asarray(x, float)
    V = embed(z, m, tau)
    n = len(V)
    G = V @ V.T
    sq = np.diag(G)
    dist = np.sqrt(np.maximum(sq[:, None] + sq[None, :] - 2 * G, 0.0))
    iu = np.triu_indices(n, k=theiler + 1)
    eps = float(np.quantile(dist[iu], rr))
    det, rec = _det_from_matrix(dist, eps, lmin, theiler)
    return float(det), eps


# ============================================================================= D6 GHKSS
def ghkss(x, m=7, q=2, k_min=50, iterations=3, mineps=1e-3, euclid=False):
    """Locally projective noise reduction (Grassberger et al. 1993), following TISEAN 3.0.1 ghkss.c
    (Hegger, Kantz & Schreiber 1999, Sect. V B) step by step:
    - the series is rescaled to [0, 1];
    - delay vectors of dimension m, unit delay;
    - for each vector the neighbourhood radius is the smallest eps = mineps * sqrt(2)^l with at
      least k_min vectors (itself included) within eps (max norm);
    - metric weights w_i = 1e3 for the first and last coordinate, 1 otherwise (euclid: all 1);
      Gamma_ij = (cov_ij of the neighbourhood) w_i w_j;
    - correction of vector n: c_i = sum over the m - q eigenvectors a of Gamma with the smallest
      eigenvalues of sum_k (s_k - mean_k) a_k a_i w_k / w_i;
    - curvature correction: c - <c>_neighbourhood; each sample receives
      sum over the vectors containing it of (c_i - <c_i>) / (trace w_i), trace = sum_i 1 / w_i,
      and is decreased by it;
    - if any vector found its neighbours at the first radius, mineps is divided by sqrt(2) for
      the next iteration."""
    x = np.asarray(x, float)
    lo, rng_ = float(x.min()), float(x.max() - x.min())
    if rng_ <= 0:
        return x.copy()
    s_ = (x - lo) / rng_
    w = np.ones(m) if euclid else np.array([1e3] + [1.0] * (m - 2) + [1e3])
    trace = float(np.sum(1.0 / w))
    fac = math.sqrt(2.0)
    for _ in range(iterations):
        V = embed(s_, m, 1)
        n = len(V)
        tree = cKDTree(V)
        level = np.full(n, -1)
        eps = mineps
        nbs = [None] * n
        lev = 0
        while (level < 0).any() and eps < 2.0:
            todo = np.nonzero(level < 0)[0]
            res = tree.query_ball_point(V[todo], eps, p=np.inf)
            for t, idx in zip(todo, res):
                if len(idx) >= k_min:
                    level[t] = lev
                    nbs[t] = np.asarray(idx)
            eps *= fac
            lev += 1
        corr = np.zeros((n, m))
        okv = level >= 0
        for i in np.nonzero(okv)[0]:
            U = V[nbs[i]]
            av = U.mean(0)
            G = (np.cov(U.T, bias=True) * w[:, None]) * w[None, :]
            ev, a = np.linalg.eigh(G)                       # ascending eigenvalues
            A = a[:, :m - q]
            corr[i] = (A @ (A.T @ ((V[i] - av) * w))) / w
        delta = np.zeros(len(s_))
        for i in np.nonzero(okv)[0]:
            c = corr[i] - corr[nbs[i]].mean(0)
            delta[i:i + m] += c / (trace * w)
        s_ = s_ - delta
        if (level == 0).any():
            mineps /= fac
    return s_ * rng_ + lo
