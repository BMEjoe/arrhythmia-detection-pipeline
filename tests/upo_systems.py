"""Deterministic synthetic systems and cached So-method runs for the UPO tests.

Systems and published reference values come from the primary sources:
  PRL = So, Ott, Schiff, Kaplan, Sauer & Grebogi, PRL 76, 4705 (1996)
  PRE = So, Ott, Sauer, Gluckman, Grebogi & Schiff, PRE 55, 5398 (1997)
"""
import sys
import pathlib
from dataclasses import replace
from functools import lru_cache

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import final_pipeline as fp  # noqa: E402

# ---------------------------------------------------------------------------
# Skewed Henon map, PRE Eq. (31): u' = b1 exp(b2 u)(1.4 - u^2) + 0.3 v, v' = u
# ---------------------------------------------------------------------------
HENON_B1 = 0.965
HENON_B2 = 0.25

# PRE Table I "exact numerical values" (Newton's method, +-0.0001)
PRE_HENON_FIXED_POINTS = (0.9302, -1.9030)
PRE_HENON_PERIOD2 = (1.3014, -0.5604)

# Source-style operational parameters for the skewed Henon runs (PRE Sec. III C):
# d = 2, M = 2 spatial neighbours, K = 1 (clusters of 2), R = 100, kappa = 3,
# N = 1024, one-sample map (tau = 1).
HENON_N = 1024


def henon_f(u, v, b1=HENON_B1, b2=HENON_B2):
    return b1 * np.exp(b2 * u) * (1.4 - u * u) + 0.3 * v


def henon_series(n=HENON_N, burn=1000, seed=0, noise=0.0, noise_seed=11):
    rng = np.random.default_rng(seed)
    u, v = rng.uniform(0.0, 0.5, 2)
    out = np.empty(n)
    for i in range(n + burn):
        u, v = henon_f(u, v), u
        if i >= burn:
            out[i - burn] = u
    if noise:
        out = out + noise * np.random.default_rng(noise_seed).uniform(-1, 1, n)
    return out


def henon_jacobian(u, b1=HENON_B1, b2=HENON_B2):
    """Analytic Jacobian of (u, v) -> (f(u, v), u)."""
    fu = b1 * np.exp(b2 * u) * (b2 * (1.4 - u * u) - 2.0 * u)
    return np.array([[fu, 0.3], [1.0, 0.0]])


def henon_fixed_points():
    """Solve u = f(u, u) by Newton's method from both published starting values."""
    roots = []
    for u in PRE_HENON_FIXED_POINTS:
        for _ in range(60):
            g = henon_f(u, u) - u
            dg = henon_jacobian(u)[0, 0] + 0.3 - 1.0
            u -= g / dg
        roots.append(float(u))
    return tuple(roots)


def henon_period2_orbit():
    """Solve the 2-cycle a -> b -> a of u_{n+1} = f(u_n, u_{n-1})."""
    a, b = PRE_HENON_PERIOD2
    for _ in range(80):
        F = np.array([henon_f(a, b) - b, henon_f(b, a) - a])
        ja, jb = henon_jacobian(a)[0, 0], henon_jacobian(b)[0, 0]
        Jm = np.array([[ja, 0.3 - 1.0], [0.3 - 1.0, jb]])
        a, b = np.array([a, b]) - np.linalg.solve(Jm, F)
    return float(a), float(b)


def henon_config(randomization="prl_norm", **kw):
    base = dict(so_jacobian_neighbors=2, so_neighbors_K=1, so_random_R=100,
                so_random_R_period_p=100, so_kappa=3.0, so_max_backbones=None,
                so_randomization=randomization, upo_map_mode="one_sample")
    base.update(kw)
    return replace(fp.CFG, **base)


@lru_cache(maxsize=None)
def henon_embedded(n=HENON_N, seed=0):
    return fp._embed_backward(henon_series(n, seed=seed), 1, 2)


@lru_cache(maxsize=None)
def henon_period1(randomization="prl_norm", rng_seed=1, n=HENON_N, kappa=3.0, seed=0):
    X = henon_embedded(n, seed)
    return fp.detect_so_fixed_points(X, 1, henon_config(randomization, so_kappa=kappa),
                                     np.random.default_rng(rng_seed))


@lru_cache(maxsize=None)
def henon_period2(randomization="prl_norm", rng_seed=2):
    return fp.detect_so_period_p(henon_embedded(), 2, 1, henon_config(randomization),
                                 np.random.default_rng(rng_seed))


@lru_cache(maxsize=None)
def henon_period2_verified(randomization="prl_norm"):
    """Period-2 detection followed by the default PROJECT verification gates."""
    import copy
    r = copy.deepcopy(henon_period2(randomization))
    return fp.apply_verification_gates(r, henon_embedded(), henon_config(randomization))


@lru_cache(maxsize=None)
def henon_significance(n=512, surrogates=12):
    x = henon_series(n)
    cfg = henon_config(so_random_R=50, so_surrogate_count=surrogates)
    r = fp.detect_so_fixed_points(fp._embed_backward(x, 1, 2), 1, cfg, np.random.default_rng(1))
    return fp.assess_so_significance(r, x, cfg, np.random.default_rng(2))


@lru_cache(maxsize=None)
def noise_significance(n=512, surrogates=12):
    x = np.random.default_rng(5).standard_normal(n)
    cfg = henon_config(so_random_R=50, so_surrogate_count=surrogates)
    r = fp.detect_so_fixed_points(fp._embed_backward(x, 1, 2), 1, cfg, np.random.default_rng(1))
    return fp.assess_so_significance(r, x, cfg, np.random.default_rng(2))


def nearest(values, target):
    values = np.asarray(values, dtype=float)
    return float(values[np.argmin(np.abs(values - target))])


def period1_locations(result):
    return [c["scalar_location"] for c in result["source_peak_candidates"]]


def orbit_distance(orbit, target):
    """Distance between two cyclic sequences, minimised over rotations."""
    orbit = np.asarray(orbit, dtype=float)
    target = np.asarray(target, dtype=float)
    return min(float(np.max(np.abs(np.roll(orbit, k) - target))) for k in range(len(orbit)))


# ---------------------------------------------------------------------------
# One-dimensional maps
# ---------------------------------------------------------------------------

def iterate_1d(f, x0, n, burn=500):
    x = x0
    out = np.empty(n)
    for i in range(n + burn):
        x = f(x)
        if i >= burn:
            out[i - burn] = x
    return out


def logistic_series(n, r=3.92, x0=0.3, burn=500):
    return iterate_1d(lambda z: r * z * (1 - z), x0, n, burn)


SKEWED_LOGISTIC_MU = 6.1


def skewed_logistic_series(n=256, x0=0.3):
    """PRE Sec. II F: f(z) = mu exp(-z) z (1 - z), mu = 6.1."""
    return iterate_1d(lambda z: SKEWED_LOGISTIC_MU * np.exp(-z) * z * (1 - z), x0, n)


def skewed_logistic_fixed_point():
    z = 0.67
    for _ in range(60):
        g = SKEWED_LOGISTIC_MU * np.exp(-z) * (1 - z) - 1.0
        dg = -SKEWED_LOGISTIC_MU * np.exp(-z) * (2 - z)
        z -= g / dg
    return float(z)


# Logistic r = 4 period-3 orbits: x_k = sin^2(pi 2^k / 7) and sin^2(pi 2^k / 9)
LOGISTIC4_PERIOD3 = (
    tuple(np.sin(np.pi * np.array([4, 1, 2]) / 7.0) ** 2),
    tuple(np.sin(np.pi * np.array([4, 1, 2]) / 9.0) ** 2),
)


@lru_cache(maxsize=None)
def logistic4_period3():
    x = logistic_series(512, r=4.0, x0=0.2345)
    cfg = replace(fp.CFG, so_jacobian_neighbors=1, so_random_R_period_p=20, so_kappa=3.0,
                  so_max_backbones=None, so_randomization="pre_tensor")
    return fp.detect_so_period_p(x[:, None], 3, 1, cfg, np.random.default_rng(0))


# ---------------------------------------------------------------------------
# Lozi map (piecewise linear, exact local Jacobians)
# ---------------------------------------------------------------------------
LOZI_A = 1.7
LOZI_B = 0.5


def lozi_series(n=2048, burn=500):
    x, y = 0.1, 0.1
    out = np.empty(n)
    for i in range(n + burn):
        x, y = 1.0 - LOZI_A * abs(x) + LOZI_B * y, x
        if i >= burn:
            out[i - burn] = x
    return out


def lozi_fixed_point():
    # x = 1 - a x + b x for x > 0
    return 1.0 / (1.0 + LOZI_A - LOZI_B)


# ---------------------------------------------------------------------------
# Ikeda map (PRE Sec. IV): a = 0.75, b = 9.0, observed v; fixed point v* = 0.537
# ---------------------------------------------------------------------------
IKEDA_FIXED_POINT_V = 0.537


def ikeda_series(n=1024, a=0.75, b=9.0, burn=500):
    u, v = 0.1, 0.1
    out = np.empty(n)
    for i in range(n + burn):
        t = 0.4 - b / (1 + u * u + v * v)
        u, v = 1.0 + a * (u * np.cos(t) - v * np.sin(t)), a * (u * np.sin(t) + v * np.cos(t))
        if i >= burn:
            out[i - burn] = v
    return out


# ---------------------------------------------------------------------------
# RR-like series for end-to-end tests (synthetic; no MIT-BIH)
# ---------------------------------------------------------------------------

def chaotic_rr(n=256, seed=0):
    return 0.8 + 0.05 * henon_series(n, seed=seed)


def periodic_rr(n=256, seed=0):
    t = np.arange(n)
    return 0.8 + 0.05 * np.sin(2 * np.pi * t / 7.3) + \
        0.001 * np.random.default_rng(seed).standard_normal(n)


def fake_upo(status, n=200, source=0, significant=None, verified=None, rJ=np.nan,
             assessed=None, map_mode="one_sample"):
    """Minimal run_upo_analysis-shaped dict for feature-semantics tests."""
    def peaks(k):
        return [{"orbit_points": np.zeros((1, 2))} for _ in range(k)]
    sig = None if significant is None else peaks(significant)
    ver = None if verified is None else peaks(verified)
    return {
        "status": status, "n_embedded_points": n, "map_mode": map_mode,
        "significance_assessed": bool(sig is not None) if assessed is None else assessed,
        "source_peak_candidates": peaks(source),
        "significant_uop_candidates": sig,
        "verified_unstable": ver,
        "coverage": {"source": {"coverage": 0.5} if source else None,
                     "significant": {"coverage": 0.25} if significant else None,
                     "verified_unstable": {"coverage": 0.125} if verified else None},
        "source_rJ": rJ,
    }
