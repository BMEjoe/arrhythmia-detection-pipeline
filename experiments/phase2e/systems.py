"""
Phase 2E synthetic systems.

Every window is generated from its own trajectory / realization (no window
is cut from a trajectory shared with another nominally independent window).
Paired designs (noise levels, map modes, surrogate counts) reuse the SAME
clean window on purpose; this dependence is documented in the README.
"""
from __future__ import annotations

import numpy as np

from . import config as C


def seed_sequence(system, window_length, seed):
    return np.random.SeedSequence([C.PHASE2E_ENTROPY, C.SYSTEM_CODE[system],
                                   int(window_length), int(seed)])


def data_seed(system, window_length, seed):
    """Integer recorded in every result row identifying the data stream."""
    return int(seed_sequence(system, window_length, seed).generate_state(1)[0])


def _rngs(system, n, seed):
    data_ss, noise_ss = seed_sequence(system, n, seed).spawn(2)
    return np.random.default_rng(data_ss), noise_ss


# ---------------------------------------------------------------------------
# Generators.  Each returns the observable passed to the pipeline.
# ---------------------------------------------------------------------------

def constant(n, rng):
    return np.full(n, rng.uniform(*C.CONSTANT_RANGE))


def white_noise(n, rng):
    return rng.standard_normal(n)


def ar1(n, rng, phi=C.AR1_PHI):
    """x_n = phi x_{n-1} + e_n, e_n ~ N(0,1); x_0 ~ N(0, 1/(1-phi^2)); TRANSIENT discarded."""
    total = n + C.TRANSIENT
    e = rng.standard_normal(total)
    x = np.empty(total)
    x[0] = rng.standard_normal() / np.sqrt(1.0 - phi * phi)
    for i in range(1, total):
        x[i] = phi * x[i - 1] + e[i]
    return x[C.TRANSIENT:]


def sinusoid(n, rng, period=C.SINUSOID_PERIOD):
    """x_n = sin(2 pi n / period + phi0), phi0 ~ U(0, 2 pi)."""
    phi0 = rng.uniform(0.0, 2.0 * np.pi)
    return np.sin(2.0 * np.pi * np.arange(n) / period + phi0)


def _logistic(n, rng, r):
    x = rng.uniform(*C.LOGISTIC_X0_RANGE)
    out = np.empty(n)
    for i in range(n + C.TRANSIENT):
        x = r * x * (1.0 - x)
        if i >= C.TRANSIENT:
            out[i - C.TRANSIENT] = x
    return out


def logistic(n, rng):
    """x_{n+1} = 4 x_n (1 - x_n), x_0 ~ U(0.05, 0.95), TRANSIENT discarded."""
    return _logistic(n, rng, C.LOGISTIC_R)


def logistic_p4(n, rng):
    """Logistic r = 3.5: stable period-4 cycle (non-chaotic periodic control)."""
    return _logistic(n, rng, C.LOGISTIC_P4_R)


def henon_xy(n, rng, a=C.HENON_A, b=C.HENON_B):
    x, y = rng.uniform(*C.HENON_IC_RANGE, size=2)
    xs = np.empty(n)
    for i in range(n + C.TRANSIENT):
        x, y = 1.0 - a * x * x + y, b * x
        if i >= C.TRANSIENT:
            xs[i - C.TRANSIENT] = x
    return xs


def henon(n, rng):
    """Standard Henon (a=1.4, b=0.3); observable x_n; (x0,y0) ~ U(-0.1,0.1)^2."""
    return henon_xy(n, rng)


# Skewed Henon, PRE Eq. (31) -- the reconstruction benchmark (tests/upo_systems.py)
SKEWED_B1, SKEWED_B2 = 0.965, 0.25


def skewed_henon(n, rng):
    u, v = rng.uniform(0.0, 0.5, 2)
    out = np.empty(n)
    for i in range(n + C.TRANSIENT):
        u, v = SKEWED_B1 * np.exp(SKEWED_B2 * u) * (1.4 - u * u) + 0.3 * v, u
        if i >= C.TRANSIENT:
            out[i - C.TRANSIENT] = u
    return out


GENERATORS = {
    "constant": constant, "white_noise": white_noise, "ar1": ar1,
    "sinusoid": sinusoid, "logistic_p4": logistic_p4, "logistic": logistic,
    "henon": henon, "skewed_henon": skewed_henon,
}


def generate(system, n, seed, snr_db=None):
    """
    Clean window for (system, n, seed); if snr_db is given, add white Gaussian
    observation noise with sigma^2 = var(clean) / 10^(snr_db/10), drawn from
    a stream that depends on (system, n, seed, snr_db) only.
    """
    rng, noise_ss = _rngs(system, n, seed)
    x = GENERATORS[system](int(n), rng)
    if not np.all(np.isfinite(x)):
        raise FloatingPointError(f"{system} seed {seed}: non-finite trajectory")
    if snr_db is None:
        return x
    sigma = np.sqrt(np.var(x) / 10.0 ** (float(snr_db) / 10.0))
    noise_rng = np.random.default_rng(
        np.random.SeedSequence([int(noise_ss.generate_state(1)[0]), int(round(snr_db * 100)) + 100000]))
    return x + sigma * noise_rng.standard_normal(len(x))


# ---------------------------------------------------------------------------
# Analytic references
# ---------------------------------------------------------------------------

def henon_fixed_points(a=C.HENON_A, b=C.HENON_B):
    """x* = [-(1-b) +- sqrt((1-b)^2 + 4a)] / (2a)."""
    s = np.sqrt((1 - b) ** 2 + 4 * a)
    return ((-(1 - b) + s) / (2 * a), (-(1 - b) - s) / (2 * a))


def henon_period2(a=C.HENON_A, b=C.HENON_B):
    """2-cycle x = [(1-b) +- sqrt(4a - 3(1-b)^2)] / (2a)."""
    s = np.sqrt(4 * a - 3 * (1 - b) ** 2)
    return ((1 - b + s) / (2 * a), (1 - b - s) / (2 * a))


def henon_delay_jacobian(x, a=C.HENON_A, b=C.HENON_B):
    """Jacobian of (x_n, x_{n-1}) -> (x_{n+1}, x_n), x_{n+1} = 1 - a x_n^2 + b x_{n-1}."""
    return np.array([[-2.0 * a * x, b], [1.0, 0.0]])


def skewed_henon_f(u, v):
    return SKEWED_B1 * np.exp(SKEWED_B2 * u) * (1.4 - u * u) + 0.3 * v


def skewed_henon_delay_jacobian(u):
    fu = SKEWED_B1 * np.exp(SKEWED_B2 * u) * (SKEWED_B2 * (1.4 - u * u) - 2.0 * u)
    return np.array([[fu, 0.3], [1.0, 0.0]])


def skewed_henon_fixed_points():
    roots = []
    for u in (0.9302, -1.9030):          # PRE Table I starting values
        for _ in range(60):
            g = skewed_henon_f(u, u) - u
            u -= g / (skewed_henon_delay_jacobian(u)[0, 0] + 0.3 - 1.0)
        roots.append(float(u))
    return tuple(roots)


def skewed_henon_period2():
    a, b = 1.3014, -0.5604               # PRE Table I starting values
    for _ in range(80):
        F = np.array([skewed_henon_f(a, b) - b, skewed_henon_f(b, a) - a])
        ja = skewed_henon_delay_jacobian(a)[0, 0]
        jb = skewed_henon_delay_jacobian(b)[0, 0]
        a, b = np.array([a, b]) - np.linalg.solve(np.array([[ja, -0.7], [-0.7, jb]]), F)
    return float(a), float(b)


REFERENCES = {
    "henon": {"fixed_points": henon_fixed_points(), "period2": henon_period2(),
              "jacobian": henon_delay_jacobian},
    "skewed_henon": {"fixed_points": skewed_henon_fixed_points(),
                     "period2": skewed_henon_period2(),
                     "jacobian": skewed_henon_delay_jacobian},
}
