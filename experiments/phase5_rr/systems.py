"""
Phase 5 RR-interval-like generators (parameters in config.py; sources in README.md).

generate(condition, seed) -> RR series in seconds, length config.N.
parts(condition, seed)    -> dict with the RR scale, the drawn parameters and
                             the unquantized series (used by realism.py).
"""
from __future__ import annotations

import numpy as np

from experiments.phase2e import systems as S2
from experiments.phase5_rr import config as C


def _ss(code, seed):
    return np.random.SeedSequence([C.ENTROPY, int(code), C.N, int(seed)])


def _standardize(z):
    z = np.asarray(z, dtype=float)
    return (z - z.mean()) / z.std()


def _segment(x, n, rng):
    start = int(rng.integers(0, len(x) - n + 1))
    return x[start:start + n]


# ---------------------------------------------------------------------------
# Base processes.  Each returns (z, params); z is later standardized and
# scaled to the window's RR mean and SD.  `mean_rr` converts Hz to cycles/beat.
# ---------------------------------------------------------------------------

def linear_rr(n, rng, mean_rr):
    """McSharry et al. (2003) RR spectrum: two Gaussian bumps (LF 0.1 Hz, HF 0.25 Hz,
    SD 0.01 Hz), amplitudes sqrt(S(f)), phases U(0, 2pi); beat index k at time k * mean_rr."""
    lfhf = rng.uniform(*C.LFHF_RANGE)
    L = C.LINEAR_SYNTH_LENGTH
    f_beat = np.fft.rfftfreq(L)                      # cycles / beat
    f = f_beat / mean_rr                             # Hz
    s = (lfhf * np.exp(-0.5 * ((f - C.LF_CENTER_HZ) / C.LF_WIDTH_HZ) ** 2) / C.LF_WIDTH_HZ
         + np.exp(-0.5 * ((f - C.HF_CENTER_HZ) / C.HF_WIDTH_HZ) ** 2) / C.HF_WIDTH_HZ)
    s[0] = 0.0
    phases = rng.uniform(0.0, 2.0 * np.pi, len(f))
    x = np.fft.irfft(np.sqrt(s) * np.exp(1j * phases), n=L)
    return _segment(x, n, rng), {"lf_hf": float(lfhf)}


def power_law(n, rng, mean_rr):
    """Gaussian 1/f^beta noise: complex Gaussian Fourier coefficients with
    |a(f)|^2 ~ f^-beta over 4096 beats; a random 256-beat segment is kept."""
    beta = rng.uniform(*C.BETA_RANGE)
    L = C.POWER_LAW_SYNTH_LENGTH
    f = np.fft.rfftfreq(L)
    amp = np.zeros_like(f)
    amp[1:] = f[1:] ** (-beta / 2.0)
    coef = amp * (rng.standard_normal(len(f)) + 1j * rng.standard_normal(len(f)))
    x = np.fft.irfft(coef, n=L)
    return _segment(x, n, rng), {"beta": float(beta)}


def noisy_rsa(n, rng, mean_rr):
    """Respiratory sinusoid (f_resp ~ U(0.2, 0.3) Hz, random phase) plus independent
    white Gaussian noise; the sinusoid carries RSA_VARIANCE_SHARE of the variance."""
    fr = rng.uniform(*C.RESP_FREQ_RANGE_HZ)
    ph = rng.uniform(0.0, 2.0 * np.pi)
    k = np.arange(n)
    s = np.sqrt(2.0) * np.sin(2.0 * np.pi * fr * mean_rr * k + ph)          # unit variance
    w = rng.standard_normal(n)
    a = C.RSA_VARIANCE_SHARE
    return np.sqrt(a) * s + np.sqrt(1.0 - a) * w, {"resp_hz": float(fr)}


def setar(n, rng, mean_rr):
    """SETAR(2;1,1): x_t = 0.7 x_{t-1} + e_t if x_{t-1} <= 0, else -0.5 x_{t-1} + e_t;
    e_t ~ N(0, 1).  Deterministic skeleton has a stable fixed point at 0."""
    total = n + C.TRANSIENT
    e = rng.standard_normal(total)
    x = np.zeros(total)
    for t in range(1, total):
        phi = C.SETAR_PHI_LOW if x[t - 1] <= 0.0 else C.SETAR_PHI_HIGH
        x[t] = phi * x[t - 1] + e[t]
    return x[C.TRANSIENT:], {}


def henon(n, rng, mean_rr):
    return S2.henon(n, rng), {}


def logistic(n, rng, mean_rr):
    return S2.logistic(n, rng), {}


def _rk4(f, y, dt):
    k1 = f(y)
    k2 = f(y + 0.5 * dt * k1)
    k3 = f(y + 0.5 * dt * k2)
    k4 = f(y + dt * k3)
    return y + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)


def lorenz_maxima(n, rng, mean_rr):
    """Successive local maxima of Lorenz z (sigma 10, rho 28, beta 8/3), RK4 dt 0.005,
    100 time units discarded; each maximum refined by a parabola through 3 samples."""
    p = C.LORENZ
    s, r, b, dt = p["sigma"], p["rho"], p["beta"], p["dt"]

    def f(v):
        return np.array([s * (v[1] - v[0]), v[0] * (r - v[2]) - v[1], v[0] * v[1] - b * v[2]])

    v = np.array([rng.uniform(-10, 10), rng.uniform(-10, 10), rng.uniform(10, 40)])
    for _ in range(int(round(p["burn_in"] / dt))):
        v = _rk4(f, v, dt)
    out, z2, z1 = [], None, None
    while len(out) < n:
        v = _rk4(f, v, dt)
        z0 = v[2]
        if z2 is not None and z1 > z2 and z1 >= z0:
            den = z2 - 2 * z1 + z0
            out.append(z1 - 0.125 * (z0 - z2) ** 2 / den if den != 0 else z1)
        z2, z1 = z1, z0
    return np.asarray(out), {}


def rossler_flow(n, rng, mean_rr):
    """Rossler x (a = b = 0.2, c = 5.7), RK4 dt 0.01, 500 time units discarded,
    sampled every 1.0 time unit (about 6 samples per mean orbital period)."""
    p = C.ROSSLER
    a, b, c, dt = p["a"], p["b"], p["c"], p["dt"]

    def f(v):
        return np.array([-v[1] - v[2], v[0] + a * v[1], b + v[2] * (v[0] - c)])

    v = np.array([rng.uniform(-5, 5), rng.uniform(-5, 5), rng.uniform(0, 1)])
    for _ in range(int(round(p["burn_in"] / dt))):
        v = _rk4(f, v, dt)
    every = int(round(p["sample_step"] / dt))
    out = np.empty(n)
    for i in range(n):
        for _ in range(every):
            v = _rk4(f, v, dt)
        out[i] = v[0]
    return out, {}


def mackey_glass(n, rng, mean_rr):
    """Mackey-Glass dx/dt = 0.2 x(t-17) / (1 + x(t-17)^10) - 0.1 x, RK4 dt 0.1 with the
    delayed term linearly interpolated at half steps; constant initial history
    ~ U(0.5, 1.3); 1000 time units discarded; sampled every 6.0 time units."""
    p = C.MACKEY_GLASS
    beta, gamma, pw, dt = p["beta"], p["gamma"], p["n"], p["dt"]
    lag = int(round(p["tau"] / dt))
    every = int(round(p["sample_step"] / dt))
    steps = int(round(p["burn_in"] / dt)) + n * every
    x = np.empty(lag + steps + 1)
    x[:lag + 1] = rng.uniform(0.5, 1.3)

    def g(xd):
        return beta * xd / (1.0 + xd ** pw)

    for i in range(lag, lag + steps):
        xd0, xd1 = x[i - lag], x[i - lag + 1]
        xdh = 0.5 * (xd0 + xd1)
        xi = x[i]
        k1 = g(xd0) - gamma * xi
        k2 = g(xdh) - gamma * (xi + 0.5 * dt * k1)
        k3 = g(xdh) - gamma * (xi + 0.5 * dt * k2)
        k4 = g(xd1) - gamma * (xi + dt * k3)
        x[i + 1] = xi + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
    start = lag + int(round(p["burn_in"] / dt))
    return x[start + every: start + every * (n + 1): every][:n], {}


BASES = {"linear_rr": linear_rr, "power_law": power_law, "noisy_rsa": noisy_rsa, "setar": setar,
         "henon": henon, "logistic": logistic, "lorenz_maxima": lorenz_maxima,
         "rossler_flow": rossler_flow, "mackey_glass": mackey_glass}


# ---------------------------------------------------------------------------
# Modifiers (applied to the scaled RR series, in the listed order)
# ---------------------------------------------------------------------------

def add_trend(x, rng, mean_rr):
    frac = rng.uniform(*C.TREND_FRACTION_RANGE) * rng.choice((-1.0, 1.0))
    k = np.arange(len(x))
    return x + frac * mean_rr * (k / (len(x) - 1) - 0.5), {"trend_fraction": float(frac)}


def add_step(x, rng, mean_rr):
    frac = rng.uniform(*C.STEP_FRACTION_RANGE) * rng.choice((-1.0, 1.0))
    k0 = int(rng.integers(int(C.STEP_POSITION_RANGE[0] * len(x)), int(C.STEP_POSITION_RANGE[1] * len(x)) + 1))
    step = np.where(np.arange(len(x)) >= k0, 1.0, 0.0)
    return x + frac * mean_rr * (step - step.mean()), {"step_fraction": float(frac), "step_index": k0}


def add_ectopics(x, rng, rate):
    """Isolated premature beats with a full compensatory pause: at an ectopic index i,
    RR_i <- c RR_i and RR_{i+1} <- RR_{i+1} + (1 - c) RR_i, c ~ U(0.6, 0.8), so every
    later sinus beat keeps its time.  round(rate * N) ectopics at indices 1..N-2,
    at least ECTOPIC_MIN_SPACING apart."""
    x = np.array(x, dtype=float)
    k = int(round(rate * len(x)))
    chosen = []
    while len(chosen) < k:
        i = int(rng.integers(1, len(x) - 1))
        if all(abs(i - j) >= C.ECTOPIC_MIN_SPACING for j in chosen):
            chosen.append(i)
    chosen.sort()
    for i in chosen:
        c = rng.uniform(*C.ECTOPIC_PREMATURITY_RANGE)
        short = c * x[i]
        x[i + 1] += x[i] - short
        x[i] = short
    return x, {"ectopic_indices": chosen}


def add_noise(x, rng, snr_db, clean_var):
    sigma = np.sqrt(clean_var / 10.0 ** (float(snr_db) / 10.0))
    return x + sigma * rng.standard_normal(len(x)), {"snr_db": float(snr_db)}


def quantize(x):
    return np.round(np.asarray(x, dtype=float) / C.QUANTUM_S) * C.QUANTUM_S


def parts(condition, seed):
    grp, base, mods, quantized = C.CONDITIONS[condition]
    rng = np.random.default_rng(_ss(C.BASE_CODE[base], seed))
    mean_rr = float(rng.uniform(*C.MEAN_RR_RANGE_S))
    sd_rr = float(rng.uniform(*C.SD_RR_RANGE_S))
    z, params = BASES[base](C.N, rng, mean_rr)
    z = _standardize(z)
    mrng = np.random.default_rng(_ss(C.COND_CODE[C.MODIFIER_TWIN.get(condition, condition)], seed))
    for mod in mods:
        if mod[0] == "warp":
            z = _standardize(np.exp(C.WARP_A * z))
    x = mean_rr + sd_rr * z
    for mod in mods:
        if mod[0] == "trend":
            x, p = add_trend(x, mrng, mean_rr)
        elif mod[0] == "step":
            x, p = add_step(x, mrng, mean_rr)
        elif mod[0] == "ectopic":
            x, p = add_ectopics(x, mrng, mod[1])
        elif mod[0] == "noise":
            x, p = add_noise(x, mrng, mod[1], sd_rr ** 2)
        else:
            p = {"warp_a": C.WARP_A}
        params.update(p)
    if not np.all(np.isfinite(x)) or np.min(x) <= 0:
        raise FloatingPointError(f"{condition} seed {seed}: invalid RR series")
    return {"condition": condition, "seed": int(seed), "mean_rr": mean_rr, "sd_rr": sd_rr,
            "params": params, "unquantized": x, "rr": quantize(x) if quantized else x,
            "quantized": bool(quantized)}


def generate(condition, seed):
    return parts(condition, seed)["rr"]


def data_seed(condition, seed):
    base = C.CONDITIONS[condition][1]
    return int(_ss(C.BASE_CODE[base], seed).generate_state(1)[0])
