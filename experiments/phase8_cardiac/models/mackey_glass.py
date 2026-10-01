"""
Mackey-Glass delay-differential equation (physiological delayed feedback).

SOURCE: Mackey MC, Glass L. Oscillation and chaos in physiological control systems.
Science 197:287-289, 1977 (not accessible here); equation as given by the original
authors in the open article Glass L, Mackey MC, "Mackey-Glass equation", Scholarpedia
5(3):6908, 2010, Eq. (1):
    dx/dt = beta * x_tau / (1 + x_tau^n) - gamma * x,     x_tau = x(t - tau)
with the standard chaotic parameter set beta = 0.2, gamma = 0.1, n = 10 (as Phase 5 G3,
which cites Farmer 1982).  Published result used for verification: the largest
Lyapunov exponent at tau = 17 is 0.0086 per time unit (Farmer JD, Physica D 4:366,
1982, Table; as quoted in Table 3 of arXiv:1509.06057).

Integration: RK4 with step dt on a uniform grid, delayed values from the grid (tau/dt
integer; half-step delayed values by linear interpolation, as Phase 5).
Ground-truth LE: the linearized DDE  d(dx)/dt = -gamma dx + g'(x_tau) dx_tau  is
integrated alongside with the same scheme; the tangent state is the whole history
segment on [t - tau, t] (the phase space of a DDE), renormalized every `renorm`
time units (Farmer's method for delay equations).

Beat events (predeclared; the model has no heartbeat): successive local maxima of x(t)
(interpolated by a parabola through the grid maximum and its neighbours).  RR = the
intervals between successive maxima, then rescaled (the model lacks physiological
units) to mean 0.8 s and the window's own SD fixed at 0.05 s by a linear map.
"""
from __future__ import annotations

import numba
import numpy as np

BETA, GAMMA, NEXP = 0.2, 0.1, 10.0


@numba.njit(cache=True)
def _g(x, beta, n):
    return beta * x / (1.0 + x ** n)


@numba.njit(cache=True)
def _dg(x, beta, n):
    xn = x ** n
    return beta * (1.0 + xn - n * xn) / (1.0 + xn) ** 2


@numba.njit(cache=True)
def integrate(x0hist, n_steps, dt, lag, beta, gamma, n, noise_sd, noise):
    """x0hist: initial history (lag+1 values).  Returns x on the grid (len lag+1+n_steps).
    noise: pre-drawn N(0,1) increments; dynamical noise enters as sqrt(dt)*noise_sd*xi
    (Euler-Maruyama additive term after each RK4 step)."""
    x = np.empty(lag + 1 + n_steps)
    x[:lag + 1] = x0hist
    for i in range(lag, lag + n_steps):
        xd0 = x[i - lag]
        xd1 = x[i - lag + 1]
        xdh = 0.5 * (xd0 + xd1)
        xi = x[i]
        k1 = _g(xd0, beta, n) - gamma * xi
        k2 = _g(xdh, beta, n) - gamma * (xi + 0.5 * dt * k1)
        k3 = _g(xdh, beta, n) - gamma * (xi + 0.5 * dt * k2)
        k4 = _g(xd1, beta, n) - gamma * (xi + dt * k3)
        x[i + 1] = xi + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4) + np.sqrt(dt) * noise_sd * noise[i - lag]
    return x


@numba.njit(cache=True)
def _lyap(x0hist, n_steps, dt, lag, beta, gamma, n, renorm_steps, v0):
    x = np.empty(lag + 1 + n_steps)
    v = np.empty(lag + 1 + n_steps)
    x[:lag + 1] = x0hist
    v[:lag + 1] = v0
    n_ren = n_steps // renorm_steps
    logs = np.empty(n_ren)
    r = 0
    for i in range(lag, lag + n_steps):
        xd0 = x[i - lag]; xd1 = x[i - lag + 1]; xdh = 0.5 * (xd0 + xd1); xi = x[i]
        k1 = _g(xd0, beta, n) - gamma * xi
        k2 = _g(xdh, beta, n) - gamma * (xi + 0.5 * dt * k1)
        k3 = _g(xdh, beta, n) - gamma * (xi + 0.5 * dt * k2)
        k4 = _g(xd1, beta, n) - gamma * (xi + dt * k3)
        x[i + 1] = xi + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        vd0 = v[i - lag]; vd1 = v[i - lag + 1]; vdh = 0.5 * (vd0 + vd1); vi = v[i]
        a0 = _dg(xd0, beta, n); a1 = _dg(xd1, beta, n); ah = _dg(xdh, beta, n)
        l1 = a0 * vd0 - gamma * vi
        l2 = ah * vdh - gamma * (vi + 0.5 * dt * l1)
        l3 = ah * vdh - gamma * (vi + 0.5 * dt * l2)
        l4 = a1 * vd1 - gamma * (vi + dt * l3)
        v[i + 1] = vi + dt / 6.0 * (l1 + 2 * l2 + 2 * l3 + l4)
        if (i - lag + 1) % renorm_steps == 0:
            s = 0.0
            for j in range(i + 1 - lag, i + 2):
                s += v[j] * v[j]
            nrm = np.sqrt(s * dt)
            logs[r] = np.log(nrm)
            r += 1
            for j in range(i + 1 - lag, i + 2):
                v[j] /= nrm
    return logs


def lyapunov(tau, total_time=2.0e5, transient=2000.0, dt=0.05, renorm=10.0, seed=0, n_blocks=20,
             beta=BETA, gamma=GAMMA, n=NEXP):
    """Largest LE per time unit; returns (lambda, SE from batch means, convergence curve)."""
    lag = int(round(tau / dt))
    rng = np.random.default_rng(seed)
    x0 = np.full(lag + 1, rng.uniform(0.5, 1.3))
    xw = integrate(x0, int(transient / dt), dt, lag, beta, gamma, n, 0.0, np.zeros(int(transient / dt)))
    hist = xw[-(lag + 1):].copy()
    v0 = rng.normal(size=lag + 1)
    v0 /= np.sqrt(np.sum(v0 ** 2) * dt)
    rs = int(round(renorm / dt))
    logs = _lyap(hist, int(total_time / dt), dt, lag, beta, gamma, n, rs, v0) / renorm
    m = len(logs)
    blocks = logs[: m - m % n_blocks].reshape(n_blocks, -1).mean(1)
    conv = {int(k * renorm): float(logs[:k].mean()) for k in (m // 16, m // 8, m // 4, m // 2, m)}
    return float(logs.mean()), float(blocks.std(ddof=1) / np.sqrt(n_blocks)), conv


def maxima_times(x, dt):
    i = np.flatnonzero((x[1:-1] > x[:-2]) & (x[1:-1] >= x[2:])) + 1
    y0, y1, y2 = x[i - 1], x[i], x[i + 1]
    den = y0 - 2 * y1 + y2
    off = np.where(den != 0, 0.5 * (y0 - y2) / np.where(den != 0, den, 1.0), 0.0)
    return (i + off) * dt


def simulate_rr(tau, n_beats, rng, dt=0.05, transient=2000.0, noise_sd=0.0, mean_rr=0.8, sd_rr=0.05,
                beta=BETA, gamma=GAMMA, n=NEXP, raw=False):
    lag = int(round(tau / dt))
    x0 = np.full(lag + 1, rng.uniform(0.5, 1.3))
    span = transient + (n_beats + 50) * (tau + 30.0)
    steps = int(span / dt)
    noise = rng.normal(size=steps) if noise_sd > 0 else np.zeros(steps)
    x = integrate(x0, steps, dt, lag, beta, gamma, n, noise_sd, noise)
    x = x[lag + 1 + int(transient / dt):]
    tm = maxima_times(x, dt)
    iv = np.diff(tm)
    if len(iv) < n_beats:
        raise RuntimeError("not enough maxima")
    iv = iv[:n_beats]
    if raw:
        return iv
    sd = iv.std()
    z = (iv - iv.mean()) / sd if sd > 0 else iv - iv.mean()
    return mean_rr + sd_rr * z
