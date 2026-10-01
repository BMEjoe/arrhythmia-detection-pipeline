"""
Three coupled modified van der Pol (Grudzinski-Zebrowski) oscillators: SA node, AV node,
His-Purkinje complex (HP) with time-delayed couplings (Gois & Savi 2009 lineage, with the
coupling form and parameters of Cheffer et al.).

SOURCE (primary for this formulation, open access): da Silva Lima G, Savi MA, Bessa WM,
"Adaptive control of cardiac rhythms", Sci Rep 14:23446, 2024 (PMC11458860):
  Eq. (1)  u'' + alpha u' (u - nu1)(u - nu2) + u (u + d)(u + e)/(d e) = F(t)
  Eqs. (2)-(4) the three nodes; with the unidirectional couplings used "for all
  simulations" (normal conduction SA -> AV -> HP; all other coupling parameters null):
    u_SA'' = F_SA(t) - a_SA u_SA'(u_SA - n_SA1)(u_SA - n_SA2) - u_SA(u_SA + d_SA)(u_SA + e_SA)/(d_SA e_SA)
    u_AV'' = - a_AV u_AV'(u_AV - n_AV1)(u_AV - n_AV2) - u_AV(u_AV + d_AV)(u_AV + e_AV)/(d_AV e_AV)
             - k_SA-AV u_AV + kt_SA-AV u_SA(t - tau_SA-AV)
    u_HP'' = - a_HP u_HP'(u_HP - n_HP1)(u_HP - n_HP2) - u_HP(u_HP + d_HP)(u_HP + e_HP)/(d_HP e_HP)
             - k_AV-HP u_HP + kt_AV-HP u_AV(t - tau_AV-HP)
  Eq. (5)  ECG = beta0 + beta1 u_SA + beta2 u_AV + beta3 u_HP, beta = (1, 0.06, 0.1, 0.3)
  F_SA(t) = rho_SA sin(omega_SA t) (external stimulus, "pathological" rhythms)
  Fig. 1b parameter table (read from the figure image); RK4 with 1 kHz sampling;
  initial conditions u0 = (-0.1, -0.6, -3.3), u0' = (0.025, 0.1, 2/3); the dimensional time
  is scaled so that the normal ECG has ~90 bpm, "beta_t = 0.1048".
Published results used for verification: the unforced model gives a regular (normal)
rhythm at ~90 bpm under the stated time scaling; the four stimulus settings
(rho_SA, omega_SA) = (5.45, 5.6), (8.625, 2.1), (8, 3.3), (9.6, 2.1) give high-rate,
"highly dispersed" (irregular) pathological rhythms (their Figs 3, 4, 8).

Time scaling: real time = model time * beta_t is assumed (beta_t = 0.1048 is the stated
scaling; the direction is checked against the ~90 bpm normal rhythm in the verification).
Beat events (predeclared): ventricular (HP) activations = upward crossings of u_HP through
the midpoint between its 5th and 95th percentile over the simulated segment (linear
interpolation between grid points).  RR in seconds after the time scaling.

Ground-truth LE: RK4 tangent (linearized) DDE integrated alongside; tangent state = current
6-D vector plus the history of du_SA and du_AV over the delays; renormalized every renorm
model time units (Farmer-type norm).  The model is non-autonomous when forced (F_SA(t)); the
LE is that of the time-dependent flow.
"""
from __future__ import annotations

import numba
import numpy as np

FIG1B = dict(a_SA=3.0, n_SA1=1.0, n_SA2=-1.9, d_SA=1.9, e_SA=0.55,
             a_AV=3.0, n_AV1=0.5, n_AV2=-0.5, d_AV=4.0, e_AV=0.67,
             a_HP=7.0, n_HP1=1.65, n_HP2=-2.0, d_HP=7.0, e_HP=0.67,
             k_SA_AV=3.0, kt_SA_AV=3.0, k_AV_HP=55.0, kt_AV_HP=55.0,
             tau_SA_AV=0.8, tau_AV_HP=0.1, rho=0.0, omega=0.0)
KEYS = list(FIG1B)
BETA_T = 0.1048
ECG_BETA = (1.0, 0.06, 0.1, 0.3)
U0 = np.array([-0.1, 0.025, -0.6, 0.1, -3.3, 2.0 / 3.0])   # (u_SA, u_SA', u_AV, u_AV', u_HP, u_HP')


def pack(P):
    return np.array([float(P[k]) for k in KEYS])


@numba.njit(cache=True)
def _f(t, y, dSA, dAV, q):
    a1, n11, n12, d1, e1, a2, n21, n22, d2, e2, a3, n31, n32, d3, e3, k12, kt12, k23, kt23, tau1, tau2, rho, om = q
    u1, w1, u2, w2, u3, w3 = y[0], y[1], y[2], y[3], y[4], y[5]
    out = np.empty(6)
    out[0] = w1
    out[1] = rho * np.sin(om * t) - a1 * w1 * (u1 - n11) * (u1 - n12) - u1 * (u1 + d1) * (u1 + e1) / (d1 * e1)
    out[2] = w2
    out[3] = -a2 * w2 * (u2 - n21) * (u2 - n22) - u2 * (u2 + d2) * (u2 + e2) / (d2 * e2) - k12 * u2 + kt12 * dSA
    out[4] = w3
    out[5] = -a3 * w3 * (u3 - n31) * (u3 - n32) - u3 * (u3 + d3) * (u3 + e3) / (d3 * e3) - k23 * u3 + kt23 * dAV
    return out


@numba.njit(cache=True)
def _fv(y, v, dvSA, dvAV, q):
    a1, n11, n12, d1, e1, a2, n21, n22, d2, e2, a3, n31, n32, d3, e3, k12, kt12, k23, kt23, tau1, tau2, rho, om = q
    u1, w1, u2, w2, u3, w3 = y[0], y[1], y[2], y[3], y[4], y[5]
    out = np.empty(6)
    out[0] = v[1]
    out[1] = (-a1 * w1 * (2 * u1 - n11 - n12) - (3 * u1 * u1 + 2 * (d1 + e1) * u1 + d1 * e1) / (d1 * e1)) * v[0] \
        - a1 * (u1 - n11) * (u1 - n12) * v[1]
    out[2] = v[3]
    out[3] = (-a2 * w2 * (2 * u2 - n21 - n22) - (3 * u2 * u2 + 2 * (d2 + e2) * u2 + d2 * e2) / (d2 * e2) - k12) * v[2] \
        - a2 * (u2 - n21) * (u2 - n22) * v[3] + kt12 * dvSA
    out[4] = v[5]
    out[5] = (-a3 * w3 * (2 * u3 - n31 - n32) - (3 * u3 * u3 + 2 * (d3 + e3) * u3 + d3 * e3) / (d3 * e3) - k23) * v[4] \
        - a3 * (u3 - n31) * (u3 - n32) * v[5] + kt23 * dvAV
    return out


@numba.njit(cache=True)
def integrate(y0, n_steps, dt, q, l1, l2, t0, tangent, v0, renorm_steps, noise_sd, noise):
    """RK4 on the grid; constant initial history y0.  l1, l2 = delays in steps (>= 1).
    Dynamical noise (optional): additive Euler-Maruyama term sqrt(dt)*noise_sd*xi on the
    three velocity equations (u'').  Returns (trajectory, log growth factors)."""
    off = max(l1, l2) + 1
    X = np.empty((off + n_steps + 1, 6))
    V = np.zeros((off + n_steps + 1, 6))
    for j in range(off + 1):
        X[j, :] = y0
        if tangent:
            V[j, :] = v0
    n_ren = n_steps // renorm_steps if tangent else 1
    logs = np.zeros(max(n_ren, 1))
    r = 0
    lmax = max(l1, l2)
    for i in range(off, off + n_steps):
        t = t0 + (i - off) * dt
        y = X[i]
        a1 = X[i - l1, 0]; b1 = X[i - l1 + 1, 0]
        a2 = X[i - l2, 2]; b2 = X[i - l2 + 1, 2]
        h1 = 0.5 * (a1 + b1); h2 = 0.5 * (a2 + b2)
        k1 = _f(t, y, a1, a2, q)
        ya = y + 0.5 * dt * k1
        k2 = _f(t + 0.5 * dt, ya, h1, h2, q)
        yb = y + 0.5 * dt * k2
        k3 = _f(t + 0.5 * dt, yb, h1, h2, q)
        yc = y + dt * k3
        k4 = _f(t + dt, yc, b1, b2, q)
        X[i + 1] = y + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        if noise_sd > 0.0:
            sq = np.sqrt(dt) * noise_sd
            X[i + 1, 1] += sq * noise[i - off, 0]
            X[i + 1, 3] += sq * noise[i - off, 1]
            X[i + 1, 5] += sq * noise[i - off, 2]
        if tangent:
            v = V[i]
            va1 = V[i - l1, 0]; vb1 = V[i - l1 + 1, 0]; va2 = V[i - l2, 2]; vb2 = V[i - l2 + 1, 2]
            vh1 = 0.5 * (va1 + vb1); vh2 = 0.5 * (va2 + vb2)
            m1 = _fv(y, v, va1, va2, q)
            m2 = _fv(ya, v + 0.5 * dt * m1, vh1, vh2, q)
            m3 = _fv(yb, v + 0.5 * dt * m2, vh1, vh2, q)
            m4 = _fv(yc, v + dt * m3, vb1, vb2, q)
            V[i + 1] = v + dt / 6.0 * (m1 + 2 * m2 + 2 * m3 + m4)
            if (i - off + 1) % renorm_steps == 0 and r < n_ren:
                s = 0.0
                for c in range(6):
                    s += V[i + 1, c] ** 2
                for j in range(i + 1 - lmax, i + 1):
                    s += (V[j, 0] ** 2 + V[j, 2] ** 2) * dt
                nrm = np.sqrt(s)
                logs[r] = np.log(nrm)
                r += 1
                for j in range(i + 1 - lmax, i + 2):
                    for c in range(6):
                        V[j, c] /= nrm
    return X[off:, :], logs[:r]


def _lags(P, dt):
    return max(1, int(round(P["tau_SA_AV"] / dt))), max(1, int(round(P["tau_AV_HP"] / dt)))


def run(P=None, model_time=500.0, transient=200.0, dt=0.001, noise_sd=0.0, rng=None, y0=None):
    """Trajectory after the transient (published initial conditions by default)."""
    P = dict(FIG1B, **(P or {}))
    q = pack(P)
    l1, l2 = _lags(P, dt)
    n = int(round((transient + model_time) / dt))
    noise = rng.normal(size=(n, 3)) if noise_sd > 0 else np.zeros((1, 3))
    X, _ = integrate(U0.copy() if y0 is None else y0, n, dt, q, l1, l2, 0.0, False, np.zeros(6), 1,
                     noise_sd, noise)
    return X[int(round(transient / dt)):], dt


def lyapunov(P=None, model_time=3000.0, transient=500.0, dt=0.001, renorm=1.0, seed=0, n_blocks=20):
    P = dict(FIG1B, **(P or {}))
    q = pack(P)
    l1, l2 = _lags(P, dt)
    rng = np.random.default_rng(seed)
    nt = int(round(transient / dt))
    Xw, _ = integrate(U0.copy(), nt, dt, q, l1, l2, 0.0, False, np.zeros(6), 1, 0.0, np.zeros((1, 3)))
    v0 = rng.normal(size=6)
    v0 /= np.linalg.norm(v0)
    rs = int(round(renorm / dt))
    # continue from the end state (history = constant end state; the first 10 renorm
    # intervals are discarded to forget it)
    X, logs = integrate(Xw[-1].copy(), int(round(model_time / dt)), dt, q, l1, l2, transient, True, v0, rs,
                        0.0, np.zeros((1, 3)))
    logs = logs[10:] / renorm
    m = len(logs)
    blocks = logs[: m - m % n_blocks].reshape(n_blocks, -1).mean(1)
    conv = {int(k * renorm): float(logs[:k].mean()) for k in (m // 8, m // 4, m // 2, m)}
    return float(logs.mean()), float(blocks.std(ddof=1) / np.sqrt(n_blocks)), conv, X


def ecg(X):
    b0, b1, b2, b3 = ECG_BETA
    return b0 + b1 * X[:, 0] + b2 * X[:, 2] + b3 * X[:, 4]


def beat_times(X, dt):
    """Upward crossings of u_HP through the midpoint of its 5th-95th percentile range."""
    u = X[:, 4]
    lo, hi = np.percentile(u, [5, 95])
    lev = 0.5 * (lo + hi)
    i = np.flatnonzero((u[:-1] < lev) & (u[1:] >= lev))
    frac = (lev - u[i]) / (u[i + 1] - u[i])
    return (i + frac) * dt


def simulate_rr(P, n_beats, rng, dt=0.001, noise_sd=0.0, transient=200.0):
    """RR (s) of ventricular events; model time scaled by BETA_T."""
    span = 50.0
    while True:
        X, _ = run(P, model_time=span, transient=transient, dt=dt, noise_sd=noise_sd, rng=rng,
                   y0=U0 + 0.001 * rng.normal(size=6))
        tb = beat_times(X, dt)
        if len(tb) > n_beats + 2:
            return np.diff(tb)[:n_beats] * BETA_T
        span *= 2.0
