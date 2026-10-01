"""
Periodically forced cardiac oscillator (phase-resetting map), Glass-Shrier lineage.

SOURCE (primary, open access): Diagne K, Bury TM, Pettebone ME, ..., Glass L et al.,
"Phase resetting in human stem cell derived cardiomyocytes explains complex cardiac
arrhythmias", PLoS Comput Biol 2026, doi:10.1371/journal.pcbi.1013935 (PMC12900431).
  - main text Eq. (1):  phi_{i+1} = f(phi_i) = phi_i + tau - g(phi_i)   (mod 1)
    phi_i = phase of the i-th stimulus in the oscillator's cycle, tau = T_stim / T
  - S1 Text Eq. (1)-(3): piecewise PRC g (normalized perturbed cycle length)
        g = 1                          0     <= phi < phi1
        g = 1 + A (phi - phi1)^4       phi1  <= phi < phi_r
        g = B (phi - phi2)^2 + C       phi_r <= phi < phi3
        g = 1 + S (phi - 1)            phi3  <= phi < 1
    phi1 = phi_r - 0.25, phi2 = phi_r + 0.04, phi3 = phi2 + S/(2B), C = 1 - S(1-phi3) - B(phi3-phi2)^2
  - S1 Table: fitted (A, phi_r, B, S) for six spheroid aggregates A-F; aggregate A is
    the PRC of Fig 1B used in the main text (intrinsic cycle 2.4 s).

Beat events (the model's own): action potentials (APs).  Stimulus i occurs at time
s_i = i*tau*T; the last AP before it is a_i = s_i - phi_i*T.  The stimulus makes the
current cycle last g(phi_i)*T, so an AP occurs at a_i + g(phi_i)*T, and the oscillator
then fires freely every T until stimulus i+1.  If stimulus i+1 arrives before that AP
(phi_i + tau - g(phi_i) < 0, the map's mod-1 wrap) the AP does not occur.  RR =
AP-to-AP intervals (simulate_rr).  By construction a_{i+1} = s_{i+1} - phi_{i+1}*T, so
the AP sequence is consistent with Eq. (1).
"""
from __future__ import annotations

import numpy as np

PRC_TABLE = {   # S1 Table: aggregate: (A, phi_r, B, S, iCL range s)
    "A": (60.0, 0.71, 50.0, 0.80, (2.07, 2.41)),
    "B": (40.0, 0.57, 10.0, 0.93, (1.99, 3.07)),
    "C": (120.0, 0.57, 0.73, 0.80, (0.61, 0.80)),
    "D": (20.0, 0.42, 10.0, 0.93, (1.13, 1.29)),
    "E": (25.0, 0.40, 10.0, 0.95, (1.32, 1.37)),
    "F": (30.0, 0.38, 20.0, 0.95, (1.36, 1.41)),
}


def prc_params(aggregate="A"):
    A, phr, B, S, _ = PRC_TABLE[aggregate]
    ph1, ph2 = phr - 0.25, phr + 0.04
    ph3 = ph2 + S / (2.0 * B)
    C = 1.0 - S * (1.0 - ph3) - B * (ph3 - ph2) ** 2
    return dict(A=A, phr=phr, B=B, S=S, ph1=ph1, ph2=ph2, ph3=ph3, C=C)


def g(phi, p):
    phi = np.asarray(phi, dtype=float)
    out = np.ones_like(phi)
    m2 = (phi >= p["ph1"]) & (phi < p["phr"])
    m3 = (phi >= p["phr"]) & (phi < p["ph3"])
    m4 = phi >= p["ph3"]
    out[m2] = 1.0 + p["A"] * (phi[m2] - p["ph1"]) ** 4
    out[m3] = p["B"] * (phi[m3] - p["ph2"]) ** 2 + p["C"]
    out[m4] = 1.0 + p["S"] * (phi[m4] - 1.0)
    return out


def dg(phi, p):
    phi = np.asarray(phi, dtype=float)
    out = np.zeros_like(phi)
    m2 = (phi >= p["ph1"]) & (phi < p["phr"])
    m3 = (phi >= p["phr"]) & (phi < p["ph3"])
    m4 = phi >= p["ph3"]
    out[m2] = 4.0 * p["A"] * (phi[m2] - p["ph1"]) ** 3
    out[m3] = 2.0 * p["B"] * (phi[m3] - p["ph2"])
    out[m4] = p["S"]
    return out


def iterate(phi0, tau, n, p, noise_sd=0.0, rng=None):
    """Eq. (1).  Optional dynamical noise: Gaussian phase noise added to each iterate (mod 1)."""
    phi = np.empty(n)
    x = float(phi0)
    for i in range(n):
        phi[i] = x
        x = (x + tau - float(g(np.array([x]), p)[0])) % 1.0
        if noise_sd > 0:
            x = (x + rng.normal(0.0, noise_sd)) % 1.0
    return phi


def lyapunov(tau, p, n=200_000, transient=10_000, phi0=0.123, n_blocks=20):
    """Ground-truth LE of the 1-D map: mean log|f'(phi_i)|, f' = 1 - g'(phi) (f is
    piecewise smooth; the discontinuity at phi_r has measure zero).  Returns
    (lambda, standard error from n_blocks batch means, convergence curve)."""
    phi = iterate(phi0, tau, n + transient, p)[transient:]
    d = np.log(np.abs(1.0 - dg(phi, p)) + 1e-300)
    blocks = d[: n - n % n_blocks].reshape(n_blocks, -1).mean(1)
    lam = float(d.mean())
    conv = {int(k): float(d[:k].mean()) for k in (n // 16, n // 8, n // 4, n // 2, n)}
    return lam, float(blocks.std(ddof=1) / np.sqrt(n_blocks)), conv


def simulate_rr(tau, n_beats, p, T=1.0, phi0=0.123, transient=2000, noise_sd=0.0, rng=None):
    """RR (AP-to-AP) intervals in units of T*seconds-per-unit, from the model's own AP events.
    Iterates Eq. (1) (with optional dynamical phase noise), then converts to AP times."""
    n_stim = transient + 4 * n_beats + 100
    while True:
        phi = iterate(phi0, tau, n_stim, p, noise_sd=noise_sd, rng=rng)
        gp = g(phi, p)
        aps = []
        for i in range(transient, n_stim - 1):
            s_i = i * tau
            t = s_i - phi[i] + gp[i]
            s_next = (i + 1) * tau
            while t < s_next:
                aps.append(t)
                t += 1.0
        if len(aps) > n_beats + 1:
            aps = np.asarray(aps[: n_beats + 1]) * T
            return np.diff(aps), phi[transient:]
        n_stim *= 2
