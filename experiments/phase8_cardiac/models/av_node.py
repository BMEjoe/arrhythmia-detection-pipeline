"""
Atrioventricular (AV) nodal conduction model of Sun, Amellal, Glass & Billette (1995).

SOURCE: Sun J, Amellal F, Glass L, Billette J. Alternans and period-doubling
bifurcations in atrioventricular nodal conduction. J Theor Biol 173:79-91, 1995
(Elsevier, not accessible here).  Equations and parameters as restated in
Zhao X, Schaeffer DG, "Alternate pacing of border-collision period-doubling
bifurcations", arXiv:math/0609106v2, Eq. (51):

  R_{n+1} = R_n exp(-(A_n + H)/tau_fat) + gamma exp(-H/tau_fat)
  A_{n+1} = A_min + R_{n+1} + (201 - 0.7 A_n) exp(-H/tau_rec)     if A_n <= 130
          = A_min + R_{n+1} + (500 - 3.0 A_n) exp(-H/tau_rec)     if A_n >= 130
  A_min = 33 ms, tau_rec = 70 ms, tau_fat = 30000 ms, gamma = 0.3 ms.
  A = atrial-His (AV nodal conduction) interval, R = drift (fatigue), H = interval from
  His activation to the next atrial activation (bifurcation parameter, ms).
Published result used for verification: border-collision period doubling at
H_bif = 56.9078 ms with A = 130 at the bifurcation point.

Beat events: His-bundle (ventricular) activations.  The His-to-His interval following
beat n is H + A_{n+1} (His_n -> atrial activation after H -> His after A_{n+1}); RR =
His-to-His intervals in seconds (the model has physiological units; no rescaling).
"""
from __future__ import annotations

import numpy as np

AMIN, TAU_REC, TAU_FAT, GAMMA = 33.0, 70.0, 30000.0, 0.3


def step(A, R, H):
    Rn = R * np.exp(-(A + H) / TAU_FAT) + GAMMA * np.exp(-H / TAU_FAT)
    if A <= 130.0:
        An = AMIN + Rn + (201.0 - 0.7 * A) * np.exp(-H / TAU_REC)
    else:
        An = AMIN + Rn + (500.0 - 3.0 * A) * np.exp(-H / TAU_REC)
    return An, Rn


def jacobian(A, R, H):
    e_f = np.exp(-(A + H) / TAU_FAT)
    dRn_dA, dRn_dR = -R / TAU_FAT * e_f, e_f
    slope = -0.7 if A <= 130.0 else -3.0
    er = np.exp(-H / TAU_REC)
    return np.array([[dRn_dA + slope * er, dRn_dR], [dRn_dA, dRn_dR]])   # rows: A_{n+1}, R_{n+1}


def iterate(H, n, A0=120.0, R0=0.0, transient=5000, noise_sd=0.0, rng=None):
    """Returns A_n, R_n after the transient.  Optional dynamical noise: Gaussian noise
    (ms) added to H at every beat (beat-to-beat variability of the atrial timing)."""
    A, R = A0, R0
    out = np.empty((n, 2))
    Hs = np.full(n + transient, H, dtype=float)
    if noise_sd > 0:
        Hs = Hs + rng.normal(0.0, noise_sd, n + transient)
    for k in range(n + transient):
        A, R = step(A, R, Hs[k])
        if k >= transient:
            out[k - transient] = (A, R)
    return out, Hs[transient:]


def lyapunov(H, n=200_000, transient=20_000, n_blocks=20, A0=120.0):
    """Ground-truth largest LE (per beat) of the 2-D map by QR (Benettin) on the exact
    Jacobian; returns (lambda1, SE from batch means, convergence curve)."""
    A, R = A0, 0.0
    for _ in range(transient):
        A, R = step(A, R, H)
    Q = np.eye(2)
    logs = np.empty(n)
    for k in range(n):
        J = jacobian(A, R, H)
        Q, Rm = np.linalg.qr(J @ Q)
        logs[k] = np.log(abs(Rm[0, 0]) + 1e-300)
        A, R = step(A, R, H)
    blocks = logs[: n - n % n_blocks].reshape(n_blocks, -1).mean(1)
    conv = {int(k): float(logs[:k].mean()) for k in (n // 16, n // 8, n // 4, n // 2, n)}
    return float(logs.mean()), float(blocks.std(ddof=1) / np.sqrt(n_blocks)), conv


def simulate_rr(H, n_beats, noise_sd=0.0, rng=None, A0=120.0):
    out, Hs = iterate(H, n_beats + 1, A0=A0, noise_sd=noise_sd, rng=rng)
    return (Hs[1:n_beats + 1] + out[1:n_beats + 1, 0]) / 1000.0
