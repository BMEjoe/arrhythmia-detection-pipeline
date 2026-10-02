"""
Phase 9 DEVELOPMENT-ONLY control family: interbeat intervals of the sine circle map (added
2026-10-02 after the split was recorded, BEFORE any detector or measure was run on any Phase 9
series; it is a development control and never enters any TEST pool).

Why: the Phase 8 TEST disclosure (HANDOFF.md) says nonlinear prediction vs IAAFT failed on the
forced QUASI-PERIODIC coupled-vdP regimes; development had no quasi-periodic family, so a chaos
measure could not be checked against quasi-periodicity on development data.  The phase map of a
periodically forced (stimulated) oscillator is a circle map (Glass & Perez, Phys Rev Lett 49:1782
(1982); Guevara & Glass, J Math Biol 14:1 (1982)); the canonical form is the sine circle map
(Arnold 1965):
    theta_{n+1} = theta_n + Omega - (K / 2 pi) sin(2 pi theta_n)      (lift, not taken mod 1)
For K < 1 the map is an invertible circle diffeomorphism: no chaos; it is mode-locked (periodic,
negative Lyapunov exponent) or quasi-periodic (irrational rotation number, exponent 0).  The time
between successive beats in units of the forcing period is the lift increment
theta_{n+1} - theta_n = Omega - (K/2pi) sin(2 pi theta_n) > 0 (Omega > K/2pi), rescaled to mean 0.8 s
as the Phase 8 models (MODELS.md Sect. 6).
Lyapunov exponent: lambda = <ln |1 - K cos(2 pi theta_n)|> (per beat).

Regimes (predeclared from the scan in devmaps_scan(); labels by the Phase 8 rule with the
flow tolerance, since a quasi-periodic exponent is exactly 0):
    K = 0.5, Omega = 0.382  quasi-periodic        K = 0.9, Omega = 0.382  quasi-periodic
    K = 0.5, Omega = 0.618  quasi-periodic        K = 0.9, Omega = 0.618  quasi-periodic
    K = 0.9, Omega = 0.5    locked 1:2 (period 2) K = 0.9, Omega = 0.65   locked 2:3 (period 3)
"""
from __future__ import annotations

import numpy as np

from experiments.phase9_waveform import families as F

REGIMES = [
    {"name": "circle:K=0.5,Om=0.382", "K": 0.5, "Om": 0.382, "kind": "quasi-periodic"},
    {"name": "circle:K=0.9,Om=0.382", "K": 0.9, "Om": 0.382, "kind": "quasi-periodic"},
    {"name": "circle:K=0.5,Om=0.618", "K": 0.5, "Om": 0.618, "kind": "quasi-periodic"},
    {"name": "circle:K=0.9,Om=0.618", "K": 0.9, "Om": 0.618, "kind": "quasi-periodic"},
    {"name": "circle:K=0.9,Om=0.5", "K": 0.9, "Om": 0.5, "kind": "locked 1:2"},
    {"name": "circle:K=0.9,Om=0.65", "K": 0.9, "Om": 0.65, "kind": "locked 2:3"},
]
BY_NAME = {r["name"]: r for r in REGIMES}


def circle_rr(K, Om, n, rng, transient=1000):
    th = rng.uniform()
    for _ in range(transient):
        th = th + Om - K / (2 * np.pi) * np.sin(2 * np.pi * th)
    out = np.empty(n)
    for i in range(n):
        nxt = th + Om - K / (2 * np.pi) * np.sin(2 * np.pi * th)
        out[i] = nxt - th
        th = nxt
    return out * (F.RR_MEAN / out.mean())


def lyapunov(K, Om, n, seed, transient=1000):
    rng = np.random.default_rng(seed)
    th = rng.uniform()
    for _ in range(transient):
        th = th + Om - K / (2 * np.pi) * np.sin(2 * np.pi * th)
    v = np.empty(n)
    for i in range(n):
        v[i] = np.log(abs(1 - K * np.cos(2 * np.pi * th)))
        th = th + Om - K / (2 * np.pi) * np.sin(2 * np.pi * th)
    nb = 20
    bm = v[: n // nb * nb].reshape(nb, -1).mean(1)
    return float(v.mean()), float(bm.std(ddof=1) / np.sqrt(nb)), float(v[: n // 2].mean())


def ground_truth(n=200000, seeds=(990100, 990101, 990102)):
    """Phase 8 rule (ground_truth.py): combined lam over 3 ICs, SE as Phase 8; NON-CHAOTIC iff
    lam - 2.576 SE <= 0 (flow tolerance: quasi-periodic exponent is 0) -- no CHAOTIC regime exists
    for K < 1."""
    out = []
    for r in REGIMES:
        ls = [lyapunov(r["K"], r["Om"], n, s) for s in seeds]
        lam = float(np.mean([v[0] for v in ls]))
        se = float(np.sqrt(np.mean([v[1] ** 2 for v in ls]) / 3 + np.var([v[0] for v in ls]) / 3))
        lab = "NON-CHAOTIC" if lam - 2.576 * se <= 0 else "AMBIGUOUS"
        out.append({**r, "lam": lam, "se": se, "label": lab})
    return out


def window_types(n):
    return np.array(["N"] * (n + 1))
