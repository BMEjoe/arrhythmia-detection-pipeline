"""
Phase 4 systems: the Phase 2E generators plus the quasi-periodic `two_tone`
control, all under the Phase 2E seed scheme.
"""
from __future__ import annotations

import numpy as np

from experiments.phase2e import config as C2
from experiments.phase2e import systems as S2
from experiments.phase4_upo import config as C


def _seed_sequence(system, n, seed):
    code = C.TWO_TONE_CODE if system == "two_tone" else C2.SYSTEM_CODE[system]
    return np.random.SeedSequence([C2.PHASE2E_ENTROPY, code, int(n), int(seed)])


def data_seed(system, n, seed):
    return int(_seed_sequence(system, n, seed).generate_state(1)[0])


def two_tone(n, rng):
    """x_k = sin(2 pi f1 k + p1) + 0.6 sin(2 pi f2 k + p2), f2/f1 = golden ratio."""
    p1, p2 = rng.uniform(0.0, 2.0 * np.pi, 2)
    k = np.arange(int(n))
    return np.sin(2 * np.pi * C.TWO_TONE_F1 * k + p1) + C.TWO_TONE_A2 * np.sin(2 * np.pi * C.TWO_TONE_F2 * k + p2)


def generate(system, n, seed, snr_db=None):
    """Phase 2E generate() for existing systems; same scheme for two_tone."""
    if system != "two_tone":
        return S2.generate(system, n, seed, snr_db)
    data_ss, noise_ss = _seed_sequence(system, n, seed).spawn(2)
    x = two_tone(n, np.random.default_rng(data_ss))
    if snr_db is None:
        return x
    sigma = np.sqrt(np.var(x) / 10.0 ** (float(snr_db) / 10.0))
    noise_rng = np.random.default_rng(
        np.random.SeedSequence([int(noise_ss.generate_state(1)[0]), int(round(snr_db * 100)) + 100000]))
    return x + sigma * noise_rng.standard_normal(len(x))


def fixed_point_reference(system):
    """On-attractor analytic period-1 fixed point(s) used for localization error."""
    if system == "henon":
        return S2.henon_fixed_points()[0]
    if system == "skewed_henon":
        return S2.skewed_henon_fixed_points()[0]
    if system == "logistic":
        return 0.75
    return None
