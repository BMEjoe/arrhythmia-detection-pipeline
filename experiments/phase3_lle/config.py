"""
Phase 3 -- LLE as a statistically meaningful, noise-robust chaos detector.

Fixed experiment configuration.  Systems, generators and the seed scheme are
the Phase 2E ones (experiments/phase2e/systems.py, SeedSequence([20260926,
system_code, window_length, seed]); noise from a child stream per (seed, SNR)).

SEED SPLIT (essential)
  DEVELOPMENT seeds = the Phase 2E core seeds (white_noise / ar1: 0-99 at 256,
      0-49 at 512; deterministic systems: 0-29; noisy systems: 0-29, the Phase 2E
      noise-arm seeds).  Used freely for baseline, diagnosis and tuning.
  TEST seeds = 1000-1099 for every condition, generated with the same scheme.
      They are NOT run before PREREGISTRATION.md is committed and pushed
      (run_phase3.py enforces this) and are the only seeds reported as final
      performance.
"""
from __future__ import annotations

import math

SYSTEMS = ("white_noise", "ar1", "sinusoid", "logistic_p4", "logistic", "henon")
NOISY_SYSTEMS = ("logistic", "henon")
SNR_LEVELS_DB = (30.0, 20.0, 10.0)
WINDOW_LENGTHS = (256, 512)
PRIMARY_WINDOW = 256          # matches CLASSIFIER_WINDOW_RR in MIT-BIH

NULL_SYSTEMS = ("white_noise", "ar1")
PERIODIC_SYSTEMS = ("sinusoid", "logistic_p4")
CHAOTIC_SYSTEMS = ("logistic", "henon")

# Largest Lyapunov exponent references, natural log per map iteration (= per
# sample).  Used ONLY to measure bias, never inside an estimator.
LLE_REFERENCE = {"logistic": math.log(2.0), "henon": 0.4192}
# Known dimension of the generating map (diagnosis only, B2).
TRUE_DIMENSION = {"logistic": 1, "henon": 2}

# ---------------------------------------------------------------------------
# Seeds
# ---------------------------------------------------------------------------
DEV_SEEDS = {
    ("white_noise", 256): range(100), ("white_noise", 512): range(50),
    ("ar1", 256): range(100), ("ar1", 512): range(50),
}
DEV_SEEDS_DEFAULT = range(30)          # deterministic and noisy conditions
TEST_SEEDS = range(1000, 1100)          # every condition; see module docstring


def conditions():
    """(system, snr_db) pairs; snr_db None = clean."""
    out = [(s, None) for s in SYSTEMS]
    out += [(s, snr) for s in NOISY_SYSTEMS for snr in SNR_LEVELS_DB]
    return out


def seeds(phase, system, n, snr_db=None):
    if phase == "dev":
        if snr_db is None:
            return list(DEV_SEEDS.get((system, n), DEV_SEEDS_DEFAULT))
        return list(DEV_SEEDS_DEFAULT)
    if phase == "test":
        return list(TEST_SEEDS)
    raise ValueError(phase)


def condition_label(system, snr_db):
    return system if snr_db is None else f"{system}@{int(snr_db)}dB"
