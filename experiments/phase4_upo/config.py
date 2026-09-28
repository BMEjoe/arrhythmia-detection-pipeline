"""
Phase 4 -- UPO detector: specificity (periodic / quasi-periodic rejection) and
noise robustness.  Fixed experiment configuration.

Systems and seed scheme are the Phase 2E ones (experiments/phase2e/systems.py:
SeedSequence([20260926, system_code, window_length, seed]); observation noise
from a child stream per (seed, SNR)).  New system: `two_tone`, a quasi-periodic
control with system code 9 under the same scheme (systems.py here).

SEED SPLIT
  DEVELOPMENT = Phase 2E seeds: white_noise / ar1 0-99 at 256 and 0-49 at 512;
      every other condition 0-29.
  TEST = seeds starting at 2000 (TEST_SEED_START), same SeedSequence scheme.
      2000+ rather than 1000-1099 so that the windows are untouched by every
      earlier phase (Phase 3 evaluated lle_chaos_test on 1000-1099).  The number
      of test seeds per condition is fixed in PREREGISTRATION.md from the
      runtime budget.  run_phase4.py refuses --phase test until
      PREREGISTRATION.md is committed, unmodified and pushed.
"""
from __future__ import annotations

SYSTEMS = ("white_noise", "ar1", "sinusoid", "two_tone", "logistic_p4",
           "logistic", "henon", "skewed_henon")
NOISY_SYSTEMS = ("logistic", "henon")
SNR_LEVELS_DB = (30.0, 20.0, 10.0)
WINDOW_LENGTHS = (256, 512)
PRIMARY_WINDOW = 256

NULL_SYSTEMS = ("white_noise", "ar1")
PERIODIC_SYSTEMS = ("sinusoid", "two_tone", "logistic_p4")
CHAOTIC_SYSTEMS = ("logistic", "henon", "skewed_henon")
# PASS set of the primary rule (user specification)
SPECIFICITY_SYSTEMS = ("white_noise", "ar1", "sinusoid", "two_tone")

TWO_TONE_CODE = 9
TWO_TONE_F1 = 1.0 / 7.3                        # same period as the Phase 2E sinusoid
TWO_TONE_F2 = TWO_TONE_F1 * (1.0 + 5 ** 0.5) / 2.0   # golden-ratio multiple: incommensurate
TWO_TONE_A2 = 0.6

DEV_SEEDS = {("white_noise", 256): range(100), ("white_noise", 512): range(50),
             ("ar1", 256): range(100), ("ar1", 512): range(50)}
DEV_SEEDS_DEFAULT = range(30)
TEST_SEED_START = 2000
N_TEST_SEEDS = None          # set in PREREGISTRATION.md (budget), copied here before the test run


def conditions():
    out = [(s, None) for s in SYSTEMS]
    out += [(s, snr) for s in NOISY_SYSTEMS for snr in SNR_LEVELS_DB]
    return out


def seeds(phase, system, n, snr_db=None):
    if phase == "dev":
        if snr_db is None:
            return list(DEV_SEEDS.get((system, n), DEV_SEEDS_DEFAULT))
        return list(DEV_SEEDS_DEFAULT)
    if phase == "test":
        if N_TEST_SEEDS is None:
            raise SystemExit("N_TEST_SEEDS is not set (fixed by PREREGISTRATION.md)")
        return list(range(TEST_SEED_START, TEST_SEED_START + int(N_TEST_SEEDS)))
    raise ValueError(phase)


def condition_label(system, snr_db):
    return system if snr_db is None else f"{system}@{int(snr_db)}dB"
