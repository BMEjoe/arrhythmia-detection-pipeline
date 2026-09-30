"""
Phase 5 -- stress test of the FROZEN combined chaos detector on RR-interval-like
signals.  Fixed experiment configuration (generator parameters are frozen at
the Step 1 commit, before PREREGISTRATION.md).

Every window is N = 256 RR intervals in seconds.  Unless a condition is marked
unquantized, the series is rounded to multiples of 1/360 s (MIT-BIH sampling
rate) AFTER all other transformations.

SEED SCHEME (Phase 2E scheme, experiments/phase2e/systems.py)
  base stream  SeedSequence([20260926, BASE_CODE[base], N, seed])
               -> RR scale (mean, SD) then the base process
  modifier     SeedSequence([20260926, COND_CODE[condition], N, seed])
               (S3 conditions use their quantized twin's code, MODIFIER_TWIN)
               -> trend / step / ectopic / observation-noise draws
  Conditions built on the same base with the same seed share the base
  realization and RR scale (paired design: e.g. N1, N3, N4, S2 and S3 differ
  only in the modification or in quantization).
SEED SPLIT
  DEVELOPMENT = Phase 2E seed numbers: 0-99 for null conditions (N*, S*),
      0-29 for positive and generalization conditions (P*, G*).
  TEST = seeds starting at 3000 (TEST_SEED_START).  The number per condition is
      fixed in PREREGISTRATION.md; run_phase5.py refuses --phase test until it
      is committed, unmodified and pushed.
"""
from __future__ import annotations

ENTROPY = 20260926                  # Phase 2E entropy
N = 256                             # RR intervals per window
QUANTUM_S = 1.0 / 360.0             # MIT-BIH sampling interval

# ---- RR scale drawn per window (shared by every condition on the same base) ----
MEAN_RR_RANGE_S = (0.65, 0.95)      # target 0.6-1.0 s; margin for trend / step / ectopics
SD_RR_RANGE_S = (0.035, 0.060)      # short-term SD from Task Force LF+HF norms (README)

# ---- N1 linear_rr: McSharry et al. (2003) bimodal spectrum, random phases ----
LF_CENTER_HZ, HF_CENTER_HZ = 0.10, 0.25
LF_WIDTH_HZ, HF_WIDTH_HZ = 0.01, 0.01          # ECGSYN defaults (c1, c2)
LFHF_RANGE = (1.5, 2.0)                         # Task Force short-term LF/HF (ECGSYN default 0.5 not used)
LINEAR_SYNTH_LENGTH = 1024                      # synthesize 1024 beats, keep a random 256 segment

# ---- N2 power_law: Gaussian 1/f^beta ----
BETA_RANGE = (0.9, 1.1)
POWER_LAW_SYNTH_LENGTH = 4096

# ---- N3 / P3 trend, N4 step (fractions of the window's mean RR) ----
TREND_FRACTION_RANGE = (0.05, 0.10)   # total linear change across the window, random sign
STEP_FRACTION_RANGE = (0.05, 0.10)    # abrupt step, random sign
STEP_POSITION_RANGE = (0.3, 0.7)      # step index uniform in [0.3 N, 0.7 N]

# ---- N5 warp: z -> exp(WARP_A z), then rescaled to (mean, SD) ----
WARP_A = 0.5                          # lognormal-type right skew (skewness ~1.75 before scaling)

# ---- N6 noisy_rsa: respiratory sinusoid + independent white Gaussian noise ----
RESP_FREQ_RANGE_HZ = (0.20, 0.30)
RSA_VARIANCE_SHARE = 0.5              # sinusoid variance / total variance

# ---- S1 SETAR(2;1,1): stable fixed point at 0, noise driven ----
SETAR_PHI_LOW, SETAR_PHI_HIGH = 0.7, -0.5   # x_t = phi x_{t-1} + e_t, regime by sign of x_{t-1}
TRANSIENT = 1000

# ---- S2 ectopic beats ----
ECTOPIC_RATES = (0.02, 0.05, 0.10)    # fraction of beats; counts round(rate * N) = 5, 13, 26
ECTOPIC_PREMATURITY_RANGE = (0.6, 0.8)  # coupling interval = c * RR_i
ECTOPIC_MIN_SPACING = 3               # beats between ectopic indices

# ---- P4 observation noise ----
SNR_LEVELS_DB = (30.0, 20.0)

# ---- G1-G3 flows (RK4) ----
LORENZ = {"sigma": 10.0, "rho": 28.0, "beta": 8.0 / 3.0, "dt": 0.005, "burn_in": 100.0}
ROSSLER = {"a": 0.2, "b": 0.2, "c": 5.7, "dt": 0.01, "burn_in": 500.0, "sample_step": 1.0}
MACKEY_GLASS = {"beta": 0.2, "gamma": 0.1, "n": 10, "tau": 17.0, "dt": 0.1, "burn_in": 1000.0,
                "sample_step": 6.0}

# ---- codes ----
BASE_CODE = {"linear_rr": 51, "power_law": 52, "noisy_rsa": 53, "setar": 54, "henon": 55,
             "logistic": 56, "lorenz_maxima": 57, "rossler_flow": 58, "mackey_glass": 59}

# condition -> (group, base, modifiers, quantized)
#   modifiers: ("warp",), ("trend",), ("step",), ("ectopic", rate), ("noise", snr_db)
CONDITIONS = {
    "N1_linear_rr":        ("null", "linear_rr", (), True),
    "N2_power_law":        ("null", "power_law", (), True),
    "N3_linear_rr_trend":  ("null", "linear_rr", (("trend",),), True),
    "N4_linear_rr_step":   ("null", "linear_rr", (("step",),), True),
    "N5_linear_rr_warped": ("null", "linear_rr", (("warp",),), True),
    "N6_noisy_rsa":        ("null", "noisy_rsa", (), True),
    "S1_setar":            ("secondary", "setar", (), True),
    "S2_ectopic_2pct":     ("secondary", "linear_rr", (("ectopic", 0.02),), True),
    "S2_ectopic_5pct":     ("secondary", "linear_rr", (("ectopic", 0.05),), True),
    "S2_ectopic_10pct":    ("secondary", "linear_rr", (("ectopic", 0.10),), True),
    "S3_linear_rr_unq":       ("secondary", "linear_rr", (), False),
    "S3_linear_rr_trend_unq": ("secondary", "linear_rr", (("trend",),), False),
    "S3_linear_rr_step_unq":  ("secondary", "linear_rr", (("step",),), False),
    "P1_henon_rr":         ("positive", "henon", (), True),
    "P2_logistic_rr":      ("positive", "logistic", (), True),
    "P3_henon_rr_trend":   ("positive", "henon", (("trend",),), True),
    "P3_logistic_rr_trend": ("positive", "logistic", (("trend",),), True),
    "P4_henon_rr_30dB":    ("positive", "henon", (("noise", 30.0),), True),
    "P4_henon_rr_20dB":    ("positive", "henon", (("noise", 20.0),), True),
    "P4_logistic_rr_30dB": ("positive", "logistic", (("noise", 30.0),), True),
    "P4_logistic_rr_20dB": ("positive", "logistic", (("noise", 20.0),), True),
    "G1_lorenz_maxima":    ("generalization", "lorenz_maxima", (), True),
    "G2_rossler_flow":     ("generalization", "rossler_flow", (), True),
    "G3_mackey_glass":     ("generalization", "mackey_glass", (), True),
}
COND_CODE = {c: 100 + i for i, c in enumerate(CONDITIONS)}
# S3 conditions are their quantized twins without quantization: they use the
# twin's modifier stream, so the two differ ONLY in quantization.
MODIFIER_TWIN = {"S3_linear_rr_unq": "N1_linear_rr", "S3_linear_rr_trend_unq": "N3_linear_rr_trend",
                 "S3_linear_rr_step_unq": "N4_linear_rr_step"}
PASS_CONDITIONS = tuple(c for c, v in CONDITIONS.items() if v[0] == "null")   # N1-N6

DEV_SEEDS_NULL = range(100)
DEV_SEEDS_OTHER = range(30)
TEST_SEED_START = 3000
N_TEST_SEEDS = 200           # fixed in PREREGISTRATION.md Section 5 (runtime budget)
N_TEST_SEEDS_BY_GROUP = {"null": 300}   # N1-N6: 300 test seeds (3000-3299); others 200 (3000-3199)


def group(condition):
    return CONDITIONS[condition][0]


def seeds(phase, condition):
    if phase == "dev":
        return list(DEV_SEEDS_NULL if group(condition) in ("null", "secondary") else DEV_SEEDS_OTHER)
    if phase == "test":
        n = (N_TEST_SEEDS_BY_GROUP or {}).get(group(condition), N_TEST_SEEDS)
        if n is None:
            raise SystemExit("N_TEST_SEEDS is not set (fixed by PREREGISTRATION.md)")
        return list(range(TEST_SEED_START, TEST_SEED_START + int(n)))
    raise ValueError(phase)
