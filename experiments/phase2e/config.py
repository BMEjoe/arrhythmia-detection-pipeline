"""
Phase 2E -- PREDECLARED experiment configuration.

Every system parameter, seed rule, window length, SNR level, surrogate count
and numerical-sensitivity grid used by Phase 2E is fixed in this file BEFORE
any experiment is run.  run_phase2e.py records the SHA-256 of this file in
results/manifest.json so that the configuration that produced a result set
is auditable.  Nothing here was chosen after inspecting Phase 2E output.

The frozen pipeline configuration is final_pipeline.CFG.  Phase 2E only
changes the fields that define an experimental arm (map mode, significance
request, surrogate count); every change is listed explicitly below.
"""
from __future__ import annotations

# ---------------------------------------------------------------------------
# Global seed rule
# ---------------------------------------------------------------------------
# Data for task (system, window_length, seed) is generated from
#   np.random.SeedSequence([PHASE2E_ENTROPY, SYSTEM_CODE[system], window_length, seed])
# so every realization is independent of every other realization and of
# worker scheduling.  Additive observation noise uses a separate child
# stream (see systems.py).
PHASE2E_ENTROPY = 20260926

SYSTEM_CODE = {
    "constant": 1,
    "white_noise": 2,
    "ar1": 3,
    "sinusoid": 4,
    "logistic_p4": 5,      # optional second periodic control (stable 4-cycle)
    "logistic": 6,
    "henon": 7,
    "skewed_henon": 8,     # only in the localization study (PRE Eq. 31 benchmark)
}

# ---------------------------------------------------------------------------
# System definitions (documented here, implemented in systems.py)
# ---------------------------------------------------------------------------
TRANSIENT = 1000               # iterations discarded for every map / AR process
CONSTANT_RANGE = (0.5, 1.5)    # constant value c ~ U(0.5, 1.5)
AR1_PHI = 0.8                  # x_n = phi x_{n-1} + e_n, e_n ~ N(0, 1)
SINUSOID_PERIOD = 7.3          # samples (non-integer; same period as tests/upo_systems.periodic_rr)
LOGISTIC_R = 4.0               # chaotic logistic
LOGISTIC_X0_RANGE = (0.05, 0.95)
LOGISTIC_P4_R = 3.5            # stable period-4 cycle
HENON_A, HENON_B = 1.4, 0.3    # x' = 1 - a x^2 + y, y' = b x; observable passed on: x
HENON_IC_RANGE = (-0.1, 0.1)   # (x0, y0) ~ U(-0.1, 0.1)^2

# ---------------------------------------------------------------------------
# Window lengths (Step 4) -- not to be changed
# ---------------------------------------------------------------------------
WINDOW_LENGTHS = (128, 256, 512)
PRODUCTION_WINDOW = 256

# ---------------------------------------------------------------------------
# Surrogates (Step 8)
# ---------------------------------------------------------------------------
# final_pipeline.CFG.so_surrogate_count = 50 (PRE Sec. IV).  The main Phase
# 2E validation uses this production value unchanged.
MAIN_SURROGATES = 50
SURROGATE_SENSITIVITY = (20, 100)   # plus MAIN_SURROGATES from the core arm

# ---------------------------------------------------------------------------
# Seeds per arm
# ---------------------------------------------------------------------------
NULL_SYSTEMS = ("white_noise", "ar1")
DETERMINISTIC_SYSTEMS = ("sinusoid", "logistic_p4", "logistic", "henon")
CORE_SYSTEMS = ("constant",) + NULL_SYSTEMS + DETERMINISTIC_SYSTEMS

N_SEEDS_CORE = {
    "constant": 3,
    # Step 7: >= 50 null windows mandatory, 100 preferred at the production window.
    "white_noise": {128: 50, 256: 100, 512: 50},
    "ar1": {128: 50, 256: 100, 512: 50},
    "sinusoid": 30, "logistic_p4": 30, "logistic": 30, "henon": 30,
}

# Step 11 -- noise robustness.  PREDECLARED SNR levels (dB).  "clean" is the
# core arm (paired: identical clean windows, same seeds).
#   SNR_dB = 10 log10( var(clean window) / sigma^2 ),  additive white Gaussian
#   OBSERVATIONAL noise (the map itself is iterated noise-free).
NOISE_SYSTEMS = ("logistic", "henon")
SNR_LEVELS_DB = (30.0, 20.0, 10.0, 5.0, 0.0)
N_SEEDS_NOISE = 30

# Step 13 -- one_sample vs tau_step, paired on the core seeds at 256.
TAU_STEP_SYSTEMS = CORE_SYSTEMS
N_SEEDS_TAU_STEP = 30               # first 30 core seeds of each system (constant: 3)

# Step 8 -- surrogate-count sensitivity, paired on the core seeds at 256.
SURROGATE_SENS_SYSTEMS = ("white_noise", "ar1", "logistic", "henon")
N_SEEDS_SURROGATE_SENS = {"white_noise": 50, "ar1": 50, "logistic": 30, "henon": 30}

# Algorithm-seed sensitivity of the null (random R draws and surrogate phases
# use config.random_seed, which is the SAME for every production window).
# Paired on the first 50 core null seeds at 256, with
#   random_seed = CFG.random_seed + 1000 + seed.
ALGO_SEED_SYSTEMS = NULL_SYSTEMS
N_SEEDS_ALGO = 50
ALGO_SEED_OFFSET = 1000

# Step 14 -- DIAGNOSTIC ONLY: the detector run at the known dimension of the
# generating map (one_sample, tau = 1).  Never substituted for production
# results; used to separate embedding-selection effects from detection.
ORACLE_DIMENSION = {"logistic": 1, "henon": 2}
N_SEEDS_ORACLE = 30

# Step 15 -- floating-point / RR-style sensitivity (detection without
# significance; Level A, Cao, LLE).  RR-like scaling rr = 0.8 + 0.05 z,
# z = standardized window, then two quantization routes at 1/360 s.
FP_SYSTEMS = ("sinusoid", "logistic", "henon")
N_SEEDS_FP = 30
RR_FS = 360.0

# Step 25 -- replicability: rerun these task ids (first k per experiment).
REPLICATE_PER_EXPERIMENT = 4

# ---------------------------------------------------------------------------
# Step 16 -- Henon localization grid (PREDECLARED)
# ---------------------------------------------------------------------------
# Source-style operating point, as in tests/upo_systems.henon_config and
# PRE Sec. III C: d = 2, one_sample map (tau = 1), M = 2, K = 1, R = 100,
# kappa = 3, N = 1024, prl_norm randomization, so_max_backbones = None.
LOC_N = 1024
LOC_SYSTEMS = ("henon", "skewed_henon")
LOC_REALIZATIONS = (0, 1, 2)
LOC_BINS = (15, 20, 30, 45, 60, 90)                 # production: 30
LOC_TUBE_PERCENTILES = (2.5, 5.0, 10.0, 20.0, 40.0) # production: 10 (tube and slab)
# Range policies.  "padXX": attractor range padded by XX% of the span on each
# side (production = pad10).  "pad10_shift": production range shifted by half
# a bin width (tests bin-phase / quantization sensitivity).
LOC_RANGE_POLICIES = ("pad0", "pad10", "pad25", "pad10_shift")
# Surrogate significance (period 1 only, MAIN_SURROGATES) on the three
# one-factor-at-a-time sweeps through the production point.
LOC_SIG_REALIZATIONS = (0, 1, 2)
# Production-Jacobian arm: the same grid point but M = 7 (CFG default).
LOC_PRODUCTION_M = 7
