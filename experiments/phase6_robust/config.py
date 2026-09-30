"""
Phase 6 -- robustness of the combined chaos detector: crash fix (Part A),
heavier ectopy (Part B), detrending candidates (Part C).

Signals: the Phase 5 generators (experiments/phase5_rr: 256 RR intervals in
seconds, 1/360 s quantization, Phase 2E seed scheme, paired base / modifier
streams) plus the Part B ectopy conditions below, all built on the N1
linear_rr base (the same base realization and RR scale as N1 for a seed).

SEED SPLIT
  DEVELOPMENT = as Phase 5: 0-99 for null-type conditions (N*, S*, E*),
      0-29 for positive and generalization conditions (P*, G*).
  TEST = seeds starting at 4000 (TEST_SEED_START); counts fixed in
      PREREGISTRATION.md; run_phase6.py refuses --phase test until it is
      committed, unmodified and pushed.
"""
from __future__ import annotations

from experiments.phase5_rr import config as C5

ENTROPY = C5.ENTROPY
N = C5.N

# ---- Part B ectopy parameters (sources: README.md) ----
PREMATURE_COUPLING_RANGE = (0.6, 0.8)   # premature interval = c * underlying sinus RR (as Phase 5 S2)
RUN_CYCLE_RANGE_S = (0.40, 0.55)        # ectopic cycle length inside a run: 109-150 bpm (> 100 bpm, VT definition)
RUN_COUNT_RANGE = (1, 3)                # runs per window (inclusive)
RUN_LENGTH_RANGE = (3, 6)               # premature beats per run (inclusive)
COUPLET_RATES = (0.05, 0.10)            # fraction of beats premature; couplets = round(rate * N / 2)
ATRIAL_RATES = (0.05, 0.10)             # isolated atrial-type premature beats, round(rate * N)
ATRIAL_POST_FACTOR_RANGE = (1.0, 1.1)   # post-ectopic interval = d * sinus RR (sinus node reset: pause < compensatory)
MIN_NORMAL_GAP = 3                      # normal intervals between ectopic groups (couplets, runs, atrial)

# condition -> (group, kind, params, quantized-raw twin for edited variants)
PARTB = {
    "E1_bigeminy":            ("ectopy", "bigeminy", {}),
    "E2_trigeminy":           ("ectopy", "trigeminy", {}),
    "E3_couplets_5pct":       ("ectopy", "couplets", {"rate": 0.05}),
    "E3_couplets_10pct":      ("ectopy", "couplets", {"rate": 0.10}),
    "E4_runs":                ("ectopy", "runs", {}),
    "E5_atrial_5pct":         ("ectopy", "atrial", {"rate": 0.05}),
    "E5_atrial_10pct":        ("ectopy", "atrial", {"rate": 0.10}),
    "E6_ectopic10_trend":     ("ectopy", "s2_trend", {}),
    "S2_ectopic_5pct_edited":  ("edited", "edited", {"raw": "S2_ectopic_5pct"}),
    "S2_ectopic_10pct_edited": ("edited", "edited", {"raw": "S2_ectopic_10pct"}),
    "E1_bigeminy_edited":      ("edited", "edited", {"raw": "E1_bigeminy"}),
    "E2_trigeminy_edited":     ("edited", "edited", {"raw": "E2_trigeminy"}),
    "E3_couplets_10pct_edited": ("edited", "edited", {"raw": "E3_couplets_10pct"}),
}
PARTB_CODE = {c: 200 + i for i, c in enumerate(PARTB)}

# every condition used in Phase 6: Phase 5 conditions + Part B
GROUP = {c: v[0] for c, v in C5.CONDITIONS.items()}
GROUP.update({c: v[0] for c, v in PARTB.items()})
CONDITIONS = list(C5.CONDITIONS) + list(PARTB)
PASS_NULLS = C5.PASS_CONDITIONS                                   # N1-N6
S2_RATES = ("S2_ectopic_2pct", "S2_ectopic_5pct", "S2_ectopic_10pct")
POWER_CHECK = ("P1_henon_rr", "P2_logistic_rr", "P4_henon_rr_30dB", "P4_henon_rr_20dB",
               "P4_logistic_rr_30dB", "P4_logistic_rr_20dB", "G1_lorenz_maxima")
TREND_POSITIVES = ("P3_henon_rr_trend", "P3_logistic_rr_trend")
PARTB_CONDITIONS = tuple(PARTB)

TEST_SEED_START = 4000
N_TEST_SEEDS = None                # fixed in PREREGISTRATION.md
N_TEST_SEEDS_BY_GROUP = None       # fixed in PREREGISTRATION.md
DEV_SEEDS_NULL = range(100)
DEV_SEEDS_OTHER = range(30)


def group(condition):
    return GROUP[condition]


def seeds(phase, condition):
    if phase == "dev":
        return list(DEV_SEEDS_OTHER if group(condition) in ("positive", "generalization") else DEV_SEEDS_NULL)
    if phase == "test":
        n = (N_TEST_SEEDS_BY_GROUP or {}).get(group(condition), N_TEST_SEEDS)
        if n is None:
            raise SystemExit("N_TEST_SEEDS is not set (fixed by PREREGISTRATION.md)")
        return list(range(TEST_SEED_START, TEST_SEED_START + int(n)))
    raise ValueError(phase)
