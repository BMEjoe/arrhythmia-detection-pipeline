"""
Phase 4 methods: fixed decision rules applied to saved detector runs.

Each method names the RUN configuration it needs (detector.run arguments
(m_rule, M, n_surrogates)) and a `detect(run_out) -> list of detected peaks`
rule.  A window is a detection ("verified detection" for candidates) when that
list is non-empty.  Candidates are frozen in PREREGISTRATION.md; nothing here
may change after it is pushed.

All three candidates include the P1-a instability gate (leading multiplier
modulus >= 1 + DELTA) on Level-B period-1 peaks.  C2 and C3 use the hybrid
PRL/PRE source-stability estimate with a ROBUST aggregate of the member
Jacobians (element-wise median, 10 % trimmed mean), never the arithmetic mean.
"""
from __future__ import annotations

from experiments.phase4_upo import detector as D

N_SUR = 50                 # surrogate count, fixed in PREREGISTRATION.md (budget)
DELTA = 0.2                # P1-a threshold: leading modulus >= 1.2

PROD_RUN = ("cao", 7, N_SUR)          # production: Cao at lag 1, M = 7
M2_M15_RUN = (2, 15, N_SUR)           # fixed m = 2, M = 15 Jacobian neighbours
M2_M7_RUN = (2, 7, N_SUR)             # fixed m = 2, production M = 7


def _baseline(run_out):
    """Unchanged production detector: >= 1 Level-B (surrogate-significant) peak."""
    return D.peaks_level_b(run_out)


def _c1(run_out):
    """C1: production run + P1-a gate on the candidate-centred extension monodromy."""
    return D.gate(D.peaks_level_b(run_out), "ext_lead", DELTA)


def _c2(run_out):
    """C2: fixed m = 2, M = 15 + P1-a gate on the hybrid estimate, element-wise median."""
    return D.gate(D.peaks_level_b(run_out), "src_median", DELTA)


def _c3(run_out):
    """C3: fixed m = 2, M = 7 + P1-a gate on the hybrid estimate, 10 % trimmed mean."""
    return D.gate(D.peaks_level_b(run_out), "src_trim10", DELTA)


METHODS = {
    "baseline": {"run": PROD_RUN, "detect": _baseline, "eligible": False, "index": 0},
    "c1_prod_extgate": {"run": PROD_RUN, "detect": _c1, "eligible": True, "index": 1},
    "c2_m2M15_mediangate": {"run": M2_M15_RUN, "detect": _c2, "eligible": True, "index": 2},
    "c3_m2M7_trimgate": {"run": M2_M7_RUN, "detect": _c3, "eligible": True, "index": 3},
}
