"""
Phase 6 detector configurations.  Every method is BASELINE-K
(combined_chaos_config() with keep_upo_on_short_lle_embedding = True) with at
most the rr_detrend option changed.  PREREGISTERED holds the baseline and the
(at most 3) Part C candidates, frozen by PREREGISTRATION.md; TUNING holds the
development-only variants used to choose the candidates' parameters.
"""
from __future__ import annotations

from dataclasses import replace

import final_pipeline as fp

BASELINE_K = replace(fp.combined_chaos_config(), keep_upo_on_short_lle_embedding=True)


def detrended(method, **kw):
    return replace(BASELINE_K, rr_detrend=method, **kw)


TUNING = {
    "baseline_k": BASELINE_K,
    "linear": detrended("linear"),
    "movmed_31": detrended("moving_median", rr_detrend_window=31),
    "movmed_61": detrended("moving_median", rr_detrend_window=61),
    "movmed_101": detrended("moving_median", rr_detrend_window=101),
    "sp_10": detrended("smoothness_priors", rr_detrend_lambda=10.0),
    "sp_50": detrended("smoothness_priors", rr_detrend_lambda=50.0),
    "sp_300": detrended("smoothness_priors", rr_detrend_lambda=300.0),
}

# Filled from the development tuning (results/dev/tables/tuning.md) BEFORE the
# preregistration; frozen afterwards.
PREREGISTERED = {
    "baseline_k": {"config": BASELINE_K, "index": 0, "eligible": False},
}
