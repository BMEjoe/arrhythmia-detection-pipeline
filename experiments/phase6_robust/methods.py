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
    "linear_g0.5": detrended("linear", rr_detrend_min_trend_sd=0.5),
    "linear_g0.7": detrended("linear", rr_detrend_min_trend_sd=0.7),
    "sp_300_g0.7": detrended("smoothness_priors", rr_detrend_lambda=300.0, rr_detrend_min_trend_sd=0.7),
}

# Filled from the development tuning (results/dev/tables/tuning.md) BEFORE the
# preregistration; frozen afterwards.
PREREGISTERED = {
    "baseline_k": {"config": BASELINE_K, "index": 0, "eligible": False},
    # D1: linear least-squares detrend, applied when the fitted trend change >= 0.5 residual SD
    "d1_linear_g05": {"config": detrended("linear", rr_detrend_min_trend_sd=0.5), "index": 1, "eligible": True},
    # D2: linear least-squares detrend, applied when the fitted trend change >= 0.7 residual SD
    "d2_linear_g07": {"config": detrended("linear", rr_detrend_min_trend_sd=0.7), "index": 2, "eligible": True},
    # D3: smoothness-priors detrend (lambda 300, beat index), applied when the change >= 0.7 residual SD
    "d3_smoothprior300_g07": {"config": detrended("smoothness_priors", rr_detrend_lambda=300.0,
                                                  rr_detrend_min_trend_sd=0.7), "index": 3, "eligible": True},
}


def effective_key(cfg, trend_ratio):
    """Configurations with the same effective key give identical analyze_segment decisions on a
    window: below its gate a detrend method returns the input unchanged (fp.detrend_rr), so
    the window is analysed exactly as by BASELINE-K."""
    if cfg.rr_detrend == "none" or (cfg.rr_detrend_min_trend_sd > 0 and trend_ratio < cfg.rr_detrend_min_trend_sd):
        return ("none",)
    key = (cfg.rr_detrend,)
    if cfg.rr_detrend == "moving_median":
        key += (int(cfg.rr_detrend_window),)
    if cfg.rr_detrend == "smoothness_priors":
        key += (float(cfg.rr_detrend_lambda),)
    return key
