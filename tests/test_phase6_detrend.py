"""Phase 6 opt-in rr_detrend option (default "none" = previous behaviour)."""
import pathlib
import sys
from dataclasses import replace

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import final_pipeline as fp  # noqa: E402

RNG = np.random.default_rng(0)
X = 0.8 + 0.04 * RNG.standard_normal(256)
TREND = np.linspace(-0.04, 0.04, 256)


def test_default_is_none_and_output_unchanged():
    assert fp.CFG.rr_detrend == "none"
    assert np.array_equal(fp.detrend_rr(X, fp.CFG), X)
    out = fp.analyze_segment(X, fp.CFG)
    assert "rr_detrend" not in out and np.array_equal(out["rr_dynamics"], X)


@pytest.mark.parametrize("method,kw", [("linear", {}), ("smoothness_priors", {"rr_detrend_lambda": 300.0}),
                                       ("moving_median", {"rr_detrend_window": 51})])
def test_detrend_removes_linear_trend_and_keeps_mean(method, kw):
    cfg = replace(fp.CFG, rr_detrend=method, **kw)
    y = fp.detrend_rr(X + TREND, cfg)
    assert abs(y.mean() - (X + TREND).mean()) < 1e-12
    slope = np.polyfit(np.arange(256), y, 1)[0] * 255
    assert abs(slope) < (1e-9 if method == "linear" else 0.03)


def test_linear_is_exact_on_a_line_plus_signal():
    y = fp.detrend_rr(X + TREND, replace(fp.CFG, rr_detrend="linear"))
    ref = X - np.polyval(np.polyfit(np.arange(256), X, 1), np.arange(256)) + (X + TREND).mean()
    assert np.allclose(y, ref)


def test_invalid_options_raise():
    with pytest.raises(ValueError):
        fp.detrend_rr(X, replace(fp.CFG, rr_detrend="bogus"))
    with pytest.raises(ValueError):
        fp.detrend_rr(X, replace(fp.CFG, rr_detrend="moving_median", rr_detrend_window=50))


def test_analyze_segment_uses_detrended_series():
    cfg = replace(fp.CFG, rr_detrend="linear", lle_chaos_test=True)
    out = fp.analyze_segment(X + TREND, cfg)
    assert out["rr_detrend"]["method"] == "linear"
    assert np.allclose(out["rr_dynamics"], fp.detrend_rr(X + TREND, cfg))
    assert np.array_equal(out["rr_raw"], X + TREND)


def test_min_trend_gate():
    cfg = replace(fp.CFG, rr_detrend="linear", rr_detrend_min_trend_sd=0.5)
    assert np.array_equal(fp.detrend_rr(X, cfg), X)                    # stationary: ratio << 0.5
    assert fp.linear_trend_ratio(X + 3 * TREND) >= 0.5
    y = fp.detrend_rr(X + 3 * TREND, cfg)
    assert np.allclose(y, fp.detrend_rr(X + 3 * TREND, replace(cfg, rr_detrend_min_trend_sd=0.0)))
    out = fp.analyze_segment(X, replace(cfg, lle_chaos_test=False))
    assert out["rr_detrend"]["applied"] is False and np.array_equal(out["rr_dynamics"], X)
