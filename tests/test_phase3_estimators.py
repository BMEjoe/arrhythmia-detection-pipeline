"""Phase 3 estimator building blocks (experiments/phase3_lle/estimators.py)."""
import sys
import pathlib

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))
import final_pipeline as fp  # noqa: E402
from experiments.phase3_lle import estimators as E  # noqa: E402
import upo_systems as us  # noqa: E402


def _logistic(n, x0=0.3, burn=500):
    x, out = x0, np.empty(n)
    for i in range(n + burn):
        x = 4.0 * x * (1.0 - x)
        if i >= burn:
            out[i - burn] = x
    return out


@pytest.mark.parametrize("tau,m", [(1, 2), (1, 4), (3, 5)])
@pytest.mark.parametrize("which", ["logistic", "henon", "noise"])
def test_divergence_curve_matches_rosenstein(which, tau, m):
    x = {"logistic": _logistic(300), "henon": us.henon_series(300),
         "noise": np.random.default_rng(3).standard_normal(300)}[which]
    X = fp._embed_backward(x, tau, m)
    _, d = fp.rosenstein_lle(X, theiler=2, require_valid_fit=False)
    y, _ = E.divergence_curve(X, 2)
    ref = np.asarray(d["mean_log_divergence"])
    assert np.array_equal(np.isfinite(ref), np.isfinite(y))
    np.testing.assert_allclose(y[np.isfinite(y)], ref[np.isfinite(ref)], rtol=0, atol=1e-12)


def test_iaaft_preserves_amplitudes_and_approximately_spectrum():
    x = us.henon_series(256)
    s = E.iaaft_surrogate(x, np.random.default_rng(0))
    np.testing.assert_array_equal(np.sort(s), np.sort(x))
    a, b = np.abs(np.fft.rfft(x)), np.abs(np.fft.rfft(s))
    assert np.linalg.norm(a - b) / np.linalg.norm(a) < 0.05
    assert not np.array_equal(s, x)


def test_window_rng_is_deterministic_and_window_specific():
    x = _logistic(64)
    assert E.window_rng(x).integers(1 << 30) == E.window_rng(x).integers(1 << 30)
    assert E.window_rng(x).integers(1 << 30) != E.window_rng(x + 1e-9).integers(1 << 30)


def test_saturation_fit_recovers_linear_rise_then_plateau():
    k = np.arange(40, dtype=float)
    y = np.minimum(0.7 * k, 5.0)
    slope, k_end = E.saturation_fit(y, k_start=1, frac=0.5)
    assert slope == pytest.approx(0.7)
    assert 3 <= k_end <= 12
