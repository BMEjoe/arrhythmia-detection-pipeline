"""Robustness of source detection to seeds, data length, kappa, noise and rescaling."""
import sys
import pathlib
from dataclasses import replace

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import upo_systems as us  # noqa: E402
from upo_systems import fp  # noqa: E402

FP0 = us.PRE_HENON_FIXED_POINTS[0]


@pytest.mark.parametrize("rng_seed", [3, 4, 5])
@pytest.mark.parametrize("form", ["prl_norm", "pre_tensor"])
def test_fixed_point_stable_across_randomization_seeds(rng_seed, form):
    locs = us.period1_locations(us.henon_period1(form, rng_seed=rng_seed))
    assert abs(us.nearest(locs, FP0) - FP0) < 0.03


@pytest.mark.parametrize("n", [512, 2048])
def test_fixed_point_stable_across_data_length(n):
    locs = us.period1_locations(us.henon_period1(n=n))
    assert abs(us.nearest(locs, FP0) - FP0) < 0.03


@pytest.mark.parametrize("kappa", [1.0, 5.0])
def test_fixed_point_stable_across_kappa(kappa):
    locs = us.period1_locations(us.henon_period1(kappa=kappa))
    assert abs(us.nearest(locs, FP0) - FP0) < 0.03


@pytest.mark.parametrize("seed", [1, 2])
def test_fixed_point_stable_across_trajectories(seed):
    locs = us.period1_locations(us.henon_period1(seed=seed))
    assert abs(us.nearest(locs, FP0) - FP0) < 0.03


def test_fixed_point_survives_small_observational_noise():
    x = us.henon_series(1024, noise=0.01)
    X = fp._embed_backward(x, 1, 2)
    r = fp.detect_so_fixed_points(X, 1, us.henon_config(so_jacobian_neighbors=7), np.random.default_rng(1))
    assert abs(us.nearest(us.period1_locations(r), FP0) - FP0) < 0.05


@pytest.mark.parametrize("a,b", [(2.0, 0.0), (0.05, 0.8)])
def test_peaks_are_equivariant_under_affine_rescaling(a, b):
    """RR-like rescaling: peak locations transform as a x + b."""
    x = us.henon_series(1024)
    X = fp._embed_backward(a * x + b, 1, 2)
    r = fp.detect_so_fixed_points(X, 1, us.henon_config(), np.random.default_rng(1))
    locs = us.period1_locations(r)
    assert abs(us.nearest(locs, a * FP0 + b) - (a * FP0 + b)) < 0.03 * a


def test_ikeda_fixed_point_pre_section_IV():
    """PRE Sec. IV (noise-free Ikeda, a = 0.75, b = 9): fixed point at v* = 0.537."""
    X = fp._embed_backward(us.ikeda_series(1024), 1, 4)
    cfg = replace(fp.CFG, so_jacobian_neighbors=7, so_random_R=100, so_kappa=8.0,
                  so_randomization="pre_tensor")
    r = fp.detect_so_fixed_points(X, 1, cfg, np.random.default_rng(0))
    assert abs(us.nearest(us.period1_locations(r), us.IKEDA_FIXED_POINT_V) - us.IKEDA_FIXED_POINT_V) < 0.02


def test_noise_has_source_peaks_but_small_rJ():
    """Stochastic data: histogram peaks exist but are not surrogate-significant."""
    r = us.noise_significance()
    assert r["significance"]["rJ"] < 2.0
    assert r["significance"]["J_W"] > 0.05


def test_significance_uses_observed_histogram_edges_for_surrogates():
    r = us.henon_significance()
    sig = r["significance"]
    assert sig["surrogate_mean_histogram"].shape == r["histogram"].shape
    assert sig["n_surrogates"] == 12
