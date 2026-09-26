"""Peak coverage: a PROJECT-derived metric, not a density of periodic orbits."""
import sys
import pathlib
from dataclasses import replace

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import upo_systems as us  # noqa: E402
from upo_systems import fp  # noqa: E402


def test_no_locations_gives_zero_coverage():
    cov = fp.peak_coverage(np.empty((0, 2)), us.henon_embedded(), fp.CFG, 0)
    assert cov["coverage"] == 0.0 and cov["n_locations"] == 0


def test_uses_200_reference_points_when_available():
    cov = fp.peak_coverage(np.array([[0.9, 0.9]]), us.henon_embedded(), fp.CFG, 0)
    assert fp.CFG.coverage_n_reference == 200
    assert cov["n_reference"] == 200


def test_reference_count_capped_by_attractor_size():
    X = us.henon_embedded()[:120]
    assert fp.peak_coverage(np.array([[0.9, 0.9]]), X, fp.CFG, 0)["n_reference"] == 120


def test_coverage_is_fraction_in_unit_interval():
    cov = fp.peak_coverage(np.array([[0.9, 0.9], [-1.9, -1.9]]), us.henon_embedded(), fp.CFG, 0)
    assert 0.0 <= cov["coverage"] <= 1.0


def test_locations_at_every_attractor_point_give_full_coverage():
    X = us.henon_embedded(300)
    assert fp.peak_coverage(X, X, fp.CFG, 0)["coverage"] == 1.0


def test_far_location_gives_zero_coverage():
    assert fp.peak_coverage(np.array([[50.0, 50.0]]), us.henon_embedded(), fp.CFG, 0)["coverage"] == 0.0


def test_radius_definition_matches_documented_procedure():
    X = us.henon_embedded(400)
    locs = np.array([[0.93, 0.93]])
    cfg = fp.CFG
    cov = fp.peak_coverage(locs, X, cfg, 7)
    mins, span = X.min(0), np.maximum(X.max(0) - X.min(0), 1e-12)
    Xn = (X - mins) / span
    rng = np.random.default_rng(7)
    ref = Xn[rng.choice(len(Xn), 200, replace=False)]
    sample = Xn[rng.choice(len(Xn), min(500, len(Xn)), replace=False)]
    nn = []
    for i in range(len(sample)):
        d = np.linalg.norm(sample - sample[i], axis=1)
        d[i] = np.inf
        nn.append(np.partition(d, cfg.coverage_radius_neighbors_k)[cfg.coverage_radius_neighbors_k])
    radius = cfg.coverage_radius_multiplier * np.median(nn)
    ln = (locs - mins) / span
    expected = np.mean(np.min(np.linalg.norm(ref[:, None] - ln[None], axis=2), axis=1) <= radius)
    assert np.isclose(cov["radius"], radius)
    assert np.isclose(cov["coverage"], expected)


def test_coverage_is_deterministic_for_seed():
    X = us.henon_embedded()
    locs = np.array([[0.93, 0.93]])
    assert fp.peak_coverage(locs, X, fp.CFG, 3) == fp.peak_coverage(locs, X, fp.CFG, 3)


@pytest.mark.parametrize("a,b", [(2.0, 0.0), (0.01, 5.0)])
def test_coverage_invariant_to_affine_rescaling(a, b):
    X = us.henon_embedded()
    locs = np.array([[0.93, 0.93], [1.3, -0.56]])
    c1 = fp.peak_coverage(locs, X, fp.CFG, 0)["coverage"]
    c2 = fp.peak_coverage(a * locs + b, a * X + b, fp.CFG, 0)["coverage"]
    assert np.isclose(c1, c2)


def test_more_locations_never_reduce_coverage():
    X = us.henon_embedded()
    one = fp.peak_coverage(np.array([[0.93, 0.93]]), X, fp.CFG, 0)["coverage"]
    two = fp.peak_coverage(np.array([[0.93, 0.93], [-0.5, 1.3]]), X, fp.CFG, 0)["coverage"]
    assert two >= one


def test_coverage_definition_states_it_is_a_project_metric():
    cov = fp.peak_coverage(np.array([[0.93, 0.93]]), us.henon_embedded(), fp.CFG, 0)
    assert "PROJECT" in cov["definition"]
    assert "density" not in fp.UPO_FEATURE_DEFINITIONS["source_peak_coverage"].replace(
        "not a density", "")


def test_population_points_stack_all_orbit_points():
    peaks = [{"orbit_points": np.zeros((1, 2))}, {"orbit_points": np.ones((2, 2))}]
    assert fp._population_points(peaks, 2).shape == (3, 2)
    assert fp._population_points([], 2).shape == (0, 2)


def test_population_specific_coverage_in_upo_analysis():
    x = us.henon_series(400)
    cfg = replace(us.henon_config(so_random_R=20), cao_max_dim=12)
    upo = fp.run_upo_analysis(x, cfg)
    src = fp._population_points(upo["source_peak_candidates"], upo["embedding_dimension"])
    expected = fp.peak_coverage(src, upo["embedded"], cfg, rng_seed=cfg.random_seed)
    assert upo["coverage"]["source"]["coverage"] == expected["coverage"]
    assert upo["coverage"]["significant"] is None  # significance not assessed
    assert upo["coverage"]["verified_unstable"] is not None
