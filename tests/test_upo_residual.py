"""Level C PROJECT EXTENSION verification: residual, R^2 and support-ratio gates."""
import copy
import sys
import pathlib

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import upo_systems as us  # noqa: E402
from upo_systems import fp  # noqa: E402


def _peak(points, count=100):
    pts = np.atleast_2d(points)
    return {"orbit_points": pts, "histogram_count": count}


def test_local_affine_model_is_exact_on_linear_lozi_branch():
    X = fp._embed_backward(us.lozi_series(4096), 1, 2)
    q = np.array([0.6, 0.2])
    m = fp._local_affine_model(X, q, 10, 1)
    assert np.isclose(m["F_q"][0], 1 - us.LOZI_A * 0.6 + us.LOZI_B * 0.2, atol=1e-8)
    assert m["F_q"][1] == 0.6
    assert m["r2"] > 1 - 1e-10
    assert np.allclose(m["jacobian"][0], [-us.LOZI_A, us.LOZI_B], atol=1e-8)


def test_residual_is_small_at_exact_fixed_point_and_large_away_from_it():
    X = fp._embed_backward(us.lozi_series(4096), 1, 2)
    xs = us.lozi_fixed_point()
    h = np.array([1.0, 1.0, 100.0])
    on = fp.candidate_verification(X, _peak([xs, xs]), h, us.henon_config(), 1)
    off = fp.candidate_verification(X, _peak([xs + 0.2, xs + 0.2]), h, us.henon_config(), 1)
    assert on["residual"] < 1e-6
    assert off["residual"] > 0.3
    assert "residual" in off["failed_gates"] and "residual" not in on["failed_gates"]


def test_residual_is_normalised_by_series_scale():
    x = us.lozi_series(4096)
    xs = us.lozi_fixed_point()
    h = np.array([100.0])
    a = fp.candidate_verification(fp._embed_backward(x, 1, 2), _peak([xs + 0.1] * 2), h,
                                  us.henon_config(), 1)
    b = fp.candidate_verification(fp._embed_backward(10 * x, 1, 2), _peak([10 * (xs + 0.1)] * 2),
                                  h, us.henon_config(), 1)
    assert np.isclose(a["residual"], b["residual"], rtol=1e-6)


def test_period2_residual_uses_cycle_successor():
    """For p = 2 the residual compares F_hat(z*(k)) with z*(k+1), not z*(k)."""
    a, b = us.henon_period2_orbit()
    X = us.henon_embedded(8192)
    pts = fp.cyclic_orbit_matrix(np.array([a, b]), 2)
    ver = fp.candidate_verification(X, _peak(pts), np.array([100.0]), us.henon_config(), 1)
    assert ver["residual"] < 0.05


def test_support_ratio_is_peak_count_over_median_occupied_cell():
    h = np.array([0, 2, 4, 6, 0, 40], dtype=float)
    ver = fp.candidate_verification(us.henon_embedded(), _peak([0.93, 0.93], count=40), h,
                                    us.henon_config(), 1)
    assert np.isclose(ver["support_ratio"], 40 / 5.0)


@pytest.mark.parametrize("gate,cfg_kw", [
    ("residual", {"verify_max_residual": 0.0}),
    ("r2", {"verify_min_r2": 1.01}),
    ("support_ratio", {"verify_min_support_ratio": 1e9}),
])
def test_each_gate_can_fail_independently(gate, cfg_kw):
    r = copy.deepcopy(us.henon_period1())
    fp.apply_verification_gates(r, us.henon_embedded(), us.henon_config(**cfg_kw))
    assert r["verification_gated_candidates"] == []
    for c in r["verification_failed_peaks"]:
        assert gate in c["verification"]["failed_gates"]


def test_model_unavailable_when_too_few_points():
    X = us.henon_embedded()[:8]
    ver = fp.candidate_verification(X, _peak([0.9, 0.9]), np.array([10.0]), us.henon_config(), 1)
    assert ver["passes"] is False and ver["failed_gates"] == ["model_unavailable"]


def test_verification_partitions_source_peaks_without_modifying_them():
    r = copy.deepcopy(us.henon_period2())
    before = copy.deepcopy(r["source_peak_candidates"])
    fp.apply_verification_gates(r, us.henon_embedded(), us.henon_config())
    assert r["verification_assessed"] is True
    n_pass = len(r["verification_gated_candidates"])
    n_fail = len(r["verification_failed_peaks"])
    assert n_pass + n_fail == len(before)
    assert len(r["source_peak_candidates"]) == len(before)
    for a, b in zip(before, r["source_peak_candidates"]):
        assert a["level"] == b["level"] == "A"
        assert np.array_equal(a["orbit_coordinates"], b["orbit_coordinates"])
        assert "verification" not in b


def test_level_C_entries_are_labelled_project_extension():
    r = us.henon_period2_verified()
    for c in r["verification_gated_candidates"] + r["verification_failed_peaks"]:
        assert c["level"] == "C"
        assert c["verification"]["provenance"] == "PROJECT EXTENSION"


def test_henon_fixed_points_pass_default_gates():
    r = copy.deepcopy(us.henon_period1())
    fp.apply_verification_gates(r, us.henon_embedded(), us.henon_config())
    passed = [c["scalar_location"] for c in r["verification_gated_candidates"]]
    for fpnt in us.henon_fixed_points():
        assert abs(us.nearest(passed, fpnt) - fpnt) < 0.03


def _published_period2_peaks(result):
    return [c for c in result["source_peak_candidates"]
            if c["minimal_period"] == 2
            and us.orbit_distance(c["orbit_coordinates"], us.PRE_HENON_PERIOD2) < 0.03]


@pytest.mark.parametrize("form", ["prl_norm", "pre_tensor"])
def test_level_C_residual_gate_can_reject_published_period2_source_peak(form):
    """PROJECT EXTENSION (Level C) behaviour, moved out of the source-fidelity
    suite.  The residual gate and its threshold are project-defined; they are
    NOT a So et al. source requirement, and nothing here is evidence that the
    So et al. method accepts or rejects this orbit.

    Thresholds are set explicitly from the peaks' own measured residuals, so the
    test does not depend on the reconstructed default (0.05): with a threshold
    below every residual the gate rejects the published 2-cycle peaks, with one
    above every residual it does not -- and in both cases the peaks remain
    unchanged Level A source peaks."""
    base = us.henon_period2(form)
    published = _published_period2_peaks(base)
    assert published  # the SOURCE detection finds the published 2-cycle
    X = us.henon_embedded()
    cfg0 = us.henon_config(form)
    residuals = [fp.candidate_verification(X, c, base["histogram"], cfg0, 1)["residual"]
                 for c in published]
    assert all(np.isfinite(residuals))

    strict = copy.deepcopy(base)
    fp.apply_verification_gates(strict, X, us.henon_config(
        form, verify_max_residual=0.5 * min(residuals), verify_min_r2=-np.inf,
        verify_min_support_ratio=0.0))
    for c in _published_period2_peaks(strict):
        failed = [f for f in strict["verification_failed_peaks"]
                  if np.array_equal(f["orbit_coordinates"], c["orbit_coordinates"])]
        assert len(failed) == 1 and failed[0]["verification"]["failed_gates"] == ["residual"]
        assert c["level"] == "A" and c["provenance"] == "SOURCE"

    loose = copy.deepcopy(base)
    fp.apply_verification_gates(loose, X, us.henon_config(
        form, verify_max_residual=2.0 * max(residuals), verify_min_r2=-np.inf,
        verify_min_support_ratio=0.0))
    passed = {tuple(c["orbit_coordinates"]) for c in loose["verification_gated_candidates"]}
    for c in _published_period2_peaks(loose):
        assert tuple(c["orbit_coordinates"]) in passed

    for r in (strict, loose):
        assert len(r["source_peak_candidates"]) == len(base["source_peak_candidates"])
        for a, b in zip(base["source_peak_candidates"], r["source_peak_candidates"]):
            assert np.array_equal(a["orbit_coordinates"], b["orbit_coordinates"])
