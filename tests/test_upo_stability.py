"""Stability: SOURCE PRL period-1 method vs PROJECT EXTENSION monodromy."""
import sys
import pathlib

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import upo_systems as us  # noqa: E402
from upo_systems import fp  # noqa: E402


def _analytic_moduli(u):
    return np.sort(np.abs(np.linalg.eigvals(us.henon_jacobian(u))))[::-1]


def test_source_stability_averages_member_jacobians():
    Js = [fp._companion([2.0, 0.3]), fp._companion([-1.0, 0.1]), None, fp._companion([0.5, 0.2])]
    st = fp.source_period1_stability(Js, [0, 1, 1, 2, 3])
    expected = np.mean([Js[0], Js[1], Js[3]], axis=0)
    assert np.allclose(st["mean_S"], expected)
    assert st["n_points_averaged"] == 3
    assert np.allclose(np.sort(np.abs(st["eigenvalues"]))[::-1], st["lyapunov_numbers"])


def test_source_stability_none_without_jacobians():
    assert fp.source_period1_stability([None, None], [0, 1]) is None


@pytest.mark.parametrize("form", ["prl_norm", "pre_tensor"])
@pytest.mark.parametrize("which", [0, 1])
def test_hybrid_period1_stability_matches_analytic_henon_multipliers(form, which):
    """Validates the HYBRID PRL/PRE stability (PRE Eq. 20 Jacobians), not the
    literal PRL temporal-S_n procedure."""
    fpnt = us.henon_fixed_points()[which]
    peaks = us.henon_period1(form)["source_peak_candidates"]
    c = min(peaks, key=lambda c: abs(c["scalar_location"] - fpnt))
    st = c["source_stability"]
    exact = _analytic_moduli(fpnt)
    assert abs(st["lyapunov_numbers"][0] - exact[0]) / exact[0] < 0.15
    assert st["unstable"] and st["saddle"]


def test_source_stability_is_labelled_hybrid_prl_pre():
    for c in us.henon_period1()["source_peak_candidates"]:
        st = c["source_stability"]
        assert st["provenance"] == fp.SOURCE_STABILITY_PROVENANCE
        assert st["jacobian_method"] == "PRE Eq. 20 spatial-neighbor least squares"
        assert st["literal_prl_temporal_S_n"] is False


def test_stability_provenance_says_hybrid_pre_eq20_and_not_literal_prl():
    """Regression for the audit finding: the period-1 stability must not be
    presented as the literal PRL temporal-S_n procedure."""
    prov = fp.SOURCE_STABILITY_PROVENANCE
    assert prov.startswith("HYBRID PRL/PRE")
    assert "PRE Eq. 20" in prov and "spatial-neighbor least-squares" in prov
    assert "not the literal PRL temporal S_n construction" in prov
    assert not prov.startswith("SOURCE")
    doc = fp.source_period1_stability.__doc__
    assert "HYBRID" in doc and "NOT the literal PRL" in doc


def test_monodromy_is_product_of_candidate_jacobians():
    X = us.henon_embedded(4096)
    pts = fp.cyclic_orbit_matrix(np.array(us.henon_period2_orbit()), 2)
    st = fp.monodromy_stability(X, pts, neighborhood=30, step=1)
    J0 = fp.estimate_jacobian_at_candidate(X, pts[0], 30, 1)
    J1 = fp.estimate_jacobian_at_candidate(X, pts[1], 30, 1)
    assert np.allclose(st["monodromy"], J1 @ J0)
    assert st["provenance"].startswith("PROJECT EXTENSION")


def test_period2_monodromy_matches_analytic_orbit_multipliers():
    a, b = us.henon_period2_orbit()
    exact = np.sort(np.abs(np.linalg.eigvals(us.henon_jacobian(b) @ us.henon_jacobian(a))))[::-1]
    st = fp.monodromy_stability(us.henon_embedded(8192), fp.cyclic_orbit_matrix(np.array([a, b]), 2))
    assert abs(st["multiplier_moduli"][0] - exact[0]) / exact[0] < 0.15
    assert st["unstable"]


def test_fixed_point_stability_is_extension_and_period1_monodromy():
    X = us.henon_embedded(4096)
    z = np.full(2, us.henon_fixed_points()[0])
    a = fp.fixed_point_stability(X, z)
    b = fp.monodromy_stability(X, z[None, :])
    assert np.allclose(a["monodromy"], b["monodromy"])
    assert a["provenance"].startswith("PROJECT EXTENSION")


def test_fixed_point_stability_requires_at_least_two_dimensions():
    assert fp.fixed_point_stability(np.random.default_rng(0).standard_normal((50, 1)), [0.0]) is None


def test_stable_fixed_point_is_not_flagged_unstable():
    Js = [fp._companion([0.5, 0.1])] * 3
    st = fp.source_period1_stability(Js, [0, 1, 2])
    assert not st["unstable"] and not st["saddle"]


def test_source_and_extension_stability_are_stored_separately():
    r = us.henon_period2_verified()
    for c in r["verification_gated_candidates"] + r["verification_failed_peaks"]:
        assert "extension_stability" in c
        assert c["source_stability"] is None  # period 2: no source procedure
    r1 = us.henon_period1()
    for c in r1["source_peak_candidates"]:
        assert "extension_stability" not in c  # only Level C carries extension stability


def test_verified_unstable_requires_gate_pass_and_unstable_monodromy():
    r = us.henon_period2_verified()
    for c in r["verification_gated_candidates"]:
        assert c["verified_unstable"] == bool(c["extension_stability"] and c["extension_stability"]["unstable"])
    for c in r["verification_failed_peaks"]:
        assert c["verified_unstable"] is False
