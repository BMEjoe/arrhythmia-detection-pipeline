"""Mandatory source-fidelity checks against So et al. PRL 76, 4705 (1996) and
PRE 55, 5398 (1997).  Numbered comments refer to the Phase 2D requirement list."""
import copy
import sys
import pathlib
from dataclasses import replace

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import upo_systems as us  # noqa: E402
from upo_systems import fp  # noqa: E402

FORMS = ["prl_norm", "pre_tensor"]


def test_citations_name_the_actual_authors():
    assert "So, Ott, Schiff, Kaplan, Sauer & Grebogi" in fp.SO_PRL_CITATION
    assert "76, 4705 (1996)" in fp.SO_PRL_CITATION
    assert "So, Ott, Sauer, Gluckman, Grebogi & Schiff" in fp.SO_PRE_CITATION
    assert "55, 5398 (1997)" in fp.SO_PRE_CITATION


# 1. PRL norm-form transformation
def test_01_prl_norm_form_transformation():
    rng = np.random.default_rng(0)
    z, Fz = rng.standard_normal(3), rng.standard_normal(3)
    J = fp._companion(rng.standard_normal(3))
    R = rng.uniform(-1, 1, (3, 3))
    kappa = 3.0
    S = J + kappa * R * np.sum(np.abs(Fz - z))
    expected = np.linalg.inv(np.eye(3) - S) @ (Fz - S @ z)
    assert np.allclose(fp.so_fixed_point_transform(z, Fz, J, R, kappa, "prl_norm"), expected)


# 2. PRE tensor-form transformation
def test_02_pre_tensor_form_transformation():
    rng = np.random.default_rng(1)
    z, Fz = rng.standard_normal(3), rng.standard_normal(3)
    J = fp._companion(rng.standard_normal(3))
    T = rng.uniform(-1, 1, (3, 3, 3))
    kappa = 2.0
    S = J + kappa * np.tensordot(T, Fz - z, axes=([2], [0]))
    expected = np.linalg.solve(np.eye(3) - S, Fz - S @ z)
    assert np.allclose(fp.so_fixed_point_transform_tensor(z, Fz, J, T, kappa), expected)


# 3. Local Jacobian (PRE Eq. 20)
def test_03_local_jacobian_pre_eq_20():
    X = us.henon_embedded()
    idx = 200
    nbrs = fp._spatial_neighbors(X, idx, K=2, exclude=1, valid_limit=len(X) - 1)
    J = fp._delay_map_jacobian_from_local_fit(X, idx, nbrs, step=1)
    # Eq. (20): w_1^k(n+1) - z_1(n+1) = grad f . (w^k - z(n)); M = d = 2 is exactly determined.
    lhs = X[nbrs + 1, 0] - X[idx + 1, 0]
    assert np.allclose((X[nbrs] - X[idx]) @ J[0], lhs)
    assert np.array_equal(J[1], [1.0, 0.0])


# 4. Correct cyclic indexing
def test_04_cyclic_indexing_k_minus_j():
    s = np.array([1.0, 2.0, 3.0])
    M = fp.cyclic_orbit_matrix(s, 3)
    for k in range(3):
        for j in range(3):
            assert M[k, j] == s[(k - j) % 3]


# 5. Period-p block system
def test_05_period_p_block_system_solution_satisfies_eq_29():
    rng = np.random.default_rng(2)
    p, d = 3, 2
    Z = [rng.standard_normal(d) for _ in range(p)]
    FZ = [rng.standard_normal(d) for _ in range(p)]
    J = [fp._companion(rng.standard_normal(d)) for _ in range(p)]
    R = [rng.uniform(-1, 1, (d, d)) for _ in range(p)]
    Zh = fp._so_period_transform(Z, FZ, J, R, 1.0)
    for k in range(p):
        S = J[k] + 1.0 * R[k] * np.sum(np.abs(FZ[k] - Z[(k + 1) % p]))
        # linearised Eq. (29): F(z_k) = zhat_{k+1} + S_k (z_k - zhat_k)
        assert np.allclose(FZ[k], Zh[(k + 1) % p] + S @ (Z[k] - Zh[k]))


def _obs(hist, peak_cell=None):
    r = {"histogram": np.asarray(hist, float), "edges": [np.arange(len(hist) + 1, dtype=float)],
         "source_peak_candidates": []}
    if peak_cell is not None:
        r["source_peak_candidates"] = [{"peak_cell": (peak_cell,)}]
    return r


def _sur(points):
    return {"reduced": np.asarray(points, float)[:, None] + 0.5}


# 6. Signed W
def test_06_W_is_signed_not_absolute():
    # Every surrogate histogram is [2, 2, 2, 2]; the observed one is [0, 0, 1, 0],
    # so obs - mean = [-2, -2, -1, -2]: signed W = -1, whereas max|.| would be 2.
    obs = _obs([0, 0, 1, 0])
    surs = [_sur([0, 0, 1, 1, 2, 2, 3, 3]), _sur([0, 0, 1, 1, 2, 2, 3, 3])]
    sig = fp.so_surrogate_significance(obs, surs)
    assert np.allclose(sig["surrogate_mean_histogram"], [2, 2, 2, 2])
    assert sig["W"] == -1.0
    assert np.allclose(sig["surrogate_W"], [0.0, 0.0])


# 7. W0 median behaviour
def test_07_W0_is_median_of_surrogate_W():
    obs = _obs([5, 0, 0, 0])
    rng = np.random.default_rng(3)
    surs = [_sur(rng.integers(0, 4, 40)) for _ in range(9)]
    sig = fp.so_surrogate_significance(obs, surs)
    assert np.isclose(sig["W0"], np.median(sig["surrogate_W"]))
    assert np.isclose(sig["rJ"], sig["W"] / sig["W0"])


# 8. J(W0) = 0.5
def test_08_J_of_W0_is_one_half():
    Ws = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    assert fp.so_J(Ws, np.median(Ws)) == 0.5
    # Surrogate i puts c_i = i+1 points in bin i, so W_i = 0.9 c_i are distinct.
    surs = [_sur([i] * (i + 1)) for i in range(10)]
    sig = fp.so_surrogate_significance(_obs([9] + [0] * 9), surs)
    assert np.allclose(np.sort(sig["surrogate_W"]), 0.9 * np.arange(1, 11))
    assert np.isclose(sig["W0"], 0.9 * 5.5)
    assert sig["J_W0"] == 0.5


def test_08b_J_is_fraction_strictly_exceeding_and_p_value_separate():
    Ws = np.array([1.0, 2.0, 2.0, 5.0])
    assert fp.so_J(Ws, 2.0) == 0.25
    sig = fp.so_surrogate_significance(_obs([3, 0, 0], peak_cell=0),
                                       [_sur([0, 1, 2]), _sur([0, 0, 1, 2]), _sur([1, 1, 2])])
    assert "p_value_finite" in sig and "J_W" in sig
    assert sig["p_value_finite"] >= 1.0 / (sig["n_surrogates"] + 1)


# 9. Gaussian-scaled phase-shuffle surrogates
def test_09_gaussian_scaled_phase_shuffle_properties():
    x = us.henon_series(1024)
    s = fp.gaussian_scaled_phase_shuffle(x, rng=np.random.default_rng(0))
    assert np.array_equal(np.sort(s), np.sort(x))  # identical amplitude distribution
    assert not np.array_equal(s, x)
    # the deterministic one-step structure is destroyed
    assert abs(np.corrcoef(x[:-1] ** 2, x[1:])[0, 1]) > 0.5 > abs(np.corrcoef(s[:-1] ** 2, s[1:])[0, 1])


def test_09b_phase_shuffle_approximately_preserves_power_spectrum():
    """AAFT (Theiler et al. 1992): approximately the same power spectrum."""
    rng = np.random.default_rng(1)
    e = rng.standard_normal(4096)
    x = np.zeros(4096)
    for t in range(2, 4096):
        x[t] = 1.6 * x[t - 1] - 0.9 * x[t - 2] + e[t]  # sharp spectral peak
    s = fp.gaussian_scaled_phase_shuffle(x, rng=np.random.default_rng(2))
    band = lambda p: np.array([b.sum() for b in np.array_split(p[1:], 32)])  # noqa: E731
    px = band(np.abs(np.fft.rfft(x - x.mean())) ** 2)
    ps = band(np.abs(np.fft.rfft(s - s.mean())) ** 2)
    assert np.argmax(px) == np.argmax(ps)
    assert np.corrcoef(np.log(px), np.log(ps))[0, 1] > 0.95


# 10. Source peak remains a source peak even if verification rejects it.
# The verification gates are PROJECT EXTENSIONS (Level C), not So et al.
# requirements; the Level C behaviour itself is tested in test_upo_residual.py.
def test_10_source_peaks_are_unaffected_by_project_verification_gates():
    """Source invariant: Level C (project) gates cannot remove or relabel a
    Level A source peak.  All peaks are forced to fail via an explicit strict
    project gate; nothing here is a source-derived accept/reject result."""
    r = copy.deepcopy(us.henon_period1())
    before = copy.deepcopy(r["source_peak_candidates"])
    fp.apply_verification_gates(r, us.henon_embedded(), us.henon_config(verify_max_residual=0.0))
    assert r["verification_gated_candidates"] == []
    assert len(r["source_peak_candidates"]) == len(r["verification_failed_peaks"]) == len(before)
    for a, b in zip(before, r["source_peak_candidates"]):
        assert b["level"] == "A" and b["provenance"] == "SOURCE"
        assert np.array_equal(a["orbit_coordinates"], b["orbit_coordinates"])


# 11. Published / skewed Henon fixed points (PRE Table I)
def test_11_newton_reference_matches_pre_table_I():
    for mine, published in zip(us.henon_fixed_points(), us.PRE_HENON_FIXED_POINTS):
        assert abs(mine - published) < 1e-4
    assert np.allclose(us.henon_period2_orbit(), us.PRE_HENON_PERIOD2, atol=1e-4)


@pytest.mark.parametrize("form", FORMS)
def test_11_henon_fixed_points_detected(form):
    """SOURCE: the published PRE Table I values.  IMPLEMENTATION-DEFINED: the
    0.03 tolerance (limited by the project's 30-bin grid and bin-mean location
    estimator; PRE reports +-0.002) and the exact peak count of 2 (depends on
    the project's padded range and median+MAD peak threshold)."""
    locs = us.period1_locations(us.henon_period1(form))
    for published in us.PRE_HENON_FIXED_POINTS:
        assert abs(us.nearest(locs, published) - published) < 0.03
    assert len(locs) == 2


# 12. Published / skewed Henon period-2 structure
@pytest.mark.parametrize("form", FORMS)
def test_12_henon_period2_orbit_detected(form):
    peaks = us.henon_period2(form)["source_peak_candidates"]
    two = [c for c in peaks if c["minimal_period"] == 2]
    # SOURCE: the published 2-cycle (PRE Table I).  IMPLEMENTATION-DEFINED: the
    # 0.03 tolerance (30-bin grid) and "strongest peak" ranking.
    # the strongest period-2 peak is the published 2-cycle (PRE Table I)
    strongest = max(two, key=lambda c: c["histogram_count"])
    assert us.orbit_distance(strongest["orbit_coordinates"], us.PRE_HENON_PERIOD2) < 0.03
    # PROJECT convention: canonical_orbit starts at the largest element; for
    # p = 2 this coincides with the Table I ordering z*(1) = (1.3014, -0.5604),
    # but Table I does not use this convention in general (e.g. its p = 6 row).
    assert np.max(np.abs(strongest["canonical_orbit"] - np.array(us.PRE_HENON_PERIOD2))) < 0.03


# 13. Source vs extension stability distinction.  The period-1 "source_stability"
# is a HYBRID PRL/PRE implementation (PRL averaging, PRE Eq. 20 Jacobians), not
# the literal PRL temporal S_n; the monodromy is a PROJECT EXTENSION.
def test_13_source_vs_extension_stability_are_distinct_methods():
    r = copy.deepcopy(us.henon_period1())
    fp.apply_verification_gates(r, us.henon_embedded(), us.henon_config())
    for c in r["verification_gated_candidates"]:
        src, ext = c["source_stability"], c["extension_stability"]
        assert src["provenance"].startswith("HYBRID PRL/PRE")
        assert src["literal_prl_temporal_S_n"] is False
        assert ext["provenance"].startswith("PROJECT EXTENSION")
        assert not np.allclose(src["mean_S"], ext["monodromy"])


# 14. one_sample vs tau_step metadata
def test_14_one_sample_metadata_and_tau_rejection():
    info = fp.resolve_upo_map("one_sample", 1)
    assert info["source_faithful_map"] is True and info["step"] == 1
    with pytest.raises(ValueError):
        fp.resolve_upo_map("one_sample", 3)
    with pytest.raises(ValueError):
        fp.detect_so_fixed_points(us.henon_embedded(), tau=2, config=fp.CFG, map_mode="one_sample")


def test_14_tau_step_is_labelled_project_extension():
    info = fp.resolve_upo_map("tau_step", 4)
    assert info["source_faithful_map"] is False and info["step"] == 4
    assert "PROJECT EXTENSION" in info["map_label"]
    X = fp._embed_backward(us.henon_series(600), 3, 2)
    r = fp.detect_so_fixed_points(X, 3, us.henon_config(so_random_R=5), map_mode="tau_step")
    assert r["source_faithful_map"] is False and "PROJECT EXTENSION" in r["detection_label"]


# 15. Feature-mode contract
def test_15_feature_mode_contract():
    assert set(fp.UPO_FEATURE_CONTRACT) == {"source", "significant_source", "extended"}
    assert "source_peak_count" in fp.UPO_FEATURE_CONTRACT["source"]
    assert "source_rJ" in fp.UPO_FEATURE_CONTRACT["significant_source"]
    assert "verified_unstable_count" in fp.UPO_FEATURE_CONTRACT["extended"]


# 16. Significance not assessed vs assessed-empty
def test_16_not_assessed_is_none_and_assessed_empty_is_list():
    """SOURCE: None-vs-[] is the API distinction.  IMPLEMENTATION-DEFINED: that
    no noise peak is significant uses the project rule J < 0.05 with 12
    surrogates (PRL used 20, PRE 50)."""
    assert us.henon_period1()["significant_uop_candidates"] is None
    noise = us.noise_significance()
    assert noise["significance_assessed"] is True
    assert noise["significant_uop_candidates"] == []
    assert noise["source_peak_candidates"]  # noise still has source peaks (Level A != Level B)


def test_16b_deterministic_structure_is_significant():
    """SOURCE quantities (W, W0, J, r_J).  IMPLEMENTATION-DEFINED: 12 surrogates
    (PRL used 20, PRE 50), the test threshold r_J > 3, the 0.03 location
    tolerance, and Level B membership via the project rule J < 0.05."""
    r = us.henon_significance()
    sig = r["significance"]
    assert sig["rJ"] > 3 and sig["J_W"] == 0.0
    locs = [c["scalar_location"] for c in r["significant_uop_candidates"]]
    assert abs(us.nearest(locs, us.PRE_HENON_FIXED_POINTS[0]) - us.PRE_HENON_FIXED_POINTS[0]) < 0.03
    assert all(c["level"] == "B" for c in r["significant_uop_candidates"])
