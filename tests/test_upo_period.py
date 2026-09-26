"""Period-p So transform (PRE Sec. III): block system, cyclic reduction, short scheme."""
import sys
import pathlib

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import upo_systems as us  # noqa: E402
from upo_systems import fp  # noqa: E402

FORMS = ["prl_norm", "pre_tensor"]


def test_cyclic_orbit_matrix_uses_k_minus_j_indexing():
    s = np.array([10.0, 20.0, 30.0])
    M = fp.cyclic_orbit_matrix(s, 4)
    expected = np.array([[10, 30, 20, 10],
                         [20, 10, 30, 20],
                         [30, 20, 10, 30]], dtype=float)
    assert np.array_equal(M, expected)


def test_cyclic_orbit_matrix_is_consistent_with_delay_dynamics():
    """z*(k+1)_{j+1} = z*(k)_j: each point is the delay-shift of its predecessor."""
    s = np.array([0.3, -1.2, 0.8, 2.0])
    M = fp.cyclic_orbit_matrix(s, 3)
    for k in range(4):
        assert np.array_equal(M[(k + 1) % 4, 1:], M[k, :-1])


def test_old_k_plus_j_indexing_is_not_the_delay_orbit_for_p3():
    s = np.array([0.1, 0.5, 0.9])
    good = fp.cyclic_orbit_matrix(s, 2)
    k = np.arange(3)[:, None]
    j = np.arange(2)[None, :]
    old = s[(k + j) % 3]
    assert not np.allclose(old, good)
    assert fp._cyclic_slab_reduction(old)[1] > 0.1


def test_slab_reduction_returns_first_components_and_zero_residual_on_orbit():
    s = np.array([0.2, 0.7, -0.4])
    s_out, resid = fp._cyclic_slab_reduction(fp.cyclic_orbit_matrix(s, 5))
    assert np.array_equal(s_out, s) and resid == 0.0


def test_batched_cyclic_residual_matches_scalar():
    Zh = np.random.default_rng(0).standard_normal((7, 3, 2))
    batch = fp._batched_cyclic_residual(Zh)
    assert np.allclose(batch, [fp._cyclic_slab_reduction(z)[1] for z in Zh])


def test_block_system_layout_is_pre_eq_30():
    p, d = 3, 2
    rng = np.random.default_rng(1)
    Z = rng.standard_normal((p, d))
    FZ = rng.standard_normal((p, d))
    S = [rng.standard_normal((d, d)) for _ in range(p)]
    B, rhs = fp._so_period_block_system(Z, FZ, S)
    for k in range(p):
        rows = slice(k * d, (k + 1) * d)
        assert np.array_equal(B[rows, k * d:(k + 1) * d], -S[k])
        nxt = (k + 1) % p
        assert np.array_equal(B[rows, nxt * d:(nxt + 1) * d], np.eye(d))
        others = [c for c in range(p) if c not in (k, nxt)]
        for c in others:
            assert np.all(B[rows, c * d:(c + 1) * d] == 0)
        assert np.allclose(rhs[rows], FZ[k] - S[k] @ Z[k])


@pytest.mark.parametrize("form", FORMS)
@pytest.mark.parametrize("seed", [0, 1])
def test_exact_period2_orbit_is_invariant_under_transform(form, seed):
    """G(Z*, R) = Z* independent of R, with S(z_k, z_{k+1}, R_k) (PRE Sec. III A)."""
    a, b = us.henon_period2_orbit()
    Z = fp.cyclic_orbit_matrix(np.array([a, b]), 2)
    FZ = [np.array([us.henon_f(z[0], z[1]), z[0]]) for z in Z]
    assert np.allclose(FZ[0], Z[1]) and np.allclose(FZ[1], Z[0])
    J = [us.henon_jacobian(z[0]) for z in Z]
    rng = np.random.default_rng(seed)
    R = [fp._draw_R(rng, 1, 2, form)[0] for _ in range(2)]
    out = fp._so_period_transform(list(Z), FZ, J, R, 3.0, form)
    assert np.allclose(out, Z, atol=1e-10)


def test_period_transform_with_p1_equals_fixed_point_transform():
    rng = np.random.default_rng(3)
    z, Fz = rng.standard_normal(2), rng.standard_normal(2)
    J, R = fp._companion(rng.standard_normal(2)), rng.uniform(-1, 1, (2, 2))
    a = fp._so_period_transform([z], [Fz], [J], [R], 2.0)
    b = fp.so_fixed_point_transform(z, Fz, J, R, 2.0)
    assert np.allclose(a[0], b)


@pytest.mark.parametrize("form", FORMS)
def test_batched_period_transform_matches_scalar(form):
    rng = np.random.default_rng(4)
    p, d = 2, 3
    Z = [rng.standard_normal(d) for _ in range(p)]
    FZ = [rng.standard_normal(d) for _ in range(p)]
    J = [fp._companion(rng.standard_normal(d)) for _ in range(p)]
    Rs = [fp._draw_R(np.random.default_rng(5 + k), 4, d, form) for k in range(p)]
    batch = fp._batched_period_transform(Z, FZ, J, Rs, 1.5, form)
    for r in range(4):
        single = fp._so_period_transform(Z, FZ, J, [Rs[k][r] for k in range(p)], 1.5, form)
        assert np.allclose(batch[r], single)


def test_canonical_rotation_and_minimal_period():
    assert np.array_equal(fp.canonical_cyclic_rotation([1.0, 3.0, 2.0]), [3.0, 2.0, 1.0])
    assert fp.minimal_cyclic_period([0.93, 0.931], tol=0.01) == 1
    assert fp.minimal_cyclic_period([1.3, -0.56], tol=0.01) == 2
    assert fp.minimal_cyclic_period([0.5, 0.2, 0.5, 0.2], tol=1e-9) == 2


def test_period_p_rejects_invalid_periods():
    X = us.henon_embedded(200)
    with pytest.raises(ValueError):
        fp.detect_so_period_p(X, 1)
    with pytest.raises(ValueError):
        fp.detect_so_period_p(X, 4)


def test_period2_histogram_contains_fixed_points_with_minimal_period_1():
    """PRE Sec. III C: a complete period-2 set contains both fixed points."""
    peaks = us.henon_period2()["source_peak_candidates"]
    ones = [c for c in peaks if c["minimal_period"] == 1]
    for fpnt in us.henon_fixed_points():
        assert min(abs(np.mean(c["orbit_coordinates"]) - fpnt) for c in ones) < 0.05


def test_period2_orbit_points_follow_cyclic_structure():
    for c in us.henon_period2()["source_peak_candidates"]:
        assert np.allclose(c["orbit_points"], fp.cyclic_orbit_matrix(c["orbit_coordinates"], 2))
        assert c["source_stability"] is None  # no source procedure for p > 1


def test_logistic_r4_period3_orbits_recovered_in_dynamical_order():
    peaks = us.logistic4_period3()["source_peak_candidates"]
    true3 = [c for c in peaks if c["minimal_period"] == 3]
    for orbit in us.LOGISTIC4_PERIOD3:
        best = min(np.max(np.abs(c["canonical_orbit"] - np.asarray(orbit))) for c in true3)
        assert best < 0.02


def test_max_backbones_caps_work():
    X = us.henon_embedded(400)
    few = fp.detect_so_period_p(X, 2, 1, us.henon_config(so_random_R_period_p=5, so_max_backbones=20),
                                np.random.default_rng(0))
    many = fp.detect_so_period_p(X, 2, 1, us.henon_config(so_random_R_period_p=5, so_max_backbones=None),
                                 np.random.default_rng(0))
    assert many["histogram"].sum() > 10 * few["histogram"].sum()


def test_period_p_slab_keeps_requested_percentile():
    X = us.henon_embedded(300)
    cfg = us.henon_config(so_random_R_period_p=10, so_slab_percentile=25.0, so_hist_range_pad=1e6)
    r = fp.detect_so_period_p(X, 2, 1, cfg, np.random.default_rng(0))
    n_total = (len(X) - 2) * 4 * 10  # backbones x (K+1)^p combos x R
    assert abs(len(r["reduced"]) / n_total - 0.25) < 0.02
