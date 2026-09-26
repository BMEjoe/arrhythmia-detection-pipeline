"""So fixed-point transform G(z, R) (PRE Eqs. 3-8; PRL) and period-1 detection."""
import sys
import pathlib
from dataclasses import replace

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import upo_systems as us  # noqa: E402
from upo_systems import fp  # noqa: E402

FORMS = ["prl_norm", "pre_tensor"]


def _R(form, d, seed=0):
    rng = np.random.default_rng(seed)
    return fp._random_R_matrix(rng, d, 1.0) if form == "prl_norm" else fp._random_R_tensor(rng, d, 1.0)


@pytest.mark.parametrize("form", FORMS)
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_fixed_point_is_invariant_for_any_R(form, seed):
    """Property (i): F(z*) = z*  =>  G(z*, R) = z* independent of R."""
    zs = us.henon_fixed_points()[0]
    z = np.array([zs, zs])
    J = us.henon_jacobian(zs)
    out = fp.so_fixed_point_transform(z, z.copy(), J, _R(form, 2, seed), 3.0, form)
    assert np.allclose(out, z, atol=1e-12)


def test_linear_map_with_R_zero_maps_every_point_to_fixed_point():
    """PRE Sec. II A: F(z) = z* + a (z - z*), R = 0 => G(z, 0) = z*."""
    a = np.array([[1.8, -0.4], [1.0, 0.0]])
    zstar = np.array([0.7, 0.7])
    for z in np.random.default_rng(0).standard_normal((20, 2)):
        Fz = zstar + a @ (z - zstar)
        out = fp.so_fixed_point_transform(z, Fz, a, np.zeros((2, 2)), 0.0)
        assert np.allclose(out, zstar)


def test_R_zero_delay_map_transforms_onto_diagonal():
    """PRE ref. [16]: with R = 0 the delay structure forces G_i = G_{i-1}."""
    X = us.henon_embedded(400)
    Js = fp._local_jacobians(X, 1, 2, 1)
    for i in range(10, 60):
        out = fp.so_fixed_point_transform(X[i], X[i + 1], Js[i], np.zeros((2, 2)), 0.0)
        assert abs(out[0] - out[1]) < 1e-8


@pytest.mark.parametrize("form", FORMS)
def test_fixed_point_is_stationary_point_of_G(form):
    """Property (ii): grad G(z*, R) = 0 (PRE Eq. 8), using the exact map."""
    zs = us.henon_fixed_points()[0]
    R = _R(form, 2, 5)

    def G(z):
        Fz = np.array([us.henon_f(z[0], z[1]), z[0]])
        return fp.so_fixed_point_transform(z, Fz, us.henon_jacobian(z[0]), R, 3.0, form)

    z0 = np.array([zs, zs])
    h = 1e-5
    for k in range(2):
        e = np.zeros(2)
        e[k] = h
        assert np.max(np.abs((G(z0 + e) - G(z0 - e)) / (2 * h))) < 1e-4


def test_random_R_matrix_is_uniform_on_minus_one_one():
    R = fp._random_R_matrix(np.random.default_rng(0), 4, 3.0)
    assert R.shape == (4, 4)
    big = np.stack([fp._random_R_matrix(np.random.default_rng(s), 3, 1.0) for s in range(2000)])
    assert big.min() >= -1 and big.max() <= 1
    assert abs(big.mean()) < 0.02
    assert abs(big.var() - 1.0 / 3.0) < 0.01


def test_random_R_tensor_shape_and_bounds():
    T = fp._random_R_tensor(np.random.default_rng(0), 3, 1.0)
    assert T.shape == (3, 3, 3)
    assert T.min() >= -1 and T.max() <= 1


def test_prl_perturbation_uses_l1_norm():
    R = np.array([[0.5, -0.25], [1.0, 0.0]])
    delta = np.array([3.0, -4.0])  # L1 = 7, L2 = 5
    assert np.allclose(fp.so_perturbation(R, delta, 2.0, "prl_norm"), 2.0 * 7.0 * R)


def test_pre_tensor_perturbation_contracts_last_index():
    T = np.random.default_rng(1).uniform(-1, 1, (2, 2, 2))
    delta = np.array([0.3, -0.7])
    expected = np.array([[sum(T[i, j, k] * delta[k] for k in range(2)) for j in range(2)]
                         for i in range(2)])
    assert np.allclose(fp.so_perturbation(T, delta, 1.5, "pre_tensor"), 1.5 * expected)


def test_unknown_randomization_raises():
    with pytest.raises(ValueError):
        fp.so_perturbation(np.eye(2), np.ones(2), 1.0, "abs")
    with pytest.raises(ValueError):
        fp._draw_R(np.random.default_rng(0), 3, 2, "gaussian")


@pytest.mark.parametrize("form", FORMS)
def test_batched_transform_matches_scalar_transform(form):
    rng = np.random.default_rng(3)
    z, Fz = rng.standard_normal(3), rng.standard_normal(3)
    J = fp._companion(rng.standard_normal(3))
    Rs = fp._draw_R(np.random.default_rng(4), 6, 3, form)
    batch = fp._batched_fixed_point_transform(z, Fz, J, Rs, 2.0, form)
    for r in range(6):
        assert np.allclose(batch[r], fp.so_fixed_point_transform(z, Fz, J, Rs[r], 2.0, form))


def test_singular_system_returns_none_and_batch_gives_nan():
    z, Fz = np.zeros(2), np.zeros(2)
    J = np.eye(2)  # I - S = 0 with kappa * ||delta|| = 0
    assert fp.so_fixed_point_transform(z, Fz, J, np.zeros((2, 2)), 0.0) is None
    out = fp._batched_fixed_point_transform(z, Fz, J, np.zeros((3, 2, 2)), 0.0, "prl_norm")
    assert np.all(np.isnan(out))


def test_kappa_zero_removes_randomization():
    rng = np.random.default_rng(6)
    z, Fz, J = rng.standard_normal(2), rng.standard_normal(2), fp._companion(rng.standard_normal(2))
    a = fp.so_fixed_point_transform(z, Fz, J, np.ones((2, 2)), 0.0)
    b = fp.so_fixed_point_transform(z, Fz, J, np.zeros((2, 2)), 0.0)
    assert np.allclose(a, b)


def test_one_dimensional_tensor_form_is_pre_eqs_21_22():
    """d = 1: S = slope + R [z(n+1) - z(n)] with signed difference (PRE Eq. 22)."""
    zn, zn1, slope, Rv, k = 0.31, 0.74, -1.3, 0.4, 10.0
    S = slope + k * Rv * (zn1 - zn)
    expected = (zn1 - S * zn) / (1 - S)  # PRE Eq. (21)
    got = fp.so_fixed_point_transform_tensor(np.array([zn]), np.array([zn1]),
                                             np.array([[slope]]), np.array([[[Rv]]]), k)
    assert np.isclose(got[0], expected)


@pytest.mark.parametrize("form", FORMS)
def test_skewed_logistic_fixed_point_pre_section_IIF(form):
    """PRE Sec. II F / Fig. 4(b): 256 iterates, 32 random R, k = 10."""
    x = us.skewed_logistic_series(256)
    cfg = replace(fp.CFG, so_jacobian_neighbors=1, so_random_R=32, so_kappa=10.0,
                  so_randomization=form)
    r = fp.detect_so_fixed_points(x[:, None], 1, cfg, np.random.default_rng(0))
    assert r["status"] == "ok"
    assert abs(us.nearest(us.period1_locations(r), us.skewed_logistic_fixed_point())
               - us.skewed_logistic_fixed_point()) < 0.01


@pytest.mark.parametrize("form", FORMS)
def test_logistic_fixed_point_prl_figure_2(form):
    """PRL Fig. 2(b): r = 3.92, 100 data points, 500 random k, kappa = 5."""
    x = us.logistic_series(100)
    cfg = replace(fp.CFG, so_jacobian_neighbors=1, so_random_R=500, so_kappa=5.0,
                  so_randomization=form)
    r = fp.detect_so_fixed_points(x[:, None], 1, cfg, np.random.default_rng(0))
    locs = us.period1_locations(r)
    assert abs(us.nearest(locs, (3.92 - 1) / 3.92) - (3.92 - 1) / 3.92) < 0.03
    # The strongest peak is the fixed point on the attractor.
    assert abs(locs[0] - (3.92 - 1) / 3.92) < 0.03
