"""Local-linear delay-map Jacobian (PRE Eq. 20) and neighbour selection."""
import sys
import pathlib

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import upo_systems as us  # noqa: E402
from upo_systems import fp  # noqa: E402


@pytest.mark.parametrize("d", [1, 2, 3, 5])
def test_companion_structure_rows_are_exact_delay_shift(d):
    grad = np.arange(1, d + 1, dtype=float)
    J = fp._companion(grad)
    assert np.array_equal(J[0], grad)
    if d > 1:
        assert np.array_equal(J[1:, :-1], np.eye(d - 1))
        assert np.all(J[1:, -1] == 0)


def test_response_is_first_coordinate_difference_of_images():
    rng = np.random.default_rng(0)
    X = rng.standard_normal((40, 3))
    idx, nbrs = 5, np.array([10, 17, 22, 30, 33])
    J = fp._delay_map_jacobian_from_local_fit(X, idx, nbrs, step=1)
    A = X[nbrs] - X[idx]
    b = X[nbrs + 1, 0] - X[idx + 1, 0]
    grad, *_ = np.linalg.lstsq(A, b, rcond=None)
    assert np.allclose(J[0], grad)


def test_step_greater_than_one_uses_images_step_ahead():
    rng = np.random.default_rng(1)
    X = rng.standard_normal((50, 2))
    idx, nbrs, s = 4, np.array([12, 20, 25, 31]), 3
    J = fp._delay_map_jacobian_from_local_fit(X, idx, nbrs, step=s)
    b = X[nbrs + s, 0] - X[idx + s, 0]
    grad, *_ = np.linalg.lstsq(X[nbrs] - X[idx], b, rcond=None)
    assert np.allclose(J[0], grad)


@pytest.mark.parametrize("M", [2, 4, 7])
def test_lozi_exact_gradient_recovered(M):
    x = us.lozi_series(4096)
    X = fp._embed_backward(x, 1, 2)
    checked = 0
    for idx in range(100, 400):
        neigh = fp._spatial_neighbors(X, idx, K=M, exclude=1, valid_limit=len(X) - 1)
        if np.all(np.sign(X[neigh, 0]) == np.sign(X[idx, 0])) and abs(X[idx, 0]) > 0.05:
            J = fp._delay_map_jacobian_from_local_fit(X, idx, neigh, step=1)
            assert np.allclose(J[0], [-us.LOZI_A * np.sign(X[idx, 0]), us.LOZI_B], atol=1e-8)
            checked += 1
    assert checked > 100


def test_henon_gradient_matches_analytic_derivative():
    X = us.henon_embedded(4096)
    errs = []
    for idx in range(50, 450):
        neigh = fp._spatial_neighbors(X, idx, K=7, exclude=1, valid_limit=len(X) - 1)
        J = fp._delay_map_jacobian_from_local_fit(X, idx, neigh, step=1)
        errs.append(np.abs(J[0] - us.henon_jacobian(X[idx, 0])[0]))
    assert np.median(np.asarray(errs)[:, 0]) < 0.05
    assert np.median(np.asarray(errs)[:, 1]) < 0.05


def test_M_equal_d_is_allowed_source_henon_setting():
    X = us.henon_embedded()
    neigh = fp._spatial_neighbors(X, 100, K=2, exclude=1, valid_limit=len(X) - 1)
    assert len(neigh) == 2
    assert fp._delay_map_jacobian_from_local_fit(X, 100, neigh, step=1) is not None


def test_fewer_than_d_neighbours_returns_none():
    X = np.random.default_rng(2).standard_normal((30, 3))
    assert fp._delay_map_jacobian_from_local_fit(X, 3, [7, 9], step=1) is None


def test_neighbours_without_image_are_dropped():
    X = np.random.default_rng(3).standard_normal((20, 2))
    assert fp._delay_map_jacobian_from_local_fit(X, 3, [19, 18], step=2) is None
    assert fp._delay_map_jacobian_from_local_fit(X, 19, [3, 5, 7], step=1) is None


def test_rank_deficient_fit_returns_none():
    X = np.zeros((20, 2))
    X[:, 0] = np.arange(20.0)
    X[:, 1] = np.arange(20.0)  # all displacements collinear
    assert fp._delay_map_jacobian_from_local_fit(X, 2, [5, 8, 11], step=1) is None


def test_invalid_step_raises():
    X = np.random.default_rng(4).standard_normal((20, 2))
    with pytest.raises(ValueError):
        fp._delay_map_jacobian_from_local_fit(X, 2, [5, 8, 11], step=0)


def test_spatial_neighbours_exclude_self_and_temporal_window():
    X = np.random.default_rng(5).standard_normal((60, 2))
    for exclude in [1, 3]:
        nb = fp._spatial_neighbors(X, 30, K=10, exclude=exclude)
        assert len(nb) == 10
        assert np.all(np.abs(nb - 30) > exclude)


def test_spatial_neighbours_are_sorted_nearest_first():
    X = np.random.default_rng(6).standard_normal((80, 3))
    nb = fp._spatial_neighbors(X, 10, K=8, exclude=1)
    d = np.linalg.norm(X[nb] - X[10], axis=1)
    assert np.all(np.diff(d) >= 0)


def test_spatial_neighbours_respect_valid_limit():
    X = np.random.default_rng(7).standard_normal((50, 2))
    nb = fp._spatial_neighbors(X, 5, K=20, exclude=1, valid_limit=30)
    assert np.all(nb < 30)


def test_local_jacobians_missing_only_for_points_without_image():
    X = us.henon_embedded(300)
    Js = fp._local_jacobians(X, step=2, M=4, exclude=2)
    assert Js[-1] is None and Js[-2] is None
    assert sum(J is not None for J in Js[:-2]) == len(X) - 2


@pytest.mark.parametrize("shift,scale", [(5.0, 1.0), (0.0, 3.0), (-2.0, 0.1)])
def test_gradient_invariant_to_affine_rescaling_of_series(shift, scale):
    x = us.henon_series(600)
    X1 = fp._embed_backward(x, 1, 2)
    X2 = fp._embed_backward(shift + scale * x, 1, 2)
    nb = fp._spatial_neighbors(X1, 50, K=5, exclude=1, valid_limit=len(X1) - 1)
    J1 = fp._delay_map_jacobian_from_local_fit(X1, 50, nb, step=1)
    J2 = fp._delay_map_jacobian_from_local_fit(X2, 50, nb, step=1)
    assert np.allclose(J1, J2, atol=1e-8)
