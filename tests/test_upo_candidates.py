"""Level A source-peak API: explicit fields, no generic "candidates", statuses."""
import sys
import pathlib
from dataclasses import replace

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import upo_systems as us  # noqa: E402
from upo_systems import fp  # noqa: E402

REQUIRED_KEYS = [
    "detection_label", "period", "map_mode", "source_faithful_map", "tau", "step",
    "randomization", "status", "levels", "source_peak_candidates",
    "significance_assessed", "significant_uop_candidates",
    "verification_assessed", "verification_gated_candidates", "verification_failed_peaks",
]


def _all_keys(obj):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield k
            yield from _all_keys(v)
    elif isinstance(obj, (list, tuple)):
        for v in obj:
            yield from _all_keys(v)


@pytest.mark.parametrize("key", REQUIRED_KEYS)
def test_detection_results_have_explicit_field(key):
    assert key in us.henon_period1()
    assert key in us.henon_period2()


@pytest.mark.parametrize("getter", [us.henon_period1, us.henon_period2,
                                    us.henon_period2_verified, us.henon_significance])
def test_no_generic_candidates_field_anywhere(getter):
    assert "candidates" not in set(_all_keys(getter()))


def test_level_field_names_are_exact():
    assert fp.LEVEL_A_FIELD == "source_peak_candidates"
    assert fp.LEVEL_B_FIELD == "significant_uop_candidates"
    assert fp.LEVEL_C_PASS_FIELD == "verification_gated_candidates"
    assert fp.LEVEL_C_FAIL_FIELD == "verification_failed_peaks"
    assert fp.UPO_LEVELS["A"]["provenance"] == "SOURCE"
    assert fp.UPO_LEVELS["C"]["provenance"] == "PROJECT EXTENSION"


def test_levels_B_and_C_are_none_before_assessment():
    r = us.henon_period1()
    assert r["significance_assessed"] is False
    assert r["significant_uop_candidates"] is None
    assert r["verification_assessed"] is False
    assert r["verification_gated_candidates"] is None
    assert r["verification_failed_peaks"] is None


def test_source_peaks_are_labelled_source_level_A():
    for c in us.henon_period1()["source_peak_candidates"]:
        assert c["provenance"] == "SOURCE"
        assert c["level"] == "A"
        assert c["period"] == 1


def test_period1_peak_location_lies_on_diagonal():
    for c in us.henon_period1()["source_peak_candidates"]:
        assert np.allclose(c["location"], c["scalar_location"])
        assert c["orbit_points"].shape == (1, 2)


def test_peaks_are_histogram_local_maxima_above_threshold():
    r = us.henon_period1()
    h = r["histogram"]
    for c in r["source_peak_candidates"]:
        i = c["peak_cell"][0]
        assert c["histogram_count"] == h[i] >= r["peak_threshold"]
        assert h[i] >= h[max(i - 1, 0)] and h[i] >= h[min(i + 1, len(h) - 1)]


def test_peaks_sorted_by_histogram_count():
    counts = [c["histogram_count"] for c in us.henon_period2()["source_peak_candidates"]]
    assert counts == sorted(counts, reverse=True)


def test_source_member_indices_are_distinct_data_points():
    for c in us.henon_period1()["source_peak_candidates"]:
        m = c["source_member_indices"]
        assert len(m) == len(np.unique(m)) > 0
        assert np.all(m < len(us.henon_embedded()) - 1)


def test_detection_label_names_method_form_and_map():
    lab = us.henon_period1()["detection_label"]
    assert "So et al." in lab and "PRL norm" in lab and "SOURCE" in lab
    lab2 = us.henon_period1("pre_tensor")["detection_label"]
    assert "PRE tensor" in lab2


def test_detection_is_deterministic_for_fixed_rng_seed():
    X = us.henon_embedded(400)
    cfg = us.henon_config(so_random_R=20)
    a = fp.detect_so_fixed_points(X, 1, cfg, np.random.default_rng(9))
    b = fp.detect_so_fixed_points(X, 1, cfg, np.random.default_rng(9))
    assert np.array_equal(a["histogram"], b["histogram"])
    assert us.period1_locations(a) == us.period1_locations(b)


def test_status_ok_when_peaks_found():
    assert us.henon_period1()["status"] == "ok"


def test_status_no_peaks_when_threshold_exceeds_histogram():
    cfg = us.henon_config(so_random_R=10, so_peak_sigma=1e9)
    r = fp.detect_so_fixed_points(us.henon_embedded(300), 1, cfg, np.random.default_rng(0))
    assert r["status"] == "no_peaks" and r["source_peak_candidates"] == []
    assert r["histogram"] is not None


def test_status_tube_empty_when_nothing_in_histogram_range():
    cfg = us.henon_config(so_random_R=10)
    r = fp.detect_so_fixed_points(us.henon_embedded(300), 1, cfg, np.random.default_rng(0),
                                  hist_range=(100.0, 101.0))
    assert r["status"] == "tube_empty" and r["source_peak_candidates"] == []


def test_status_too_few_in_tube():
    cfg = us.henon_config(so_random_R=10)
    X = us.henon_embedded(300)
    wide = fp.detect_so_fixed_points(X, 1, cfg, np.random.default_rng(0), hist_range=(-1e9, 1e9))
    s = np.sort(wide["scalar"])
    k = len(s) // 2
    narrow = (s[k], np.nextafter(s[k + 2], np.inf))  # exactly 3 tube points
    r = fp.detect_so_fixed_points(X, 1, cfg, np.random.default_rng(0), hist_range=narrow)
    assert len(r["scalar"]) == 3 < fp.CFG.so_peak_min_count
    assert r["status"] == "too_few_in_tube"


def test_status_too_few_points():
    r = fp.detect_so_fixed_points(us.henon_embedded()[:10], 1, us.henon_config(), None)
    assert r["status"] == "too_few_points" and r["source_peak_candidates"] == []


def test_histogram_range_is_padded_attractor_range():
    r = us.henon_period1()
    x = us.henon_embedded()[:, 0]
    pad = fp.CFG.so_hist_range_pad * (x.max() - x.min())
    assert np.isclose(r["edges"][0][0], x.min() - pad)
    assert np.isclose(r["edges"][0][-1], x.max() + pad)
    assert len(r["edges"][0]) == fp.CFG.so_hist_bins + 1


def test_tube_keeps_requested_percentile_closest_to_diagonal():
    Z = np.random.default_rng(0).standard_normal((1000, 3))
    scalar, keep = fp._project_fixed_point_tube(Z, 10.0)
    assert abs(keep.mean() - 0.10) < 0.005
    dist = np.linalg.norm(Z - Z.mean(1, keepdims=True), axis=1)
    assert dist[keep].max() <= dist[~keep].min()
    assert np.allclose(scalar, Z[keep].mean(1))


def test_histogram_peaks_one_per_plateau():
    h = np.array([0, 5, 9, 9, 3, 0, 7, 0], dtype=float)
    cells = fp._histogram_peaks(h, 5)
    assert sorted(c[0] for c in cells) == [2, 6]


def test_raw_histogram_resolves_adjacent_peaks():
    """Smoothing is used for the threshold only; close peaks stay separate."""
    h = np.zeros(30)
    h[10], h[12] = 50, 40
    cfg = replace(fp.CFG, so_peak_min_count=5)
    cells = fp._histogram_peaks(h, fp._histogram_peak_threshold(h, cfg))
    assert sorted(c[0] for c in cells) == [10, 12]
