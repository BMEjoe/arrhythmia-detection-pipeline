"""Failure semantics: NO PEAKS -> 0, PIPELINE FAILURE -> NaN."""
import sys
import pathlib
from dataclasses import replace

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import upo_systems as us  # noqa: E402
from upo_systems import fp  # noqa: E402

UPO_COLS = [c for mode in fp.UPO_FEATURE_CONTRACT.values() for c in mode if c != "lle_per_beat"]


def test_failure_and_no_peak_status_sets_are_disjoint_and_complete():
    assert not (fp.UPO_FAILURE_STATUSES & fp.UPO_NO_PEAK_STATUSES)
    assert fp.UPO_STATUS_OK not in fp.UPO_FAILURE_STATUSES | fp.UPO_NO_PEAK_STATUSES
    assert {"constant_data", "nonfinite_input", "too_few_points",
            "embedding_not_saturated"} <= fp.UPO_FAILURE_STATUSES
    assert fp.UPO_NO_PEAK_STATUSES == {"no_peaks", "tube_empty", "too_few_in_tube"}


@pytest.mark.parametrize("status", sorted(fp.UPO_FAILURE_STATUSES))
def test_pipeline_failure_gives_nan_features(status):
    feats = fp.upo_summary_features(us.fake_upo(status))
    assert all(np.isnan(feats[c]) for c in UPO_COLS)


@pytest.mark.parametrize("status", sorted(fp.UPO_NO_PEAK_STATUSES))
def test_no_peaks_gives_zero_features(status):
    feats = fp.upo_summary_features(us.fake_upo(status, significant=0, verified=0, rJ=0.9))
    for c in ["source_peak_count", "source_peaks_per_point", "source_peak_coverage",
              "significant_peak_count", "significant_peaks_per_point", "significant_peak_coverage",
              "verified_unstable_count", "verified_unstable_per_point", "verified_unstable_coverage"]:
        assert feats[c] == 0.0
    assert feats["source_rJ"] == 0.9  # the significance test still ran


def test_ok_status_counts_rates_and_coverage():
    feats = fp.upo_summary_features(us.fake_upo("ok", n=200, source=4, significant=2, verified=1, rJ=3.0))
    assert feats["source_peak_count"] == 4 and feats["source_peaks_per_point"] == 4 / 200
    assert feats["source_peak_coverage"] == 0.5
    assert feats["significant_peak_count"] == 2 and feats["significant_peak_coverage"] == 0.25
    assert feats["verified_unstable_count"] == 1 and feats["verified_unstable_coverage"] == 0.125
    assert feats["source_rJ"] == 3.0


def test_unassessed_significance_gives_nan_significant_features_not_zero():
    feats = fp.upo_summary_features(us.fake_upo("ok", source=3, significant=None))
    assert np.isnan(feats["significant_peak_count"]) and np.isnan(feats["source_rJ"])


def test_nonfinite_input_status():
    x = us.henon_series(300)
    x[10] = np.nan
    assert fp.run_upo_analysis(x, fp.CFG)["status"] == "nonfinite_input"
    X = us.henon_embedded(300).copy()
    X[5, 0] = np.inf
    assert fp.detect_so_fixed_points(X, 1, fp.CFG)["status"] == "nonfinite_input"


def test_constant_data_status():
    assert fp.run_upo_analysis(np.full(300, 0.8), fp.CFG)["status"] == "constant_data"
    assert fp.detect_so_fixed_points(np.full((100, 2), 0.8), 1, fp.CFG)["status"] == "constant_data"


def test_too_few_points_status():
    upo = fp.run_upo_analysis(us.henon_series(300), replace(fp.CFG, so_min_embedded_points=10_000))
    assert upo["status"] == "too_few_points"
    assert fp.detect_so_period_p(us.henon_embedded()[:12], 2, 1, fp.CFG)["status"] == "too_few_points"


def test_very_short_series_fails_as_embedding_not_saturated():
    """Cao cannot select a dimension from 25 points, so no UPO analysis is run."""
    upo = fp.run_upo_analysis(us.henon_series(25), replace(fp.CFG, cao_max_dim=3))
    assert upo["status"] == "embedding_not_saturated"


def test_embedding_not_saturated_status(monkeypatch):
    monkeypatch.setattr(fp, "cao_method", lambda x, tau, max_dim, tol, theiler: (max_dim + 1, None, None))
    upo = fp.run_upo_analysis(us.henon_series(300), fp.CFG)
    assert upo["status"] == "embedding_not_saturated"
    assert upo["source_peak_candidates"] == [] and upo["significance_assessed"] is False
    assert all(np.isnan(fp.upo_summary_features(upo)[c]) for c in UPO_COLS)


def test_no_valid_transforms_status(monkeypatch):
    monkeypatch.setattr(fp, "_batched_fixed_point_transform",
                        lambda z, Fz, J, Rs, k, r: np.full((len(Rs), len(z)), np.nan))
    r = fp.detect_so_fixed_points(us.henon_embedded(300), 1, us.henon_config(so_random_R=3))
    assert r["status"] == "no_valid_transforms"


def test_significance_not_assessed_for_failure_status():
    r = fp.detect_so_fixed_points(np.full((100, 2), 0.8), 1, fp.CFG)
    fp.assess_so_significance(r, np.full(101, 0.8), fp.CFG)
    assert r["significance_assessed"] is False and r["significant_uop_candidates"] is None
    assert "not_assessed_reason" in r["significance"]


def test_linalg_failure_inside_upo_becomes_analysis_error_and_lle_is_kept(monkeypatch):
    def boom(*a, **k):
        raise np.linalg.LinAlgError("forced")
    monkeypatch.setattr(fp, "run_upo_analysis", boom)
    res = fp.analyze_segment(us.chaotic_rr(), fp.CFG)
    assert res["upo_status"] == "analysis_error"
    assert np.isfinite(res["lle_per_beat"])
    feats = fp.extract_classifier_features(res, "source")
    assert np.isfinite(feats["lle_per_beat"]) and np.isnan(feats["source_peak_count"])


def test_short_rr_still_raises_in_analyze_segment():
    with pytest.raises(ValueError):
        fp.analyze_segment(us.chaotic_rr(80), fp.CFG)


# ---------------------------------------------------------------------------
# Significance requested but not assessable (e.g. an insufficient surrogate
# count) must never be read as "zero significant peaks".
# ---------------------------------------------------------------------------
SIG_COLS = ["significant_peak_count", "significant_peaks_per_point",
            "significant_peak_coverage", "source_rJ"]


def _sig_cfg(surrogates):
    return replace(fp.CFG, upo_feature_mode="significant_source", so_assess_significance=True,
                   so_surrogate_count=surrogates)


def test_insufficient_surrogates_gives_explicit_not_assessed_status():
    res = fp.analyze_segment(us.chaotic_rr(), _sig_cfg(1))
    upo = res["upo"]
    assert upo["status"] == "ok"  # the detector itself did not fail
    assert upo["significance_requested"] is True
    assert upo["significance_assessed"] is False
    assert upo["significance_status"] == fp.SIGNIFICANCE_NOT_ASSESSED == "significance_not_assessed"
    assert res["upo_significance_status"] == "significance_not_assessed"
    assert upo["significant_uop_candidates"] is None
    assert "two surrogate" in upo["significance_not_assessed_reasons"][1]


def test_insufficient_surrogates_gives_nan_significance_features_not_zero():
    feats = fp.extract_classifier_features(fp.analyze_segment(us.chaotic_rr(), _sig_cfg(1)),
                                           "significant_source")
    assert all(np.isnan(feats[c]) for c in SIG_COLS)
    assert np.isfinite(feats["lle_per_beat"])


def test_source_peak_features_remain_available_when_significance_not_assessed():
    res = fp.analyze_segment(us.chaotic_rr(), _sig_cfg(1))
    src = fp.extract_classifier_features(res, "source")
    assert src["source_peak_count"] == len(res["upo"]["source_peak_candidates"]) >= 1
    assert np.isfinite(src["source_peaks_per_point"]) and np.isfinite(src["source_peak_coverage"])


def test_assessed_with_no_significant_peaks_stays_empty_list_and_zero():
    noise = 0.8 + 0.05 * np.random.default_rng(5).standard_normal(256)
    res = fp.analyze_segment(noise, _sig_cfg(6))
    upo = res["upo"]
    assert upo["significance_status"] == fp.SIGNIFICANCE_ASSESSED
    assert upo["significant_uop_candidates"] == []
    assert upo["source_peak_candidates"]  # source peaks exist, none significant
    feats = fp.extract_classifier_features(res, "significant_source")
    assert feats["significant_peak_count"] == 0.0
    assert feats["significant_peaks_per_point"] == 0.0
    assert feats["significant_peak_coverage"] == 0.0
    assert np.isfinite(feats["source_rJ"])


def test_significance_status_distinguishes_detector_failure_and_not_requested():
    assert fp.run_upo_analysis(np.full(300, 0.8), _sig_cfg(6))["significance_status"] == \
        fp.SIGNIFICANCE_DETECTOR_FAILURE
    assert fp.run_upo_analysis(np.full(300, 0.8), fp.CFG)["significance_status"] == \
        fp.SIGNIFICANCE_NOT_REQUESTED
    assert fp.run_upo_analysis(us.henon_series(300), fp.CFG)["significance_status"] == \
        fp.SIGNIFICANCE_NOT_REQUESTED
    assert len(set(fp.SIGNIFICANCE_STATUSES)) == 4


def test_significant_source_without_requesting_significance_is_still_a_config_error():
    res = fp.analyze_segment(us.chaotic_rr(), fp.CFG)
    with pytest.raises(ValueError):
        fp.extract_classifier_features(res, "significant_source")
