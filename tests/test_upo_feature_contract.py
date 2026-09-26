"""Feature contract, explicit classifier plumbing, and leakage protection."""
import inspect
import sys
import pathlib
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import upo_systems as us  # noqa: E402
from upo_systems import fp  # noqa: E402

EXPECTED = {
    "source": ["lle_per_beat", "source_peak_count", "source_peaks_per_point",
               "source_peak_coverage"],
    "significant_source": ["lle_per_beat", "significant_peak_count",
                           "significant_peaks_per_point", "significant_peak_coverage",
                           "source_rJ"],
    "extended": ["lle_per_beat", "verified_unstable_count", "verified_unstable_per_point",
                 "verified_unstable_coverage"],
}


@pytest.mark.parametrize("mode", list(EXPECTED))
def test_contract_columns_are_exact(mode):
    assert list(fp.UPO_FEATURE_CONTRACT[mode]) == EXPECTED[mode]


@pytest.mark.parametrize("mode", list(EXPECTED))
def test_classifier_feature_columns_follow_config(mode):
    cfg = replace(fp.CFG, upo_feature_mode=mode, so_assess_significance=(mode == "significant_source"))
    assert fp.classifier_feature_columns(cfg) == EXPECTED[mode]


@pytest.mark.parametrize("old", ["n_unstable_upos", "upo_coverage", "upo_candidate_density"])
def test_old_ambiguous_feature_names_are_gone(old):
    for cols in fp.UPO_FEATURE_CONTRACT.values():
        assert old not in cols
    assert not hasattr(fp, "CLASSIFIER_FEATURE_COLUMNS")


def test_every_contract_feature_has_a_definition():
    for cols in fp.UPO_FEATURE_CONTRACT.values():
        for c in cols:
            assert c in fp.UPO_FEATURE_DEFINITIONS
    assert "histogram peaks" in fp.UPO_FEATURE_DEFINITIONS["source_peak_count"]
    assert "not a density" in fp.UPO_FEATURE_DEFINITIONS["source_peaks_per_point"]


def test_significant_source_requires_significance_assessment():
    with pytest.raises(ValueError):
        fp.classifier_feature_columns(replace(fp.CFG, upo_feature_mode="significant_source"))


@pytest.mark.parametrize("field,value", [("upo_feature_mode", "all"),
                                         ("upo_map_mode", "lag"),
                                         ("so_randomization", "abs")])
def test_invalid_modes_are_rejected(field, value):
    with pytest.raises(ValueError):
        fp.validate_upo_config(replace(fp.CFG, **{field: value}))


def test_default_feature_and_map_modes():
    assert fp.CFG.upo_feature_mode == "source"
    assert fp.CFG.upo_map_mode == "one_sample"
    assert fp.CFG.so_assess_significance is False


@pytest.fixture(scope="module")
def chaotic_result():
    cfg = replace(fp.CFG, so_assess_significance=True, so_surrogate_count=6)
    return fp.analyze_segment(us.chaotic_rr(), cfg)


@pytest.mark.parametrize("mode", list(EXPECTED))
def test_extract_returns_exactly_contract_columns(chaotic_result, mode):
    assert list(fp.extract_classifier_features(chaotic_result, mode)) == EXPECTED[mode]


def test_extract_rejects_unknown_mode(chaotic_result):
    with pytest.raises(ValueError):
        fp.extract_classifier_features(chaotic_result, "legacy")


def test_significant_source_on_unassessed_result_raises():
    res = fp.analyze_segment(us.chaotic_rr(), fp.CFG)
    with pytest.raises(ValueError):
        fp.extract_classifier_features(res, "significant_source")


def test_row_metadata_fields(chaotic_result):
    meta = fp.upo_row_metadata(chaotic_result, "extended")
    assert set(meta) == set(fp.UPO_ROW_METADATA_COLUMNS)
    assert meta == {"upo_feature_mode": "extended", "upo_map_mode": "one_sample",
                    "upo_status": chaotic_result["upo_status"],
                    "upo_significance_assessed": True,
                    "upo_significance_status": "assessed"}


def test_grouped_logistic_regression_requires_explicit_feature_columns():
    sig = inspect.signature(fp.run_grouped_logistic_regression)
    assert sig.parameters["feature_columns"].default is inspect.Parameter.empty
    with pytest.raises(TypeError):
        fp.run_grouped_logistic_regression(pd.DataFrame())


@pytest.mark.parametrize("bad", [[], ["lle_per_beat", "label"], ["record"], ["abnormal_fraction"]])
def test_grouped_logistic_regression_rejects_empty_or_leaking_columns(bad):
    with pytest.raises(ValueError):
        fp.run_grouped_logistic_regression(pd.DataFrame({"label": [0, 1], "record": ["a", "b"]}), bad)


def test_feature_extraction_signature_takes_no_label_information():
    params = list(inspect.signature(fp.extract_classifier_features).parameters)
    assert params == ["result", "feature_mode"]
    assert list(inspect.signature(fp.analyze_segment).parameters) == ["rr_intervals", "config"]


WIN = 200


def _synthetic_beats(fs):
    rr = us.chaotic_rr(2 * WIN + 20, seed=1)
    return np.concatenate([[1000], 1000 + np.cumsum(np.round(rr * fs))]).astype(int)


def _patched_dataset(monkeypatch, symbols, record="r1", config=None):
    """Synthetic record through the MIT-BIH plumbing (no download, no MIT-BIH data)."""
    fs = 360.0
    beats = _synthetic_beats(fs)
    monkeypatch.setattr(fp, "load_mitbih_record_with_annotations",
                        lambda name, dl_dir=None, db_name=None: (np.zeros(beats[-1] + 1000), fs,
                                                                 beats, np.array(symbols(len(beats)))))
    monkeypatch.setattr(fp, "detect_r_peaks", lambda sig, fs_, cfg: beats)
    cfg = config or fp.CFG
    return fp.build_mitbih_classification_dataset(records=[record], config=cfg, window_size=WIN,
                                                  step_size=WIN)


def test_dataset_rows_record_feature_mode_and_upo_metadata(monkeypatch):
    df, meta = _patched_dataset(monkeypatch, lambda n: ["N"] * n)
    assert len(df) == 2
    assert list(df["analysis_error"]) == ["", ""]
    for col in EXPECTED["source"] + list(fp.UPO_ROW_METADATA_COLUMNS):
        assert col in df.columns
    assert set(df["upo_feature_mode"]) == {"source"}
    assert set(df["upo_map_mode"]) == {"one_sample"}
    assert not df["upo_significance_assessed"].any()


def test_features_do_not_depend_on_labels_or_record_identity(monkeypatch):
    df_n, _ = _patched_dataset(monkeypatch, lambda n: ["N"] * n, record="r1")
    df_v, _ = _patched_dataset(monkeypatch, lambda n: ["V"] * n, record="zz9")
    assert list(df_n["label"]) == [0, 0] and list(df_v["label"]) == [1, 1]
    for col in EXPECTED["source"]:
        assert np.array_equal(df_n[col].to_numpy(), df_v[col].to_numpy(), equal_nan=True)


def test_window_features_do_not_depend_on_other_windows(monkeypatch):
    df, _ = _patched_dataset(monkeypatch, lambda n: ["N"] * n)
    rr_det, _, _ = fp.extract_rr_intervals(_synthetic_beats(360.0), 360.0)  # same arithmetic as the pipeline
    alone = fp.extract_classifier_features(fp.analyze_segment(rr_det[WIN:2 * WIN], fp.CFG), "source")
    assert df.loc[1, "analysis_error"] == ""
    for col in EXPECTED["source"]:
        assert np.isclose(df.loc[1, col], alone[col], equal_nan=True)


def test_grouped_cv_runs_on_synthetic_frame_with_explicit_columns():
    rng = np.random.default_rng(0)
    n = 40
    df = pd.DataFrame({"record": np.repeat([f"r{i}" for i in range(8)], 5),
                       "label": np.tile([0, 1, 0, 1, 1], 8),
                       "f1": rng.standard_normal(n), "f2": rng.standard_normal(n)})
    out = fp.run_grouped_logistic_regression(df, ["f1", "f2"], n_splits=4)
    for _, fold in out["fold_results"].iterrows():
        assert fold["n_train_records"] + fold["n_test_records"] == 8


def test_insufficient_surrogates_do_not_discard_windows_or_record(monkeypatch):
    """Regression: significance requested but not assessable (1 surrogate) used
    to raise outside the per-window error handling and drop the whole record."""
    cfg = replace(fp.CFG, upo_feature_mode="significant_source", so_assess_significance=True,
                  so_surrogate_count=1)
    df, meta = _patched_dataset(monkeypatch, lambda n: ["N"] * n, config=cfg)
    assert list(meta["status"]) == ["ok"]
    assert len(df) == 2
    assert list(df["analysis_error"]) == ["", ""]
    assert set(df["upo_significance_status"]) == {"significance_not_assessed"}
    assert not df["upo_significance_assessed"].any()
    for col in ["significant_peak_count", "significant_peaks_per_point",
                "significant_peak_coverage", "source_rJ"]:
        assert df[col].isna().all()
    assert df["lle_per_beat"].notna().all()
