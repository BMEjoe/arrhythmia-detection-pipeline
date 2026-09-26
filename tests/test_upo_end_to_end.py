"""analyze_segment / run_upo_analysis integration on synthetic RR-like series."""
import sys
import pathlib
from dataclasses import replace

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import upo_systems as us  # noqa: E402
from upo_systems import fp  # noqa: E402

META_KEYS = ["upo_map_mode", "upo_source_faithful_map", "upo_tau", "upo_embedding_dimension",
             "upo_n_embedded_points", "upo_status", "upo_significance_assessed"]


@pytest.fixture(scope="module")
def results():
    rr = us.chaotic_rr()
    out = {}
    for mode in ["one_sample", "tau_step"]:
        for sig in [False, True]:
            cfg = replace(fp.CFG, upo_map_mode=mode, so_assess_significance=sig, so_surrogate_count=6)
            out[(mode, sig)] = fp.analyze_segment(rr, cfg)
    return out


@pytest.mark.parametrize("mode", ["one_sample", "tau_step"])
def test_analyze_segment_records_upo_metadata(results, mode):
    res = results[(mode, False)]
    for k in META_KEYS:
        assert k in res
    assert res["upo_map_mode"] == mode
    assert res["upo_source_faithful_map"] is (mode == "one_sample")
    assert res["upo_n_embedded_points"] == len(res["upo"]["embedded"])


def test_one_sample_uses_lag1_embedding_with_cao_at_lag1(results):
    res = results[("one_sample", False)]
    m1, _, _ = fp.cao_method(res["rr_dynamics"], 1, fp.CFG.cao_max_dim, fp.CFG.cao_tol, fp.CFG.cao_theiler)
    assert res["upo_tau"] == 1 and res["upo_embedding_dimension"] == m1
    assert np.array_equal(res["upo"]["embedded"], fp._embed_backward(res["rr_dynamics"], 1, m1))
    assert res["tau"] > 1  # the LLE keeps its TDMI tau


def test_tau_step_reuses_tdmi_tau_and_lle_embedding(results):
    res = results[("tau_step", False)]
    assert res["upo_tau"] == res["tau"]
    assert res["upo_embedding_dimension"] == res["embedding_dimension"]
    assert np.array_equal(res["upo"]["embedded"], res["embedded"])
    assert "PROJECT EXTENSION" in res["upo"]["map_label"]


def test_lle_is_independent_of_upo_settings(results):
    lles = {k: r["lle_per_beat"] for k, r in results.items()}
    assert len(set(lles.values())) == 1
    embs = [r["embedded"] for r in results.values()]
    assert all(np.array_equal(embs[0], e) for e in embs[1:])


def test_upo_significance_flag_is_independent_of_lle_surrogates(results):
    for (mode, sig), res in results.items():
        assert res["upo_significance_assessed"] is sig
        assert res["surrogate_analysis"] is None  # compute_surrogates (LLE) stays off


def test_lle_surrogate_analysis_has_no_upo_component():
    cfg = replace(fp.CFG, compute_surrogates=True, surrogate_count=3)
    res = fp.analyze_segment(us.chaotic_rr(), cfg)
    assert "upo" not in res["surrogate_analysis"]
    assert res["upo_significance_assessed"] is False


def test_significance_assessed_gives_list_else_none(results):
    assert results[("one_sample", False)]["upo"]["significant_uop_candidates"] is None
    assert isinstance(results[("one_sample", True)]["upo"]["significant_uop_candidates"], list)
    assert np.isnan(results[("one_sample", False)]["upo"]["source_rJ"])
    assert np.isfinite(results[("one_sample", True)]["upo"]["source_rJ"])


def test_aggregate_levels_concatenate_period_results():
    cfg = replace(us.henon_config(so_random_R=20, so_random_R_period_p=5, so_max_backbones=60),
                  so_periods=(1, 2))
    upo = fp.run_upo_analysis(us.henon_series(400), cfg)
    assert set(upo["periods"]) == {1, 2}
    n = sum(len(r["source_peak_candidates"]) for r in upo["periods"].values())
    assert len(upo["source_peak_candidates"]) == n
    n_c = sum(len(r["verification_gated_candidates"]) + len(r["verification_failed_peaks"])
              for r in upo["periods"].values())
    assert len(upo["verification_gated_candidates"]) + len(upo["verification_failed_peaks"]) == n_c


def test_periodic_orbits_reported_when_enabled():
    cfg = replace(fp.CFG, so_periods=(1, 2), so_random_R_period_p=3, so_max_backbones=30)
    res = fp.analyze_segment(us.chaotic_rr(), cfg)
    assert 2 in res["so_periodic_orbits"] and res["so_periodic_orbits"][2]["period"] == 2
    assert res["so_fixed_points"]["period"] == 1


def test_analyze_segment_is_deterministic():
    a = fp.analyze_segment(us.chaotic_rr(), fp.CFG)
    b = fp.analyze_segment(us.chaotic_rr(), fp.CFG)
    fa, fb = fp.extract_classifier_features(a, "extended"), fp.extract_classifier_features(b, "extended")
    assert fa == fb


def _keys(obj):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield k
            yield from _keys(v)
    elif isinstance(obj, (list, tuple)):
        for v in obj:
            yield from _keys(v)


def test_no_generic_candidates_field_in_segment_result(results):
    for res in results.values():
        assert "candidates" not in set(_keys(res))


def test_periodic_series_has_no_verified_unstable_peaks():
    res = fp.analyze_segment(us.periodic_rr(), fp.CFG)
    feats = fp.extract_classifier_features(res, "extended")
    assert feats["verified_unstable_count"] == 0.0
    assert res["upo_status"] in {"ok"} | fp.UPO_NO_PEAK_STATUSES


def test_run_upo_analysis_rejects_invalid_config():
    with pytest.raises(ValueError):
        fp.run_upo_analysis(us.chaotic_rr(), replace(fp.CFG, upo_map_mode="lag"))


def test_twin_summary_uses_source_peak_coverage():
    s = fp.summarize_twin_results([fp.analyze_segment(us.chaotic_rr(), fp.CFG)],
                                  [fp.analyze_segment(us.periodic_rr(), fp.CFG)])
    assert "source_peak_coverage" in s and "upo_coverage" not in s
