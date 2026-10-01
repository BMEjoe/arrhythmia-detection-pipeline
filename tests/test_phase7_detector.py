"""Phase 7 detector wrapper (experiments/phase7_mitbih/detector.py) on synthetic windows."""
from dataclasses import replace

import numpy as np
import pytest

import final_pipeline as fp
from experiments.phase5_rr import systems as S5
from experiments.phase7_mitbih import detector as DT


def test_detector_is_the_phase6_recommended_configuration():
    assert DT.DETECTOR == fp.combined_chaos_config(replace(fp.CFG, keep_upo_on_short_lle_embedding=True),
                                                   detrend=True)
    assert DT.DETECTOR.rr_detrend == "linear" and DT.DETECTOR.rr_detrend_min_trend_sd == 0.7
    assert DT.DETECTOR.keep_upo_on_short_lle_embedding and DT.DETECTOR.lle_chaos_test
    assert DT.DETECTOR.upo_fixed_dimension == 2 and DT.DETECTOR.upo_instability_gate


@pytest.mark.parametrize("cond,seed", [("N1_linear_rr", 0), ("S2_ectopic_10pct", 1), ("P1_henon_rr", 0)])
def test_record_matches_direct_pipeline_call(cond, seed):
    rr = S5.generate(cond, seed)
    rec = DT.evaluate(rr)
    out = fp.analyze_segment(rr, DT.DETECTOR)
    assert rec["error"] is None
    assert rec["and_detected"] == fp.combined_chaos_detected(out)
    assert rec["lle_p"] == out["lle_chaos_test"]["p"]
    assert rec["upo_detected"] == out["upo"]["instability_gate_detected"]
    assert rec["upo_score_consistent"]
    sur = out["lle_chaos_test"]["surrogate_lle"]
    assert rec["lle_z"] == pytest.approx((out["lle_chaos_test"]["lle"] - np.mean(sur)) / np.std(sur, ddof=1))
    if rec["upo_detected"]:
        # a significant gate-passing peak exists, so the score exceeds 1 (J < 0.05 => dev > median W)
        assert rec["upo_score"] > 1.0


def test_error_counts_as_not_detected():
    rec = DT.evaluate(np.ones(50))                     # < 100 intervals: analyze_segment raises
    assert rec["error"] and rec["and_detected"] is False and rec["upo_score"] is None
