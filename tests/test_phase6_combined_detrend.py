"""Phase 6 C5: combined_chaos_config(detrend=True) is the preregistered winner D2
(experiments/phase6_robust/PREREGISTRATION.md); the default call is unchanged."""
import json
import pathlib
import sys
from dataclasses import replace

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import final_pipeline as fp  # noqa: E402
from experiments.phase6_robust import methods as MM  # noqa: E402
from experiments.phase6_robust import systems as S  # noqa: E402

STORED = ROOT / "experiments/phase6_robust/results/test/test.jsonl"
TESTED = replace(fp.CFG, keep_upo_on_short_lle_embedding=True)


def test_default_call_unchanged():
    c = fp.combined_chaos_config()
    assert c == fp.phase4_upo_config(replace(fp.CFG, lle_chaos_test=True))
    assert c.rr_detrend == "none" and c.rr_detrend_min_trend_sd == 0.0
    assert fp.combined_chaos_config(detrend=False) == c


def test_detrend_is_the_preregistered_winner():
    c = fp.combined_chaos_config(TESTED, detrend=True)
    assert c == MM.PREREGISTERED["d2_linear_g07"]["config"]
    assert (c.rr_detrend, c.rr_detrend_min_trend_sd) == ("linear", 0.7)
    assert fp.combined_chaos_config(TESTED) == MM.BASELINE_K


@pytest.mark.parametrize("condition,seed", [("P3_henon_rr_trend", 4000), ("P3_logistic_rr_trend", 4001),
                                            ("N4_linear_rr_step", 4000), ("G1_lorenz_maxima", 4000),
                                            ("E6_ectopic10_trend", 4000)])
def test_reproduces_stored_d2_decisions(condition, seed):
    ref = None
    for line in open(STORED):
        r = json.loads(line)
        if r["task"]["condition"] == condition and r["task"]["seed"] == seed:
            ref = r["methods"]["d2_linear_g07"]
            break
    out = fp.analyze_segment(S.generate(condition, seed), fp.combined_chaos_config(TESTED, detrend=True))
    assert out["lle_chaos_test"]["detected"] == ref["lle"]
    assert out["upo"].get("instability_gate_detected", False) == ref["upo"]
    assert fp.combined_chaos_detected(out) == ref["and"]
