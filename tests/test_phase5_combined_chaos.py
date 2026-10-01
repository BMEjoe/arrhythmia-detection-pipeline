"""Phase 5 opt-in helpers: combined_chaos_config() is the configuration tested in
experiments/phase5_rr/PREREGISTRATION.md and combined_chaos_detected() is its AND
decision on an analyze_segment result."""
import json
import pathlib
import sys
from dataclasses import replace

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import final_pipeline as fp  # noqa: E402
from experiments.phase5_rr import systems as S  # noqa: E402

STORED = ROOT / "experiments/phase5_rr/results/test/test.jsonl"


def _stored(conditions, seeds):
    out = {}
    for line in open(STORED):
        r = json.loads(line)
        t = r["task"]
        if t["condition"] in conditions and t["seed"] in seeds:
            out[(t["condition"], t["seed"])] = r["m2"]
    return out


def test_config_is_the_tested_one_and_defaults_unchanged():
    c = fp.combined_chaos_config()
    assert c == fp.phase4_upo_config(replace(fp.CFG, lle_chaos_test=True))
    assert c.lle_chaos_test and c.upo_instability_gate and c.upo_fixed_dimension == 2
    assert not c.use_corrected_rr_for_dynamics
    assert fp.CFG.lle_chaos_test is False and fp.CFG.upo_instability_gate is False
    assert fp.CFG.upo_fixed_dimension is None
    base = replace(fp.CFG, random_seed=7)
    assert fp.combined_chaos_config(base).random_seed == 7


def test_missing_components_raise():
    rr = S.generate("N1_linear_rr", 0)
    with pytest.raises(ValueError):
        fp.combined_chaos_detected(fp.analyze_segment(rr, fp.CFG))
    with pytest.raises(ValueError):
        fp.combined_chaos_detected({"lle_chaos_test": {"detected": True}, "upo": {}})
    with pytest.raises(ValueError):
        fp.combined_chaos_detected({"lle_chaos_test": {"detected": True}, "upo": {"status": "ok"}})


def test_upo_failure_status_counts_as_not_detected():
    """Phase 6 fix: a constant window ends the UPO analysis with 'constant_data' (no gate fields)."""
    import numpy as np
    out = fp.analyze_segment(np.full(256, 0.8), fp.combined_chaos_config())
    assert out["upo"]["status"] in fp.UPO_FAILURE_STATUSES
    assert "instability_gate_detected" not in out["upo"]
    assert fp.combined_chaos_detected(out) is False
    assert fp.combined_chaos_detected({"lle_chaos_test": {"detected": True},
                                       "upo": {"status": "no_valid_transforms"}}) is False


@pytest.mark.parametrize("condition,seed", [("N1_linear_rr", 3000), ("S2_ectopic_10pct", 3000),
                                            ("P1_henon_rr", 3000), ("P3_henon_rr_trend", 3000),
                                            ("G3_mackey_glass", 3000)])
def test_reproduces_stored_phase5_decisions(condition, seed):
    ref = _stored({condition}, {seed})[(condition, seed)]
    out = fp.analyze_segment(S.generate(condition, seed), fp.combined_chaos_config())
    assert out["lle_chaos_test"]["detected"] == ref["lle"]
    assert out["upo"]["instability_gate_detected"] == ref["upo"]
    assert fp.combined_chaos_detected(out) == (ref["lle"] and ref["upo"])
