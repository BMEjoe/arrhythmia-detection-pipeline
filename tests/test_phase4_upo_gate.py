"""Phase 4 opt-in UPO options (default off): upo_fixed_dimension and the P1-a
instability gate with a robust hybrid-stability aggregate; phase4_upo_config()
is the preregistered winner C2 (experiments/phase4_upo/PREREGISTRATION.md)."""
import sys
import pathlib
from dataclasses import replace

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import final_pipeline as fp  # noqa: E402
from experiments.phase4_upo import detector as D  # noqa: E402
from experiments.phase4_upo import methods as MM  # noqa: E402
from experiments.phase4_upo import systems as S  # noqa: E402

NEW_KEYS = ("instability_gate", "instability_gated_uop_candidates", "instability_gate_detected")
C2 = fp.phase4_upo_config()


def test_defaults_are_off_and_default_output_has_no_new_fields():
    c = fp.CFG
    assert c.upo_fixed_dimension is None and c.upo_instability_gate is False
    assert (c.upo_instability_gate_delta, c.upo_instability_gate_aggregate) == (0.2, "median")
    out = fp.run_upo_analysis(S.generate("henon", 256, 0), config=c)
    assert not any(k in out for k in NEW_KEYS)


def test_phase4_config_is_the_preregistered_c2():
    assert (C2.upo_fixed_dimension, C2.so_jacobian_neighbors, C2.so_surrogate_count) == (2, 15, 50)
    assert C2.so_assess_significance and C2.upo_instability_gate
    assert (C2.upo_instability_gate_delta, C2.upo_instability_gate_aggregate) == (0.2, "median")
    assert MM.METHODS["c2_m2M15_mediangate"]["run"] == (2, 15, 50) and MM.DELTA == 0.2


def test_robust_aggregate_resists_an_ill_conditioned_member():
    good = [np.array([[-1.9, 0.3], [1.0, 0.0]]) + 0.01 * k for k in range(20)]
    bad = np.array([[4000.0, 0.3], [1.0, 0.0]])
    Js = good + [bad]
    members = list(range(len(Js)))
    med = fp.robust_source_stability(Js, members, "median")["leading_modulus"]
    trim = fp.robust_source_stability(Js, members, "trim10")["leading_modulus"]
    mean = fp.source_period1_stability(Js, members)["lyapunov_numbers"][0]
    assert abs(med - 2.05) < 0.2 and abs(trim - 2.05) < 0.3
    assert mean > 100                                     # the arithmetic mean is dominated
    assert fp.robust_source_stability([None, None], [0, 1]) is None
    with pytest.raises(ValueError):
        fp.robust_source_stability(Js, members, "mean")


@pytest.mark.parametrize("system,snr", [("henon", None), ("henon", 30.0), ("sinusoid", None),
                                        ("two_tone", None), ("white_noise", None)])
def test_pipeline_reproduces_phase4_c2_exactly(system, snr):
    x = S.generate(system, 256, 3, snr)
    out = fp.run_upo_analysis(x, config=C2)
    assert out["embedding_dimension"] == 2
    ref = MM.METHODS["c2_m2M15_mediangate"]["detect"](D.run(x, 2, 15, 50))
    new = out["instability_gated_uop_candidates"] or []
    assert bool(ref) == out["instability_gate_detected"]
    assert sorted(p["loc"] for p in ref) == sorted(float(p["scalar_location"]) for p in new)


def test_decisions_on_clear_cases():
    h = fp.run_upo_analysis(S.generate("henon", 256, 0), config=C2)
    assert h["instability_gate_detected"]
    locs = [p["scalar_location"] for p in h["instability_gated_uop_candidates"]]
    assert min(abs(v - S.fixed_point_reference("henon")) for v in locs) < 0.05
    assert not fp.run_upo_analysis(S.generate("sinusoid", 256, 0), config=C2)["instability_gate_detected"]


def test_fixed_dimension_alone_leaves_levels_unchanged_in_meaning():
    x = S.generate("henon", 256, 1)
    cfg = replace(fp.CFG, upo_fixed_dimension=2, so_assess_significance=True)
    out = fp.run_upo_analysis(x, config=cfg)
    assert out["embedding_dimension"] == 2 and not any(k in out for k in NEW_KEYS)
    assert out["significant_uop_candidates"] is not None
