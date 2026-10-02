"""Phase 9 opt-in masked_growth_chaos_test (winner k3_mnlp_growth_ann); analyze_segment unchanged."""
import numpy as np
import pytest

import final_pipeline as fp


def logistic_rr(n=512, seed=0):
    rng = np.random.default_rng(seed)
    v = rng.uniform(0.1, 0.9)
    for _ in range(500):
        v = 4.0 * v * (1 - v)
    out = np.empty(n)
    for i in range(n):
        v = 4.0 * v * (1 - v)
        out[i] = v
    return 0.6 + 0.4 * out


def circle_rr(n=512, K=0.9, Om=0.382, seed=0):
    rng = np.random.default_rng(seed)
    th = rng.uniform()
    out = np.empty(n)
    for _ in range(1000):
        th = th + Om - K / (2 * np.pi) * np.sin(2 * np.pi * th)
    for i in range(n):
        nxt = th + Om - K / (2 * np.pi) * np.sin(2 * np.pi * th)
        out[i] = nxt - th
        th = nxt
    return out * (0.8 / out.mean())


def test_fields_fixed_at_preregistered_values():
    assert fp.CFG.masked_growth_z == 9.4521 and fp.CFG.masked_growth_g == 0.25
    assert fp.CFG.masked_growth_surrogates == 39


def test_default_output_unchanged():
    rr = 0.8 + 0.05 * np.random.default_rng(1).standard_normal(300)
    out = fp.analyze_segment(rr)
    assert "masked_growth_test" not in out


def test_requires_labels():
    r = fp.masked_growth_chaos_test(logistic_rr(), None)
    assert r["detected"] is False and r["analysable"] is False
    with pytest.raises(ValueError):
        fp.masked_growth_chaos_test(logistic_rr(), ["N"] * 10)


def test_chaotic_map_detected_quasiperiodic_not():
    lab = ["N"] * 513
    assert fp.masked_growth_chaos_test(logistic_rr(), lab)["detected"] is True
    r = fp.masked_growth_chaos_test(circle_rr(), lab)
    assert r["detected"] is False


def test_mask_excludes_ectopy_and_bigeminy_unanalysable():
    x = logistic_rr()
    lab = np.array(["N"] * 513)
    lab[1::2] = "V"                                     # bigeminy: every interval ectopy-related
    r = fp.masked_growth_chaos_test(x, lab)
    assert r["analysable"] is False and r["detected"] is False and r["n_masked"] == 512


def test_matches_preregistered_candidate():
    from experiments.phase9_waveform import candidates9 as C
    x = circle_rr(seed=3) + 0.004 * np.random.default_rng(3).standard_normal(512)
    lab = np.array(["N"] * 513)
    lab[[40, 41, 200, 333]] = "V"
    a = fp.masked_growth_chaos_test(x, lab)
    b = C.run_all(x, lab, names=("k3_mnlp_growth_ann",))["k3_mnlp_growth_ann"]
    assert a["zmax"] == b["zmax"] and a["G"] == b["G"] and a["detected"] == b["detected"]


def test_analyze_segment_signature_unchanged():
    import inspect
    assert list(inspect.signature(fp.analyze_segment).parameters) == ["rr_intervals", "config"]
