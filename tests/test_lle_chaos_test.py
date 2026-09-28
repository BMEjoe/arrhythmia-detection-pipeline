"""PipelineConfig.lle_chaos_test (PROJECT, Phase 3 B5; default off).

The pipeline implementation must reproduce the preregistered Phase 3 winner
(experiments/phase3_lle/estimators.py::c1_rosenstein_m2_iaaft) exactly.
"""
import sys
import pathlib
from dataclasses import replace

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))
import final_pipeline as fp  # noqa: E402
from experiments.phase3_lle import estimators as E  # noqa: E402
import upo_systems as us  # noqa: E402

ON = replace(fp.CFG, lle_chaos_test=True)


def _logistic(n, r=4.0, x0=0.3, burn=500):
    x, out = x0, np.empty(n)
    for i in range(n + burn):
        x = r * x * (1.0 - x)
        if i >= burn:
            out[i - burn] = x
    return out


def _series():
    rng = np.random.default_rng(11)
    ar = np.zeros(256)
    for i in range(1, 256):
        ar[i] = 0.8 * ar[i - 1] + rng.standard_normal()
    return {"logistic": _logistic(256), "henon": us.henon_series(256),
            "noisy_logistic": _logistic(256) + 0.1 * np.std(_logistic(256)) * rng.standard_normal(256),
            "white_noise": rng.standard_normal(256), "ar1": ar,
            "sinusoid": np.sin(2 * np.pi * np.arange(256) / 7.3 + 0.4),
            "period4": _logistic(256, r=3.5)}


def test_defaults_are_the_preregistered_values_and_off():
    c = fp.CFG
    assert c.lle_chaos_test is False
    assert (c.lle_chaos_test_m, c.lle_chaos_test_fit_first, c.lle_chaos_test_fit_last) == (2, 1, 5)
    assert (c.lle_chaos_test_surrogates, c.lle_chaos_test_alpha) == (99, 0.05)
    assert c.lle_chaos_test_seed_entropy == E.SURROGATE_ENTROPY


@pytest.mark.parametrize("name", list(_series()))
def test_pipeline_reproduces_phase3_estimator_exactly(name):
    x = _series()[name]
    a = fp.lle_chaos_test(x, config=ON)
    b = E.c1_rosenstein_m2_iaaft(x)
    assert a["detected"] == b["detected"]
    if b["lle"] is None:
        assert np.isnan(a["lle"]) and b["p"] is None
    else:
        assert a["lle"] == b["lle"]
        assert a["p"] == b["p"]
        assert a["surrogate_median"] == b["surrogate_median"]


def test_decisions_on_clear_cases():
    s = _series()
    assert fp.lle_chaos_test(s["logistic"], ON)["detected"]
    assert fp.lle_chaos_test(s["henon"], ON)["detected"]
    assert abs(fp.lle_chaos_test(s["logistic"], ON)["lle"] - np.log(2)) < 0.05
    assert not fp.lle_chaos_test(s["sinusoid"], ON)["detected"]
    r = fp.lle_chaos_test(s["period4"], ON)
    assert not r["detected"] and r["status"] != "ok"


def test_iaaft_surrogate_properties():
    x = us.henon_series(256)
    sur = fp.iaaft_surrogate(x, np.random.default_rng(0))
    np.testing.assert_array_equal(np.sort(sur), np.sort(x))
    np.testing.assert_array_equal(sur, E.iaaft_surrogate(x, np.random.default_rng(0)))


def test_analyze_segment_default_unchanged_and_opt_in_adds_result():
    x = us.henon_series(512)
    base = fp.analyze_segment(x, config=fp.CFG)
    on = fp.analyze_segment(x, config=ON)
    assert "lle_chaos_test" not in base
    assert on["lle_chaos_test"]["detected"] is True
    assert base["lle_per_beat"] == on["lle_per_beat"]            # production LLE untouched
    assert base["upo_status"] == on["upo_status"]
