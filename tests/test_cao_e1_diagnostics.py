"""Cao E1 diagnostics (PROJECT, Phase 3 A5).

`embedding_not_saturated` is returned both when Cao's E1 never plateaus and when
E1 is undefined (exact duplicate delay vectors).  The status is NOT changed; the
opt-in fields below tell the two cases apart.
"""
import sys
import pathlib
from dataclasses import replace

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from upo_systems import fp  # noqa: E402

KEYS = ("cao_e1_undefined", "cao_e1_undefined_duplicate_vectors", "cao_e1_undefined_too_few_vectors",
        "cao_not_saturated", "cao_diagnostics")


def _logistic(r, n, x0=0.3, burn=1000):
    x, out = x0, np.empty(n)
    for i in range(n + burn):
        x = r * x * (1.0 - x)
        if i >= burn:
            out[i - burn] = x
    return out


def _cao(x, tau=1):
    c = fp.CFG
    return fp.cao_method(x, tau, c.cao_max_dim, c.cao_tol, c.cao_theiler, return_diagnostics=True)


def test_default_cao_return_is_unchanged():
    x = _logistic(4.0, 256)
    c = fp.CFG
    plain = fp.cao_method(x, 1, c.cao_max_dim, c.cao_tol, c.cao_theiler)
    diag = _cao(x)
    assert len(plain) == 3 and len(diag) == 4
    assert plain[0] == diag[0]
    np.testing.assert_array_equal(plain[1], diag[1])
    np.testing.assert_array_equal(plain[2], diag[2])


def test_duplicate_vectors_flagged_and_status_unchanged():
    x = _logistic(3.5, 256)                      # stable 4-cycle: 4 distinct values
    m, _, _, d = _cao(x)
    assert m == fp.CFG.cao_max_dim + 1
    assert d["e1_undefined"] and d["e1_undefined_duplicate_vectors"]
    assert not d["e1_undefined_too_few_vectors"]
    assert d["e_undefined_reason_by_m"][1] == "all_neighbours_duplicate"
    cfg = replace(fp.CFG, upo_report_cao_diagnostics=True)
    res = fp.run_upo_analysis(x, config=cfg)
    assert res["status"] == "embedding_not_saturated"          # status unchanged
    assert res["cao_e1_undefined"] and res["cao_e1_undefined_duplicate_vectors"]
    assert res["cao_not_saturated"]


def test_genuine_non_saturation_is_not_flagged_undefined():
    x = np.random.default_rng(4).standard_normal(256)
    m, _, _, d = _cao(x)
    assert not d["e1_undefined"]
    assert d["not_saturated"] == (m > fp.CFG.cao_max_dim)
    assert all(v == 0 for v in d["duplicate_neighbour_count_by_m"].values())


def test_chaotic_series_defined_and_saturated():
    m, _, _, d = _cao(_logistic(4.0, 256))
    assert m <= fp.CFG.cao_max_dim
    assert not d["e1_undefined"] and not d["not_saturated"]


def test_too_few_vectors_flagged():
    _, _, _, d = _cao(_logistic(4.0, 40))       # n < max(30, 3m) for large m
    assert d["e1_undefined"] and d["e1_undefined_too_few_vectors"]
    assert not d["e1_undefined_duplicate_vectors"]


@pytest.mark.parametrize("mode", ["one_sample", "tau_step"])
def test_fields_absent_by_default(mode):
    x = _logistic(3.5, 256)
    cfg = replace(fp.CFG, upo_map_mode=mode)
    assert fp.CFG.upo_report_cao_diagnostics is False
    res = fp.run_upo_analysis(x, config=cfg)
    assert not any(k in res for k in KEYS)
    res_on = fp.run_upo_analysis(x, config=replace(cfg, upo_report_cao_diagnostics=True))
    assert all(k in res_on for k in KEYS)
    assert res_on["status"] == res["status"]
