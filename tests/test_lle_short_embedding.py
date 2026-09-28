"""keep_upo_on_short_lle_embedding (PROJECT option, default off).

When the LLE delay embedding has fewer than 50 points, analyze_segment raises
by default (previous behaviour).  With the option on, the UPO result is kept and
only the LLE is marked as failed.
"""
import sys
import pathlib
from dataclasses import replace

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import upo_systems as us  # noqa: E402
from upo_systems import fp  # noqa: E402


def _short_series():
    # 128 samples of the skewed Henon map: TDMI tau = 8, Cao m = 13 -> 32 LLE
    # embedded points (< 50), while the lag-1 UPO embedding is long enough.
    x = us.henon_series(128, seed=0)
    emb = fp.takens_embed(x, config=fp.CFG)[0]
    assert len(emb) < 50
    return x


def test_default_is_off_and_still_raises():
    assert fp.CFG.keep_upo_on_short_lle_embedding is False
    with pytest.raises(ValueError, match="Insufficient embedded points after reconstruction"):
        fp.analyze_segment(_short_series(), config=fp.CFG)


def test_option_keeps_upo_and_fails_only_lle():
    x = _short_series()
    cfg = replace(fp.CFG, keep_upo_on_short_lle_embedding=True)
    res = fp.analyze_segment(x, config=cfg)
    assert res["lle_embedding_too_short"] is True
    assert np.isnan(res["lle_per_beat"])
    assert res["lle_diagnostics"]["reason"] == "insufficient embedded points after reconstruction"
    assert res["lle_diagnostics"]["valid_fit"] is False
    assert res["lle_diagnostics"]["n_embedded_points"] < 50
    # the UPO result is exactly what run_upo_analysis returns on its own
    direct = fp.run_upo_analysis(x, config=cfg, lle_tau=res["tau"], lle_m=res["embedding_dimension"],
                                 rng=np.random.default_rng(cfg.random_seed))
    assert res["upo_status"] == direct["status"]
    assert res["upo_status"] != "analysis_error"
    assert res["upo_embedding_dimension"] == direct["embedding_dimension"]
    assert res["upo_n_embedded_points"] == direct["n_embedded_points"]
    feats = fp.upo_summary_features(res["upo"])
    feats_direct = fp.upo_summary_features(direct)
    assert feats.keys() == feats_direct.keys()
    for k in feats:
        a, b = feats[k], feats_direct[k]
        assert (a == b) or (isinstance(a, float) and isinstance(b, float) and np.isnan(a) and np.isnan(b)), k


def test_option_does_not_change_long_windows():
    x = us.henon_series(512, seed=0)
    base = fp.analyze_segment(x, config=fp.CFG)
    opt = fp.analyze_segment(x, config=replace(fp.CFG, keep_upo_on_short_lle_embedding=True))
    assert opt["lle_embedding_too_short"] is False
    assert "lle_embedding_too_short" not in base
    assert base["lle_per_beat"] == opt["lle_per_beat"]
    assert base["upo_status"] == opt["upo_status"]
    assert fp.upo_summary_features(base["upo"]).keys() == fp.upo_summary_features(opt["upo"]).keys()
