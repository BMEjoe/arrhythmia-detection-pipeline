"""
Phase 5 detector under test (FROZEN; not tuned in this phase).

    config = fp.phase4_upo_config(replace(fp.CFG, lle_chaos_test=True))
    out = fp.analyze_segment(rr, config)       # raw RR (use_corrected_rr_for_dynamics False)
    LLE component  out["lle_chaos_test"]["detected"]
    UPO component  out["upo"]["instability_gate_detected"]
    PRIMARY        LLE AND UPO

`decide(x, m)` returns both component decisions.  m = 2 is the primary path
through analyze_segment.  The m-sensitivity arms (m = 3, 4; report only) set
BOTH lle_chaos_test_m and upo_fixed_dimension to m and call the two components
directly (lle_chaos_test(rr, cfg) and run_upo_analysis(rr, cfg,
rng=default_rng(cfg.random_seed))), which is exactly what analyze_segment does
for them, without the production LLE / TDMI / Cao work that does not affect
either decision.  Checked against analyze_segment on development windows
(results/dev/m_path_check.json).

If analyze_segment raises (for example ValueError "Insufficient embedded points
after reconstruction" from the production LLE embedding), the error is recorded
and the window counts as NOT detected by every rule; the component decisions
from the direct calls are also recorded so the effect can be reported.
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np

import final_pipeline as fp

CONFIG = fp.phase4_upo_config(replace(fp.CFG, lle_chaos_test=True))


def config_m(m):
    return replace(CONFIG, lle_chaos_test_m=int(m), upo_fixed_dimension=int(m))


def _f(v):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    return v if np.isfinite(v) else None


def _summ(upo, lc):
    return {"upo": bool(upo.get("instability_gate_detected", False)),
            "upo_status": upo.get("status"),
            "upo_n_level_b": (None if upo.get("significant_uop_candidates") is None
                              else len(upo["significant_uop_candidates"])),
            "upo_n_gated": len(upo.get("instability_gated_uop_candidates") or []),
            "lle": bool(lc["detected"]), "lle_p": _f(lc["p"]), "lle_value": _f(lc["lle"]),
            "lle_status": lc["status"]}


def direct(x, cfg):
    x = np.asarray(x, dtype=float)
    try:
        upo = fp.run_upo_analysis(x, config=cfg, rng=np.random.default_rng(cfg.random_seed))
    except np.linalg.LinAlgError as exc:
        upo = {"status": f"analysis_error {exc!r}"}
    return _summ(upo, fp.lle_chaos_test(x, config=cfg))


def decide(x, m=2):
    """Decisions for one window.  m = 2: analyze_segment (primary path)."""
    x = np.asarray(x, dtype=float)
    if m != 2:
        d = direct(x, config_m(m))
        d.update(error=None, path="direct")
        return d
    try:
        out = fp.analyze_segment(x, CONFIG)
    except Exception as exc:                                    # recorded; counts as not detected
        d = direct(x, CONFIG)
        d.update(error=f"{type(exc).__name__}: {exc}", path="direct_after_error",
                 upo_if_no_error=d["upo"], lle_if_no_error=d["lle"], upo=False, lle=False)
        return d
    d = _summ(out["upo"], out["lle_chaos_test"])
    d.update(error=None, path="analyze_segment", tdmi_tau=int(out["tau"]),
             cao_m=int(out["embedding_dimension"]),
             lle_production=_f(out["lle_per_beat"]))
    return d
