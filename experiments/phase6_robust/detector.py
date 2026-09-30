"""Phase 6: one window through fp.analyze_segment with a method's configuration.
A raised exception is recorded and counts as not detected by every rule."""
from __future__ import annotations

import numpy as np

import final_pipeline as fp


def _f(v):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    return v if np.isfinite(v) else None


def decide(rr, cfg):
    try:
        out = fp.analyze_segment(np.asarray(rr, dtype=float), cfg)
    except Exception as exc:                                     # recorded; not detected
        return {"lle": False, "upo": False, "error": f"{type(exc).__name__}: {exc}"}
    lc, upo = out["lle_chaos_test"], out["upo"]
    return {"lle": bool(lc["detected"]), "lle_p": _f(lc["p"]), "lle_status": lc["status"],
            "upo": bool(upo.get("instability_gate_detected", False)), "upo_status": upo.get("status"),
            "upo_n_level_b": (None if upo.get("significant_uop_candidates") is None
                              else len(upo["significant_uop_candidates"])),
            "upo_n_gated": len(upo.get("instability_gated_uop_candidates") or []),
            "and": fp.combined_chaos_detected(out),
            "lle_embedding_too_short": bool(out.get("lle_embedding_too_short", False)), "error": None}
