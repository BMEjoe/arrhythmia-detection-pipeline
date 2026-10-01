"""
Phase 7 frozen detector and per-window record (PREREGISTRATION.md Section A).

DETECTOR = combined_chaos_config(replace(CFG, keep_upo_on_short_lle_embedding=True), detrend=True)
decision = combined_chaos_detected(analyze_segment(rr, DETECTOR))

Continuous scores (read from the existing analyze_segment output; nothing in the
detector is changed or re-run):

  lle_z      = (S - mean(S_sur)) / SD(S_sur, ddof=1), with S = lle_chaos_test["lle"] (the
               Phase 3 statistic) and S_sur its finite IAAFT surrogate values
               (lle_chaos_test["surrogate_lle"]).  NaN when the statistic is undefined.

  upo_score  = max(0, max over period-1 Level-A peaks that pass the Phase 4 instability
               gate of the peak's rJ = deviation / W0), with deviation and W0 from the
               period-1 So surrogate test (upo.periods[1].significance.per_peak, aligned
               with upo.periods[1].source_peak_candidates) and the gate = robust
               (element-wise median) source-stability leading modulus >= 1 + delta
               (fp.robust_source_stability on upo.periods[1].jacobians, exactly as
               fp._attach_instability_gate does for Level-B peaks).  0 when no peak
               passes the gate or the UPO analysis ended in a failure status / without
               significance.  The UPO decision is "some gate-passing peak has J < 0.05";
               J is a decreasing function of the deviation and W0 is fixed within a
               window, so the decision is a within-window threshold on this score.
"""
from __future__ import annotations

import time
import traceback
from dataclasses import replace

import numpy as np

import final_pipeline as fp

DETECTOR = fp.combined_chaos_config(replace(fp.CFG, keep_upo_on_short_lle_embedding=True), detrend=True)
METHOD = "phase6_recommended_and"


def lle_z(lc):
    if lc is None or not np.isfinite(lc.get("lle", np.nan)):
        return np.nan
    sur = np.asarray(lc.get("surrogate_lle", []), dtype=float)
    sur = sur[np.isfinite(sur)]
    if len(sur) < 2:
        return np.nan
    sd = float(np.std(sur, ddof=1))
    return float((lc["lle"] - np.mean(sur)) / sd) if sd > 0 else np.nan


def upo_peak_table(upo, config=DETECTOR):
    """Per period-1 Level-A peak: rJ, J, deviation, robust leading modulus, gate pass."""
    if upo.get("status") in fp.UPO_FAILURE_STATUSES:
        return []
    p1 = (upo.get("periods") or {}).get(1)
    if p1 is None or not p1.get("significance_assessed"):
        return []
    sig = p1["significance"]
    rows = []
    for peak, info in zip(p1[fp.LEVEL_A_FIELD], sig["per_peak"]):
        st = None
        if p1.get("jacobians") is not None:
            st = fp.robust_source_stability(p1["jacobians"], peak["source_member_indices"],
                                            config.upo_instability_gate_aggregate)
        lead = np.nan if st is None else st["leading_modulus"]
        gate = bool(np.isfinite(lead) and lead >= 1.0 + config.upo_instability_gate_delta)
        rows.append({"rJ": float(info["rJ"]), "J": float(info["J"]), "deviation": float(info["deviation"]),
                     "leading_modulus": float(lead), "gate": gate, "significant": bool(info["significant"]),
                     "location": _scalar(peak.get("scalar_location"))})
    return rows


def _scalar(v):
    try:
        return float(np.asarray(v, dtype=float).reshape(-1)[0])
    except (TypeError, ValueError, IndexError):
        return float("nan")


def upo_score(peaks):
    vals = [p["rJ"] for p in peaks if p["gate"] and np.isfinite(p["rJ"])]
    return float(max(0.0, max(vals))) if vals else 0.0


def evaluate(rr, config=DETECTOR):
    """Run the frozen detector on one window; return a JSON-serialisable record.
    Errors (analyze_segment raising) count as not detected and are recorded."""
    t = time.time()
    rec = {"method": METHOD}
    try:
        out = fp.analyze_segment(np.asarray(rr, dtype=float), config)
    except Exception as exc:                                # noqa: BLE001 - recorded, counted not detected
        rec.update({"error": f"{type(exc).__name__}: {exc}", "traceback": traceback.format_exc(limit=3),
                    "and_detected": False, "lle_detected": False, "upo_detected": False,
                    "lle_stat": None, "lle_p": None, "lle_z": None, "upo_score": None,
                    "runtime_s": time.time() - t})
        return rec
    lc = out["lle_chaos_test"]
    upo = out["upo"]
    peaks = upo_peak_table(upo, config)
    gated_b = upo.get("instability_gated_uop_candidates") or []
    upo_det = bool(upo.get("instability_gate_detected", False))
    # consistency: the UPO decision equals "a gate-passing peak is significant"
    consistent = upo_det == any(p["gate"] and p["significant"] for p in peaks)
    sur = np.asarray(lc.get("surrogate_lle", []), dtype=float)
    fin = sur[np.isfinite(sur)] if sur.size else sur
    det = out.get("rr_detrend") or {}
    ratio = det.get("linear_trend_ratio")
    rec.update({
        "error": None,
        "and_detected": bool(fp.combined_chaos_detected(out)),
        "lle_detected": bool(lc["detected"]), "upo_detected": upo_det,
        "lle_stat": None if not np.isfinite(lc.get("lle", np.nan)) else float(lc["lle"]),
        "lle_p": None if not np.isfinite(lc.get("p", np.nan)) else float(lc["p"]),
        "lle_status": lc.get("status"),
        "lle_sur_mean": float(np.mean(fin)) if fin.size else None,
        "lle_sur_sd": float(np.std(fin, ddof=1)) if fin.size > 1 else None,
        "lle_n_sur_undefined": int(sur.size - fin.size),
        "lle_z": None if not np.isfinite(lle_z(lc)) else lle_z(lc),
        "upo_score": upo_score(peaks),
        "upo_status": upo.get("status"), "upo_significance_status": upo.get("significance_status"),
        "upo_n_level_a": int(len(upo.get(fp.LEVEL_A_FIELD) or [])),
        "upo_n_level_b": None if upo.get(fp.LEVEL_B_FIELD) is None else int(len(upo[fp.LEVEL_B_FIELD])),
        "upo_n_gated": int(len(gated_b)),
        "upo_max_gated_modulus": (float(max(p["robust_source_stability"]["leading_modulus"] for p in gated_b))
                                  if gated_b else None),
        "upo_source_rJ": (float(upo["source_rJ"]) if upo.get("source_rJ") is not None
                          and np.isfinite(upo.get("source_rJ")) else None),
        "upo_peaks": peaks, "upo_score_consistent": bool(consistent),
        "lle_embedding_too_short": bool(out.get("lle_embedding_too_short", False)),
        "detrend_ratio": None if ratio is None or not np.isfinite(ratio) else float(ratio),
        "detrend_applied": bool(det.get("applied", False)),
        "runtime_s": time.time() - t,
    })
    return rec
