"""
Phase 2E result contract.

Rows are extracted from the ACTUAL outputs of final_pipeline.analyze_segment /
run_upo_analysis (upo_summary_features is used for every Level A/B/C
population so that Phase 2E counts are exactly the classifier's features).
Failed windows are never dropped: they keep their status and NaN fields.
"""
from __future__ import annotations

import math

import numpy as np

import final_pipeline as fp

NAN = float("nan")


def _f(v):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return NAN
    return v


def _len_or_nan(v):
    return NAN if v is None else float(len(v))


def _lle_fields(result):
    diag = result.get("lle_diagnostics") or {}
    lle = _f(result.get("lle_per_beat"))
    m = result.get("embedding_dimension")
    return {
        "lle_per_beat": lle,
        "lle_valid": bool(math.isfinite(lle)),
        "lle_r2": _f(diag.get("r2")),
        "lle_valid_fit_quality": bool(diag.get("valid_fit_quality", False)),
        "lle_reason": diag.get("reason", "") if not math.isfinite(lle) else "",
        "lle_theiler_beats": result.get("lle_theiler_beats"),
        "lle_tau": result.get("tau"),
        "lle_m": m,
        "lle_embedding_not_saturated": bool(m is not None and m > fp.CFG.cao_max_dim),
        "lle_n_embedded_points": int(len(result["embedded"])) if result.get("embedded") is not None else 0,
    }


def _significance_fields(p1, config):
    out = {"W": NAN, "W0": NAN, "J_W": NAN, "J_W0": NAN, "rJ": NAN,
           "p_value_finite": NAN, "n_surrogates_used": NAN,
           "n_surrogates_requested": int(config.so_surrogate_count) if config.so_assess_significance else 0,
           "significance_not_assessed_reason": ""}
    sig = None if p1 is None else p1.get("significance")
    if not sig:
        return out
    if "not_assessed_reason" in sig:
        out["significance_not_assessed_reason"] = str(sig["not_assessed_reason"])
        return out
    out.update({"W": _f(sig["W"]), "W0": _f(sig["W0"]), "J_W": _f(sig["J_W"]),
                "J_W0": _f(sig["J_W0"]), "rJ": _f(sig["rJ"]),
                "p_value_finite": _f(sig["p_value_finite"]),
                "n_surrogates_used": int(sig["n_surrogates"])})
    return out


def _peak_detail(peak, sig_by_cell, gated, failed):
    stab = peak.get("source_stability")
    entry = {
        "period": int(peak["period"]),
        "location": [float(v) for v in np.atleast_1d(peak["orbit_coordinates"])],
        "minimal_period": int(peak.get("minimal_period", peak["period"])),
        "histogram_count": int(peak["histogram_count"]),
        "peak_cell": [int(c) for c in peak["peak_cell"]],
        "source_stability_lyapunov_numbers": (None if stab is None else
                                              [float(v) for v in stab["lyapunov_numbers"]]),
        "source_stability_unstable": None if stab is None else bool(stab["unstable"]),
        "source_stability_n_points": None if stab is None else int(stab["n_points_averaged"]),
    }
    info = sig_by_cell.get(tuple(entry["peak_cell"]) + (entry["period"],))
    if info is not None:
        entry.update({"J": float(info["J"]), "deviation": float(info["deviation"]),
                      "peak_rJ": _f(info["rJ"]), "significant": bool(info["significant"])})
    key = tuple(entry["peak_cell"]) + (entry["period"],)
    lvl_c = gated.get(key) or failed.get(key)
    if lvl_c is not None:
        ver = lvl_c["verification"]
        ext = lvl_c.get("extension_stability")
        entry.update({
            "levelC_passes": bool(ver["passes"]),
            "levelC_failed_gates": list(ver["failed_gates"]),
            "levelC_residual": _f(ver["residual"]), "levelC_r2": _f(ver["r2"]),
            "levelC_support_ratio": _f(ver["support_ratio"]),
            "extension_multiplier_moduli": (None if ext is None else
                                            [float(v) for v in ext["multiplier_moduli"]]),
            "verified_unstable": bool(lvl_c.get("verified_unstable", False)),
        })
    return entry


def upo_fields(upo, config):
    """Row fields + detail for a run_upo_analysis output."""
    feats = fp.upo_summary_features(upo)
    p1 = (upo.get("periods") or {}).get(1)
    status = upo["status"]
    row = {
        "map_mode": upo["map_mode"],
        "source_faithful_map": bool(upo["source_faithful_map"]),
        "upo_tau": upo["tau"],
        "embedding_dimension": upo["embedding_dimension"],
        "n_embedded_points": int(upo["n_embedded_points"]),
        "embedding_not_saturated": status == "embedding_not_saturated",
        "upo_status": status,
        "upo_failure": status in fp.UPO_FAILURE_STATUSES,
        "upo_no_peak_status": status in fp.UPO_NO_PEAK_STATUSES,
        "period1_status": None if p1 is None else p1["status"],
        "upo_significance_status": upo["significance_status"],
        "significance_requested": bool(upo["significance_requested"]),
        "significance_assessed": bool(upo["significance_assessed"]),
        "n_in_tube": NAN if p1 is None or p1.get("scalar") is None else int(len(p1["scalar"])),
        "peak_threshold": _f(None if p1 is None else p1.get("peak_threshold")),
        **feats,
        "verification_assessed": bool(upo["verification_assessed"]),
        "verification_gated_count": (NAN if status in fp.UPO_FAILURE_STATUSES
                                     else _len_or_nan(upo.get(fp.LEVEL_C_PASS_FIELD))),
        "verification_failed_count": (NAN if status in fp.UPO_FAILURE_STATUSES
                                      else _len_or_nan(upo.get(fp.LEVEL_C_FAIL_FIELD))),
        **_significance_fields(p1, config),
    }
    # Source-stability diagnostic summary (HYBRID PRL/PRE, attached after detection)
    peaks = upo.get(fp.LEVEL_A_FIELD) or []
    stab = [c.get("source_stability") for c in peaks]
    row["source_stability_available_count"] = float(sum(s is not None for s in stab)) \
        if status not in fp.UPO_FAILURE_STATUSES else NAN
    row["source_stability_unstable_count"] = float(sum(bool(s and s["unstable"]) for s in stab)) \
        if status not in fp.UPO_FAILURE_STATUSES else NAN

    sig_by_cell, gated, failed = {}, {}, {}
    for per, r in (upo.get("periods") or {}).items():
        sig = r.get("significance") or {}
        for peak, info in zip(r.get(fp.LEVEL_A_FIELD, []), sig.get("per_peak", []) or []):
            sig_by_cell[tuple(int(c) for c in peak["peak_cell"]) + (int(per),)] = info
        for c in r.get(fp.LEVEL_C_PASS_FIELD) or []:
            gated[tuple(int(v) for v in c["peak_cell"]) + (int(per),)] = c
        for c in r.get(fp.LEVEL_C_FAIL_FIELD) or []:
            failed[tuple(int(v) for v in c["peak_cell"]) + (int(per),)] = c
    detail = {"peaks": [_peak_detail(c, sig_by_cell, gated, failed) for c in peaks]}
    if p1 is not None:
        sig = p1.get("significance") or {}
        if "surrogate_W" in sig:
            detail.update({
                "surrogate_W": [float(v) for v in sig["surrogate_W"]],
                "observed_histogram": [float(v) for v in np.ravel(sig["observed_histogram"])],
                "surrogate_mean_histogram": [float(v) for v in np.ravel(sig["surrogate_mean_histogram"])],
            })
        if p1.get("edges") is not None:
            detail["histogram_edges"] = [float(v) for v in p1["edges"][0]]
    cov = upo.get("coverage") or {}
    src_cov = cov.get("source") or {}
    row["coverage_radius"] = _f(src_cov.get("radius"))
    return row, detail


def segment_fields(result, config):
    row = _lle_fields(result)
    upo_row, detail = upo_fields(result["upo"], config)
    row.update(upo_row)
    row["analysis_error"] = ""
    return row, detail


def error_fields(exc, upo_direct, config):
    """
    analyze_segment raised.  The production MIT-BIH wrapper then records
    upo_status = "analysis_error" with NaN features; that is reproduced here.
    upo_direct_status is a Phase 2E DIAGNOSTIC: the status run_upo_analysis
    itself returns for the same window (shows what the wrapper masks).
    """
    row = {
        "lle_per_beat": NAN, "lle_valid": False, "lle_r2": NAN, "lle_valid_fit_quality": False,
        "lle_reason": "analyze_segment raised", "lle_theiler_beats": None, "lle_tau": None,
        "lle_m": None, "lle_embedding_not_saturated": False, "lle_n_embedded_points": 0,
        "map_mode": config.upo_map_mode, "source_faithful_map": config.upo_map_mode == "one_sample",
        "upo_tau": None, "embedding_dimension": None, "n_embedded_points": 0,
        "embedding_not_saturated": False, "upo_status": "analysis_error", "upo_failure": True,
        "upo_no_peak_status": False, "period1_status": None,
        "upo_significance_status": (fp.SIGNIFICANCE_DETECTOR_FAILURE if config.so_assess_significance
                                    else fp.SIGNIFICANCE_NOT_REQUESTED),
        "significance_requested": bool(config.so_assess_significance),
        "significance_assessed": False,
        "analysis_error": f"{type(exc).__name__}: {exc}",
    }
    for k in fp.upo_summary_features({"status": "analysis_error"}):
        row[k] = NAN
    row.update({"verification_assessed": False, "verification_gated_count": NAN,
                "verification_failed_count": NAN, "n_in_tube": NAN, "peak_threshold": NAN,
                "source_stability_available_count": NAN, "source_stability_unstable_count": NAN,
                "coverage_radius": NAN, **_significance_fields(None, config)})
    if upo_direct is not None:
        row["upo_direct_status"] = upo_direct["status"]
        row["upo_direct_embedding_dimension"] = upo_direct["embedding_dimension"]
    return row
