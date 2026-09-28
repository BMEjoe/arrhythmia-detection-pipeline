"""
Phase 4 UPO detector runs and per-peak features.

A RUN is one execution of the So et al. period-1 detector with surrogate
significance and Level-C verification, under one configuration:
    m_rule   "cao"  -> Cao at lag 1 (production), via fp.run_upo_analysis unchanged
             int m  -> fixed lag-1 embedding dimension m (same pipeline steps)
    M        so_jacobian_neighbors (production 7)
    n_sur    so_surrogate_count
A METHOD (baseline / candidate) is a fixed decision rule applied to the saved
per-peak features of one run (gates.py-style functions below), so methods that
share a run configuration share the computation.

Per-peak features: scalar location, surrogate J and significance (Level B),
Level-C verification (residual, R2, support, passes), candidate-centred
extension monodromy moduli, and the hybrid PRL/PRE source-stability leading
modulus computed from the member Jacobians with three aggregates:
arithmetic mean (production, outlier-sensitive; Phase 3 A2), element-wise
median, and 10 % trimmed mean.  Pipeline functions are imported, not copied.
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np
from scipy.stats import trim_mean

import final_pipeline as fp

CFG = fp.CFG


def _f(v):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    return v if np.isfinite(v) else None


def run_key(m_rule, M, n_sur):
    return f"m{m_rule}_M{int(M)}_S{int(n_sur)}"


def _config(M, n_sur):
    return replace(CFG, so_assess_significance=True, so_surrogate_count=int(n_sur),
                   so_jacobian_neighbors=int(M))


def _lead_modulus(S):
    if S is None or not np.all(np.isfinite(S)):
        return None
    return _f(np.max(np.abs(np.linalg.eigvals(S))))


def _source_moduli(jacobians, members):
    Js = [jacobians[i] for i in np.unique(np.asarray(members, dtype=int)) if jacobians[i] is not None]
    if not Js:
        return {"n_members": 0, "src_mean": None, "src_median": None, "src_trim10": None}
    A = np.asarray(Js, dtype=float)
    return {"n_members": int(len(Js)),
            "src_mean": _lead_modulus(A.mean(axis=0)),
            "src_median": _lead_modulus(np.median(A, axis=0)),
            "src_trim10": _lead_modulus(trim_mean(A, 0.1, axis=0))}


def _period1(x, m_rule, cfg):
    """Returns (period-1 result dict or None, status, m, n_embedded)."""
    x = np.asarray(x, dtype=float)
    if m_rule == "cao":
        upo = fp.run_upo_analysis(x, config=cfg)            # production path, unchanged
        return upo["periods"].get(1), upo["status"], upo["embedding_dimension"], upo["n_embedded_points"]
    m = int(m_rule)
    if not np.all(np.isfinite(x)) or float(np.ptp(x)) <= 1e-12:
        return None, "constant_or_nonfinite", m, 0
    emb = fp._embed_backward(x, 1, m)
    r = fp.detect_so_fixed_points(emb, tau=1, config=cfg, rng=np.random.default_rng(cfg.random_seed))
    fp.assess_so_significance(r, x, config=cfg, rng=np.random.default_rng(cfg.random_seed + 104729),
                              embedding_dimension=m)
    if cfg.so_verify_peaks and r["status"] not in fp.UPO_FAILURE_STATUSES:
        fp.apply_verification_gates(r, emb, config=cfg)
    return r, r["status"], m, len(emb)


def run(x, m_rule="cao", M=7, n_sur=50):
    cfg = _config(M, n_sur)
    r, status, m, ne = _period1(x, m_rule, cfg)
    out = {"status": status, "m": None if m is None else int(m), "n_embedded": int(ne or 0),
           "significance_assessed": bool(r and r.get("significance_assessed")), "peaks": []}
    if r is None:
        return out
    sig = r.get("significance") or {}
    out["J_W"] = _f(sig.get("J_W"))
    out["rJ"] = _f(sig.get("rJ"))
    per_peak = sig.get("per_peak") or []
    levelc = {}
    for e in (r.get(fp.LEVEL_C_PASS_FIELD) or []) + (r.get(fp.LEVEL_C_FAIL_FIELD) or []):
        levelc[tuple(np.atleast_1d(e["peak_cell"]).tolist())] = e
    jac = r.get("jacobians")
    for i, p in enumerate(r.get(fp.LEVEL_A_FIELD) or []):
        cell = tuple(np.atleast_1d(p["peak_cell"]).tolist())
        c = levelc.get(cell, {})
        ver = c.get("verification") or {}
        stab = c.get("extension_stability")
        info = per_peak[i] if i < len(per_peak) else {}
        rec = {"loc": _f(p["scalar_location"]), "count": int(p["histogram_count"]),
               "J": _f(info.get("J")), "significant": bool(info.get("significant", False)),
               "peak_rJ": _f(info.get("rJ")),
               "levelC_passes": bool(ver.get("passes", False)), "residual": _f(ver.get("residual")),
               "r2": _f(ver.get("r2")), "support": _f(ver.get("support_ratio")),
               "ext_lead": None if stab is None else _f(stab["multiplier_moduli"][0]),
               "ext_moduli": None if stab is None else [_f(v) for v in stab["multiplier_moduli"]],
               "verified_unstable": bool(c.get("verified_unstable", False))}
        rec.update(_source_moduli(jac, p["source_member_indices"]) if jac is not None
                   else {"n_members": 0, "src_mean": None, "src_median": None, "src_trim10": None})
        out["peaks"].append(rec)
    return out


# -----------------------------------------------------------------------------
# Decision rules on a run's per-peak features
# -----------------------------------------------------------------------------

def peaks_level_b(run_out):
    return [p for p in run_out.get("peaks", []) if p["significant"]]


def gate(peaks, stat, delta):
    """Keep peaks whose `stat` leading modulus is >= 1 + delta (undefined -> removed)."""
    return [p for p in peaks if p.get(stat) is not None and p[stat] >= 1.0 + delta]


def loc_error(peaks, ref):
    if ref is None or not peaks:
        return None
    return float(min(abs(p["loc"] - ref) for p in peaks if p["loc"] is not None))
