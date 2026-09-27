"""
Phase 2E -- POST-HOC DIAGNOSTICS (NOT part of the predeclared grid).

These diagnostics were run during the read-only audit of 2026-09-27, AFTER
part of the predeclared Phase 2E output had been inspected.  They reproduce
exactly the audit diagnostics quoted in docs/PHASE2E_SYNTHETIC_VALIDATION.md
(Part B) and nothing else.  They call the frozen final_pipeline functions
unchanged; nothing here is substituted into any predeclared result.

    python -m experiments.phase2e.posthoc_diagnostics

Writes results/posthoc/{lle,tdmi,levelc}.json.
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import json
import math
import pathlib
import sys
import warnings

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import final_pipeline as fp  # noqa: E402
from experiments.phase2e import systems as S  # noqa: E402

CFG = fp.CFG
OUT = pathlib.Path(__file__).resolve().parent / "results" / "posthoc"
HENON_LLE = 0.4192   # literature value for a=1.4, b=0.3 (per iteration, natural log)
TRUE_FP = {"logistic": 0.75, "henon": S.henon_fixed_points()[0]}


def _f(v):
    v = float(v)
    return v if math.isfinite(v) else None


def _lle(emb, theiler):
    return fp.rosenstein_lle(emb, theiler=theiler, max_iter=CFG.lle_max_iter,
                             min_fit_points=CFG.lle_min_fit_points,
                             max_fit_fraction=CFG.lle_max_fit_fraction, min_r2=CFG.lle_min_r2)


# -----------------------------------------------------------------------------
# D1 -- LLE: production (TDMI tau, Cao m at tau) vs map-compatible tau = 1
# -----------------------------------------------------------------------------

def lle_diagnostics():
    out = {"note": "lle in natural log per sample; one sample = one map iteration",
           "references": {"logistic": math.log(2.0), "henon": HENON_LLE},
           "production_vs_tau1_N512_seeds0_9": {}, "divergence_curves_seed0_N512": {},
           "grid_median_seeds0_4": {}, "theiler_fit_range_N2048_seed0": {}}
    for system in ("logistic", "henon", "white_noise"):
        rows = []
        for s in range(10):
            x = S.generate(system, 512, s)
            emb, tau, m, *_ = fp.takens_embed(x, config=CFG)
            th = fp._mean_period_beats(x)
            v, d = _lle(emb, th)
            rows.append({"seed": s, "tau": tau, "m": m, "theiler": th, "lle": _f(v), "r2": _f(d.get("r2", np.nan))})
        out["production_vs_tau1_N512_seeds0_9"][system] = rows

        x = S.generate(system, 512, 0)
        emb, tau, m, *_ = fp.takens_embed(x, config=CFG)
        th = fp._mean_period_beats(x)
        curves = {}
        v, d = _lle(emb, th)
        curves[f"production tau={tau} m={m}"] = {"slope": _f(v), "r2": _f(d["r2"]),
                                                 "y": [_f(t) for t in d["mean_log_divergence"][:12]]}
        for mm in ((1, 2) if system != "white_noise" else (2,)):
            v1, d1 = _lle(fp._embed_backward(x, 1, mm), th)
            curves[f"tau=1 m={mm}"] = {"slope": _f(v1), "r2": _f(d1["r2"]),
                                       "y": [_f(t) for t in d1["mean_log_divergence"][:12]]}
        curves["ln_attractor_span"] = _f(np.log(np.ptp(x)))
        out["divergence_curves_seed0_N512"][system] = curves

    for system, mtrue in (("logistic", 1), ("henon", 2), ("white_noise", 2)):
        for N in (256, 512, 1024, 2048):
            res = {}
            for s in range(5):
                x = S.generate(system, N, s)
                th = fp._mean_period_beats(x)
                tau_p, m_p = fp.takens_embed(x, config=CFG)[1:3]
                m1 = fp.cao_method(x, 1, CFG.cao_max_dim, CFG.cao_tol, CFG.cao_theiler)[0]
                for lab, (t, m) in {"production(tau*,m*)": (tau_p, m_p), "tau1,m_cao_lag1": (1, min(m1, 12)),
                                    f"tau1,m={mtrue}": (1, mtrue), f"tau*,m={mtrue}": (tau_p, mtrue)}.items():
                    res.setdefault(lab, []).append(_lle(fp._embed_backward(x, t, m), th)[0])
            out["grid_median_seeds0_4"][f"{system}/N={N}"] = {k: _f(np.nanmedian(v)) for k, v in res.items()}

    for system, mtrue in (("logistic", 1), ("henon", 2)):
        x = S.generate(system, 2048, 0)
        emb = fp._embed_backward(x, 1, mtrue)
        for th in (0, 1, 5, 20):
            v, d = _lle(emb, th)
            y = d["mean_log_divergence"]
            out["theiler_fit_range_N2048_seed0"][f"{system}/theiler={th}"] = {
                "production_slope": _f(v), "r2": _f(d["r2"]),
                "slope_first_K": {K: _f(np.polyfit(np.arange(K), y[:K], 1)[0]) for K in (2, 3, 4, 6, 10)}}
    return out


# -----------------------------------------------------------------------------
# D2 -- TDMI first-minimum behaviour
# -----------------------------------------------------------------------------

def tdmi_diagnostics():
    out = {"config": {"tdmi_max_tau": CFG.tdmi_max_tau, "tdmi_bins": CFG.tdmi_bins,
                      "tdmi_smooth_window": CFG.tdmi_smooth_window}, "curves": {}, "tau_core_seeds_N512": {}}
    for system in ("logistic", "henon", "white_noise", "ar1", "sinusoid"):
        for N in (512, 4096):
            x = S.generate(system, N, 0)
            mi = fp.time_delayed_mutual_information(x, CFG.tdmi_max_tau, CFG.tdmi_bins)
            rng = np.random.default_rng(0)
            floor = np.mean([fp.time_delayed_mutual_information(rng.permutation(x), 5, CFG.tdmi_bins).mean()
                             for _ in range(5)])
            out["curves"][f"{system}/N={N}"] = {
                "tau_selected": int(fp.find_optimal_tau(mi, CFG.tdmi_smooth_window)),
                "shuffled_floor": _f(floor), "mi_1_12": [_f(v) for v in mi[:12]],
                "smoothed_1_12": [_f(v) for v in np.convolve(mi, np.ones(3) / 3, mode="same")[:12]]}
    for system in ("logistic", "henon"):
        out["tau_core_seeds_N512"][system] = sorted(
            int(fp.find_optimal_tau(fp.time_delayed_mutual_information(S.generate(system, 512, s),
                                                                        CFG.tdmi_max_tau, CFG.tdmi_bins),
                                    CFG.tdmi_smooth_window)) for s in range(30))
    return out


# -----------------------------------------------------------------------------
# D3 -- Level C geometry at the period-1 peak nearest the analytic fixed point
# -----------------------------------------------------------------------------

def _map_next(system, q):
    return 4 * q[0] * (1 - q[0]) if system == "logistic" else 1 - 1.4 * q[0] ** 2 + 0.3 * q[1]


def _fit_geometry(X, q, k, step=1):
    n, d = X.shape
    dist = np.linalg.norm(X[:n - step] - q, axis=1)
    idx = np.argsort(dist, kind="stable")[:max(k, d + 2)]
    A = np.column_stack([np.ones(len(idx)), X[idx] - q])
    sv = np.linalg.svd(A, compute_uv=False)
    C = X[idx] - X[idx].mean(0)
    _, _, Vt = np.linalg.svd(C, full_matrices=False)
    P = Vt[:min(2, d)]
    r = q - X[idx].mean(0)
    return {"rank": int(np.linalg.matrix_rank(A)), "cols": int(A.shape[1]), "cond": _f(sv[0] / sv[-1]),
            "n_neighbors": int(len(idx)), "nn_dist_median": _f(np.median(dist[idx])),
            "nn_dist_max": _f(dist[idx].max()), "q_off_local_2plane": _f(np.linalg.norm(r - P.T @ (P @ r)))}


def levelc_diagnostics():
    out = {"config": {"verify_neighbors": CFG.verify_neighbors, "verify_max_residual": CFG.verify_max_residual,
                      "verify_min_r2": CFG.verify_min_r2,
                      "verify_min_support_ratio": CFG.verify_min_support_ratio},
           "rows": []}
    for system in ("henon", "logistic"):
        for seed in range(3):
            for m_forced in (None, 2 if system == "henon" else 1, 3):
                x = S.generate(system, 512, seed)
                m = (fp.cao_method(x, 1, CFG.cao_max_dim, CFG.cao_tol, CFG.cao_theiler)[0]
                     if m_forced is None else m_forced)
                X = fp._embed_backward(x, 1, m)
                r = fp.detect_so_fixed_points(X, 1, CFG, np.random.default_rng(CFG.random_seed))
                fp.apply_verification_gates(r, X, CFG)
                sd = float(np.std(X[:, 0]))
                xs = TRUE_FP[system]
                for c in (r[fp.LEVEL_C_PASS_FIELD] or []) + (r[fp.LEVEL_C_FAIL_FIELD] or []):
                    q = c["orbit_points"][0]
                    if abs(q[0] - xs) > 0.15:
                        continue
                    ver, st = c["verification"], c["extension_stability"]
                    model = fp._local_affine_model(X, q, CFG.verify_neighbors, 1)
                    model_t = fp._local_affine_model(X, np.full(m, xs), CFG.verify_neighbors, 1)
                    out["rows"].append({
                        "system": system, "seed": seed, "m": int(m), "m_source": "Cao lag 1" if m_forced is None
                        else "forced (diagnostic)", "true_fixed_point": xs, "candidate": _f(q[0]),
                        "candidate_error": _f(q[0] - xs), "pipeline_residual": _f(ver["residual"]),
                        "r2": _f(ver["r2"]), "support_ratio": _f(ver["support_ratio"]),
                        "failed_gates": ver["failed_gates"], "fhat_at_candidate": _f(model["F_q"][0]),
                        "exact_map_at_candidate": _f(_map_next(system, q)),
                        "exact_map_residual_at_candidate": _f(abs(_map_next(system, q) - q[0]) / sd),
                        "fhat_residual_at_true_fixed_point": _f(abs(model_t["F_q"][0] - xs) / sd),
                        "extension_multiplier_moduli": None if st is None else [_f(v) for v in st["multiplier_moduli"]],
                        **_fit_geometry(X, q, CFG.verify_neighbors)})
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for name, fn in (("tdmi", tdmi_diagnostics), ("levelc", levelc_diagnostics), ("lle", lle_diagnostics)):
            (OUT / f"{name}.json").write_text(json.dumps(fn(), indent=1))
            print("wrote", OUT / f"{name}.json", flush=True)


if __name__ == "__main__":
    main()
