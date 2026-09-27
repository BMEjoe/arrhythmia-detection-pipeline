"""
Phase 2E -- synthetic end-to-end validation of the FROZEN So et al. UPO pipeline.

Reproduce everything (from the repository root):

    python -m experiments.phase2e.run_phase2e --experiments all --workers 2
    python -m experiments.phase2e.analysis

Raw per-window results: experiments/phase2e/results/<experiment>.jsonl
(one JSON object per task: task spec, result row, per-peak detail) and a
flat <experiment>.csv.  No window is ever dropped; failures keep their
status and NaN fields.  Nothing in final_pipeline.py is modified; the only
runtime substitutions (stability disabled, analytic Jacobian) are scoped
diagnostic context managers in the "stability" / "henon_loc_oracleJ" arms.
"""
from __future__ import annotations

import os

# Single-threaded BLAS: bitwise reproducibility across runs and workers.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import contextlib
import copy
import datetime as _dt
import hashlib
import json
import multiprocessing as mp
import pathlib
import platform
import subprocess
import sys
import time
import warnings
from dataclasses import replace

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import final_pipeline as fp  # noqa: E402
from experiments.phase2e import config as C  # noqa: E402
from experiments.phase2e import metrics as M  # noqa: E402
from experiments.phase2e import systems as S  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
RESULTS = HERE / "results"


# =============================================================================
# Task construction
# =============================================================================

def _n_seeds(spec, n):
    return spec[n] if isinstance(spec, dict) else spec


def _window_task(experiment, system, n, seed, map_mode="one_sample",
                 surrogates=C.MAIN_SURROGATES, snr_db=None, algo_seed=None, **extra):
    t = {"experiment": experiment, "kind": "window", "system": system,
         "window_length": int(n), "seed": int(seed), "map_mode": map_mode,
         "surrogates": int(surrogates), "snr_db": snr_db, "algo_seed": algo_seed}
    t.update(extra)
    return t


def tasks_core():
    out = []
    for n in C.WINDOW_LENGTHS:
        for system in C.CORE_SYSTEMS:
            for s in range(_n_seeds(C.N_SEEDS_CORE[system], n)):
                out.append(_window_task("core", system, n, s))
    return out


def tasks_noise():
    return [_window_task("noise", system, C.PRODUCTION_WINDOW, s, snr_db=snr)
            for system in C.NOISE_SYSTEMS for snr in C.SNR_LEVELS_DB
            for s in range(C.N_SEEDS_NOISE)]


def tasks_tau_step():
    out = []
    for system in C.TAU_STEP_SYSTEMS:
        k = min(C.N_SEEDS_TAU_STEP, _n_seeds(C.N_SEEDS_CORE[system], C.PRODUCTION_WINDOW))
        out += [_window_task("tau_step", system, C.PRODUCTION_WINDOW, s, map_mode="tau_step")
                for s in range(k)]
    return out


def tasks_surrogate_sens():
    return [_window_task("surrogate_sens", system, C.PRODUCTION_WINDOW, s, surrogates=ns)
            for ns in C.SURROGATE_SENSITIVITY for system in C.SURROGATE_SENS_SYSTEMS
            for s in range(C.N_SEEDS_SURROGATE_SENS[system])]


def tasks_algo_seed():
    return [_window_task("algo_seed", system, C.PRODUCTION_WINDOW, s,
                         algo_seed=fp.CFG.random_seed + C.ALGO_SEED_OFFSET + s)
            for system in C.ALGO_SEED_SYSTEMS for s in range(C.N_SEEDS_ALGO)]


def tasks_oracle_dim():
    return [{"experiment": "oracle_dim", "kind": "oracle", "system": system,
             "window_length": C.PRODUCTION_WINDOW, "seed": s, "dimension": d,
             "surrogates": C.MAIN_SURROGATES}
            for system, d in C.ORACLE_DIMENSION.items() for s in range(C.N_SEEDS_ORACLE)]


def tasks_fp_sens():
    return [{"experiment": "fp_sens", "kind": "fp", "system": system,
             "window_length": C.PRODUCTION_WINDOW, "seed": s, "variant": v, "surrogates": 0}
            for system in C.FP_SYSTEMS for s in range(C.N_SEEDS_FP)
            for v in ("original", "nextafter", "rr_route_a", "rr_route_b")]


def _loc_pad(policy):
    return {"pad0": 0.0, "pad10": 0.10, "pad25": 0.25, "pad10_shift": 0.10}[policy]


def tasks_henon_loc():
    out = []
    prod = (30, 10.0, "pad10")
    for system in C.LOC_SYSTEMS:
        for real in C.LOC_REALIZATIONS:
            for bins in C.LOC_BINS:
                for tube in C.LOC_TUBE_PERCENTILES:
                    for pol in C.LOC_RANGE_POLICIES:
                        on_sweep = sum([bins != prod[0], tube != prod[1], pol != prod[2]]) <= 1
                        for period in (1, 2):
                            out.append({"experiment": "henon_loc", "kind": "loc",
                                        "system": system, "realization": real, "seed": real,
                                        "window_length": C.LOC_N, "bins": bins, "tube": tube,
                                        "range_policy": pol, "period": period,
                                        "M": 2, "jacobian": "estimated",
                                        "significance": bool(period == 1 and on_sweep
                                                             and real in C.LOC_SIG_REALIZATIONS)})
            # production-Jacobian arm at the production point
            for period in (1, 2):
                out.append({"experiment": "henon_loc", "kind": "loc", "system": system,
                            "realization": real, "seed": real, "window_length": C.LOC_N,
                            "bins": 30, "tube": 10.0, "range_policy": "pad10", "period": period,
                            "M": C.LOC_PRODUCTION_M, "jacobian": "estimated",
                            "significance": period == 1})
            # analytic-Jacobian diagnostic along the bins sweep
            for bins in C.LOC_BINS:
                for period in (1, 2):
                    out.append({"experiment": "henon_loc_oracleJ", "kind": "loc",
                                "system": system, "realization": real, "seed": real,
                                "window_length": C.LOC_N, "bins": bins, "tube": 10.0,
                                "range_policy": "pad10", "period": period, "M": 2,
                                "jacobian": "analytic", "significance": False})
    return out


def tasks_stability():
    out = [{"experiment": "stability", "kind": "stability", "system": system,
            "realization": real, "seed": real, "window_length": C.LOC_N, "M": Mn}
           for system in C.LOC_SYSTEMS for real in C.LOC_REALIZATIONS for Mn in (2, C.LOC_PRODUCTION_M)]
    out += [{"experiment": "stability", "kind": "stability_pipeline", "system": system,
             "window_length": C.PRODUCTION_WINDOW, "seed": s, "surrogates": C.MAIN_SURROGATES}
            for system in ("henon", "logistic", "white_noise") for s in range(3)]
    return out


def tasks_flag_independence():
    return [{"experiment": "flag_independence", "kind": "flag", "system": system,
             "window_length": C.PRODUCTION_WINDOW, "seed": s, "surrogates": C.MAIN_SURROGATES}
            for system in ("henon", "logistic", "white_noise") for s in range(2)]


EXPERIMENTS = {
    "core": tasks_core, "noise": tasks_noise, "tau_step": tasks_tau_step,
    "surrogate_sens": tasks_surrogate_sens, "algo_seed": tasks_algo_seed,
    "oracle_dim": tasks_oracle_dim, "fp_sens": tasks_fp_sens,
    "henon_loc": tasks_henon_loc, "stability": tasks_stability,
    "flag_independence": tasks_flag_independence,
}


def all_tasks(names):
    """Task ids are numbered within each experiment, independent of selection."""
    out = []
    for name in names:
        counters = {}
        for t in EXPERIMENTS[name]():   # henon_loc also yields the henon_loc_oracleJ arm
            k = counters.get(t["experiment"], 0)
            counters[t["experiment"]] = k + 1
            t["task_id"] = f"{t['experiment']}:{k:05d}"
            out.append(t)
    return out


# =============================================================================
# Execution
# =============================================================================

def pipeline_config(task):
    ns = int(task.get("surrogates", 0))
    kw = {"upo_map_mode": task.get("map_mode", "one_sample"),
          "so_assess_significance": ns > 0}
    if ns > 0:
        kw["so_surrogate_count"] = ns
    if task.get("algo_seed") is not None:
        kw["random_seed"] = int(task["algo_seed"])
    return replace(fp.CFG, **kw)


def _base_row(task, x=None):
    row = {k: task.get(k) for k in ("task_id", "experiment", "system", "seed", "window_length",
                                    "snr_db", "surrogates", "variant")}
    row["data_seed"] = S.data_seed(task["system"], task["window_length"], task["seed"])
    row["algo_random_seed"] = task.get("algo_seed") or fp.CFG.random_seed
    if x is not None:
        row["signal_mean"] = float(np.mean(x))
        row["signal_std"] = float(np.std(x))
    return row


def _run_segment(x, cfg):
    try:
        res = fp.analyze_segment(x, config=cfg)
    except (ValueError, np.linalg.LinAlgError) as exc:
        try:
            direct = fp.run_upo_analysis(x, config=cfg)
        except Exception as exc2:  # diagnostic only
            direct = {"status": f"raised {type(exc2).__name__}", "embedding_dimension": None}
        return M.error_fields(exc, direct, cfg), {"peaks": []}, None
    row, detail = M.segment_fields(res, cfg)
    return row, detail, res


def run_window(task):
    x = S.generate(task["system"], task["window_length"], task["seed"], task.get("snr_db"))
    cfg = pipeline_config(task)
    row = _base_row(task, x)
    upo_row, detail, _ = _run_segment(x, cfg)
    row.update(upo_row)
    if task.get("map_mode", "one_sample") == "one_sample" and task["system"] != "constant":
        m, e1, _ = fp.cao_method(x, 1, cfg.cao_max_dim, cfg.cao_tol, cfg.cao_theiler)
        detail["cao_lag1_E1"] = [None if not np.isfinite(v) else float(v) for v in e1]
    return row, detail


def run_oracle(task):
    """DIAGNOSTIC: run_upo_analysis steps at a fixed known dimension (one_sample, tau=1)."""
    x = S.generate(task["system"], task["window_length"], task["seed"])
    cfg = pipeline_config(task)
    d = int(task["dimension"])
    emb = fp._embed_backward(x, 1, d)
    r = fp.detect_so_fixed_points(emb, tau=1, config=cfg, rng=np.random.default_rng(cfg.random_seed))
    fp.assess_so_significance(r, x, config=cfg, rng=np.random.default_rng(cfg.random_seed + 104729),
                              embedding_dimension=d)
    if r["status"] not in fp.UPO_FAILURE_STATUSES:
        fp.apply_verification_gates(r, emb, config=cfg)
    # Assemble the same summary structure run_upo_analysis would.
    upo = fp._empty_upo_analysis(fp.resolve_upo_map("one_sample", 1), cfg, r["status"], 1, d,
                                 n_points=len(emb), embedded=emb)
    upo["periods"] = {1: r}
    upo[fp.LEVEL_A_FIELD] = r[fp.LEVEL_A_FIELD]
    if r["status"] not in fp.UPO_FAILURE_STATUSES:
        upo["significance_assessed"] = r["significance_assessed"]
        upo["significance_status"] = (fp.SIGNIFICANCE_ASSESSED if r["significance_assessed"]
                                      else fp.SIGNIFICANCE_NOT_ASSESSED)
        if r["significance_assessed"]:
            upo[fp.LEVEL_B_FIELD] = r[fp.LEVEL_B_FIELD]
            upo["source_rJ"] = r["significance"]["rJ"]
        upo["verification_assessed"] = r["verification_assessed"]
        upo[fp.LEVEL_C_PASS_FIELD] = r[fp.LEVEL_C_PASS_FIELD]
        upo[fp.LEVEL_C_FAIL_FIELD] = r[fp.LEVEL_C_FAIL_FIELD]
        upo["verified_unstable"] = [c for c in r[fp.LEVEL_C_PASS_FIELD] or [] if c["verified_unstable"]]
        seed = cfg.random_seed
        upo["coverage"] = {
            "source": fp.peak_coverage(fp._population_points(upo[fp.LEVEL_A_FIELD], d), emb, cfg, seed),
            "significant": (fp.peak_coverage(fp._population_points(upo[fp.LEVEL_B_FIELD], d), emb, cfg, seed)
                            if upo[fp.LEVEL_B_FIELD] is not None else None),
            "verified_unstable": fp.peak_coverage(fp._population_points(upo["verified_unstable"], d),
                                                  emb, cfg, seed),
        }
    row = _base_row(task, x)
    row["oracle_dimension"] = d
    upo_row, detail = M.upo_fields(upo, cfg)
    row.update(upo_row)
    return row, detail


def run_fp(task):
    x = S.generate(task["system"], task["window_length"], task["seed"])
    v = task["variant"]
    if v == "nextafter":
        x = np.nextafter(x, np.inf)
    elif v.startswith("rr_route"):
        rr = 0.8 + 0.05 * (x - np.mean(x)) / np.std(x)
        x = (np.round(rr * C.RR_FS) / C.RR_FS if v == "rr_route_a"
             else np.round(rr * C.RR_FS) * (1.0 / C.RR_FS))
    cfg = pipeline_config(task)
    row = _base_row(task, x)
    row["variant"] = v
    upo_row, detail, _ = _run_segment(x, cfg)
    row.update(upo_row)
    row["n_distinct_values"] = int(len(np.unique(x)))
    return row, detail


# -------------------------------------------------------------------------
# Henon localization
# -------------------------------------------------------------------------

def loc_config(task):
    return replace(fp.CFG, so_jacobian_neighbors=int(task["M"]), so_neighbors_K=1, so_random_R=100,
                   so_random_R_period_p=100, so_kappa=3.0, so_max_backbones=None,
                   so_randomization="prl_norm", upo_map_mode="one_sample",
                   so_hist_bins=int(task["bins"]), so_diagonal_percentile=float(task["tube"]),
                   so_slab_percentile=float(task["tube"]),
                   so_hist_range_pad=_loc_pad(task["range_policy"]),
                   so_surrogate_count=C.MAIN_SURROGATES)


@contextlib.contextmanager
def analytic_jacobians(system):
    """DIAGNOSTIC: replace the estimated local Jacobians by the exact map Jacobian."""
    jac = S.REFERENCES[system]["jacobian"]
    original = fp._local_jacobians

    def exact(X, step, Mn, exclude):
        n = len(X)
        return [jac(X[i, 0]) if i < n - step else None for i in range(n)]

    fp._local_jacobians = exact
    try:
        yield
    finally:
        fp._local_jacobians = original


def _orbit_distance(orbit, target):
    orbit = np.asarray(orbit, float)
    target = np.asarray(target, float)
    return min(float(np.max(np.abs(np.roll(orbit, k) - target))) for k in range(len(orbit)))


def _kde_mode(values, center, half_width, bandwidth):
    v = values[np.abs(values - center) <= half_width]
    if len(v) < 5:
        return float("nan")
    grid = np.linspace(center - half_width, center + half_width, 801)
    dens = np.exp(-0.5 * ((grid[:, None] - v[None, :]) / bandwidth) ** 2).sum(axis=1)
    return float(grid[int(np.argmax(dens))])


def run_loc(task):
    system = task["system"]
    x = S.generate(system, task["window_length"], task["seed"])
    X = fp._embed_backward(x, 1, 2)
    cfg = loc_config(task)
    ref = S.REFERENCES[system]
    hist_range = None
    if task["range_policy"] == "pad10_shift":
        lo, hi = fp._so_histogram_range(X, 0.10)
        shift = 0.5 * (hi - lo) / cfg.so_hist_bins
        hist_range = (lo + shift, hi + shift)
    ctx = analytic_jacobians(system) if task["jacobian"] == "analytic" else contextlib.nullcontext()
    period = int(task["period"])
    with ctx:
        if period == 1:
            r = fp.detect_so_fixed_points(X, 1, cfg, np.random.default_rng(1), hist_range=hist_range)
        else:
            r = fp.detect_so_period_p(X, 2, 1, cfg, np.random.default_rng(2), hist_range=hist_range)
        if task["significance"]:
            fp.assess_so_significance(r, x, cfg, np.random.default_rng(cfg.random_seed + 104729),
                                      embedding_dimension=2)
    rng_range = r.get("histogram_range") or (np.nan, np.nan)
    bw = (rng_range[1] - rng_range[0]) / cfg.so_hist_bins
    row = {k: task[k] for k in ("task_id", "experiment", "system", "realization", "seed",
                                "window_length", "bins", "tube", "range_policy", "period", "M",
                                "jacobian", "significance")}
    row.update({"data_seed": S.data_seed(system, task["window_length"], task["seed"]),
                "status": r["status"], "bin_width": float(bw),
                "source_peak_count": len(r[fp.LEVEL_A_FIELD]),
                "significance_assessed": bool(r["significance_assessed"]),
                "significant_peak_count": (len(r[fp.LEVEL_B_FIELD]) if r[fp.LEVEL_B_FIELD] is not None
                                           else float("nan")),
                "rJ": (float(r["significance"]["rJ"]) if r["significance_assessed"] else float("nan")),
                "J_W": (float(r["significance"]["J_W"]) if r["significance_assessed"] else float("nan"))})
    peaks = r[fp.LEVEL_A_FIELD]
    sig_info = (r["significance"] or {}).get("per_peak") if r["significance_assessed"] else None
    detail = {"peaks": []}
    for i, c in enumerate(peaks):
        detail["peaks"].append({"orbit": [float(v) for v in np.atleast_1d(c["orbit_coordinates"])],
                                "minimal_period": int(c.get("minimal_period", period)),
                                "count": int(c["histogram_count"]),
                                "J": None if sig_info is None else float(sig_info[i]["J"])})
    if period == 1:
        locs = np.array([c["scalar_location"] for c in peaks])
        counts = np.array([c["histogram_count"] for c in peaks])
        scal = np.asarray(r.get("scalar", np.empty(0)))
        for j, fx in enumerate(ref["fixed_points"]):
            tag = f"fp{j + 1}"
            row[f"{tag}_true"] = float(fx)
            if len(locs):
                k = int(np.argmin(np.abs(locs - fx)))
                row[f"{tag}_nearest"] = float(locs[k])
                row[f"{tag}_error"] = float(locs[k] - fx)
                row[f"{tag}_abs_error"] = float(abs(locs[k] - fx))
                row[f"{tag}_rank"] = int(np.sum(counts > counts[k]) + 1)
                row[f"{tag}_within_bin"] = bool(abs(locs[k] - fx) <= bw)
                row[f"{tag}_significant"] = (None if sig_info is None else bool(sig_info[k]["significant"]))
                row[f"{tag}_J"] = (None if sig_info is None else float(sig_info[k]["J"]))
                # Histogram-free estimates from the tube scalars around the detected peak
                cell_lo = r["edges"][0][peaks[k]["peak_cell"][0]]
                row[f"{tag}_cell_center_error"] = float(cell_lo + 0.5 * bw - fx)
                row[f"{tag}_kde_mode_error"] = _kde_mode(scal, locs[k], bw, 0.002) - fx
                near = scal[np.abs(scal - fx) <= 0.05]
                row[f"{tag}_median_within_0.05_error"] = (float(np.median(near) - fx) if len(near)
                                                          else float("nan"))
            else:
                for key in ("nearest", "error", "abs_error", "rank", "within_bin", "significant", "J",
                            "cell_center_error", "kde_mode_error", "median_within_0.05_error"):
                    row[f"{tag}_{key}"] = float("nan")
    else:
        target = ref["period2"]
        two = [c for c in peaks if c["minimal_period"] == 2]
        one = [c for c in peaks if c["minimal_period"] == 1]
        row["p2_true_a"], row["p2_true_b"] = float(target[0]), float(target[1])
        row["n_minimal_period2"] = len(two)
        row["n_minimal_period1"] = len(one)
        if two:
            d = [_orbit_distance(c["orbit_coordinates"], target) for c in two]
            strongest = max(range(len(two)), key=lambda i: two[i]["histogram_count"])
            row["p2_nearest_error"] = float(min(d))
            row["p2_strongest_error"] = float(d[strongest])
            row["p2_nearest_rank"] = int(sorted(d).index(min(d)) + 1)
            row["p2_within_bin"] = bool(min(d) <= bw)
        else:
            row.update({"p2_nearest_error": float("nan"), "p2_strongest_error": float("nan"),
                        "p2_nearest_rank": float("nan"), "p2_within_bin": False})
        fx = ref["fixed_points"][0]
        row["p2_fp1_error"] = (float(min(abs(np.mean(c["orbit_coordinates"]) - fx) for c in one))
                               if one else float("nan"))
    return row, detail


# -------------------------------------------------------------------------
# Stability (Step 17) and flag independence (Step 21)
# -------------------------------------------------------------------------

@contextlib.contextmanager
def stability_disabled():
    original = fp.source_period1_stability
    fp.source_period1_stability = lambda jacobians, member_indices: None
    try:
        yield
    finally:
        fp.source_period1_stability = original


def _peak_signature(r):
    return [(tuple(c["peak_cell"]), int(c["histogram_count"]),
             float(np.round(c["scalar_location"], 12)), int(c["transformed_points_in_peak"]))
            for c in r[fp.LEVEL_A_FIELD]]


def run_stability(task):
    system = task["system"]
    x = S.generate(system, task["window_length"], task["seed"])
    X = fp._embed_backward(x, 1, 2)
    cfg = loc_config({"M": task["M"], "bins": 30, "tube": 10.0, "range_policy": "pad10"})
    r_on = fp.detect_so_fixed_points(X, 1, cfg, np.random.default_rng(1))
    with stability_disabled():
        r_off = fp.detect_so_fixed_points(X, 1, cfg, np.random.default_rng(1))
    r2 = fp.detect_so_period_p(X, 2, 1, cfg, np.random.default_rng(2))
    fp.apply_verification_gates(r_on, X, cfg)
    ref = S.REFERENCES[system]
    row = {"task_id": task["task_id"], "experiment": "stability", "kind": "detector",
           "system": system, "realization": task["realization"], "seed": task["seed"],
           "window_length": task["window_length"], "M": task["M"],
           "data_seed": S.data_seed(system, task["window_length"], task["seed"]),
           "status_with_stability": r_on["status"], "status_without_stability": r_off["status"],
           "peaks_identical_without_stability": _peak_signature(r_on) == _peak_signature(r_off),
           "all_off_stability_none": all(c["source_stability"] is None for c in r_off[fp.LEVEL_A_FIELD]),
           "period2_source_stability_all_none": all(c["source_stability"] is None
                                                    for c in r2[fp.LEVEL_A_FIELD]),
           "n_peaks": len(r_on[fp.LEVEL_A_FIELD])}
    ext = {tuple(c["peak_cell"]): c.get("extension_stability")
           for c in (r_on[fp.LEVEL_C_PASS_FIELD] or []) + (r_on[fp.LEVEL_C_FAIL_FIELD] or [])}
    for j, fx in enumerate(ref["fixed_points"]):
        tag = f"fp{j + 1}"
        true = np.sort(np.abs(np.linalg.eigvals(ref["jacobian"](fx))))[::-1]
        row[f"{tag}_analytic_moduli"] = [float(v) for v in true]
        peaks = r_on[fp.LEVEL_A_FIELD]
        if not peaks:
            continue
        k = int(np.argmin([abs(c["scalar_location"] - fx) for c in peaks]))
        c = peaks[k]
        row[f"{tag}_peak_abs_error"] = float(abs(c["scalar_location"] - fx))
        st = c["source_stability"]
        row[f"{tag}_hybrid_moduli"] = None if st is None else [float(v) for v in st["lyapunov_numbers"]]
        row[f"{tag}_hybrid_unstable"] = None if st is None else bool(st["unstable"])
        row[f"{tag}_hybrid_n_points"] = None if st is None else int(st["n_points_averaged"])
        e = ext.get(tuple(c["peak_cell"]))
        row[f"{tag}_extension_moduli"] = None if e is None else [float(v) for v in e["multiplier_moduli"]]
        # Analytic Jacobian at the DETECTED location (separates location error from J error)
        at_loc = np.sort(np.abs(np.linalg.eigvals(ref["jacobian"](c["scalar_location"]))))[::-1]
        row[f"{tag}_analytic_moduli_at_detected"] = [float(v) for v in at_loc]
    return row, {}


def run_stability_pipeline(task):
    """Full production analysis with and without the stability diagnostic."""
    x = S.generate(task["system"], task["window_length"], task["seed"])
    cfg = pipeline_config(task)
    row_on, det_on, _ = _run_segment(x, cfg)
    with stability_disabled():
        row_off, det_off, _ = _run_segment(x, cfg)
    ignore = {"source_stability_available_count", "source_stability_unstable_count"}
    diff = sorted(k for k in row_on if k not in ignore and
                  json.dumps(row_on[k], default=str) != json.dumps(row_off.get(k), default=str))
    strip = lambda d: [{k: v for k, v in p.items() if not k.startswith("source_stability")}  # noqa: E731
                       for p in d["peaks"]]
    row = {"task_id": task["task_id"], "experiment": "stability", "kind": "pipeline",
           "system": task["system"], "seed": task["seed"], "window_length": task["window_length"],
           "data_seed": S.data_seed(task["system"], task["window_length"], task["seed"]),
           "upo_status": row_on["upo_status"], "significance_status": row_on["upo_significance_status"],
           "source_peak_count": row_on["source_peak_count"],
           "significant_peak_count": row_on["significant_peak_count"],
           "fields_differing_without_stability": diff,
           "peaks_identical_without_stability": strip(det_on) == strip(det_off)}
    return row, {}


def run_flag(task):
    x = S.generate(task["system"], task["window_length"], task["seed"])
    base = pipeline_config(task)
    rows = {}
    for flag in (False, True):
        cfg = replace(base, compute_surrogates=flag)
        res = fp.analyze_segment(x, config=cfg)
        r, d = M.segment_fields(res, cfg)
        rows[flag] = (r, d, res.get("surrogate_analysis"))
    keys = ["upo_status", "upo_significance_status", "significance_assessed", "W", "W0", "J_W",
            "rJ", "source_peak_count", "significant_peak_count", "lle_per_beat"]
    same = {k: json.dumps(rows[False][0][k]) == json.dumps(rows[True][0][k]) for k in keys}
    row = {"task_id": task["task_id"], "experiment": "flag_independence", "system": task["system"],
           "seed": task["seed"], "window_length": task["window_length"],
           "data_seed": S.data_seed(task["system"], task["window_length"], task["seed"]),
           "all_upo_fields_identical": all(same.values()),
           "surrogate_W_identical": rows[False][1].get("surrogate_W") == rows[True][1].get("surrogate_W"),
           "lle_surrogates_ran_when_flag_true": rows[True][2] is not None,
           "lle_surrogates_absent_when_flag_false": rows[False][2] is None,
           **{f"same_{k}": v for k, v in same.items()}}
    return row, {}


RUNNERS = {"window": run_window, "oracle": run_oracle, "fp": run_fp, "loc": run_loc,
           "stability": run_stability, "stability_pipeline": run_stability_pipeline,
           "flag": run_flag}


def execute(task):
    t0 = time.perf_counter()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        try:
            row, detail = RUNNERS[task["kind"]](task)
        except Exception as exc:  # never drop a task: record the failure
            row, detail = {"task_id": task["task_id"], "experiment": task["experiment"],
                           "system": task.get("system"), "seed": task.get("seed"),
                           "harness_error": f"{type(exc).__name__}: {exc}"}, {}
    row["n_runtime_warnings"] = len(w)
    runtime = time.perf_counter() - t0
    return {"task": task, "row": row, "detail": detail, "runtime_s": runtime}


# =============================================================================
# I/O
# =============================================================================

def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, tuple):
        return list(o)
    raise TypeError(type(o))


def dumps(obj):
    return json.dumps(obj, sort_keys=True, default=_json_default)


def git_state():
    def g(*a):
        return subprocess.run(["git", *a], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    return {"head": g("rev-parse", "HEAD"), "status_porcelain": g("status", "--porcelain")}


def sha256(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def write_csv(experiment, records):
    import pandas as pd
    rows = []
    for rec in sorted(records, key=lambda r: r["task"]["task_id"]):
        row = dict(rec["row"])
        for k, v in list(row.items()):
            if isinstance(v, (list, dict, tuple)):
                row[k] = dumps(v)
        row["runtime_s"] = rec["runtime_s"]
        rows.append(row)
    pd.DataFrame(rows).to_csv(RESULTS / f"{experiment}.csv", index=False)


def load_records(experiment):
    path = RESULTS / f"{experiment}.jsonl"
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def run(names, workers, resume=True):
    RESULTS.mkdir(parents=True, exist_ok=True)
    tasks = all_tasks(names)
    by_exp = {}
    for t in tasks:
        by_exp.setdefault(t["experiment"], []).append(t)
    done = {}
    for exp in by_exp:
        recs = load_records(exp) if resume else []
        if not resume and (RESULTS / f"{exp}.jsonl").exists():
            (RESULTS / f"{exp}.jsonl").unlink()
        done.update({r["task"]["task_id"]: r for r in recs})
    todo = [t for t in tasks if t["task_id"] not in done]
    # expensive first for load balance
    cost = lambda t: (t.get("window_length", 256) ** 2) * (1 + t.get("surrogates", 0) / 10)  # noqa: E731
    todo.sort(key=cost, reverse=True)
    print(f"{len(tasks)} tasks, {len(done)} already done, {len(todo)} to run, {workers} workers",
          flush=True)
    handles = {exp: open(RESULTS / f"{exp}.jsonl", "a") for exp in by_exp}
    t0 = time.time()
    ctx = mp.get_context("fork")
    with ctx.Pool(workers, maxtasksperchild=50) as pool:
        for i, rec in enumerate(pool.imap_unordered(execute, todo, chunksize=1), 1):
            handles[rec["task"]["experiment"]].write(dumps(rec) + "\n")
            handles[rec["task"]["experiment"]].flush()
            done[rec["task"]["task_id"]] = rec
            if i % 25 == 0 or i == len(todo):
                print(f"  {i}/{len(todo)} done, {time.time() - t0:.0f}s elapsed", flush=True)
    for h in handles.values():
        h.close()
    for exp in by_exp:
        write_csv(exp, load_records(exp))


def replicate(names, per_experiment, workers):
    """Step 25: rerun the first k tasks of every experiment and compare exactly."""
    tasks = all_tasks(names)
    chosen, seen = [], {}
    for t in tasks:
        k = seen.get(t["experiment"], 0)
        if k < per_experiment:
            chosen.append(t)
            seen[t["experiment"]] = k + 1
    stored = {}
    for exp in seen:
        stored.update({r["task"]["task_id"]: r for r in load_records(exp)})
    ctx = mp.get_context("fork")
    with ctx.Pool(workers) as pool:
        reruns = pool.map(execute, chosen, chunksize=1)
    out = []
    for rec in reruns:
        tid = rec["task"]["task_id"]
        old = stored.get(tid)
        strip = lambda r: {k: v for k, v in r.items() if k != "n_runtime_warnings"}  # noqa: E731
        out.append({"task_id": tid, "experiment": rec["task"]["experiment"],
                    "stored_found": old is not None,
                    "row_identical": old is not None and dumps(strip(old["row"])) == dumps(strip(rec["row"])),
                    "detail_identical": old is not None and dumps(old["detail"]) == dumps(rec["detail"])})
    (RESULTS / "replicability.json").write_text(json.dumps(out, indent=1))
    n_ok = sum(o["row_identical"] and o["detail_identical"] for o in out)
    print(f"replicability: {n_ok}/{len(out)} tasks bitwise identical")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--experiments", default="all",
                    help="comma list of " + ",".join(EXPERIMENTS) + " or 'all'")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--no-resume", action="store_true")
    ap.add_argument("--replicate", action="store_true", help="run the Step 25 replicability check")
    args = ap.parse_args(argv)
    names = list(EXPERIMENTS) if args.experiments == "all" else args.experiments.split(",")
    RESULTS.mkdir(parents=True, exist_ok=True)
    if args.replicate:
        replicate(names, C.REPLICATE_PER_EXPERIMENT, args.workers)
        return
    manifest_path = RESULTS / "manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {"runs": []}
    manifest["runs"].append({
        "started_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "experiments": names, "workers": args.workers,
        "git": git_state(),
        "sha256": {"final_pipeline.py": sha256(ROOT / "final_pipeline.py"),
                   "config.py": sha256(HERE / "config.py"), "systems.py": sha256(HERE / "systems.py"),
                   "metrics.py": sha256(HERE / "metrics.py"),
                   "run_phase2e.py": sha256(HERE / "run_phase2e.py")},
        "python": sys.version, "platform": platform.platform(),
        "numpy": np.__version__, "blas_threads": os.environ.get("OMP_NUM_THREADS"),
        "n_tasks": len(all_tasks(names)),
    })
    manifest_path.write_text(json.dumps(manifest, indent=1))
    run(names, args.workers, resume=not args.no_resume)
    manifest = json.loads(manifest_path.read_text())
    manifest["runs"][-1]["finished_utc"] = _dt.datetime.now(_dt.timezone.utc).isoformat()
    manifest_path.write_text(json.dumps(manifest, indent=1))


if __name__ == "__main__":
    main()
