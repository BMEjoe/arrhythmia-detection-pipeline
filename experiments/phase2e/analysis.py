"""
Phase 2E analysis: tables A-I, audits, and plots from the raw results.

    python -m experiments.phase2e.analysis

Reads results/*.jsonl (raw) and writes results/tables/*.csv|.md,
results/audits.json and plots/*.png.  Every table reports N, valid N,
NaN N and failure N; nothing is dropped before counting.
"""
from __future__ import annotations

import json
import math
import pathlib
import sys

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import final_pipeline as fp  # noqa: E402
from experiments.phase2e import config as C  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
RESULTS = HERE / "results"
TABLES = RESULTS / "tables"
PLOTS = HERE / "plots"

# Validated reference categorical palette (dataviz skill, light mode), fixed order.
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
SYSTEM_ORDER = ["constant", "white_noise", "ar1", "sinusoid", "logistic_p4", "logistic", "henon"]
SYSTEM_COLOR = {s: PALETTE[i] for i, s in enumerate(SYSTEM_ORDER)}

# Step 18 -- predeclared interpretation of every status (A = evidence of no
# detectable structure, B = insufficient information / numerical failure).
STATUS_SEMANTICS = {
    "ok": ("detector completed; >= 1 Level-A peak", "neither", "not a failure"),
    "no_peaks": ("detector completed; histogram had no peak above threshold", "A (weak)",
                 "evidence of no detectable period-1 peak at this resolution, not absence of UPOs"),
    "tube_empty": ("no transformed point inside histogram range after tube", "B",
                   "pipeline counts it as 0 peaks (UPO_NO_PEAK_STATUSES) -- conflation"),
    "too_few_in_tube": ("< so_peak_min_count points in tube", "B",
                        "pipeline counts it as 0 peaks (UPO_NO_PEAK_STATUSES) -- conflation"),
    "constant_data": ("input has zero range", "B", "degenerate input; NaN features"),
    "nonfinite_input": ("NaN/inf in input", "B", "NaN features"),
    "too_few_points": ("embedded points below detector minimum", "B", "NaN features"),
    "embedding_not_saturated": ("Cao found no E1 plateau up to cao_max_dim (or E1 undefined)", "B",
                                "also triggered by exact duplicate vectors, where E1 is undefined"),
    "no_valid_transforms": ("every So transform singular / no Jacobian", "B", "NaN features"),
    "analysis_error": ("analyze_segment raised (production wrapper)", "B",
                       "masks the underlying UPO status; NaN features"),
    "significance_not_assessed": ("significance requested, < 2 usable surrogates", "B",
                                  "Level B None, features NaN"),
    "not_assessed_detector_failure": ("significance requested, detector failed", "B", "NaN"),
    "assessed": ("surrogate test completed (Level B is a list, possibly empty)", "neither", ""),
    "not_requested": ("significance not requested", "neither", ""),
}


# =============================================================================
# Loading
# =============================================================================

def load(experiment):
    path = RESULTS / f"{experiment}.jsonl"
    if not path.exists():
        return pd.DataFrame(), []
    recs = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    recs.sort(key=lambda r: r["task"]["task_id"])
    rows = []
    for r in recs:
        row = dict(r["row"])
        row["runtime_s"] = r["runtime_s"]
        rows.append(row)
    return pd.DataFrame(rows), recs


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (max(0.0, c - h), min(1.0, c + h))


def desc(s):
    """N / valid / NaN / mean / median / SD / IQR of a numeric series."""
    s = pd.to_numeric(s, errors="coerce")
    v = s.dropna()
    return {"N": int(len(s)), "valid_N": int(len(v)), "NaN_N": int(s.isna().sum()),
            "mean": float(v.mean()) if len(v) else np.nan,
            "median": float(v.median()) if len(v) else np.nan,
            "SD": float(v.std(ddof=1)) if len(v) > 1 else np.nan,
            "Q1": float(v.quantile(0.25)) if len(v) else np.nan,
            "Q3": float(v.quantile(0.75)) if len(v) else np.nan}


def fmt_desc(s, digits=3):
    d = desc(s)
    if d["valid_N"] == 0:
        return f"NA (valid 0/{d['N']})"
    return (f"{d['median']:.{digits}g} [{d['Q1']:.{digits}g}, {d['Q3']:.{digits}g}]; "
            f"mean {d['mean']:.{digits}g}; valid {d['valid_N']}/{d['N']}")


def frac(mask_num, denom):
    n = int(denom)
    k = int(mask_num)
    lo, hi = wilson(k, n)
    return f"{k}/{n} = {k / n:.3f} [{lo:.3f}, {hi:.3f}]" if n else "NA"


def counts_str(s):
    vc = s.fillna("NA").astype(str).value_counts()
    return "; ".join(f"{k}: {v}" for k, v in vc.items())


def save(df, name, title, note=""):
    TABLES.mkdir(parents=True, exist_ok=True)
    df.to_csv(TABLES / f"{name}.csv", index=False)
    with open(TABLES / f"{name}.md", "w") as fh:
        fh.write(f"### {title}\n\n")
        if note:
            fh.write(note.strip() + "\n\n")
        fh.write(df.to_markdown(index=False))
        fh.write("\n")
    return df


# =============================================================================
# Tables
# =============================================================================

def window_summary(g):
    n = len(g)
    fail = g["upo_failure"].fillna(True).astype(bool)
    ok = g[~fail]
    sig_assessed = g["significance_assessed"].fillna(False).astype(bool)
    return {
        "N": n,
        "failure_N": int(fail.sum()),
        "upo_status": counts_str(g["upo_status"]),
        "cao_dim (UPO)": fmt_desc(g["embedding_dimension"]),
        "n_embedded": fmt_desc(g["n_embedded_points"]),
        "LLE valid rate": (frac(g["lle_valid"].astype(bool).sum(), n) if "lle_valid" in g
                           else "not computed (detector-only arm)"),
        "LLE (valid)": fmt_desc(g["lle_per_beat"]) if "lle_per_beat" in g else "not computed (detector-only arm)",
        ">=1 Level-A peak": frac((ok["source_peak_count"] >= 1).sum(), n),
        "source_peak_count": fmt_desc(g["source_peak_count"]),
        "sig assessed N": int(sig_assessed.sum()),
        ">=1 Level-B peak": frac((g["significant_peak_count"] >= 1).sum(), n),
        "significant_peak_count": fmt_desc(g["significant_peak_count"]),
        "rJ": fmt_desc(g["source_rJ"]),
        "source_coverage": fmt_desc(g["source_peak_coverage"]),
        "verified_unstable_count": fmt_desc(g["verified_unstable_count"]),
    }


def table_A(core):
    rows = []
    for (system, n), g in core.groupby(["system", "window_length"], sort=False):
        rows.append({"system": system, "window": n, **window_summary(g)})
    df = pd.DataFrame(rows)
    df["_o"] = df["system"].map(SYSTEM_ORDER.index)
    df = df.sort_values(["_o", "window"]).drop(columns="_o")
    return save(df, "table_A_core", "Table A -- core system x window length (one_sample, 50 surrogates)",
                "Cells: median [Q1, Q3]; mean; valid/N.  Fractions: k/N = p [Wilson 95% CI], "
                "denominator = ALL windows (failures count as 'no peak').")


def null_table(frames, system, letter):
    rows = []
    for label, df in frames:
        for n, g in df[df["system"] == system].groupby("window_length"):
            N = len(g)
            ksig = int((g["significant_peak_count"] >= 1).sum())
            assessed = int(g["significance_assessed"].astype(bool).sum())
            ns = int(g["surrogates"].iloc[0])
            # nominal size of "J(W) < 0.05" under exact exchangeability: #exceed <= ceil(0.05 n) - 1
            kmax = math.ceil(0.05 * ns) - 1
            rows.append({
                "arm": label, "window": n, "surrogates": ns, "N": N,
                "assessed_N": assessed,
                "failure_N": int(g["upo_failure"].astype(bool).sum()),
                ">=1 Level-A peak": frac((g["source_peak_count"] >= 1).sum(), N),
                ">=1 Level-B peak (empirical FPR)": frac(ksig, N),
                "nominal size (k+1)/(n+1)": f"{(kmax + 1) / (ns + 1):.4f}",
                "J_W < 0.05 (window-level)": frac((g["J_W"] < 0.05).sum(), N),
                "source_peak_count": fmt_desc(g["source_peak_count"]),
                "significant_peak_count": fmt_desc(g["significant_peak_count"]),
                "J_W": fmt_desc(g["J_W"]),
                "rJ": fmt_desc(g["source_rJ"]),
                "source_coverage": fmt_desc(g["source_peak_coverage"]),
                "significance_status": counts_str(g["upo_significance_status"]),
                "upo_status": counts_str(g["upo_status"]),
            })
    df = pd.DataFrame(rows)
    name = {"white_noise": "white noise", "ar1": f"AR(1), phi={C.AR1_PHI}"}[system]
    return save(df, f"table_{letter}_{system}_null", f"Table {letter} -- {name} empirical null",
                "FPR = fraction of windows with >= 1 Level-B (surrogate-significant) peak.  "
                "Nominal size: under exact exchangeability of observed and surrogate W, "
                "P(#{W_i > W} <= ceil(0.05 n) - 1) = ceil(0.05 n)/(n+1) (3/51 = 0.0588 for n = 50).")


def table_D(core, noise):
    clean = core[(core["window_length"] == C.PRODUCTION_WINDOW) & core["system"].isin(C.NOISE_SYSTEMS)
                 & (core["seed"] < C.N_SEEDS_NOISE)].copy()
    clean["snr_db"] = np.inf
    df = pd.concat([clean, noise], ignore_index=True)
    rows = []
    for (system, snr), g in df.groupby(["system", "snr_db"]):
        rows.append({"system": system, "SNR_dB": "clean" if np.isinf(snr) else snr, **window_summary(g),
                     "significant_coverage": fmt_desc(g["significant_peak_coverage"])})
    out = pd.DataFrame(rows)
    out["_k"] = out["SNR_dB"].map(lambda v: 1e9 if v == "clean" else v)
    out = out.sort_values(["system", "_k"], ascending=[True, False]).drop(columns="_k")
    return save(out, "table_D_noise", "Table D -- noise robustness (256 samples, paired seeds 0-29)",
                "Additive white Gaussian observation noise; SNR = 10 log10(var(clean)/sigma^2); "
                "'clean' rows are the matching core windows.")


def table_E(core, tau):
    rows = []
    base = core[core["window_length"] == C.PRODUCTION_WINDOW]
    for system in SYSTEM_ORDER:
        t = tau[tau["system"] == system]
        o = base[(base["system"] == system) & base["seed"].isin(t["seed"])]
        for label, g in (("one_sample (SOURCE)", o), ("tau_step (PROJECT EXT.)", t)):
            s = window_summary(g)
            rows.append({"system": system, "map_mode": label, "N": s["N"],
                         "upo_tau": fmt_desc(g["upo_tau"]),
                         "cao_dim": s["cao_dim (UPO)"], "n_embedded": s["n_embedded"],
                         "failure_N": s["failure_N"], "upo_status": s["upo_status"],
                         ">=1 Level-A": s[">=1 Level-A peak"], "source_count": s["source_peak_count"],
                         ">=1 Level-B": s[">=1 Level-B peak"], "sig_count": s["significant_peak_count"],
                         "rJ": s["rJ"], "coverage": s["source_coverage"]})
        m = o.merge(t, on="seed", suffixes=("_o", "_t"))
        if len(m):
            rows.append({"system": system, "map_mode": "paired: same Level-B count",
                         "N": len(m), "upo_tau": "",
                         "cao_dim": f"dim equal in {int((m['embedding_dimension_o'] == m['embedding_dimension_t']).sum())}/{len(m)}",
                         "n_embedded": "", "failure_N": "", "upo_status": "",
                         ">=1 Level-A": "", "source_count": "",
                         ">=1 Level-B": f"{int((m['significant_peak_count_o'].fillna(-1) == m['significant_peak_count_t'].fillna(-1)).sum())}/{len(m)}",
                         "sig_count": "", "rJ": "", "coverage": ""})
    return save(pd.DataFrame(rows), "table_E_map_modes",
                "Table E -- one_sample (SOURCE) vs tau_step (PROJECT EXTENSION), 256 samples, paired seeds")


def table_F(loc, locJ):
    L = pd.concat([loc, locJ], ignore_index=True)
    rows = []
    keys = ["system", "jacobian", "M", "bins", "tube", "range_policy"]
    for key, g in L.groupby(keys):
        p1, p2 = g[g["period"] == 1], g[g["period"] == 2]
        fp2_det = (p1["fp2_abs_error"] <= 2 * p1["bin_width"])
        sig = p1[p1["significance"].astype(bool)]
        rows.append({
            **dict(zip(keys, key)),
            "n_real": len(p1),
            "bin_width": round(float(p1["bin_width"].mean()), 4),
            "fp1 |err| (bin-mean) mean": float(p1["fp1_abs_error"].mean()),
            "fp1 err signed mean": float(p1["fp1_error"].mean()),
            "fp1 err / bin_width": float((p1["fp1_abs_error"] / p1["bin_width"]).mean()),
            "fp1 KDE-mode |err| mean": float(p1["fp1_kde_mode_error"].abs().mean()),
            "fp1 cell-center |err| mean": float(p1["fp1_cell_center_error"].abs().mean()),
            "fp1 strongest (rank 1)": f"{int((p1['fp1_rank'] == 1).sum())}/{len(p1)}",
            "fp2 (off-attractor) detected <=2 bins": f"{int(fp2_det.sum())}/{len(p1)}",
            "p1 source peaks mean": float(p1["source_peak_count"].mean()),
            "p1 sig peaks (N assessed)": (f"{sig['significant_peak_count'].mean():.2f} ({len(sig)})"
                                          if len(sig) else "not assessed"),
            "fp1 significant": (f"{int(sig['fp1_significant'].fillna(False).astype(bool).sum())}/{len(sig)}"
                                if len(sig) else ""),
            "p2 nearest orbit err mean": float(p2["p2_nearest_error"].mean()) if len(p2) else np.nan,
            "p2 strongest-is-true (err<=bin)": (f"{int((p2['p2_strongest_error'] <= p2['bin_width']).sum())}/{len(p2)}"
                                                if len(p2) else ""),
            "p2 n min-period-2 peaks mean": float(p2["n_minimal_period2"].mean()) if len(p2) else np.nan,
            "p2 source peaks mean": float(p2["source_peak_count"].mean()) if len(p2) else np.nan,
        })
    df = pd.DataFrame(rows)
    return save(df, "table_F_henon_localization",
                "Table F -- Henon localization vs numerical parameters (N=1024, d=2, one_sample)",
                "fp1 = on-attractor fixed point, fp2 = off-attractor fixed point.  'bin-mean' is the "
                "pipeline's location estimator (mean of tube points in the peak cell).  KDE-mode: "
                "histogram-free diagnostic (Gaussian KDE, bw 0.002, of tube scalars within one bin of "
                "the peak).  jacobian=analytic is a DIAGNOSTIC arm (exact map Jacobian substituted).")


def table_G(core, tau, noise):
    rows = []
    frames = [("one_sample", core), ("tau_step", tau), ("one_sample+noise", noise)]
    for mode, df in frames:
        keys = ["system", "window_length"] + (["snr_db"] if mode.endswith("noise") else [])
        for key, g in df.groupby(keys):
            key = key if isinstance(key, tuple) else (key,)
            system = key[0]
            if mode.endswith("noise"):
                system = f"noisy_{system} ({key[2]:g} dB)"
            for status, gg in g.groupby(g["upo_status"].fillna("NA")):
                sem = STATUS_SEMANTICS.get(status, ("", "?", ""))
                sig = counts_str(gg["upo_significance_status"])
                rows.append({"system": system, "window": key[1], "map_mode": mode,
                             "expected information": EXPECTED_INFO.get(key[0], ""),
                             "actual status": status, "N": len(gg), "of": len(g),
                             "significance_status": sig,
                             "no-structure evidence?": "yes (weak)" if sem[1].startswith("A") else "no",
                             "insufficient information?": "yes" if sem[1] == "B" else "no",
                             "notes": sem[2] + (" | underlying: " + counts_str(gg["upo_direct_status"])
                                                if "upo_direct_status" in gg and status == "analysis_error"
                                                else "")})
    df = pd.DataFrame(rows)
    return save(df, "table_G_status_matrix", "Table G -- status / failure-mode matrix")


EXPECTED_INFO = {
    "constant": "none (degenerate)",
    "white_noise": "stochastic, no deterministic structure",
    "ar1": "stochastic, linear correlation only",
    "sinusoid": "deterministic periodic (linear map, neutral fixed point at the centre)",
    "logistic_p4": "deterministic stable 4-cycle (4 distinct values)",
    "logistic": "deterministic chaos, 1-D, fixed point 0.75 unstable",
    "henon": "deterministic chaos, 2-D, fixed points 0.631 / -1.131",
}


def table_H(core, noise, tau):
    rows = []
    for label, df, keys in (("core", core, ["system", "window_length"]),
                            ("noise", noise, ["system", "snr_db"]),
                            ("tau_step (same LLE path)", tau, ["system", "window_length"])):
        for key, g in df.groupby(keys):
            v = g["lle_per_beat"]
            valid = v.notna()
            rows.append({"arm": label, "system": key[0], "window/SNR": key[1], "N": len(g),
                         "LLE valid_N": int(valid.sum()), "LLE NaN_N": int((~valid).sum()),
                         "valid-fit rate": frac(valid.sum(), len(g)),
                         "R2>=0.9 rate (of valid)": frac(g.loc[valid, "lle_valid_fit_quality"].astype(bool).sum(),
                                                         int(valid.sum())),
                         "LLE among valid": fmt_desc(v),
                         "SD": desc(v)["SD"],
                         "fraction > 0 (of valid)": frac((v > 0).sum(), int(valid.sum())),
                         "LLE tau (TDMI)": fmt_desc(g["lle_tau"]), "LLE m (Cao)": fmt_desc(g["lle_m"]),
                         "LLE m not saturated": int(g["lle_embedding_not_saturated"].astype(bool).sum()),
                         "NaN reasons": counts_str(g.loc[~valid, "lle_reason"]) if (~valid).any() else ""})
    return save(pd.DataFrame(rows), "table_H_lle",
                "Table H -- Rosenstein LLE validity and distributions (frozen estimator, per sample)",
                "Theoretical references (NOT used by the pipeline): logistic r=4 ln 2 = 0.693, "
                "Henon ~0.42, periodic <= 0, stochastic: undefined.")


def table_I(core, recs):
    peaks = []
    for r in recs:
        row = r["row"]
        for p in r["detail"].get("peaks", []):
            peaks.append({"system": row["system"], "window_length": row["window_length"], **p})
    P = pd.DataFrame(peaks)
    rows = []
    for (system, n), g in core.groupby(["system", "window_length"]):
        pp = P[(P["system"] == system) & (P["window_length"] == n)] if len(P) else pd.DataFrame()
        ok = g[~g["upo_failure"].astype(bool)]
        nA = int(ok["source_peak_count"].sum())
        nB = int(ok["significant_peak_count"].fillna(0).sum())
        nC = int(ok["verification_gated_count"].fillna(0).sum())
        nV = int(ok["verified_unstable_count"].fillna(0).sum())
        gates = {}
        if len(pp) and "levelC_failed_gates" in pp:
            for fg in pp["levelC_failed_gates"].dropna():
                for x in fg:
                    gates[x] = gates.get(x, 0) + 1
        rows.append({
            "system": system, "window": n, "windows": len(g), "non-failed": len(ok),
            "Level A peaks": nA, "Level B peaks": nB,
            "B/A": f"{nB / nA:.3f}" if nA else "NA",
            "Level C gated": nC, "C/A": f"{nC / nA:.3f}" if nA else "NA",
            "verified unstable": nV,
            "C failed gates": "; ".join(f"{k}: {v}" for k, v in sorted(gates.items())),
            "C residual (per peak)": fmt_desc(pp["levelC_residual"]) if len(pp) and "levelC_residual" in pp else "",
            "C R2 (per peak)": fmt_desc(pp["levelC_r2"]) if len(pp) and "levelC_r2" in pp else "",
            "C support ratio": fmt_desc(pp["levelC_support_ratio"]) if len(pp) and "levelC_support_ratio" in pp else "",
            "hybrid source_stability unstable": (f"{int(pp['source_stability_unstable'].fillna(False).astype(bool).sum())}"
                                                 f"/{int(pp['source_stability_unstable'].notna().sum())}"
                                                 if len(pp) and "source_stability_unstable" in pp else ""),
        })
    df = pd.DataFrame(rows)
    df["_o"] = df["system"].map(SYSTEM_ORDER.index)
    df = df.sort_values(["_o", "window"]).drop(columns="_o")
    return save(df, "table_I_levels", "Table I -- Level A vs Level B vs Level C (core arm)",
                "Level C = PROJECT EXTENSION gates (residual <= 0.05, R2 >= 0.9, support >= 2) and "
                "candidate-centred monodromy; never part of source detection."), P


# =============================================================================
# Audits
# =============================================================================

def significance_audit(recs):
    """Step 21: recompute W, W0, J(W0), J(W), rJ and per-peak J from stored histograms."""
    out = {"n_windows": 0, "W_mismatch": 0, "W0_mismatch": 0, "rJ_mismatch": 0, "J_W_mismatch": 0,
           "per_peak_J_mismatch": 0, "n_peaks": 0, "J_W0_values": {}, "Ws_recompute_note":
           "surrogate W_i are stored as computed by the pipeline; W, W0, J, rJ re-derived here"}
    for r in recs:
        d, row = r["detail"], r["row"]
        if "surrogate_W" not in d:
            continue
        out["n_windows"] += 1
        Ws = np.array(d["surrogate_W"])
        obs = np.array(d["observed_histogram"])
        mean = np.array(d["surrogate_mean_histogram"])
        W = float(np.max(obs - mean))
        W0 = float(np.median(Ws))
        J = lambda w: float(np.mean(Ws > w))  # noqa: E731
        out["W_mismatch"] += int(not np.isclose(W, row["W"], rtol=0, atol=1e-9))
        out["W0_mismatch"] += int(not np.isclose(W0, row["W0"], rtol=0, atol=1e-12))
        out["J_W_mismatch"] += int(J(W) != row["J_W"])
        rj = W / W0 if W0 > 0 else np.nan
        out["rJ_mismatch"] += int(not (np.isclose(rj, row["rJ"]) or (np.isnan(rj) and np.isnan(row["rJ"]))))
        k = f"{J(W0):.3f}"
        out["J_W0_values"][k] = out["J_W0_values"].get(k, 0) + 1
        for p in d["peaks"]:
            if p.get("period") != 1 or "J" not in p:
                continue
            out["n_peaks"] += 1
            dev = obs[p["peak_cell"][0]] - mean[p["peak_cell"][0]]
            out["per_peak_J_mismatch"] += int(J(dev) != p["J"] or
                                              p["significant"] != (J(dev) < row.get("alpha", 0.05)))
    # Hand-constructed check of fp.so_surrogate_significance against the formula.
    rng = np.random.default_rng(0)
    edges = [np.linspace(0, 1, 11)]
    sur = [{"reduced": rng.uniform(0, 1, (200, 1))} for _ in range(50)]
    obs_pts = np.concatenate([rng.uniform(0, 1, (190, 1)), np.full((30, 1), 0.55)])
    obs_hist = fp._histogram_on_edges(obs_pts, edges)
    res = {"edges": edges, "histogram": obs_hist,
           fp.LEVEL_A_FIELD: [{"peak_cell": (5,)}]}
    sig = fp.so_surrogate_significance(res, sur, alpha=0.05)
    H = np.array([fp._histogram_on_edges(s["reduced"], edges) for s in sur])
    m = H.mean(axis=0)
    Ws = (H - m).max(axis=1)
    W = (obs_hist - m).max()
    out["hand_check"] = {
        "W": bool(np.isclose(sig["W"], W)), "W0": bool(np.isclose(sig["W0"], np.median(Ws))),
        "J_W": bool(sig["J_W"] == np.mean(Ws > W)), "J_W0": float(sig["J_W0"]),
        "rJ": bool(np.isclose(sig["rJ"], W / np.median(Ws))),
        "signed_not_abs": bool(np.all(np.isclose(sig["deviation"], obs_hist - m))),
        "peak_J": bool(sig["per_peak"][0]["J"] == np.mean(Ws > (obs_hist - m)[5])),
    }
    # Exchangeability asymmetry: surrogate i is included in the mean it is compared with.
    out["mean_includes_self_note"] = (
        "W_i uses a mean histogram that contains surrogate i (deviation shrunk by (n-1)/n), while "
        "the observed W uses a mean that does not contain the observation; with n = 50 the scale "
        "ratio is sqrt((n+1)/(n-1)) = %.4f -- a small anti-conservative asymmetry." % math.sqrt(51 / 49))
    return out


def status_consistency(df):
    """Level B None-vs-[] and significance-status semantics on every window."""
    bad = []
    for _, r in df.iterrows():
        st, ss = r["upo_status"], r["upo_significance_status"]
        if r["upo_failure"] and ss not in ("not_assessed_detector_failure", "not_requested"):
            bad.append((r["task_id"], "failure but significance status " + str(ss)))
        if ss == "assessed" and (pd.isna(r["significant_peak_count"]) or pd.isna(r["J_W"])):
            bad.append((r["task_id"], "assessed but Level B NaN"))
        if ss != "assessed" and not pd.isna(r["significant_peak_count"]):
            bad.append((r["task_id"], "not assessed but Level B count present"))
        if st in ("tube_empty", "too_few_in_tube", "no_peaks") and r["source_peak_count"] != 0:
            bad.append((r["task_id"], f"{st} but source count {r['source_peak_count']}"))
        if r["upo_failure"] and not pd.isna(r["source_peak_count"]):
            bad.append((r["task_id"], "failure but source count not NaN"))
    return bad


def coverage_audit(core):
    from scipy.stats import spearmanr
    rows = []
    ok = core[~core["upo_failure"].astype(bool)]
    for system, g in [("ALL (non-failed)", ok)] + list(ok.groupby("system")):
        g = g[["source_peak_count", "source_peak_coverage"]].dropna()
        if len(g) < 5 or g["source_peak_count"].nunique() < 2:
            rows.append({"system": system, "N": len(g), "spearman_rho": np.nan, "p": np.nan,
                         "coverage_SD_within_count (median over counts)": np.nan,
                         "coverage==0 when count==0": ""})
            continue
        rho, p = spearmanr(g["source_peak_count"], g["source_peak_coverage"])
        within = g.groupby("source_peak_count")["source_peak_coverage"].std().dropna()
        zero = g[g["source_peak_count"] == 0]
        rows.append({"system": system, "N": len(g), "spearman_rho": rho, "p": p,
                     "coverage_SD_within_count (median over counts)": float(within.median()) if len(within) else np.nan,
                     "coverage==0 when count==0": f"{int((zero['source_peak_coverage'] == 0).sum())}/{len(zero)}"})
    return save(pd.DataFrame(rows), "table_coverage_vs_count",
                "Step 23 -- source_peak_coverage vs source_peak_count (core, non-failed windows)")


def fp_audit(fpdf):
    rows = []
    keys = ["embedding_dimension", "upo_status", "source_peak_count", "lle_per_beat", "lle_m", "lle_tau"]
    for system, g in fpdf.groupby("system"):
        piv = {v: gg.set_index("seed") for v, gg in g.groupby("variant")}
        for a, b in (("original", "nextafter"), ("rr_route_a", "rr_route_b")):
            A, B = piv[a], piv[b]
            common = A.index.intersection(B.index)
            rec = {"system": system, "comparison": f"{a} vs {b}", "N": len(common)}
            for k in keys:
                x, y = A.loc[common, k], B.loc[common, k]
                same = [(xx == yy) or (pd.isna(xx) and pd.isna(yy)) for xx, yy in zip(x, y)]
                rec[f"{k} identical"] = f"{sum(same)}/{len(common)}"
            rec["distinct values (route a) median"] = float(A.loc[common, "n_distinct_values"].median())
            rec["status (b)"] = counts_str(B.loc[common, "upo_status"])
            rows.append(rec)
        for v in ("rr_route_a",):
            A, O = piv[v], piv["original"]
            common = A.index.intersection(O.index)
            rows.append({"system": system, "comparison": f"original vs {v} (quantized to 1/{C.RR_FS:g} s)",
                         "N": len(common),
                         "embedding_dimension identical": f"{int((A.loc[common,'embedding_dimension'].fillna(-1) == O.loc[common,'embedding_dimension'].fillna(-1)).sum())}/{len(common)}",
                         "upo_status identical": f"{int((A.loc[common,'upo_status'] == O.loc[common,'upo_status']).sum())}/{len(common)}",
                         "status (b)": counts_str(A.loc[common, "upo_status"]),
                         "distinct values (route a) median": float(A.loc[common, "n_distinct_values"].median())})
    return save(pd.DataFrame(rows), "table_fp_sensitivity",
                "Step 15 -- floating-point / RR-quantization sensitivity (256 samples, no significance)")


def rr_route_bitdiff():
    """How many samples differ between the two RR quantization routes (same inputs)."""
    from experiments.phase2e import systems as S
    out = {}
    for system in C.FP_SYSTEMS:
        diffs = []
        for s in range(C.N_SEEDS_FP):
            x = S.generate(system, C.PRODUCTION_WINDOW, s)
            rr = 0.8 + 0.05 * (x - np.mean(x)) / np.std(x)
            a = np.round(rr * C.RR_FS) / C.RR_FS
            b = np.round(rr * C.RR_FS) * (1.0 / C.RR_FS)
            diffs.append(int(np.sum(a != b)))
        out[system] = {"median_differing_samples": float(np.median(diffs)), "max": int(max(diffs)),
                       "max_abs_difference": 2.2e-16}
    return out


def surrogate_count_table(core, sens, algo):
    rows = []
    base = core[core["window_length"] == C.PRODUCTION_WINDOW]
    for system in C.SURROGATE_SENS_SYSTEMS:
        seeds = range(C.N_SEEDS_SURROGATE_SENS[system])
        for ns in sorted(set(C.SURROGATE_SENSITIVITY) | {C.MAIN_SURROGATES}):
            g = (base[(base["system"] == system) & base["seed"].isin(seeds)] if ns == C.MAIN_SURROGATES
                 else sens[(sens["system"] == system) & (sens["surrogates"] == ns)])
            kmax = math.ceil(0.05 * ns) - 1
            rows.append({"system": system, "surrogates": ns, "N": len(g),
                         ">=1 Level-B": frac((g["significant_peak_count"] >= 1).sum(), len(g)),
                         "nominal size": f"{(kmax + 1) / (ns + 1):.4f}",
                         "significant_peak_count": fmt_desc(g["significant_peak_count"]),
                         "rJ": fmt_desc(g["source_rJ"]), "J_W": fmt_desc(g["J_W"]),
                         "J resolution": f"1/{ns}"})
    for system in C.ALGO_SEED_SYSTEMS:
        g = algo[algo["system"] == system]
        rows.append({"system": system, "surrogates": f"{C.MAIN_SURROGATES} (algo seed varied)", "N": len(g),
                     ">=1 Level-B": frac((g["significant_peak_count"] >= 1).sum(), len(g)),
                     "nominal size": f"{3 / 51:.4f}",
                     "significant_peak_count": fmt_desc(g["significant_peak_count"]),
                     "rJ": fmt_desc(g["source_rJ"]), "J_W": fmt_desc(g["J_W"]), "J resolution": "1/50"})
    return save(pd.DataFrame(rows), "table_surrogate_count",
                "Step 8 -- surrogate-count and algorithm-seed sensitivity (256 samples, paired seeds)")


def oracle_table(core, oracle):
    rows = []
    base = core[core["window_length"] == C.PRODUCTION_WINDOW]
    for system, d in C.ORACLE_DIMENSION.items():
        o = oracle[oracle["system"] == system]
        c = base[(base["system"] == system) & base["seed"].isin(o["seed"])]
        for label, g in ((f"production (Cao, lag 1)", c), (f"DIAGNOSTIC fixed d={d}", o)):
            s = window_summary(g)
            rows.append({"system": system, "arm": label, "N": s["N"], "dim": s["cao_dim (UPO)"],
                         ">=1 Level-A": s[">=1 Level-A peak"], "source_count": s["source_peak_count"],
                         ">=1 Level-B": s[">=1 Level-B peak"], "sig_count": s["significant_peak_count"],
                         "rJ": s["rJ"], "coverage": s["source_coverage"],
                         "verified_unstable": s["verified_unstable_count"]})
    return save(pd.DataFrame(rows), "table_oracle_dimension",
                "Step 14 -- production Cao dimension vs DIAGNOSTIC known dimension (256 samples)")


def peak_location_summary(P):
    """Where Level-A / Level-B period-1 peaks sit for the chaotic controls (256, core)."""
    rows = []
    refs = {"logistic": [0.75, 0.0], "henon": [0.6313544770895047, -1.1313544770895048],
            "sinusoid": [0.0]}
    for (system, n), g in P[P["system"].isin(refs)].groupby(["system", "window_length"]):
        loc = g["location"].map(lambda v: v[0])
        near = [min(abs(l - r) for r in refs[system]) for l in loc]
        g = g.assign(near=near)
        sig = g[g.get("significant", pd.Series(False, index=g.index)).fillna(False).astype(bool)]
        rows.append({"system": system, "window": n, "Level-A peaks": len(g),
                     "A within 0.05 of a known fixed point": f"{int((g['near'] <= 0.05).sum())}/{len(g)}",
                     "Level-B peaks": len(sig),
                     "B within 0.05 of a known fixed point": f"{int((sig['near'] <= 0.05).sum())}/{len(sig)}",
                     "A location": fmt_desc(loc)})
    return save(pd.DataFrame(rows), "table_peak_locations",
                "Level-A / Level-B peak locations vs analytic fixed points (core; logistic 0.75, "
                "Henon 0.631/-1.131, sinusoid centre 0)")


# =============================================================================
# Plots
# =============================================================================

def _style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 300, "font.size": 9, "axes.titlesize": 10,
        "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
        "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
        "figure.facecolor": "white", "axes.facecolor": "white", "lines.linewidth": 2,
    })
    return plt


def _savefig(fig, name):
    PLOTS.mkdir(parents=True, exist_ok=True)
    fig.savefig(PLOTS / f"{name}.png", bbox_inches="tight")


def _strip(ax, groups, labels, colors, ylabel, jitter=0.12, seed=0):
    rng = np.random.default_rng(seed)
    for i, (vals, c) in enumerate(zip(groups, colors)):
        v = pd.to_numeric(pd.Series(vals), errors="coerce").dropna().to_numpy()
        if len(v) == 0:
            ax.text(i, 0, "all NaN", ha="center", va="bottom", color=INK2, fontsize=7)
            continue
        ax.scatter(i + rng.uniform(-jitter, jitter, len(v)), v, s=10, color=c, alpha=0.55,
                   edgecolor="none")
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        ax.plot([i - 0.25, i + 0.25], [med, med], color=INK, lw=2)
        ax.plot([i, i], [q1, q3], color=INK, lw=1)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel(ylabel)


def plots(core, noise, tau, loc, locJ, sens, algo):
    plt = _style()
    systems = [s for s in SYSTEM_ORDER if s != "constant"]

    # 1-3, 6: distributions by system at each window length
    for col, name, ylabel in (("lle_per_beat", "01_lle_by_system", "Rosenstein LLE (per sample)"),
                              ("source_peak_count", "02_source_peaks_by_system", "Level-A source peak count"),
                              ("significant_peak_count", "03_significant_peaks_by_system",
                               "Level-B significant peak count"),
                              ("source_rJ", "06_rJ_by_system", "r_J = W / W0")):
        fig, axes = plt.subplots(1, 3, figsize=(11, 3.4), sharey=True)
        for ax, n in zip(axes, C.WINDOW_LENGTHS):
            g = core[core["window_length"] == n]
            _strip(ax, [g.loc[g["system"] == s, col] for s in systems], systems,
                   [SYSTEM_COLOR[s] for s in systems], ylabel if n == 128 else "")
            nan_n = [int(g.loc[g["system"] == s, col].isna().sum()) for s in systems]
            ax.set_title(f"{n} samples  (NaN per system: {', '.join(map(str, nan_n))})", fontsize=8)
            if col == "lle_per_beat":
                ax.axhline(0, color=INK2, lw=0.8)
        fig.suptitle({"lle_per_beat": "Rosenstein LLE by system (constant: all NaN, not shown)",
                      "source_peak_count": "Level-A source peaks by system (failed windows = NaN, not shown)",
                      "significant_peak_count": "Level-B significant peaks by system",
                      "source_rJ": "So r_J by system"}[col], fontsize=10)
        _savefig(fig, name)
        plt.close(fig)

    # 4-5: null distributions
    for system, name in (("white_noise", "04_white_noise_null"), ("ar1", "05_ar1_null")):
        g = core[(core["system"] == system) & (core["window_length"] == C.PRODUCTION_WINDOW)]
        fig, axes = plt.subplots(1, 4, figsize=(12, 3))
        for ax, col, lab in zip(axes, ["source_peak_count", "significant_peak_count", "J_W", "source_rJ"],
                                ["Level-A count", "Level-B count", "J(W)", "r_J"]):
            v = g[col].dropna()
            if col in ("source_peak_count", "significant_peak_count"):
                vc = v.value_counts().sort_index()
                ax.bar(vc.index, vc.values, width=0.8, color=PALETTE[0])
                ax.set_xticks(vc.index)
            else:
                ax.hist(v, bins=20, color=PALETTE[0], edgecolor="white", linewidth=0.5)
            if col == "J_W":
                ax.axvline(0.05, color=PALETTE[7], lw=1)
                ax.text(0.06, ax.get_ylim()[1] * 0.9, "0.05", color=INK2, fontsize=7)
            ax.set_xlabel(lab)
            ax.set_ylabel("windows")
        k = int((g["significant_peak_count"] >= 1).sum())
        fig.suptitle(f"{system} empirical null, 256 samples, 50 surrogates: "
                     f"{k}/{len(g)} windows with >= 1 Level-B peak", fontsize=10)
        fig.tight_layout()
        _savefig(fig, name)
        plt.close(fig)

    # 7: noise robustness curves
    clean = core[(core["window_length"] == 256) & core["system"].isin(C.NOISE_SYSTEMS) & (core["seed"] < 30)].copy()
    clean["snr_db"] = 40.0  # plotted position for 'clean'
    df = pd.concat([clean, noise])
    fig, axes = plt.subplots(1, 4, figsize=(13, 3.2))
    metrics = [("lle_per_beat", "median LLE (valid)", np.nanmedian),
               ("source_peak_count", "mean Level-A count", np.nanmean),
               ("sig_any", "fraction >= 1 Level-B peak", np.nanmean),
               ("source_rJ", "median r_J", np.nanmedian)]
    df["sig_any"] = (df["significant_peak_count"] >= 1).astype(float)
    for ax, (col, lab, fn) in zip(axes, metrics):
        for i, system in enumerate(C.NOISE_SYSTEMS):
            g = df[df["system"] == system].groupby("snr_db")[col].agg(fn).sort_index(ascending=False)
            ax.plot(g.index, g.values, marker="o", ms=5, color=PALETTE[i], label=system)
        ax.set_xticks([40, 30, 20, 10, 5, 0])
        ax.set_xticklabels(["clean", "30", "20", "10", "5", "0"])
        ax.invert_xaxis()
        ax.set_xlabel("SNR (dB)")
        ax.set_ylabel(lab)
    axes[0].legend()
    fig.suptitle("Noise robustness, 256 samples, 30 paired seeds per point", fontsize=10)
    fig.tight_layout()
    _savefig(fig, "07_noise_robustness")
    plt.close(fig)

    # 8: window-length comparison
    fig, axes = plt.subplots(1, 4, figsize=(13, 3.2))
    core2 = core.assign(a_any=(core["source_peak_count"] >= 1).astype(float),
                        b_any=(core["significant_peak_count"] >= 1).astype(float),
                        lle_ok=core["lle_valid"].astype(float))
    for ax, (col, lab) in zip(axes, [("embedding_dimension", "median Cao dimension (lag 1)"),
                                     ("a_any", "fraction >= 1 Level-A"),
                                     ("b_any", "fraction >= 1 Level-B"),
                                     ("lle_ok", "LLE valid-fit rate")]):
        for system in systems:
            g = core2[core2["system"] == system].groupby("window_length")[col]
            g = g.median() if col == "embedding_dimension" else g.mean()
            ax.plot(g.index, g.values, marker="o", ms=5, color=SYSTEM_COLOR[system], label=system)
        ax.set_xticks(C.WINDOW_LENGTHS)
        ax.set_xlabel("window length")
        ax.set_ylabel(lab)
    axes[0].legend(fontsize=7)
    fig.suptitle("Window-length behaviour (one_sample, core arm)", fontsize=10)
    fig.tight_layout()
    _savefig(fig, "08_window_length")
    plt.close(fig)

    # 9: one_sample vs tau_step
    base = core[core["window_length"] == 256]
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.4))
    for ax, (col, lab) in zip(axes, [("embedding_dimension", "embedding dimension"),
                                     ("source_peak_count", "Level-A count"),
                                     ("significant_peak_count", "Level-B count")]):
        for j, (label, d) in enumerate((("one_sample", base), ("tau_step", tau))):
            vals = [d.loc[(d["system"] == s) & (d["seed"] < 30), col].mean() for s in systems]
            ax.bar(np.arange(len(systems)) + (j - 0.5) * 0.38, vals, width=0.36, color=PALETTE[j], label=label)
        ax.set_xticks(range(len(systems)))
        ax.set_xticklabels(systems, rotation=30, ha="right")
        ax.set_ylabel(f"mean {lab} (non-NaN)")
    axes[0].legend()
    fig.suptitle("one_sample (SOURCE) vs tau_step (PROJECT EXTENSION), 256 samples, paired seeds", fontsize=10)
    fig.tight_layout()
    _savefig(fig, "09_map_modes")
    plt.close(fig)

    # 10-11: Henon localization
    L = loc[(loc["M"] == 2)]
    for period, name, col, lab in ((1, "10_henon_fixed_point_localization", "fp1_abs_error",
                                    "|bin-mean location - x*|"),
                                   (2, "11_henon_period2_localization", "p2_nearest_error",
                                    "orbit distance to true 2-cycle")):
        fig, axes = plt.subplots(1, 3, figsize=(13, 3.4), sharey=True)
        for ax, system in zip(axes[:2], C.LOC_SYSTEMS):
            g = L[(L["period"] == period) & (L["system"] == system) & (L["range_policy"] == "pad10")]
            for i, tube in enumerate(C.LOC_TUBE_PERCENTILES):
                gg = g[g["tube"] == tube].groupby("bins")[col].mean()
                ax.plot(gg.index, gg.values, marker="o", ms=4, color=PALETTE[i], label=f"tube {tube:g}%")
            bw = g.groupby("bins")["bin_width"].mean()
            ax.plot(bw.index, bw.values / 2, ls="--", color=INK2, lw=1, label="half bin width")
            if period == 1:
                gk = g[g["tube"] == 10.0].groupby("bins")["fp1_kde_mode_error"].apply(lambda v: v.abs().mean())
                ax.plot(gk.index, gk.values, marker="s", ms=4, color=INK, lw=1, label="KDE mode (tube 10%)")
            ax.set_xscale("log")
            ax.set_xticks(C.LOC_BINS)
            ax.set_xticklabels(C.LOC_BINS)
            ax.set_xlabel("histogram bins (range pad10)")
            ax.set_ylabel(lab)
            ax.set_title(system)
        ax = axes[2]
        g = L[(L["period"] == period) & (L["tube"] == 10.0)]
        for i, pol in enumerate(C.LOC_RANGE_POLICIES):
            for j, system in enumerate(C.LOC_SYSTEMS):
                gg = g[(g["range_policy"] == pol) & (g["system"] == system)].groupby("bins")[col].mean()
                ax.plot(gg.index, gg.values, marker="o" if j == 0 else "^", ms=4, color=PALETTE[i],
                        ls="-" if j == 0 else ":", label=f"{pol} ({'std' if j == 0 else 'skewed'})")
        ax.set_xscale("log")
        ax.set_xticks(C.LOC_BINS)
        ax.set_xticklabels(C.LOC_BINS)
        ax.set_xlabel("histogram bins (tube 10%)")
        ax.set_title("range policy")
        for a in axes:
            a.legend(fontsize=6)
        fig.suptitle(f"Henon period-{period} localization (N=1024, d=2, M=2; mean over 3 realizations)",
                     fontsize=10)
        fig.tight_layout()
        _savefig(fig, name)
        plt.close(fig)

    # 12: status frequencies
    frames = [("one_sample", core), ("tau_step", tau)]
    statuses = sorted(set(core["upo_status"].dropna()) | set(tau["upo_status"].dropna()))
    fig, axes = plt.subplots(1, 4, figsize=(14, 3.6), sharey=True)
    panels = [(m, n) for n in C.WINDOW_LENGTHS for m in ("one_sample",)] + [("tau_step", 256)]
    for ax, (mode, n) in zip(axes, panels):
        d = dict(frames)[mode]
        d = d[d["window_length"] == n]
        bottom = np.zeros(len(SYSTEM_ORDER))
        for i, st in enumerate(statuses):
            vals = np.array([(d.loc[d["system"] == s, "upo_status"] == st).mean() if (d["system"] == s).any()
                             else 0 for s in SYSTEM_ORDER])
            ax.bar(range(len(SYSTEM_ORDER)), vals, bottom=bottom, color=PALETTE[i % 8], label=st,
                   edgecolor="white", linewidth=1)
            bottom += vals
        ax.set_xticks(range(len(SYSTEM_ORDER)))
        ax.set_xticklabels(SYSTEM_ORDER, rotation=35, ha="right")
        ax.set_title(f"{mode}, {n} samples")
    axes[0].set_ylabel("fraction of windows")
    axes[-1].legend(fontsize=7, loc="upper left", bbox_to_anchor=(1, 1))
    fig.suptitle("UPO status frequencies", fontsize=10)
    fig.tight_layout()
    _savefig(fig, "12_status_frequencies")
    plt.close(fig)

    # 13: surrogate count sensitivity
    fig, ax = plt.subplots(figsize=(5.5, 3.4))
    base = core[core["window_length"] == 256]
    for i, system in enumerate(C.SURROGATE_SENS_SYSTEMS):
        xs, ys, lo, hi = [], [], [], []
        for ns in (20, 50, 100):
            g = (base[(base["system"] == system) & (base["seed"] < C.N_SEEDS_SURROGATE_SENS[system])]
                 if ns == 50 else sens[(sens["system"] == system) & (sens["surrogates"] == ns)])
            k, N = int((g["significant_peak_count"] >= 1).sum()), len(g)
            a, b = wilson(k, N)
            xs.append(ns); ys.append(k / N); lo.append(k / N - a); hi.append(b - k / N)
        ax.errorbar(xs, ys, yerr=[lo, hi], marker="o", ms=5, color=PALETTE[i], label=system, capsize=3, lw=1.5)
    ax.axhline(0.05, color=INK2, lw=0.8, ls="--")
    ax.set_xticks([20, 50, 100])
    ax.set_xlabel("surrogates")
    ax.set_ylabel("fraction with >= 1 Level-B peak")
    ax.legend(fontsize=7)
    ax.set_title("Surrogate-count sensitivity (Wilson 95% CI)")
    fig.tight_layout()
    _savefig(fig, "13_surrogate_count")
    plt.close(fig)

    # 14: coverage vs count
    ok = core[~core["upo_failure"].astype(bool)]
    fig, ax = plt.subplots(figsize=(5.5, 3.6))
    rng = np.random.default_rng(1)
    for s in systems:
        g = ok[ok["system"] == s]
        ax.scatter(g["source_peak_count"] + rng.uniform(-0.15, 0.15, len(g)), g["source_peak_coverage"],
                   s=9, alpha=0.5, color=SYSTEM_COLOR[s], label=s, edgecolor="none")
    ax.set_xlabel("Level-A source peak count (jittered)")
    ax.set_ylabel("source_peak_coverage (PROJECT)")
    ax.legend(fontsize=7)
    ax.set_title("Coverage vs count (core, all window lengths)")
    fig.tight_layout()
    _savefig(fig, "14_coverage_vs_count")
    plt.close(fig)


# =============================================================================

def main():
    core, core_recs = load("core")
    noise, noise_recs = load("noise")
    tau, tau_recs = load("tau_step")
    sens, sens_recs = load("surrogate_sens")
    algo, algo_recs = load("algo_seed")
    oracle, _ = load("oracle_dim")
    fpdf, _ = load("fp_sens")
    loc, _ = load("henon_loc")
    locJ, _ = load("henon_loc_oracleJ")
    stab, _ = load("stability")
    flag, _ = load("flag_independence")

    table_A(core)
    frames = [("core", core), ("surrogate_sens", sens), ("algo_seed (random_seed varied)", algo)]
    null_table(frames, "white_noise", "B")
    null_table(frames, "ar1", "C")
    table_D(core, noise)
    table_E(core, tau)
    table_F(loc, locJ)
    table_G(core, tau, noise)
    table_H(core, noise, tau)
    _, P = table_I(core, core_recs)
    coverage_audit(core)
    fp_audit(fpdf)
    surrogate_count_table(core, sens, algo)
    oracle_table(core, oracle)
    peak_location_summary(P)

    all_recs = core_recs + noise_recs + tau_recs + sens_recs + algo_recs
    all_windows = pd.concat([core, noise, tau, sens, algo, fpdf], ignore_index=True)
    audits = {
        "significance_audit": significance_audit(all_recs),
        "status_consistency_violations": status_consistency(all_windows),
        "harness_errors": (all_windows["harness_error"].dropna().tolist()
                           if "harness_error" in all_windows else []),
        "n_windows_total": int(len(all_windows)),
        "runtime_warnings_total": int(all_windows["n_runtime_warnings"].sum()),
        "rr_route_bit_differences": rr_route_bitdiff(),
        "stability_detector": stab[stab["kind"] == "detector"].to_dict(orient="records"),
        "stability_pipeline": stab[stab["kind"] == "pipeline"].to_dict(orient="records"),
        "flag_independence": flag.to_dict(orient="records"),
        "total_runtime_s": float(sum(d["runtime_s"].sum() for d in (core, noise, tau, sens, algo, oracle,
                                                                    fpdf, loc, locJ, stab, flag) if len(d))),
    }
    rep = RESULTS / "replicability.json"
    if rep.exists():
        r = json.loads(rep.read_text())
        audits["replicability"] = {"n": len(r), "row_identical": sum(x["row_identical"] for x in r),
                                   "detail_identical": sum(x["detail_identical"] for x in r),
                                   "non_identical": [x for x in r if not (x["row_identical"] and x["detail_identical"])]}
    (RESULTS / "audits.json").write_text(json.dumps(audits, indent=1, default=str))
    plots(core, noise, tau, loc, locJ, sens, algo)
    print("tables:", sorted(p.name for p in TABLES.glob("*.md")))
    print(json.dumps({k: v for k, v in audits.items()
                      if k in ("significance_audit", "status_consistency_violations", "harness_errors",
                               "n_windows_total", "replicability")}, indent=1, default=str)[:4000])


if __name__ == "__main__":
    main()
