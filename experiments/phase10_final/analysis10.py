"""
Phase 10 analyses (Questions 1-4) on the runner's JSONL output; identical code for DEVELOPMENT and
CONFIRMATORY data (definitions fixed in PREREGISTRATION.md).

    python -m experiments.phase10_final.analysis10 --phase dev
    python -m experiments.phase10_final.analysis10 --phase conf

Writes results/<phase>/analysis/analysis.json, tables.md and plots/*.png.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import warnings

import numpy as np
import pandas as pd
from scipy import stats

HERE = pathlib.Path(__file__).resolve().parent
B_BOOT = 10000          # cluster bootstrap resamples (rates, differences)
B_CURVE = 2000          # cluster bootstrap resamples (detection curves)
F_FINE = np.round(np.arange(0.05, 0.9001, 0.01), 2)
THRESHOLDS = (0.05, 0.20)
SEED = 20261011
# Phase 9 TEST pooled false-positive counts in the 101 PASS conditions (PHASE9 report Section 6):
SYNTH_FP = {"K3": (0, 10100), "K1": (97, 10100)}
MIN_POS_GEE = 10
ASSUMED_START_H = 9.9   # 09:54 = median start of the 33 development long-term records (no-start-time fallback)


def _clock(h, t_start):
    if h is not None and np.isfinite(h):
        return float(h), False
    return float((ASSUMED_START_H + t_start / 3600.0) % 24.0), True


# ------------------------------------------------------------------------------- basic statistics
def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (float(max(0.0, c - h)), float(min(1.0, c + h)))


def cp_upper(k, n, alpha=0.05):
    """One-sided (1 - alpha) Clopper-Pearson upper bound."""
    return 1.0 if k >= n else float(stats.beta.ppf(1 - alpha, k + 1, n - k))


def cp_lower(k, n, alpha=0.05):
    return 0.0 if k <= 0 else float(stats.beta.ppf(alpha, k, n - k + 1))


def subject_rates(df, col, subj="subject"):
    g = df.groupby(subj)[col].mean()
    return g.values.astype(float), list(g.index)


def boot_mean(vals, B=B_BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    vals = np.asarray(vals, float)
    n = len(vals)
    if n == 0:
        return np.array([np.nan])
    return vals[rng.integers(0, n, size=(B, n))].mean(1)


def boot_diff(a, b, B=B_BOOT, seed=SEED + 1):
    rng = np.random.default_rng(seed)
    a, b = np.asarray(a, float), np.asarray(b, float)
    return a[rng.integers(0, len(a), (B, len(a)))].mean(1) - b[rng.integers(0, len(b), (B, len(b)))].mean(1)


def perm_p(a, b, B=B_BOOT, seed=SEED + 2):
    rng = np.random.default_rng(seed)
    a, b = np.asarray(a, float), np.asarray(b, float)
    d0 = a.mean() - b.mean()
    allv = np.concatenate([a, b])
    cnt = 0
    for _ in range(B):
        q = rng.permutation(allv)
        cnt += abs(q[:len(a)].mean() - q[len(a):].mean()) >= abs(d0) - 1e-15
    return float((1 + cnt) / (B + 1))


def ci(v, lo=2.5, hi=97.5):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return [float(np.percentile(v, lo)), float(np.percentile(v, hi))] if len(v) else [float("nan")] * 2


# ------------------------------------------------------------------------------- Firth logistic
def firth_logit(X, y, max_iter=100, tol=1e-8):
    """Firth (1993) bias-reduced logistic regression (Jeffreys-prior penalised likelihood), Newton
    iterations on the modified score X'(y - p + h (1/2 - p)).  Finite under separation."""
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    b = np.zeros(X.shape[1])
    for _ in range(max_iter):
        eta = np.clip(X @ b, -30, 30)
        p = 1 / (1 + np.exp(-eta))
        w = p * (1 - p)
        XtWX = X.T @ (X * w[:, None])
        try:
            inv = np.linalg.inv(XtWX)
        except np.linalg.LinAlgError:
            inv = np.linalg.pinv(XtWX)
        h = np.einsum("ij,jk,ik->i", X * np.sqrt(w)[:, None], inv, X * np.sqrt(w)[:, None])
        U = X.T @ (y - p + h * (0.5 - p))
        step = inv @ U
        # step halving for stability
        mx = np.max(np.abs(step))
        if mx > 5:
            step *= 5 / mx
        b = b + step
        if np.max(np.abs(step)) < tol:
            break
    return b


def curve(b, f):
    return 1 / (1 + np.exp(-(b[0] + b[1] * np.asarray(f))))


# ------------------------------------------------------------------------------- loading
EXCLUDE = set()          # Amendment 1 sensitivity analysis: subjects dropped from every frame


def load(phase, part):
    p = HERE / "results" / phase / f"{part}.jsonl"
    if not p.exists():
        return []
    R = [json.loads(line) for line in open(p)]
    R = [r for r in R if r.get("subject") not in EXCLUDE]
    seen, out = set(), []
    for r in R:
        if r["id"] not in seen:
            seen.add(r["id"])
            out.append(r)
    return out


def real_frame(R):
    rows = []
    for r in R:
        if r.get("error"):
            rows.append({"id": r["id"], "error": r["error"]})
            continue
        k1, k1e = r["K1"], r["K1_edited"]
        rows.append({
            "id": r["id"], "subject": r["subject"], "group": r["group"], "j": r["j"], "clock_hour": r["clock_hour"],
            "t_start": r["t_start"],
            "burden": r["burden"], "n_masked": r["n_masked"], "nn_fraction": r["nn_fraction"],
            "mean_rr": r["mean_rr"], "sdnn": r["sdnn"], "rmssd": r["rmssd"],
            "K1": bool(k1["and_detected"]), "LLE": bool(k1["lle_detected"]), "UPO": bool(k1["upo_detected"]),
            "K1_err": k1.get("error") is not None,
            "K1e_eligible": k1e.get("eligible", True), "K1e": bool(k1e.get("and_detected", False)),
            "LLEe": bool(k1e.get("lle_detected", False)), "UPOe": bool(k1e.get("upo_detected", False)),
            "K3": bool(r["K3"]["detected"]), "K3_analysable": bool(r["K3"].get("analysable")),
            "K4": bool(r["K4"]["detected"]), "K3RR": bool(r["K3RR"]["detected"]),
            "TIT": bool(r["TIT_raw"]["positive"]), "TIT_NL": r["TIT_raw"]["NL"],
            "TITm": bool(r["TIT_masked"].get("positive", False)),
            "TITm_analysable": bool(r["TIT_masked"].get("analysable", False)),
            "TITe": bool(r["TIT_edited"].get("positive", False)),
            "TITe_eligible": bool(r["TIT_edited"].get("eligible", False)),
        })
    df = pd.DataFrame(rows)
    if "error" in df:
        df = df[df["error"].isna()].drop(columns=["error"])
    df["chf"] = (df["group"] == "CHF").astype(int)
    df["x"] = np.log2(1 + df["burden"])
    ck = [_clock(h, t) for h, t in zip(df["clock_hour"], df["t_start"])]
    df["clock_eff"] = [c[0] for c in ck]
    df["clock_approx"] = [c[1] for c in ck]
    df["night"] = (df["clock_eff"] >= 0) & (df["clock_eff"] < 5)
    return df


# ------------------------------------------------------------------------------- Question 1
def spike_frame(R):
    rows = []
    for r in R:
        if r.get("error") or "designs" not in r:
            continue
        for d, v in r["designs"].items():
            rows.append({"subject": r["subject"], "group": r["group"], "family": r["family"], "param": r["param"],
                         "lam": r["lam"], "design": "replace" if d == "replace" else "add",
                         "f": v.get("f", np.nan), "K1": bool(v["K1"]["and_detected"]),
                         "LLE": bool(v["K1"]["lle_detected"]), "UPO": bool(v["K1"]["upo_detected"]),
                         "K3": bool(v["K3"]["detected"])})
    return pd.DataFrame(rows)


def fit_detection(sp, det, B=B_CURVE, seed=SEED + 10):
    """Per (family, param): Firth logistic logit p = b0 + b1 f on the additive design; subject-cluster
    bootstrap pointwise bands on F_FINE.  Empirical rates (Wilson) at each design f; replacement rate."""
    out = {}
    for (fam, par), g in sp.groupby(["family", "param"]):
        a = g[g.design == "add"]
        subs = sorted(a.subject.unique())
        X = np.column_stack([np.ones(len(a)), a.f.values])
        y = a[det].values.astype(float)
        b = firth_logit(X, y)
        rng = np.random.default_rng(seed)
        bysub = {s: a[a.subject == s] for s in subs}
        curves = []
        for _ in range(B):
            pick = rng.choice(subs, len(subs), replace=True)
            aa = pd.concat([bysub[s] for s in pick])
            bb = firth_logit(np.column_stack([np.ones(len(aa)), aa.f.values]), aa[det].values.astype(float))
            curves.append(curve(bb, F_FINE))
        curves = np.array(curves)
        emp = {}
        for f, gf in a.groupby("f"):
            k, n = int(gf[det].sum()), len(gf)
            emp[f"{f:g}"] = {"k": k, "n": n, "p": k / n if n else np.nan, "wilson": wilson(k, n)}
        rep = g[g.design == "replace"]
        out[f"{fam}|{par}"] = {
            "family": fam, "param": par, "lam": float(g.lam.iloc[0]), "beta": b.tolist(), "n_subjects": len(subs),
            "p_hat": curve(b, F_FINE).tolist(), "p_lo95_1s": np.percentile(curves, 5, axis=0).tolist(),
            "p_ci95": [np.percentile(curves, 2.5, axis=0).tolist(), np.percentile(curves, 97.5, axis=0).tolist()],
            "empirical": emp,
            "replacement": {"k": int(rep[det].sum()), "n": len(rep), "wilson": wilson(int(rep[det].sum()), len(rep))}}
    return out


def real_upper(df, det, seed=SEED + 20):
    """One-sided 95 % upper bound on the real window-level detection rate: max of the subject-cluster
    bootstrap 95th percentile of the mean subject rate and the Clopper-Pearson bound on pooled windows."""
    rates, _ = subject_rates(df, det)
    bs = boot_mean(rates, seed=seed)
    k, n = int(df[det].sum()), len(df)
    u_boot = float(np.percentile(bs, 95))
    u_cp = cp_upper(k, n)
    return {"k": k, "n_windows": n, "n_subjects": len(rates), "rate_subject_mean": float(np.mean(rates)),
            "U_boot95": u_boot, "U_cp95": u_cp, "U": max(u_boot, u_cp),
            "U_cp95_subjects": cp_upper(int(np.sum(rates > 0)), len(rates))}


def exclusion(curves, U, L=0.0):
    """PRIMARY (model-free): at each design f, pi_upper = min(1, max(0, U - L) / p_emp(f)) with p_emp = k / n
    (pi_upper = 1 if p_emp = 0); f_min(th) = smallest design f with pi_upper < th at that f and at every larger
    design f.  SECONDARY: the same on the fine grid with the Firth curve p_hat and with its one-sided 95 %
    lower band (the Firth penalty adds pseudo-detections, so these can show 'exclusion' where no detection
    was observed; they are descriptive only)."""
    out = {}
    num = max(0.0, U - L)
    for key, c in curves.items():
        res = {}
        fs = sorted(float(x) for x in c["empirical"])
        pe = [c["empirical"][f"{f:g}"]["p"] for f in fs]
        pi_e = [float(min(1.0, num / p)) if p and p > 0 else 1.0 for p in pe]
        fm = {}
        for th in THRESHOLDS:
            ok = [i for i in range(len(fs)) if all(v < th for v in pi_e[i:])]
            fm[f"{th:g}"] = fs[ok[0]] if ok else None
        res["design"] = {"f": fs, "p_emp": pe, "pi_upper": pi_e, "f_min": fm}
        pw = [c["empirical"][f"{f:g}"]["wilson"][0] for f in fs]
        pi_w = [float(min(1.0, num / p)) if p and p > 0 else 1.0 for p in pw]
        fmw = {}
        for th in THRESHOLDS:
            ok = [i for i in range(len(fs)) if all(v < th for v in pi_w[i:])]
            fmw[f"{th:g}"] = fs[ok[0]] if ok else None
        res["design_wilson_lower"] = {"pi_upper": pi_w, "f_min": fmw}
        for which in ("p_hat", "p_lo95_1s"):
            p = np.asarray(c[which])
            pi = np.minimum(1.0, num / np.maximum(p, 1e-12))
            fm = {}
            for th in THRESHOLDS:
                ok = [i for i in range(len(F_FINE)) if np.all(pi[i:] < th)]
                fm[f"{th:g}"] = float(F_FINE[ok[0]]) if ok else None
            res[which] = {"pi_upper": pi.tolist(), "f_min": fm}
        rp = c["replacement"]
        pr = rp["k"] / rp["n"] if rp["n"] else np.nan
        res["replacement_pi_upper"] = float(min(1.0, num / pr)) if pr and pr > 0 else 1.0
        out[key] = res
    return out


def q1(phase, df_real):
    sp = spike_frame(load(phase, "spike"))
    if sp.empty:
        return {}
    res = {"n_spike_rows": len(sp), "detectors": {}}
    for det, real_col in (("K3", "K3"), ("K1", "K1"), ("LLE", "LLE"), ("UPO", "UPO")):
        curves = fit_detection(sp, det, B=B_CURVE if det in ("K3", "K1") else 500)
        up = real_upper(df_real, real_col)
        r = {"real": up, "curves": curves, "exclusion_L0": exclusion(curves, up["U"])}
        if det in SYNTH_FP:
            k, n = SYNTH_FP[det]
            L = cp_lower(k, n)
            r["L_fp"] = L
            r["exclusion_Lfp"] = exclusion(curves, up["U"], L)
        for grp in ("NSR", "CHF"):
            dg = df_real[df_real.group == grp]
            if len(dg):
                ug = real_upper(dg, real_col, seed=SEED + 21)
                r[f"real_{grp}"] = ug
                r[f"exclusion_{grp}"] = {k: {"f_min_design": v["design"]["f_min"], "f_min_fit": v["p_hat"]["f_min"]}
                                         for k, v in exclusion(curves, ug["U"]).items()}
        res["detectors"][det] = r
    return res


# ------------------------------------------------------------------------------- Question 2 (segments)
def seg_frame(R):
    rows = []
    for r in R:
        if r.get("error"):
            continue
        for s in r["segments"]:
            if s["gap"]:
                rows.append({"subject": r["subject"], "group": r["group"], "seg": s["seg"], "gap": True})
                continue
            rows.append({"subject": r["subject"], "group": r["group"], "seg": s["seg"], "gap": False,
                         "clock_hour": s.get("clock_hour"), "t_start": s.get("t_start", np.nan),
                         "burden": s["burden"], "n_masked": s["n_masked"],
                         "raw": s["TIT_raw"]["positive"], "raw_NL": s["TIT_raw"]["NL"],
                         "masked": s["TIT_masked"].get("positive", False),
                         "masked_ok": bool(s["TIT_masked"].get("analysable", False)),
                         "edited": s["TIT_edited"].get("positive", False),
                         "edited_ok": bool(s["TIT_edited"].get("eligible", False)),
                         "wu": s["TIT_wu"].get("positive", False),
                         "wu_ok": bool(s["TIT_wu"].get("eligible", True) and s["TIT_wu"].get("analysable", True)),
                         "wu_NL": s["TIT_wu"].get("NL", 0.0)})
    df = pd.DataFrame(rows)
    if not df.empty:
        d = ~df["gap"]
        ck = [_clock(h, t) if g else (np.nan, False) for h, t, g in zip(df["clock_hour"], df["t_start"], d)]
        df["clock_eff"] = [c[0] for c in ck]
        df["clock_approx"] = [c[1] for c in ck]
        df["night"] = (df["clock_eff"] >= 0) & (df["clock_eff"] < 5)
    return df


def group_rates(df, col, ok=None):
    d = df if ok is None else df[df[ok]]
    out = {}
    per = {}
    for g in ("NSR", "CHF"):
        dg = d[d.group == g]
        if dg.empty:
            continue
        v, subs = subject_rates(dg, col)
        per[g] = v
        out[g] = {"n_subjects": len(v), "n_segments": len(dg), "mean_subject_rate": float(np.mean(v)),
                  "ci95": ci(boot_mean(v)), "positive_segments": int(dg[col].sum())}
    if "NSR" in per and "CHF" in per:
        out["CHF_minus_NSR"] = {"diff": float(per["CHF"].mean() - per["NSR"].mean()),
                                "ci95": ci(boot_diff(per["CHF"], per["NSR"])),
                                "perm_p": perm_p(per["CHF"], per["NSR"])}
    return out


def paired_change(df, a, b, ok_b, by_group=True):
    """Per subject rate(b) - rate(a) on the segments/windows usable in arm b; cluster bootstrap."""
    d = df[df[ok_b]]
    out = {}
    for g in (("NSR", "CHF") if by_group else ("all",)):
        dg = d if g == "all" else d[d.group == g]
        if dg.empty:
            continue
        ra, _ = subject_rates(dg, a)
        rb, _ = subject_rates(dg, b)
        delta = rb - ra
        both = int((dg[a] & dg[b]).sum())
        out[g] = {"n_subjects": len(delta), "n_units": len(dg), "rate_a": float(ra.mean()), "rate_b": float(rb.mean()),
                  "change": float(delta.mean()), "ci95": ci(boot_mean(delta)),
                  "a_pos": int(dg[a].sum()), "b_pos": int(dg[b].sum()), "both_pos": both,
                  "a_pos_removed_frac": float(1 - both / dg[a].sum()) if dg[a].sum() else None}
    return out


def q2(phase):
    sg = seg_frame(load(phase, "seg"))
    if sg.empty:
        return {}, sg
    d = sg[~sg.gap]
    res = {"n_segments": len(sg), "n_gap": int(sg.gap.sum()), "n_analysable": len(d),
           "P2_raw_by_group": group_rates(d, "raw"),
           "P2_change_masked": paired_change(d, "raw", "masked", "masked_ok"),
           "P2_change_masked_all": paired_change(d, "raw", "masked", "masked_ok", by_group=False),
           "P2_change_masked_nonanalysable_negative": paired_change(d.assign(all_ok=True), "raw", "masked", "all_ok"),
           "change_edited": paired_change(d, "raw", "edited", "edited_ok"),
           "change_wu": paired_change(d, "raw", "wu", "wu_ok"),
           "wu_by_group": group_rates(d, "wu", "wu_ok"),
           "masked_by_group": group_rates(d, "masked", "masked_ok"),
           "edited_by_group": group_rates(d, "edited", "edited_ok"),
           "masked_analysable_frac": {g: float(d[d.group == g].masked_ok.mean()) for g in ("NSR", "CHF") if (d.group == g).any()},
           "edited_eligible_frac": {g: float(d[d.group == g].edited_ok.mean()) for g in ("NSR", "CHF") if (d.group == g).any()},
           "NL_among_positive_raw": {g: float(d[(d.group == g) & d.raw].raw_NL.mean()) for g in ("NSR", "CHF")
                                     if ((d.group == g) & d.raw).any()},
           "NL_among_positive_wu": {g: float(d[(d.group == g) & d.wu].wu_NL.mean()) for g in ("NSR", "CHF")
                                    if ((d.group == g) & d.wu).any()},
           "raw_night_day": {g: {"night": float(d[(d.group == g) & d.night].raw.mean()) if ((d.group == g) & d.night).any() else None,
                                 "day": float(d[(d.group == g) & ~d.night].raw.mean()) if ((d.group == g) & ~d.night).any() else None}
                             for g in ("NSR", "CHF")}}
    res["seg_gee_burden"] = gee(d, "raw", ["x", "chf"]) if not d.empty else None
    return res, sg


# ------------------------------------------------------------------------------- Question 3
def gee(df, outcome, covars):
    import statsmodels.api as sm
    d = df.copy()
    if "x" not in d:
        d["x"] = np.log2(1 + d["burden"])
    if "chf" not in d:
        d["chf"] = (d["group"] == "CHF").astype(int)
    y = d[outcome].astype(float).values
    npos = int(y.sum())
    out = {"outcome": outcome, "covariates": covars, "n": len(d), "n_pos": npos,
           "n_subjects": int(d.subject.nunique())}
    if npos < MIN_POS_GEE or len(d) - npos < MIN_POS_GEE or d["x"].std() == 0:
        out["estimated"] = False
        out["reason"] = f"fewer than {MIN_POS_GEE} positives or negatives (or no burden variation)"
        return out
    X = sm.add_constant(d[covars].astype(float).values)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m = sm.GEE(y, X, groups=d.subject.values, family=sm.families.Binomial(),
                       cov_struct=sm.cov_struct.Exchangeable()).fit()
        names = ["const"] + list(covars)
        out["estimated"] = True
        out["coef"] = {n: {"beta": float(b), "se": float(s), "OR": float(np.exp(b)),
                           "OR_ci95": [float(np.exp(b - 1.96 * s)), float(np.exp(b + 1.96 * s))],
                           "p": float(p)} for n, b, s, p in zip(names, m.params, m.bse, m.pvalues)}
    except Exception as exc:                                    # noqa: BLE001
        out["estimated"] = False
        out["reason"] = f"{type(exc).__name__}: {exc}"
    return out


BURDEN_BINS = [(0, 0), (1, 1), (2, 4), (5, 15), (16, 10 ** 9)]


def burden_table(df, cols):
    rows = []
    for lo, hi in BURDEN_BINS:
        d = df[(df.burden >= lo) & (df.burden <= hi)]
        row = {"burden": f"{lo}" if lo == hi else (f"{lo}-{hi}" if hi < 10 ** 9 else f">={lo}"), "n": len(d),
               "n_subjects": int(d.subject.nunique())}
        for c in cols:
            k = int(d[c].sum())
            row[c] = {"k": k, "rate": k / len(d) if len(d) else None, "wilson": wilson(k, len(d))}
        rows.append(row)
    return rows


def q3(df, covars=("x", "chf")):
    cols = ("TIT", "LLE", "UPO", "K1", "K3")
    covars = list(covars)
    res = {"P3_gee": {c: gee(df, c, covars) for c in cols},
           "gee_unadjusted_burden": {c: gee(df, c, ["x"]) for c in cols},
           "gee_K3_analysable_only": gee(df[df.K3_analysable], "K3", covars),
           "Q3c_titration_group_unadjusted": gee(df, "TIT", ["chf"]) if "chf" in covars else None,
           "burden_table": burden_table(df, cols),
           "burden_ge1_vs_0": {}}
    for c in cols:
        a = df[df.burden >= 1]
        b = df[df.burden == 0]
        if len(a) and len(b):
            ra, _ = subject_rates(a, c)
            rb, _ = subject_rates(b, c)
            res["burden_ge1_vs_0"][c] = {"rate_ge1": float(a[c].mean()), "rate_0": float(b[c].mean()),
                                         "n_ge1": len(a), "n_0": len(b),
                                         "diff_subject_means": float(ra.mean() - rb.mean()),
                                         "ci95": ci(boot_diff(ra, rb))}
    # 3b paired: raw vs masked / edited (windows usable in both arms), all and n_masked >= 1
    pairs = {"TIT_masked": ("TIT", "TITm", "TITm_analysable"), "TIT_edited": ("TIT", "TITe", "TITe_eligible"),
             "K1_edited": ("K1", "K1e", "K1e_eligible"), "LLE_edited": ("LLE", "LLEe", "K1e_eligible"),
             "UPO_edited": ("UPO", "UPOe", "K1e_eligible")}
    res["Q3b_paired"] = {}
    for name, (a, b, ok) in pairs.items():
        res["Q3b_paired"][name] = {"all": paired_change(df, a, b, ok, by_group=False),
                                   "n_masked_ge1": paired_change(df[df.n_masked >= 1], a, b, ok, by_group=False)}
    return res


# ------------------------------------------------------------------------------- Question 4
def agreement(a, b):
    a, b = np.asarray(a, bool), np.asarray(b, bool)
    n11, n10, n01, n00 = int((a & b).sum()), int((a & ~b).sum()), int((~a & b).sum()), int((~a & ~b).sum())
    n = len(a)
    po = (n11 + n00) / n if n else np.nan
    pe = ((n11 + n10) * (n11 + n01) + (n00 + n01) * (n00 + n10)) / n ** 2 if n else np.nan
    kappa = (po - pe) / (1 - pe) if n and pe < 1 else None
    mc = float(stats.binomtest(min(n10, n01), n10 + n01, 0.5).pvalue) if n10 + n01 else 1.0
    return {"n": n, "both": n11, "a_only": n10, "b_only": n01, "neither": n00, "agreement": po,
            "kappa": kappa, "mcnemar_exact_p": mc}


def q4(phase, df):
    res = {}
    rb = load(phase, "robust")
    if rb:
        idx = {(r["subject"], r["j"]): r for r in rb if not r.get("error")}
        m = df.set_index(["subject", "j"])
        a499 = {}
        rows = []
        for (s, j), r in idx.items():
            if (s, j) not in m.index:
                continue
            w = m.loc[(s, j)]
            rows.append({"subject": s, "group": r["group"], "j": j,
                         "K1": w.K1, "LLE": w.LLE, "UPO": w.UPO, "K3": w.K3,
                         "K1_499": bool(r["K1_499"]["and_detected"]), "LLE_499": bool(r["K1_499"]["lle_detected"]),
                         "UPO_499": bool(r["K1_499"]["upo_detected"]), "K3_499": bool(r["K3_499"]["detected"]),
                         **{f"{d}_{L}": (bool(r[f"L{L}"][dd][key]) if r[f"L{L}"].get("available") else None)
                            for L in (256, 1024) for d, dd, key in (("K1", "K1", "and_detected"),
                                                                     ("LLE", "K1", "lle_detected"),
                                                                     ("UPO", "K1", "upo_detected"),
                                                                     ("K3", "K3", "detected"))}})
        R = pd.DataFrame(rows)
        for d in ("K1", "LLE", "UPO", "K3"):
            a499[d] = agreement(R[d], R[f"{d}_499"])
        res["Q4a_499"] = a499
        lens = {}
        for d in ("K1", "LLE", "UPO", "K3"):
            ok = R[f"{d}_256"].notna() & R[f"{d}_1024"].notna()
            Rk = R[ok]
            lens[d] = {"n_windows": int(ok.sum()), **{f"rate_{L}": float(Rk[f"{d}_{L}"].astype(bool).mean()) if len(Rk) else None
                                                     for L in (256, 1024)},
                       "rate_512": float(Rk[d].mean()) if len(Rk) else None,
                       "by_group": {g: {L: float(Rk[Rk.group == g][c].astype(bool).mean()) if (Rk.group == g).any() else None
                                        for L, c in ((256, f"{d}_256"), (512, d), (1024, f"{d}_1024"))}
                                    for g in ("NSR", "CHF")}}
        res["Q4b_lengths"] = lens
    nd = {}
    for d in ("TIT", "LLE", "UPO", "K1", "K3", "K4"):
        nd[d] = {g: {"night": {"k": int(df[(df.group == g) & df.night][d].sum()), "n": int(((df.group == g) & df.night).sum())},
                     "day": {"k": int(df[(df.group == g) & ~df.night][d].sum()), "n": int(((df.group == g) & ~df.night).sum())}}
                 for g in ("NSR", "CHF")}
    res["Q4c_night_day"] = nd
    res["Q4c_titration_gee_night"] = gee(df.assign(night_i=df.night.astype(int)), "TIT", ["x", "chf", "night_i"])
    res["Q4d_labels"] = {"rates": {d: {g: float(df[df.group == g][d].mean()) for g in ("NSR", "CHF") if (df.group == g).any()}
                                   for d in ("K3", "K4", "K3RR")},
                         "K3_vs_K4": agreement(df.K3, df.K4), "K3_vs_K3RR": agreement(df.K3, df.K3RR),
                         "by_burden": burden_table(df, ("K3", "K4", "K3RR"))}
    return res


# ------------------------------------------------------------------------------- synthetic (Q2e)
def q2e(phase):
    R = load(phase, "synth")
    rows = []
    for r in R:
        if r.get("error"):
            continue
        c = r["cond"] if r["kind"] == "null" else f"vdp{r['cond'][0]}_{r['cond'][1]}"
        rows.append({"cond": c, "n": r["n"], "n_ectopic": r["n_ectopic"],
                     **{a: bool(r[f"TIT_{a}"].get("positive", False)) for a in ("raw", "masked", "edited")},
                     "raw_q128": bool(r.get("TIT_raw_q128", {}).get("positive", False)),
                     "masked_ok": bool(r["TIT_masked"].get("analysable", False)),
                     "edited_ok": bool(r["TIT_edited"].get("eligible", False))})
    df = pd.DataFrame(rows)
    out = {}
    if df.empty:
        return out
    for (c, n), g in df.groupby(["cond", "n"]):
        out[f"{c}|{n}"] = {"n_windows": len(g), "raw": int(g.raw.sum()), "raw_q128": int(g.raw_q128.sum()),
                           "masked": int(g.masked.sum()), "masked_analysable": int(g.masked_ok.sum()),
                           "edited": int(g.edited.sum()), "edited_eligible": int(g.edited_ok.sum())}
    return out


def _json_default(o):
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


# ------------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", choices=["dev", "conf"], required=True)
    ap.add_argument("--sensitivity-flagged", action="store_true",
                    help="Amendment 1: drop the subjects flagged by the original overlap rule")
    a = ap.parse_args(argv)
    out = HERE / "results" / a.phase / "analysis"
    if a.sensitivity_flagged:
        EXCLUDE.update(json.load(open(HERE / "results" / "conf" / "exclusions.json"))["originally_flagged_subjects"])
        out = HERE / "results" / a.phase / "analysis_sensitivity_flagged"
    out.mkdir(parents=True, exist_ok=True)
    df = real_frame(load(a.phase, "real"))
    res = {"n_windows": len(df), "n_subjects": int(df.subject.nunique()) if len(df) else 0}
    main_df = df[df.group.isin(["NSR", "CHF"])] if len(df) else df
    res["Q1"] = q1(a.phase, main_df)
    res["Q2"], _ = q2(a.phase)
    res["Q2e"] = q2e(a.phase)
    res["Q3"] = q3(main_df) if len(main_df) else {}
    if a.phase == "dev" and len(df) and (df.group == "MIT").any():
        res["Q3_mitbih_exploratory"] = q3(df[df.group == "MIT"], covars=("x",))
    res["Q4"] = q4(a.phase, main_df) if len(main_df) else {}
    json.dump(res, open(out / "analysis.json", "w"), indent=1, default=_json_default)
    print(json.dumps({k: (v if k in ("n_windows", "n_subjects") else "...") for k, v in res.items()}))


if __name__ == "__main__":
    main()
