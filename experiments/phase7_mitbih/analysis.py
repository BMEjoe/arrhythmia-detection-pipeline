"""
Phase 7 preregistered analyses (PREREGISTRATION.md Section C).

    python -m experiments.phase7_mitbih.analysis              # results/run -> results/analysis
    python -m experiments.phase7_mitbih.analysis --dir results/synthetic_check --out results/synthetic_check/analysis

Every confidence interval is a SUBJECT-cluster percentile bootstrap (10,000
resamples, fixed seeds); subjects, not records or windows, are resampled.
Errors (analyze_segment raised) count as not detected in Q1/Q2 and are excluded
from Q3 (complete cases), with counts reported.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
import pandas as pd
from scipy import stats

HERE = pathlib.Path(__file__).resolve().parent
B = 10_000
SEED = 20261001
ALPHA = 0.05
Q2_BOUND = 0.064            # Phase 6 ectopy-only AND upper 95 % bound (E5 atrial 10 %)
Q2_MAX_POINT = 0.037        # Phase 6 highest ectopy-only AND rate (11/300)
LOGREG_C = 1.0
HRV = ["sdnn_ms", "rmssd_ms", "pnn50"]
MODELS = {"M0": HRV, "M1": ["lle_z"], "M2": ["upo_score"], "M3": ["lle_z", "upo_score"],
          "M4": HRV + ["lle_z", "upo_score"]}
DIFFS = [("M3", "M1"), ("M3", "M2"), ("M4", "M0")]


# ---------------------------------------------------------------- loading
def load(run_dir, arm):
    path = pathlib.Path(run_dir) / f"{arm}.jsonl"
    rows = [json.loads(l) for l in open(path) if l.strip()]
    df = pd.DataFrame(rows).drop_duplicates("task_id", keep="first")
    df["error_flag"] = df["error"].notna()
    for c in ("and_detected", "lle_detected", "upo_detected"):
        df[c] = df[c].fillna(False).astype(bool)
    df["or_detected"] = df.lle_detected | df.upo_detected
    df["record"] = df["record"].astype(str)
    return df.sort_values(["record", "window"]).reset_index(drop=True)


# ---------------------------------------------------------------- bootstrap helpers
def subject_draws(subjects, seed, b=B):
    """(b, n_subjects) multiplicity matrix of a subject-cluster bootstrap."""
    rng = np.random.default_rng(seed)
    k = len(subjects)
    idx = rng.integers(0, k, size=(b, k))
    return np.stack([np.bincount(row, minlength=k) for row in idx])


def pct_ci(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if not len(x):
        return [None, None]
    return [float(np.percentile(x, 2.5)), float(np.percentile(x, 97.5))]


def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return [None, None]
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [float(max(0.0, c - h)), float(min(1.0, c + h))]


# ---------------------------------------------------------------- Q1 / Q2
def classification(df, decision="and_detected", seed=SEED):
    subj = sorted(df.subject.unique())
    si = df.subject.map({s: i for i, s in enumerate(subj)}).to_numpy()
    y = df.label.to_numpy().astype(int)
    d = df[decision].to_numpy().astype(bool)
    k = len(subj)
    # per-subject counts: TP, FN, FP, TN
    cnt = np.zeros((k, 4))
    np.add.at(cnt, (si, 0), (y == 1) & d)
    np.add.at(cnt, (si, 1), (y == 1) & ~d)
    np.add.at(cnt, (si, 2), (y == 0) & d)
    np.add.at(cnt, (si, 3), (y == 0) & ~d)

    def metrics(c):
        tp, fn, fp, tn = (c[..., i] for i in range(4))
        with np.errstate(invalid="ignore", divide="ignore"):
            sens = tp / (tp + fn)
            fpr = fp / (fp + tn)
            m = {"rate_abnormal": sens, "rate_normal": fpr, "difference": sens - fpr,
                 "sensitivity": sens, "specificity": 1 - fpr, "PPV": tp / (tp + fp),
                 "NPV": tn / (tn + fn), "balanced_accuracy": (sens + 1 - fpr) / 2}
        return m

    point = {kk: float(v) if np.isfinite(v) else None for kk, v in metrics(cnt.sum(0)).items()}
    W = subject_draws(subj, seed)
    bc = W @ cnt                                    # (B, 4)
    valid = (bc[:, 0] + bc[:, 1] > 0) & (bc[:, 2] + bc[:, 3] > 0)
    bm = metrics(bc[valid])
    ci = {kk: pct_ci(v) for kk, v in bm.items()}
    undefined = {kk: int(np.sum(~np.isfinite(v))) for kk, v in bm.items()}
    tot = cnt.sum(0)
    return {"decision": decision, "n_windows": int(len(df)), "n_abnormal": int(tot[0] + tot[1]),
            "n_normal": int(tot[2] + tot[3]), "n_subjects": k,
            "n_subjects_with_abnormal": int(np.sum(cnt[:, 0] + cnt[:, 1] > 0)),
            "TP": int(tot[0]), "FN": int(tot[1]), "FP": int(tot[2]), "TN": int(tot[3]),
            "point": point, "ci95": ci, "bootstrap_valid": int(valid.sum()),
            "bootstrap_undefined": undefined,
            "wilson_abnormal": wilson(int(tot[0]), int(tot[0] + tot[1])),
            "wilson_normal": wilson(int(tot[2]), int(tot[2] + tot[3]))}


def q1(df):
    r = classification(df, "and_detected", SEED)
    lo, hi = r["ci95"]["difference"]
    r["Q1_supported"] = bool(lo is not None and lo > 0)
    r["Q1_opposite_direction"] = bool(hi is not None and hi < 0)
    return r


def q2(q1res):
    lo, hi = q1res["ci95"]["rate_abnormal"]
    return {"rate_abnormal": q1res["point"]["rate_abnormal"], "ci95_cluster": [lo, hi],
            "wilson95_ignoring_clustering": q1res["wilson_abnormal"],
            "phase6_ectopy_only_max_rate": Q2_MAX_POINT, "phase6_ectopy_only_upper_bound": Q2_BOUND,
            "exceeds_ectopy_only_expectation": bool(lo is not None and lo > Q2_BOUND)}


# ---------------------------------------------------------------- Q3
def loso_predictions(df, features, C=LOGREG_C):
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    X = df[features].to_numpy(dtype=float)
    y = df.label.to_numpy().astype(int)
    g = df.subject.to_numpy()
    oof = np.full(len(df), np.nan)
    for s in np.unique(g):
        te = g == s
        tr = ~te
        if len(np.unique(y[tr])) < 2:
            continue
        model = make_pipeline(StandardScaler(), LogisticRegression(C=C, solver="lbfgs",  # L2 (sklearn default)
                                                                   max_iter=10_000))
        model.fit(X[tr], y[tr])
        oof[te] = model.predict_proba(X[te])[:, 1]
    return oof


def _pair_matrix(score, y):
    pos, neg = score[y == 1], score[y == 0]
    return (pos[:, None] > neg[None, :]).astype(float) + 0.5 * (pos[:, None] == neg[None, :])


def weighted_auc(M, wp, wn):
    """AUC for bootstrap weights: wp (B, n_pos), wn (B, n_neg)."""
    num = np.einsum("bi,ij,bj->b", wp, M, wn)
    den = wp.sum(1) * wn.sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return num / den


def delong(scores, y):
    """Fast DeLong (Sun & Xu 2014): AUCs and covariance for a list of score vectors."""
    y = np.asarray(y).astype(int)
    pos = [np.asarray(s)[y == 1] for s in scores]
    neg = [np.asarray(s)[y == 0] for s in scores]
    m, n = len(pos[0]), len(neg[0])
    aucs, v10, v01 = [], [], []
    for p, q in zip(pos, neg):
        M = (p[:, None] > q[None, :]) + 0.5 * (p[:, None] == q[None, :])
        aucs.append(M.mean())
        v10.append(M.mean(1))
        v01.append(M.mean(0))
    v10, v01 = np.array(v10), np.array(v01)
    S = np.cov(v10) / m + np.cov(v01) / n
    return np.array(aucs), np.atleast_2d(S)


def q3(df, seed=SEED + 3):
    feats = sorted(set(sum(MODELS.values(), [])))
    usable = (~df.error_flag) & np.isfinite(df[feats].astype(float)).all(axis=1)
    cc = df[usable].reset_index(drop=True)
    excluded = df[~usable]
    y = cc.label.to_numpy().astype(int)
    preds = {name: loso_predictions(cc, f) for name, f in MODELS.items()}
    subj = sorted(cc.subject.unique())
    si = cc.subject.map({s: i for i, s in enumerate(subj)}).to_numpy()
    W = subject_draws(subj, seed)                   # (B, k)
    ww = W[:, si]                                   # window weights
    wp, wn = ww[:, y == 1], ww[:, y == 0]
    valid = (wp.sum(1) > 0) & (wn.sum(1) > 0)
    from sklearn.metrics import roc_auc_score
    res = {"n_windows": int(len(cc)), "n_abnormal": int(y.sum()), "n_normal": int((1 - y).sum()),
           "n_subjects": len(subj), "excluded_windows": int(len(excluded)),
           "excluded_by_label": {int(k): int(v) for k, v in excluded.label.value_counts().items()},
           "C": LOGREG_C, "bootstrap_valid": int(valid.sum()), "models": {}, "differences": {}}
    boots = {}
    for name, p in preds.items():
        ok = np.isfinite(p)
        auc = float(roc_auc_score(y[ok], p[ok]))
        M = _pair_matrix(p, y)
        bt = weighted_auc(M, wp[valid], wn[valid])
        boots[name] = bt
        a, S = delong([p], y)
        se = float(np.sqrt(S[0, 0]))
        res["models"][name] = {"features": MODELS[name], "auc": auc, "ci95_cluster": pct_ci(bt),
                               "delong_ci95": [auc - 1.959964 * se, auc + 1.959964 * se],
                               "n_oof_missing": int((~ok).sum())}
    for a_, b_ in DIFFS:
        d = boots[a_] - boots[b_]
        aucs, S = delong([preds[a_], preds[b_]], y)
        var = S[0, 0] + S[1, 1] - 2 * S[0, 1]
        z = (aucs[0] - aucs[1]) / np.sqrt(var) if var > 0 else np.nan
        res["differences"][f"{a_}-{b_}"] = {
            "difference": float(res["models"][a_]["auc"] - res["models"][b_]["auc"]),
            "ci95_cluster_paired": pct_ci(d),
            "delong_z": None if not np.isfinite(z) else float(z),
            "delong_p_two_sided": None if not np.isfinite(z) else float(2 * stats.norm.sf(abs(z)))}
    oof = cc[["task_id", "record", "subject", "window", "label"]].copy()
    for name, p in preds.items():
        oof[f"oof_{name}"] = p
    return res, oof


# ---------------------------------------------------------------- descriptive
def rate_table(df, by, decisions=("and_detected", "lle_detected", "upo_detected", "or_detected")):
    rows = []
    for key, g in df.groupby(by, dropna=False):
        row = {by if isinstance(by, str) else "group": key, "n": len(g), "subjects": g.subject.nunique(),
               "errors": int(g.error_flag.sum())}
        for d in decisions:
            k = int(g[d].sum())
            lo, hi = wilson(k, len(g))
            row[d] = f"{k}/{len(g)} ({100 * k / len(g):.1f} %; {100 * lo:.1f}–{100 * hi:.1f})"
        rows.append(row)
    return pd.DataFrame(rows)


def per_record(df):
    rows = []
    for r, g in df.groupby("record"):
        n, a = g[g.label == 0], g[g.label == 1]
        rows.append({"record": r, "subject": g.subject.iloc[0], "windows N": len(n), "windows A": len(a),
                     "AND N": int(n.and_detected.sum()), "AND A": int(a.and_detected.sum()),
                     "LLE N": int(n.lle_detected.sum()), "LLE A": int(a.lle_detected.sum()),
                     "UPO N": int(n.upo_detected.sum()), "UPO A": int(a.upo_detected.sum()),
                     "errors": int(g.error_flag.sum()),
                     "median lle_z": float(np.nanmedian(g.lle_z.astype(float))) if g.lle_z.notna().any() else None,
                     "median upo_score": float(np.nanmedian(g.upo_score.astype(float))) if g.upo_score.notna().any() else None})
    return pd.DataFrame(rows)


def descriptives(df, arm):
    out = {}
    out["by_label"] = rate_table(df, "label")
    ab = df[df.label == 1]
    if len(ab):
        out["abnormal_by_beat_type"] = rate_table(ab, "beat_type")
    af_records = {"201", "202", "203", "210", "217", "219", "221", "222"}
    df = df.copy()
    df["af_record"] = df.record.isin(af_records)
    df["af_overlap"] = df.af_fraction.fillna(0) > 0
    df["label_af_record"] = df.label.astype(str) + "/" + np.where(df.af_record, "AF record", "other")
    out["by_label_and_af_record"] = rate_table(df, "label_af_record")
    df["label_af_overlap"] = df.label.astype(str) + "/" + np.where(df.af_overlap, "AF overlap", "no AF")
    out["by_label_and_af_overlap"] = rate_table(df, "label_af_overlap")
    if arm != "annotation":
        df["label_missed"] = df.label.astype(str) + "/" + np.where(df.n_missed_beats_in_span.fillna(0) > 0,
                                                                   "missed beat", "no missed beat")
        out["by_label_and_missed_beat"] = rate_table(df, "label_missed")
    out["per_record"] = per_record(df)
    return out


def fmt_ci(p, ci, pct=True):
    if p is None:
        return "undefined"
    f = (lambda v: f"{100 * v:.1f} %") if pct else (lambda v: f"{v:.3f}")
    lo, hi = ci
    return f"{f(p)} [{f(lo) if lo is not None else 'NA'}, {f(hi) if hi is not None else 'NA'}]"


def q1_markdown(r, title):
    p, c = r["point"], r["ci95"]
    lines = [f"### {title}", "",
             f"Windows {r['n_windows']} ({r['n_normal']} normal, {r['n_abnormal']} abnormal) from "
             f"{r['n_subjects']} subjects ({r['n_subjects_with_abnormal']} with abnormal windows). "
             f"TP {r['TP']}, FN {r['FN']}, FP {r['FP']}, TN {r['TN']}. "
             f"Valid bootstrap resamples {r['bootstrap_valid']}/{B}.", "",
             "| quantity | estimate [95 % subject-cluster bootstrap CI] |", "|---|---|"]
    for k in ("rate_abnormal", "rate_normal", "difference", "sensitivity", "specificity", "PPV", "NPV",
              "balanced_accuracy"):
        note = f" (undefined in {r['bootstrap_undefined'][k]} resamples)" if r["bootstrap_undefined"][k] else ""
        lines.append(f"| {k} | {fmt_ci(p[k], c[k])}{note} |")
    if "Q1_supported" in r:
        lines += ["", f"**Q1 supported (lower bound of the difference > 0): {r['Q1_supported']}.**"
                  + (" The CI lies entirely below 0." if r["Q1_opposite_direction"] else "")]
    return "\n".join(lines)


def q3_markdown(r, title):
    lines = [f"### {title}", "",
             f"Complete-case windows {r['n_windows']} ({r['n_normal']} normal, {r['n_abnormal']} abnormal), "
             f"{r['n_subjects']} subjects; excluded {r['excluded_windows']} {r['excluded_by_label']}. "
             f"L2 logistic regression, C = {r['C']}, standardization fit in each training fold, "
             f"leave-one-subject-out. Valid bootstrap resamples {r['bootstrap_valid']}/{B}.", "",
             "| model | features | out-of-fold AUC [95 % subject-cluster bootstrap] | DeLong 95 % CI (secondary) |",
             "|---|---|---|---|"]
    for name, m in r["models"].items():
        lo, hi = m["ci95_cluster"]
        dl = m["delong_ci95"]
        lines.append(f"| {name} | {', '.join(m['features'])} | {m['auc']:.3f} [{lo:.3f}, {hi:.3f}] | "
                     f"[{dl[0]:.3f}, {dl[1]:.3f}] |")
    lines += ["", "| difference | ΔAUC [95 % paired subject-cluster bootstrap] | DeLong z, p (secondary) |",
              "|---|---|---|"]
    for k, d in r["differences"].items():
        lo, hi = d["ci95_cluster_paired"]
        z = "NA" if d["delong_z"] is None else f"{d['delong_z']:.2f}, {d['delong_p_two_sided']:.3g}"
        lines.append(f"| {k} | {d['difference']:+.3f} [{lo:+.3f}, {hi:+.3f}] | {z} |")
    return "\n".join(lines)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=str(HERE / "results" / "run"))
    ap.add_argument("--out", default=str(HERE / "results" / "analysis"))
    a = ap.parse_args(argv)
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    allres, md = {}, ["# Phase 7 preregistered analyses (generated by analysis.py)", ""]
    arms = [("raw", "PRIMARY: raw RR (detected peaks, 150 ms labels)"),
            ("annotation", "SENSITIVITY: annotation beat times"),
            ("edited", "SECONDARY: annotation-edited NN (NN fraction >= 0.80)")]
    for arm, title in arms:
        df = load(a.dir, arm)
        res = {"n_errors": int(df.error_flag.sum()),
               "errors_by_label": {int(k): int(v) for k, v in df[df.error_flag].label.value_counts().items()},
               "error_messages": df[df.error_flag].error.value_counts().to_dict(),
               "upo_score_inconsistent": int((df.get("upo_score_consistent", pd.Series(True)) == False).sum())}
        res["Q1"] = q1(df)
        md += [f"## {title}", "", f"Errors: {res['n_errors']} {res['errors_by_label']}.", "",
               q1_markdown(res["Q1"], "Q1: AND detection, abnormal vs normal windows"), ""]
        if arm in ("raw", "annotation"):
            res["Q2"] = q2(res["Q1"])
            r2 = res["Q2"]
            lo, hi = r2["ci95_cluster"]
            md += ["### Q2: abnormal-window AND rate vs the Phase 6 ectopy-only expectation", "",
                   f"Rate {fmt_ci(r2['rate_abnormal'], r2['ci95_cluster'])} (Wilson, ignoring clustering: "
                   f"{100 * r2['wilson95_ignoring_clustering'][0]:.1f}–{100 * r2['wilson95_ignoring_clustering'][1]:.1f} %). "
                   f"Exceeds the ectopy-only expectation (lower bound > 6.4 %): "
                   f"**{r2['exceeds_ectopy_only_expectation']}**.", ""]
        for dec in ("lle_detected", "upo_detected", "or_detected"):
            res[f"classification_{dec}"] = classification(df, dec, SEED + 7)
        r3, oof = q3(df)
        res["Q3"] = r3
        oof.to_csv(out / f"{arm}_oof_predictions.csv", index=False)
        md += [q3_markdown(r3, "Q3: leave-one-subject-out logistic regression"), ""]
        desc = descriptives(df, arm)
        for name, tab in desc.items():
            tab.to_csv(out / f"{arm}_{name}.csv", index=False)
            md += [f"### Descriptive: {name}", "", tab.to_markdown(index=False), ""]
        if arm == "raw":
            sub = df[df.retained_tol75.fillna(False).astype(bool)].copy()
            sub["label"] = sub.label_tol75.astype(int)
            res["tol75_Q1_descriptive"] = classification(sub, "and_detected", SEED + 11)
            md += [q1_markdown(res["tol75_Q1_descriptive"], "Descriptive: 75 ms labels (same detector runs)"), ""]
            nm = df[df.n_missed_beats_in_span.fillna(0) == 0]
            res["no_missed_beat_Q1_descriptive"] = classification(nm, "and_detected", SEED + 13)
            md += [q1_markdown(res["no_missed_beat_Q1_descriptive"],
                               "Descriptive: windows without a missed beat"), ""]
        allres[arm] = res
    json.dump(allres, open(out / "results.json", "w"), indent=1, default=lambda o: None)
    open(out / "results.md", "w").write("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
