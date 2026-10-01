"""
Phase 7 EXPLORATORY analyses.  Everything here was chosen AFTER the preregistered
detector outputs were seen (results/analysis/).  Nothing here is confirmatory.

    python -m experiments.phase7_mitbih.exploratory

E1  Q3 out-of-fold AUCs below 0.5 (M2 raw 0.19): is pooled leave-one-subject-out
    prediction biased when labels cluster within subjects?  Controls: intercept-only
    and pure-noise feature models through the same analysis.loso_predictions; and
    fit-free AUCs of the single scores (no model, so no fold effect).
E2  The 6 AND detections of the annotation-time arm, and the raw-RR windows that
    cover the same time span.
E3  Raw vs annotation windows paired by time overlap (>= 90 % of the raw window).
"""
from __future__ import annotations

import json
import pathlib

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from experiments.phase7_mitbih import analysis as A

HERE = pathlib.Path(__file__).resolve().parent
RUN = HERE / "results" / "run"
OUT = HERE / "results" / "exploratory"


def direct_auc(df, col, seed):
    """Fit-free AUC of a score (higher = abnormal) with a subject-cluster bootstrap."""
    d = df[np.isfinite(df[col].astype(float))]
    y = d.label.to_numpy().astype(int)
    s = d[col].to_numpy(dtype=float)
    subj = sorted(d.subject.unique())
    si = d.subject.map({k: i for i, k in enumerate(subj)}).to_numpy()
    W = A.subject_draws(subj, seed)[:, si]
    wp, wn = W[:, y == 1], W[:, y == 0]
    ok = (wp.sum(1) > 0) & (wn.sum(1) > 0)
    bt = A.weighted_auc(A._pair_matrix(s, y), wp[ok], wn[ok])
    return {"auc": float(roc_auc_score(y, s)), "ci95_cluster": A.pct_ci(bt), "n": int(len(d))}


def e1(df, arm):
    res = {}
    rng = np.random.default_rng(7)
    noise_aucs = []
    d = df.copy()
    for k in range(20):
        d["noise"] = rng.normal(size=len(d))
        p = A.loso_predictions(d, ["noise"])
        noise_aucs.append(roc_auc_score(d.label, p))
    res["noise_feature_loso_auc"] = {"median": float(np.median(noise_aucs)), "min": float(np.min(noise_aucs)),
                                     "max": float(np.max(noise_aucs)), "n_repeats": 20}
    d["const"] = 0.0
    p = A.loso_predictions(d, ["const"])                # intercept only (scaler maps 0 -> 0)
    res["intercept_only_loso_auc"] = float(roc_auc_score(d.label, p))
    # held-out subject's abnormal share vs the training prevalence (the fitted intercept)
    sub = d.groupby("subject").label.agg(["mean", "size"])
    n1, n = d.label.sum(), len(d)
    prev_tr = (n1 - sub["mean"] * sub["size"]) / (n - sub["size"])
    res["corr_subject_abnormal_share_vs_training_prevalence"] = float(np.corrcoef(sub["mean"], prev_tr)[0, 1])
    for col in ("lle_z", "upo_score", "sdnn_ms", "rmssd_ms", "pnn50", "lle_stat", "upo_source_rJ"):
        if col in d:
            res[f"direct_auc_{col}"] = direct_auc(d, col, 20261101)
    # UPO score distribution
    res["upo_score_zero_fraction_by_label"] = {int(k): float(v) for k, v in
                                               (d.upo_score == 0).groupby(d.label).mean().items()}
    res["upo_score_median_by_label"] = {int(k): float(v) for k, v in d.upo_score.groupby(d.label).median().items()}
    res["lle_z_median_by_label"] = {int(k): float(v) for k, v in d.lle_z.groupby(d.label).median().items()}
    return res


def overlap_pairs(raw, ann):
    rows = []
    for r, g in raw.groupby("record"):
        a = ann[ann.record == r]
        for _, x in g.iterrows():
            ov = np.minimum(x.t1, a.t1) - np.maximum(x.t0, a.t0)
            if not len(a):
                continue
            j = int(np.argmax(ov.to_numpy()))
            frac = float(ov.iloc[j]) / float(x.t1 - x.t0)
            if frac >= 0.90:
                y = a.iloc[j]
                rows.append({"record": r, "subject": x.subject, "raw_window": int(x.window),
                             "ann_window": int(y.window), "overlap": frac, "label_raw": int(x.label),
                             "label_ann": int(y.label),
                             **{f"{k}_raw": x[k] for k in ("and_detected", "lle_detected", "upo_detected",
                                                            "lle_z", "upo_score", "n_missed_beats_in_span")},
                             **{f"{k}_ann": y[k] for k in ("and_detected", "lle_detected", "upo_detected",
                                                            "lle_z", "upo_score")}})
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    out = {}
    md = ["# Phase 7 EXPLORATORY analyses (chosen after the preregistered outputs were seen)", ""]
    arms = {a: A.load(RUN, a) for a in ("raw", "annotation", "edited")}
    for arm in ("raw", "annotation", "edited"):
        out[f"E1_{arm}"] = e1(arms[arm], arm)
    md += ["## E1: out-of-fold AUC artifact and fit-free score AUCs", "", "```",
           json.dumps({k: v for k, v in out.items() if k.startswith("E1")}, indent=1), "```", ""]
    raw, ann = arms["raw"], arms["annotation"]
    det = ann[ann.and_detected]
    cols = ["record", "subject", "window", "label", "abnormal_fraction", "beat_type", "nn_fraction", "lle_p",
            "lle_z", "upo_score", "upo_n_gated", "upo_max_gated_modulus", "detrend_applied", "af_fraction"]
    det[cols].to_csv(OUT / "E2_annotation_detections.csv", index=False)
    pairs = overlap_pairs(raw, ann)
    pairs.to_csv(OUT / "E3_raw_vs_annotation_pairs.csv", index=False)
    e2 = []
    for _, x in det.iterrows():
        p = pairs[(pairs.record == x.record) & (pairs.ann_window == x.window)]
        e2.append({"record": x.record, "ann_window": int(x.window),
                   "raw_pair": None if not len(p) else p.iloc[0][["raw_window", "overlap", "and_detected_raw",
                                                                 "lle_detected_raw", "upo_detected_raw",
                                                                 "upo_score_raw", "upo_score_ann"]].to_dict()})
    out["E2"] = e2
    md += ["## E2: annotation-time AND detections", "", det[cols].to_markdown(index=False), "",
           "Raw-RR windows covering the same span (>= 90 % overlap):", "", "```",
           json.dumps(e2, indent=1, default=str), "```", ""]
    agree = {}
    for k in ("and_detected", "lle_detected", "upo_detected"):
        agree[k] = pd.crosstab(pairs[f"{k}_raw"], pairs[f"{k}_ann"]).to_dict()
    out["E3"] = {"n_pairs": int(len(pairs)), "label_agreement": float((pairs.label_raw == pairs.label_ann).mean()),
                 "crosstabs_raw_vs_ann": {k: {str(a): {str(b): int(c) for b, c in v.items()} for a, v in t.items()}
                                          for k, t in agree.items()},
                 "corr_lle_z": float(np.corrcoef(pairs.lle_z_raw, pairs.lle_z_ann)[0, 1]),
                 "corr_upo_score": float(np.corrcoef(pairs.upo_score_raw, pairs.upo_score_ann)[0, 1])}
    for lab in (0, 1):
        g = pairs[pairs.label_raw == lab]
        out["E3"][f"label{lab}"] = {"n": int(len(g)),
                                    "upo_raw": int(g.upo_detected_raw.sum()), "upo_ann": int(g.upo_detected_ann.sum()),
                                    "lle_raw": int(g.lle_detected_raw.sum()), "lle_ann": int(g.lle_detected_ann.sum()),
                                    "and_raw": int(g.and_detected_raw.sum()), "and_ann": int(g.and_detected_ann.sum()),
                                    "median_upo_score_raw": float(g.upo_score_raw.median()),
                                    "median_upo_score_ann": float(g.upo_score_ann.median())}
    md += ["## E3: raw vs annotation windows paired by time", "", "```", json.dumps(out["E3"], indent=1), "```", ""]
    json.dump(out, open(OUT / "exploratory.json", "w"), indent=1, default=str)
    open(OUT / "exploratory.md", "w").write("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
