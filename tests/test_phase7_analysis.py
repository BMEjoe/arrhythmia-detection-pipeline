"""Phase 7 analysis helpers (experiments/phase7_mitbih/analysis.py); synthetic inputs only."""
import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import roc_auc_score

from experiments.phase7_mitbih import analysis as A


def _toy(seed=0, n_subj=12):
    rng = np.random.default_rng(seed)
    rows = []
    for s in range(n_subj):
        for w in range(int(rng.integers(3, 8))):
            lab = int(rng.random() < 0.4)
            rows.append({"subject": f"S{s}", "label": lab, "and_detected": bool(rng.random() < (0.5 if lab else 0.1)),
                         "x": rng.normal(lab, 1.0), "z": rng.normal(0, 1.0)})
    return pd.DataFrame(rows)


def test_weighted_auc_unit_weights_equals_sklearn():
    df = _toy()
    y, s = df.label.to_numpy(), df.x.to_numpy()
    M = A._pair_matrix(s, y)
    auc = A.weighted_auc(M, np.ones((1, y.sum())), np.ones((1, (1 - y).sum())))[0]
    assert auc == pytest.approx(roc_auc_score(y, s))


def test_weighted_auc_integer_weights_equals_duplicated_data():
    df = _toy(1)
    y, s = df.label.to_numpy(), df.x.to_numpy()
    w = np.random.default_rng(3).integers(0, 4, len(y))
    M = A._pair_matrix(s, y)
    auc = A.weighted_auc(M, w[y == 1][None, :].astype(float), w[y == 0][None, :].astype(float))[0]
    assert auc == pytest.approx(roc_auc_score(np.repeat(y, w), np.repeat(s, w)))


def test_delong_auc_matches_and_variance_positive():
    df = _toy(2)
    y = df.label.to_numpy()
    aucs, S = A.delong([df.x.to_numpy(), df.z.to_numpy()], y)
    assert aucs[0] == pytest.approx(roc_auc_score(y, df.x))
    assert aucs[1] == pytest.approx(roc_auc_score(y, df.z))
    assert S.shape == (2, 2) and S[0, 0] > 0 and S[1, 1] > 0


def test_classification_point_counts():
    df = _toy(4)
    r = A.classification(df)
    tp = int(((df.label == 1) & df.and_detected).sum())
    fp = int(((df.label == 0) & df.and_detected).sum())
    assert (r["TP"], r["FP"]) == (tp, fp)
    assert r["point"]["difference"] == pytest.approx(tp / (df.label == 1).sum() - fp / (df.label == 0).sum())
    lo, hi = r["ci95"]["difference"]
    assert lo <= r["point"]["difference"] <= hi


def test_subject_draws_resample_subjects():
    W = A.subject_draws(["a", "b", "c"], seed=1, b=50)
    assert W.shape == (50, 3) and np.all(W.sum(1) == 3)


def test_loso_never_trains_on_the_test_subject():
    df = _toy(5)
    df["x"] = np.where(df.subject == "S0", 100.0, df.x)        # S0 is an outlier
    oof = A.loso_predictions(df, ["x"])
    assert np.all(np.isfinite(oof))


def test_wilson_bounds():
    lo, hi = A.wilson(0, 300)
    assert lo == pytest.approx(0.0, abs=1e-12) and hi == pytest.approx(0.01264, abs=1e-4)
