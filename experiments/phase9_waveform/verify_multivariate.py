"""
Phase 9 Part C2: verification of the multivariate IAAFT (Schreiber & Schmitz 2000, Sect. 4.6, Fig. 9)
on synthetic data with known structure (no Phase 9 series).

    python -m experiments.phase9_waveform.verify_multivariate

Published behaviour: the multivariate surrogate keeps the individual distributions exactly and the
individual spectra and the cross-correlation function (approximately), whereas surrogates made for
each channel separately destroy the cross-correlation (Fig. 9).
Checks (predeclared):
 V1 bivariate VAR(1) x_t = 0.6 x_{t-1} + 0.3 y_{t-1} + e1, y_t = 0.5 y_{t-1} + 0.4 x_{t-1} + e2
    (N = 512), channel 2 transformed by exp(): distributions identical (max |diff| <= 1e-12, machine
    precision after standardize / restore; first run used == 0 and failed at 1.8e-15); cross-
    correlation at lags -10..10: mean max |CCF_sur - CCF_data| over 20 surrogates <= 0.1, and at least
    3 x smaller than for per-channel IAAFT (fp.iaaft_surrogate on each channel).
 V2 size of the multivariate prediction test (mv_nlp_error, m = 2, target channel 1, 19 surrogates,
    one-sided 5 %) on 100 VAR(1) realizations with the exp() transform: rejection rate <= 10 %.
"""
from __future__ import annotations

import json
import pathlib

import numpy as np

from experiments.phase9_waveform import multivariate as MV

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "results" / "verification"


def var1(n, seed):
    rng = np.random.default_rng(seed)
    x = np.zeros((n + 200, 2))
    for t in range(1, n + 200):
        x[t, 0] = 0.6 * x[t - 1, 0] + 0.3 * x[t - 1, 1] + rng.standard_normal()
        x[t, 1] = 0.5 * x[t - 1, 1] + 0.4 * x[t - 1, 0] + rng.standard_normal()
    x = x[200:]
    x[:, 1] = np.exp(0.5 * x[:, 1])
    return x


def ccf(a, b, lags=10):
    a = (a - a.mean()) / a.std()
    b = (b - b.mean()) / b.std()
    n = len(a)
    return np.array([np.mean(a[max(0, -l):n - max(0, l)] * b[max(0, l):n - max(0, -l)]) for l in range(-lags, lags + 1)])


def main():
    import final_pipeline as fp
    res = {}
    X = var1(512, 1)
    rng = np.random.default_rng(2)
    c0 = ccf(X[:, 0], X[:, 1])
    dmv, duni, dist = [], [], []
    for _ in range(20):
        S = MV.mv_iaaft(X, rng)
        dist.append(float(np.max(np.abs(np.sort(S, 0) - np.sort(X, 0)))))
        dmv.append(float(np.max(np.abs(ccf(S[:, 0], S[:, 1]) - c0))))
        U = np.stack([fp.iaaft_surrogate(X[:, m], rng) for m in range(2)], axis=1)
        duni.append(float(np.max(np.abs(ccf(U[:, 0], U[:, 1]) - c0))))
    res["V1"] = {"max_distribution_diff": max(dist), "ccf_maxdiff_mv_mean": float(np.mean(dmv)),
                 "ccf_maxdiff_univariate_mean": float(np.mean(duni)), "ccf_data_lag0": float(c0[10])}
    rej = 0
    for s in range(100):
        Y = var1(256, 100 + s)
        r2 = np.random.default_rng(1000 + s)
        obs = MV.mv_nlp_error(Y, m=2)
        sur = [MV.mv_nlp_error(MV.mv_iaaft(Y, r2), m=2) for _ in range(19)]
        rej += (1 + np.sum(np.array(sur) <= obs)) / 20 <= 0.05
    res["V2_rejections_of_100"] = int(rej)
    res["VERIFIED"] = bool(max(dist) <= 1e-12 and np.mean(dmv) <= 0.1 and np.mean(duni) >= 3 * np.mean(dmv) and rej <= 10)
    OUT.mkdir(parents=True, exist_ok=True)
    json.dump(res, open(OUT / "multivariate_iaaft.json", "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
