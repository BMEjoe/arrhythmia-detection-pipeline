"""
Phase 9 Part B2: verification of the waveform surrogates on systems with published behaviour (no
Phase 9 series).  Results: results/verification/surrogates.json.

    python -m experiments.phase9_waveform.verify_surrogates --workers 4

PPS and CS (Small et al. 2001; Luo, Nakamura & Small, arXiv:nlin/0404054, Sect. III): Rossler
x' = -y - z, y' = x + a y, z' = b + z (x - c), b = 2, c = 4, sampling 0.1, 1,000 transient points
discarded:
  P-white  a = 0.39095 (period-6 limit cycle) + 5 % white Gaussian observational noise: the PPS null
           holds -> rejection rate near the nominal level;
  P-AR1    same with 5 % AR(1) observational noise xi_{i+1} = 0.8 xi_i + e_i (the paper's coloured
           case): PPS "may incorrectly reject the null hypothesis if the intercycles ... have a linear
           stochastic dependence induced by colored additive observational noise" (Luo et al., citing
           Small et al.) -> rejections expected (published limitation);
  C        chaotic Rossler a = 0.398 (LE computed here) + 5 % white noise -> rejection.
  Statistic: correlation dimension (GP slope, C in [1e-3, 5e-2]), two-sided rank test, Theiler 60 samples;
  embedding by the predeclared rule; 19 surrogates, one-sided alpha = 0.05; 10 realizations of 5,000
  points.  Criterion: P-white <= 2/10 rejections and C >= 8/10 rejections, for PPS and CS (CS cut at
  the local maxima of x, one per cycle).
TS (Thiel et al. 2008, Sect. 5): two coupled Rossler oscillators
  x1,2' = -(1 +/- nu) y1,2 - z1,2 + eps (x2,1 - x1,2), y1,2' = (1 +/- nu) x1,2 + 0.15 y1,2,
  z1,2' = 0.2 + z1,2 (x1,2 - 10), nu = 0.015;  phase synchronization index R_L = |<exp(i dPhi)>| with
  Hilbert phases of x1, x2; 100 TS of the 6-D state; reject 'no PS' iff R_L >= (1 - alpha) max_i R_L,si
  with alpha = 0.01 (their rule).  Published: eps = 0.02 (no PS): 0/100 false rejections; eps = 0.045
  (PS): 100/100 rejections.  Criterion: eps = 0.02 <= 1/10 and eps = 0.045 >= 9/10 rejections.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import pathlib
from multiprocessing import Pool

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402

from experiments.phase9_waveform import surrogates_wave as SW  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "results" / "verification"


def rossler3(n, a, b=2.0, c=4.0, dt=0.1, sub=20, transient=1000, seed=0):
    rng = np.random.default_rng(seed)
    v = np.array([0.1, 0.1, 0.1]) + rng.normal(0, 0.01, 3)
    h = dt / sub

    def f(u):
        return np.array([-u[1] - u[2], u[0] + a * u[1], b + u[2] * (u[0] - c)])
    out = np.empty((n, 3))
    for k in range(transient + n):
        for _ in range(sub):
            k1 = f(v); k2 = f(v + 0.5 * h * k1); k3 = f(v + 0.5 * h * k2); k4 = f(v + h * k3)
            v = v + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        if k >= transient:
            out[k - transient] = v
    return out


def rossler_le(a, b=2.0, c=4.0, T=3000.0, h=0.01):
    v = np.array([0.1, 0.1, 0.1]); w = np.array([1.0, 0, 0])

    def f(u):
        return np.array([-u[1] - u[2], u[0] + a * u[1], b + u[2] * (u[0] - c)])

    def J(u):
        return np.array([[0, -1, -1], [1, a, 0], [u[2], 0, u[0] - c]])
    for _ in range(int(200 / h)):
        k1 = f(v); k2 = f(v + 0.5 * h * k1); k3 = f(v + 0.5 * h * k2); k4 = f(v + h * k3)
        v = v + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    acc = 0.0
    n = int(T / h)
    for k in range(n):
        a1 = f(v); b1 = J(v) @ w
        a2 = f(v + .5 * h * a1); b2 = J(v + .5 * h * a1) @ (w + .5 * h * b1)
        a3 = f(v + .5 * h * a2); b3 = J(v + .5 * h * a2) @ (w + .5 * h * b2)
        a4 = f(v + h * a3); b4 = J(v + h * a3) @ (w + h * b3)
        v = v + h / 6 * (a1 + 2 * a2 + 2 * a3 + a4)
        w = w + h / 6 * (b1 + 2 * b2 + 2 * b3 + b4)
        if k % 50 == 49:
            nw = np.linalg.norm(w); acc += math.log(nw); w /= nw
    return acc / (n * h)


def obs(kind, seed, n=5000):
    rng = np.random.default_rng(1000 + seed)
    if kind == "C":
        x = rossler3(n, 0.398, seed=seed)[:, 0]
    else:
        x = rossler3(n, 0.39095, seed=seed)[:, 0]
    if kind == "P-AR1":
        xi = np.zeros(n)
        e = rng.standard_normal(n)
        for i in range(1, n):
            xi[i] = 0.8 * xi[i - 1] + e[i]
        noise = xi / xi.std()
    else:
        noise = rng.standard_normal(n)
    return x + 0.05 * x.std() * noise


def one_ppscs(task):
    kind, seed = task
    x = obs(kind, seed)
    rng = np.random.default_rng(seed)
    tau, m = SW.choose_embedding(x)
    th = 60
    # Statistic for verification (revised after a smoke test, before any full run): the published
    # tests used the correlation dimension; one-step prediction error rejects even periodic + white
    # noise (PPS jumps and CS junctions make every surrogate rougher at one-sample scale).
    stat = lambda y: SW.corr_dim(y, m, tau, th)
    # PPS
    sd = np.std(x)
    rhos = sd * np.array([0.02, 0.04, 0.08, 0.16, 0.32, 0.64])
    rho, counts = SW.noise_radius(x, m, tau, rhos, rng)
    z = SW.embed_fwd(x, m, tau)
    obs_pps = stat(z[:, 0])
    sur = np.array([stat(SW.pps(x, m, tau, rho, rng)) for _ in range(19)])
    p_pps = (1 + min(np.sum(sur <= obs_pps), np.sum(sur >= obs_pps))) / 20 * 2   # two-sided
    # CS: cut at local maxima of x above its median (one per cycle)
    from scipy.signal import find_peaks
    pk, _ = find_peaks(x, height=np.median(x), distance=30)
    s0, dat = SW.cycle_shuffle(x, pk, rng)
    obs_cs = stat(dat)
    sc = np.array([stat(SW.cycle_shuffle(x, pk, rng)[0]) for _ in range(19)])
    p_cs = (1 + min(np.sum(sc <= obs_cs), np.sum(sc >= obs_cs))) / 20 * 2
    return {"kind": kind, "seed": seed, "tau": tau, "m": m, "rho_over_sd": rho / sd, "p_pps": float(p_pps),
            "p_cs": float(p_cs), "D2": obs_pps, "D2_pps": [float(np.min(sur)), float(np.max(sur))],
            "D2_cs_data": obs_cs, "D2_cs": [float(np.min(sc)), float(np.max(sc))]}


def coupled_rossler(n, eps, nu=0.015, dt=0.2, sub=40, transient=500, seed=0):
    rng = np.random.default_rng(seed)
    v = rng.uniform(-1, 1, 6) + np.array([1, 1, 0, -1, 1, 0])
    h = dt / sub

    def f(u):
        x1, y1, z1, x2, y2, z2 = u
        return np.array([-(1 + nu) * y1 - z1 + eps * (x2 - x1), (1 + nu) * x1 + 0.15 * y1, 0.2 + z1 * (x1 - 10),
                         -(1 - nu) * y2 - z2 + eps * (x1 - x2), (1 - nu) * x2 + 0.15 * y2, 0.2 + z2 * (x2 - 10)])
    out = np.empty((n, 6))
    for k in range(transient + n):
        for _ in range(sub):
            k1 = f(v); k2 = f(v + 0.5 * h * k1); k3 = f(v + 0.5 * h * k2); k4 = f(v + h * k3)
            v = v + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        if k >= transient:
            out[k - transient] = v
    return out


def one_ts(task):
    from scipy.signal import hilbert
    eps, seed = task
    X = coupled_rossler(3000, eps, seed=seed)
    Z = (X - X.mean(0)) / X.std(0)
    group, delta = SW.twin_sets(Z, rr=0.05)
    rng = np.random.default_rng(seed)
    ph1 = np.unwrap(np.angle(hilbert(X[:, 0] - X[:, 0].mean())))
    ph2 = np.unwrap(np.angle(hilbert(X[:, 3] - X[:, 3].mean())))
    RL = abs(np.mean(np.exp(1j * (ph1 - ph2))))
    rs = []
    for _ in range(100):
        idx = SW.twin_surrogate(Z, group, rng, length=len(Z))
        x2s = X[idx, 3]
        p2s = np.unwrap(np.angle(hilbert(x2s - x2s.mean())))
        rs.append(abs(np.mean(np.exp(1j * (ph1 - p2s)))))
    n_twins = int(sum(1 for g in group if len(g) > 1))
    return {"eps": eps, "seed": seed, "RL": float(RL), "max_RL_sur": float(max(rs)), "n_points_with_twins": n_twins,
            "reject": bool(RL >= 0.99 * max(rs))}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    res = {"rossler_LE_a0.398": rossler_le(0.398), "rossler_LE_a0.39095": rossler_le(0.39095)}
    with Pool(a.workers) as pool:
        r1 = list(pool.imap_unordered(one_ppscs, [(k, s) for k in ("P-white", "P-AR1", "C") for s in range(10)]))
        r2 = list(pool.imap_unordered(one_ts, [(e, s) for e in (0.02, 0.045) for s in range(10)]))
    summ = {}
    for k in ("P-white", "P-AR1", "C"):
        rr = [r for r in r1 if r["kind"] == k]
        summ[k] = {"PPS_rejections_of_10": int(sum(r["p_pps"] <= 0.1 for r in rr)),
                   "CS_rejections_of_10": int(sum(r["p_cs"] <= 0.1 for r in rr)),
                   "rho_over_sd": [r["rho_over_sd"] for r in rr], "m_tau": [(r["m"], r["tau"]) for r in rr]}
    res["PPS_CS"] = summ
    res["TS"] = {str(e): {"rejections_of_10": int(sum(r["reject"] for r in r2 if r["eps"] == e)),
                          "RL": [r["RL"] for r in r2 if r["eps"] == e],
                          "max_RL_sur": [r["max_RL_sur"] for r in r2 if r["eps"] == e],
                          "points_with_twins": [r["n_points_with_twins"] for r in r2 if r["eps"] == e]} for e in (0.02, 0.045)}
    res["PPS_VERIFIED"] = bool(summ["P-white"]["PPS_rejections_of_10"] <= 2 and summ["C"]["PPS_rejections_of_10"] >= 8)
    res["CS_VERIFIED"] = bool(summ["P-white"]["CS_rejections_of_10"] <= 2 and summ["C"]["CS_rejections_of_10"] >= 8)
    res["TS_VERIFIED"] = bool(res["TS"]["0.02"]["rejections_of_10"] <= 1 and res["TS"]["0.045"]["rejections_of_10"] >= 9)
    res["raw"] = {"ppscs": r1, "ts": r2}
    json.dump(res, open(OUT / "surrogates.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "raw"}, indent=1))


if __name__ == "__main__":
    main()
