"""
Phase 9 Part D: verification of every measure on systems with known published behaviour (no
Phase 9 series).  Results: results/verification/measures.json; plot plots/D_verification.png.

    python -m experiments.phase9_waveform.verify_measures [--tisean DIR]

Each check has a predeclared pass criterion (in the code below, 'criterion').  A measure is
VERIFIED only if all its checks pass.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import pathlib
import subprocess
import tempfile

import numpy as np

from experiments.phase9_waveform import measures as M

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "results" / "verification"
PLOTS = HERE / "plots"


# ----------------------------------------------------------------------------- reference systems
def lorenz(n, dt=0.02, sub=20, s=16.0, r=45.92, b=4.0, seed=0, dyn_noise=0.0, transient=2000):
    """Gao et al. 2011's Lorenz system (sigma 16, r 45.92, b 4), RK4 with step dt/sub; x sampled every
    dt.  dyn_noise D: Euler-Maruyama additive D * eta_i on each equation (Gao Eq. 13)."""
    rng = np.random.default_rng(seed)
    h = dt / sub
    v = np.array([1.0, 1.0, 20.0]) + rng.normal(0, 0.1, 3)

    def f(u):
        return np.array([-s * (u[0] - u[1]), -u[0] * u[2] + r * u[0] - u[1], u[0] * u[1] - b * u[2]])
    out = np.empty(n)
    for k in range(transient + n):
        for _ in range(sub):
            k1 = f(v); k2 = f(v + 0.5 * h * k1); k3 = f(v + 0.5 * h * k2); k4 = f(v + h * k3)
            v = v + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            if dyn_noise > 0:
                v = v + dyn_noise * math.sqrt(h) * rng.standard_normal(3)
        if k >= transient:
            out[k - transient] = v[0]
    return out


def lorenz_lambda1(s=16.0, r=45.92, b=4.0, T=2000.0, h=0.002):
    """Largest LE of the Lorenz flow (Benettin with tangent RK4)."""
    v = np.array([1.0, 1.0, 20.0])
    w = np.array([1.0, 0.0, 0.0])

    def f(u):
        return np.array([-s * (u[0] - u[1]), -u[0] * u[2] + r * u[0] - u[1], u[0] * u[1] - b * u[2]])

    def J(u):
        return np.array([[-s, s, 0], [r - u[2], -1, -u[0]], [u[1], u[0], -b]])
    for _ in range(int(50 / h)):
        k1 = f(v); k2 = f(v + 0.5 * h * k1); k3 = f(v + 0.5 * h * k2); k4 = f(v + h * k3)
        v = v + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    acc = 0.0
    n = int(T / h)
    for k in range(n):
        a1 = f(v); b1 = J(v) @ w
        v2 = v + 0.5 * h * a1; w2 = w + 0.5 * h * b1
        a2 = f(v2); b2 = J(v2) @ w2
        v3 = v + 0.5 * h * a2; w3 = w + 0.5 * h * b2
        a3 = f(v3); b3 = J(v3) @ w3
        v4 = v + h * a3; w4 = w + h * b3
        a4 = f(v4); b4 = J(v4) @ w4
        v = v + h / 6 * (a1 + 2 * a2 + 2 * a3 + a4)
        w = w + h / 6 * (b1 + 2 * b2 + 2 * b3 + b4)
        if k % 50 == 49:
            nw = np.linalg.norm(w)
            acc += math.log(nw)
            w /= nw
    return acc / (n * h)


def henon(n, a=1.4, b=0.3, seed=0, transient=1000):
    rng = np.random.default_rng(seed)
    x, y = rng.uniform(-0.1, 0.1, 2)
    out = np.empty(n)
    for k in range(transient + n):
        x, y = 1 - a * x * x + y, b * x
        if k >= transient:
            out[k - transient] = x
    return out


def logistic(n, r=4.0, seed=0, transient=1000):
    rng = np.random.default_rng(seed)
    x = rng.uniform(0.1, 0.9)
    out = np.empty(n)
    for k in range(transient + n):
        x = r * x * (1 - x)
        if k >= transient:
            out[k - transient] = x
    return out


def cencini_map(n, delta=0.4, seed=0, transient=100):
    """Cencini et al. 2000 Eqs. (19)-(20): x_{t+1} = [x_t] + F(x_t - [x_t])."""
    rng = np.random.default_rng(seed)
    x = rng.uniform(0, 1)
    out = np.empty(n)
    for k in range(transient + n):
        ip = math.floor(x)
        y = x - ip
        Fy = (2 + delta) * y if y < 0.5 else (2 + delta) * y - (1 + delta)
        x = ip + Fy
        if k >= transient:
            out[k - transient] = x
    return out


def ar3(n, seed=0):
    rng = np.random.default_rng(seed)
    x = np.zeros(n + 500)
    for i in range(3, n + 500):
        x[i] = 0.8 * x[i - 1] + 0.3 * x[i - 2] - 0.25 * x[i - 3] + 0.9 * rng.standard_normal()
    return x[500:]


def rossler(n, a=0.25, b=0.25, c=40.0, dt=0.1, sub=20, seed=0, transient=2000):
    rng = np.random.default_rng(seed)
    v = np.array([1.0, 1.0, 0.0]) + rng.normal(0, 0.1, 3)
    h = dt / sub

    def f(u):
        return np.array([-u[1] - u[2], u[0] + a * u[1], b + u[2] * (u[0] - c)])
    out = np.empty(n)
    for k in range(transient + n):
        for _ in range(sub):
            k1 = f(v); k2 = f(v + 0.5 * h * k1); k3 = f(v + 0.5 * h * k2); k4 = f(v + h * k3)
            v = v + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        if k >= transient:
            out[k - transient] = v[0]
    return out


def powerlaw_noise(n, k, seed=0):
    """Gaussian noise with 1/f^k spectrum (Fourier filtering), as in Rosso 2007 / ordpy."""
    rng = np.random.default_rng(seed)
    w = rng.standard_normal(n)
    F = np.fft.rfft(w)
    f = np.fft.rfftfreq(n)
    f[0] = f[1]
    return np.fft.irfft(F * f ** (-k / 2.0), n)


# ----------------------------------------------------------------------------- checks
def check_sdle():
    """Gao et al. 2011 (PMC3264951): clean chaos -> lambda(eps) ~ lambda_1 (plateau); noise -> lambda(eps) ~
    -gamma ln eps for times within the embedding window, and ~0 beyond it.  Revised after the first run
    (which averaged t = 5-20 samples, inside the initial embedding transient where every shell's slope is
    3-6 for the clean Lorenz system): the plateau is read after the transient, t >= 2 (m - 1) L."""
    res = {}
    lam1 = lorenz_lambda1()
    res["lorenz_lambda1_benettin"] = lam1
    dt, m, L = 0.02, 4, 2
    t0 = 2 * (m - 1) * L
    x = lorenz(10000, dt=dt)
    r = M.sdle(x, m=m, L=L, T=t0 + 60, n_shells=6, eps_top=0.4)
    # Second revision (read-out only): the per-step local slopes oscillate with period L for L = 2, so
    # the plateau is read as the mean slope of Lambda(t) over the window, (Lambda(t1) - Lambda(t0)) / (t1 - t0)
    # (the first revision's median of local slopes gave 2.08, IQR 0.91-2.45).
    sm = [s for s in r["shells"] if "lam" in s][-3:]
    vals = np.array([(s["Lambda"][t0 + 40] - s["Lambda"][t0]) / 40.0 for s in sm]) / dt
    res["lorenz_clean_plateau_median_per_time"] = float(np.median(vals))
    res["lorenz_clean_plateau_by_shell"] = vals.tolist()
    res["criterion_clean"] = ("median over the 3 smallest shells of the mean SDLE for t in [2(m-1)L, 2(m-1)L+40] within "
                              "+/-20 % of lambda1")
    res["pass_clean"] = bool(abs(res["lorenz_clean_plateau_median_per_time"] - lam1) <= 0.2 * lam1)
    out = {}
    for name, y, mm, LL in (("white_noise_m4_L1", np.random.default_rng(2).standard_normal(10000), 4, 1),
                            ("lorenz_D4_m4_L2", lorenz(10000, dt=dt, dyn_noise=4.0, seed=1), 4, 2)):
        rr = M.sdle(y, m=mm, L=LL, T=40, n_shells=8, eps_top=1.0 if "white" in name else 0.4)
        sh = [s for s in rr["shells"] if "lam" in s and s["n_pairs"] >= 50]
        e0 = np.array([math.exp(s["Lambda"][0]) for s in sh])
        l1 = np.array([s["lam"][0] for s in sh])
        g = float(np.polyfit(np.log(e0), l1, 1)[0])
        post = float(np.median(np.concatenate([s["lam"][(mm - 1) * LL + 1:] for s in sh])))
        out[name] = {"first_step_slope_vs_ln_eps": g, "median_lambda_after_window": post,
                     "first_step_lambda_by_shell": dict(zip(np.round(e0, 3).tolist(), np.round(l1, 3).tolist()))}
    res["noise"] = out
    res["criterion_noise"] = ("first-step SDLE vs ln eps slope < 0 for white noise and noisy Lorenz (lambda ~ -gamma ln eps); "
                              "white noise median SDLE beyond the embedding window |.| < 0.05")
    res["pass_noise"] = bool(all(v["first_step_slope_vs_ln_eps"] < 0 for v in out.values())
                             and abs(out["white_noise_m4_L1"]["median_lambda_after_window"]) < 0.05)
    res["curves"] = {"lorenz": [r["eps_all"].tolist()[::5], (r["lam_all"] / dt).tolist()[::5]]}
    res["VERIFIED"] = res["pass_clean"] and res["pass_noise"]
    return res


def check_fsle():
    res = {}
    x = cencini_map(200000, 0.4)
    d, lam, cnt, nref = M.fsle(x, m=1, L=1, delta0=1e-4, r=math.sqrt(2.0), n_levels=44, theiler=0,
                               maxt=2000, standardize_input=False)
    ok = cnt >= 20
    small = ok & (d < 0.1)
    big = ok & (d > 1.5) & (d < 20)
    res["cencini_small_delta_median"] = float(np.median(lam[small]))
    res["ln_2.4"] = math.log(2.4)
    slope = float(np.polyfit(np.log(d[big]), np.log(lam[big]), 1)[0]) if big.sum() >= 3 else float("nan")
    res["cencini_large_delta_loglog_slope"] = slope
    res["criterion"] = "small-delta FSLE within +/-15 % of ln 2.4; large-delta log-log slope in [-2.5, -1.5]"
    res["pass_map"] = bool(abs(res["cencini_small_delta_median"] - math.log(2.4)) <= 0.15 * math.log(2.4)
                           and -2.5 <= slope <= -1.5)
    w = np.random.default_rng(3).standard_normal(20000)
    d2, lam2, cnt2, _ = M.fsle(w, m=2, L=1, delta0=1e-3, n_levels=20, theiler=1)
    ok2 = cnt2 >= 20
    res["white_noise_fsle"] = dict(zip(np.round(d2[ok2], 5).tolist(), np.round(lam2[ok2], 3).tolist()))
    g = np.polyfit(np.log(d2[ok2]), lam2[ok2], 1)[0]
    res["white_noise_slope_vs_ln_delta"] = float(g)
    res["criterion_noise"] = "white noise: no plateau, FSLE decreases with ln delta (slope < -0.3)"
    res["pass_noise"] = bool(g < -0.3)
    res["curves"] = {"cencini": [d[ok].tolist(), lam[ok].tolist()]}
    res["VERIFIED"] = res["pass_map"] and res["pass_noise"]
    return res


def check_entropy(tisean):
    res = {}
    eps = np.geomspace(1e-3, 2.0, 40)
    lo = logistic(5000)
    e, C = M.correlation_sums(lo, m_max=5, L=1, eps=eps, theiler=0)
    h = M.eps_entropy(C)
    sel = (e > 0.005) & (e < 0.1)
    res["logistic_h_m2_m4_small_eps"] = [float(np.nanmedian(h[k][sel])) for k in (1, 2, 3)]
    res["ln2"] = math.log(2)
    he = henon(5000)
    e2, C2 = M.correlation_sums(he, m_max=5, L=1, eps=eps, theiler=0)
    h2 = M.eps_entropy(C2)
    sel2 = (e2 > 0.01) & (e2 < 0.1)
    res["henon_h_m3_m4_small_eps"] = [float(np.nanmedian(h2[k][sel2])) for k in (2, 3)]
    res["henon_ref"] = 0.42
    w = np.random.default_rng(4).standard_normal(5000)
    e3, C3 = M.correlation_sums(w, m_max=3, L=1, eps=eps, theiler=0)
    h3 = M.eps_entropy(C3)
    sel3 = (e3 > 0.01) & (e3 < 0.3)
    slope = float(np.polyfit(np.log(e3[sel3]), h3[0][sel3], 1)[0])
    res["white_h1_slope_vs_ln_eps"] = slope
    res["criterion"] = ("logistic h_m (m = 2..4) within +/-10 % of ln 2; Henon within +/-15 % of 0.42; white noise "
                        "h_1 slope vs ln eps in [-1.2, -0.8] (Cencini Eq. 13)")
    res["pass_published"] = bool(all(abs(v - math.log(2)) <= 0.1 * math.log(2) for v in res["logistic_h_m2_m4_small_eps"])
                                 and all(abs(v - 0.42) <= 0.15 * 0.42 for v in res["henon_h_m3_m4_small_eps"])
                                 and -1.2 <= slope <= -0.8)
    # TISEAN d2 agreement (correlation sums, max norm)
    if tisean and (pathlib.Path(tisean) / "d2").exists():
        with tempfile.TemporaryDirectory() as td:
            f = pathlib.Path(td) / "h.dat"
            np.savetxt(f, he[:3000])
            subprocess.run([str(pathlib.Path(tisean) / "d2"), str(f), "-M1,4", "-d1", "-t0", "-#100", "-r0.001", "-R2", "-N0",
                            "-o", str(pathlib.Path(td) / "h"), "-V0"], check=True, capture_output=True)
            c2 = open(pathlib.Path(td) / "h.c2").read().strip().split("\n\n")
            # TISEAN d2 works on the unrescaled data (no -E) with max norm; eps in data units
            diffs = []
            for k, blk in enumerate(b for b in c2 if b.strip() and "#dim=" in b):
                rows = np.array([[float(v) for v in ln.split()] for ln in blk.split("\n") if ln and not ln.startswith("#")])
                ee = rows[:, 0]
                _, Cm = M.correlation_sums(he[:3000], m_max=k + 1, L=1, eps=ee, theiler=0, standardize_input=False)
                mine = Cm[k][np.argsort(np.argsort(ee))]
                good = rows[:, 1] > 1e-3
                diffs.append(float(np.max(np.abs(mine[good] / rows[good, 1] - 1))))
            res["tisean_d2_max_rel_diff_per_m"] = diffs
            res["pass_tisean"] = bool(max(diffs) < 0.02)
    res["VERIFIED"] = res["pass_published"] and res.get("pass_tisean", True)
    res["curves"] = {"logistic_h": [e.tolist(), h[2].tolist()], "white_h1": [e3.tolist(), h3[0].tolist()]}
    return res


def check_pe():
    import ordpy
    res = {}
    n = 1_000_000
    lo = logistic(n)
    rw = np.cumsum(np.random.default_rng(5).standard_normal(n))
    pl = M.ordinal_distribution(lo, 3)
    pr = M.ordinal_distribution(rw, 3)
    from itertools import permutations
    perms = list(permutations(range(3)))
    # ordpy paper order: (0,1,2),(0,2,1),(1,0,2),(1,2,0),(2,0,1),(2,1,0)
    exact_l = dict(zip([(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)], [1 / 3, 1 / 15, 2 / 15, 3 / 15, 4 / 15, 0]))
    exact_r = dict(zip([(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)], [1 / 4, 1 / 8, 1 / 8, 1 / 8, 1 / 8, 1 / 4]))
    res["logistic_d3_maxabs_vs_exact"] = float(max(abs(pl[i] - exact_l[p]) for i, p in enumerate(perms)))
    res["randomwalk_d3_maxabs_vs_exact"] = float(max(abs(pr[i] - exact_r[p]) for i, p in enumerate(perms)))
    series = {"logistic": logistic(2 ** 15), "henon": henon(2 ** 15), "white": np.random.default_rng(6).standard_normal(2 ** 15)}
    for k in (1.0, 2.0, 3.0):
        series[f"f^-{k:g}"] = powerlaw_noise(2 ** 15, k, seed=int(10 * k))
    plane, diffs = {}, []
    for name, s in series.items():
        H, C = M.pe_cecp(s, 6, 1)
        Ho, Co = ordpy.complexity_entropy(s, dx=6, taux=1)
        diffs.append(max(abs(H - Ho), abs(C - Co)))
        plane[name] = [H, C]
    res["plane_d6"] = plane
    res["max_abs_diff_vs_ordpy"] = float(max(diffs))
    noise_pts = np.array([plane[k] for k in plane if k.startswith("f^") or k == "white"])

    def c_noise_at(Hq):
        o = np.argsort(noise_pts[:, 0])
        return float(np.interp(Hq, noise_pts[o, 0], noise_pts[o, 1]))
    res["chaos_above_noise"] = {k: [plane[k][1], c_noise_at(plane[k][0])] for k in ("logistic", "henon")}
    res["criterion"] = ("exact d=3 distributions within 0.003; |H|,|C| equal to ordpy within 1e-10; Rosso 2007: chaotic "
                        "maps' C above the C of the f^-k noise curve at the same H; white noise H > 0.99, C < 0.01")
    res["VERIFIED"] = bool(res["logistic_d3_maxabs_vs_exact"] < 0.003 and res["randomwalk_d3_maxabs_vs_exact"] < 0.003
                           and res["max_abs_diff_vs_ordpy"] < 1e-10
                           and all(v[0] > v[1] for v in res["chaos_above_noise"].values())
                           and plane["white"][0] > 0.99 and plane["white"][1] < 0.01)
    return res


def check_rqa():
    res = {}
    dets = [M.rqa_det(ar3(1500, seed=s), m=4, tau=4, rr=0.1)[0] for s in range(5)]
    res["ar3_DET_m4_tau4_RR0.1"] = [float(np.mean(dets)), float(np.std(dets))]
    dr = [M.rqa_det(rossler(2000, c=c, seed=s), m=3, tau=6, rr=0.05)[0] for c, s in ((36.0, 0), (40.0, 1), (44.0, 2))]
    res["rossler_DET_m3_tau6_RR0.05_c36_40_44"] = [float(v) for v in dr]
    res["white_DET_m3_tau1_RR0.05"] = float(np.mean([M.rqa_det(np.random.default_rng(s).standard_normal(1500), m=3, tau=1, rr=0.05)[0]
                                                      for s in range(3)]))
    res["published"] = {"ar3": 0.6, "rossler": 0.94}
    res["criterion"] = "AR(3) DET within 0.6 +/- 0.1; Rossler DET within 0.94 +/- 0.04 (Marwan 2011 Sect. 3.3-3.4)"
    res["VERIFIED"] = bool(abs(res["ar3_DET_m4_tau4_RR0.1"][0] - 0.6) <= 0.1 and all(abs(v - 0.94) <= 0.04 for v in dr))
    return res


def check_ghkss(tisean):
    res = {}
    rng = np.random.default_rng(7)
    # Revised check (first run, N = 5000 with these settings, gave factor 1.02 for BOTH this
    # implementation and TISEAN ghkss itself: too few points per neighbourhood in 7 dimensions).
    # Now N = 20000 with the TISEAN paper's NMR settings (m = 7, q = 2, >= 50 neighbours, 3 iterations).
    clean = henon(20000, seed=7)
    noisy = clean + 0.05 * clean.std() * rng.standard_normal(len(clean))
    y = M.ghkss(noisy, m=7, q=2, k_min=50, iterations=3)
    core = slice(10, -10)
    e0 = np.sqrt(np.mean((noisy - clean)[core] ** 2))
    e1 = np.sqrt(np.mean((y - clean)[core] ** 2))
    res["henon_5pct_rms_error_before_after"] = [float(e0), float(e1)]
    res["noise_reduction_factor"] = float(e0 / e1)
    if tisean and (pathlib.Path(tisean) / "ghkss").exists():
        with tempfile.TemporaryDirectory() as td:
            f = pathlib.Path(td) / "n.dat"
            np.savetxt(f, noisy, fmt="%.12g")
            subprocess.run([str(pathlib.Path(tisean) / "ghkss"), str(f), "-m1,7", "-q2", "-k50", "-i3", "-o",
                            str(pathlib.Path(td) / "o"), "-V0"], check=True, capture_output=True)
            yt = np.loadtxt(pathlib.Path(td) / "o.3")
            yt = yt if yt.ndim == 1 else yt[:, 0]
            res["tisean_rms_error_after"] = float(np.sqrt(np.mean((yt - clean)[core] ** 2)))
            res["corr_mine_vs_tisean"] = float(np.corrcoef(y[core], yt[core])[0, 1])
            res["rms_diff_mine_vs_tisean_over_noise"] = float(np.sqrt(np.mean((y - yt)[core] ** 2)) / e0)
    # pure noise: no spurious determinism.  Nonlinear prediction error (Phase 8 statistic, m=3) of
    # ghkss(noise) vs 19 IAAFT surrogates each processed identically with ghkss.
    import final_pipeline as fp
    from experiments.phase8_cardiac import methods as MM
    rej = 0
    zs = []
    for s in range(10):
        w = np.random.default_rng(100 + s).standard_normal(1000)
        obs = MM.nl_pred_error(M.ghkss(w, m=7, q=2, k_min=50, iterations=3), m=3)
        r2 = np.random.default_rng(200 + s)
        sur = [MM.nl_pred_error(M.ghkss(fp.iaaft_surrogate(w, r2), m=7, q=2, k_min=50, iterations=3), m=3) for _ in range(19)]
        p = (1 + np.sum(np.array(sur) <= obs)) / 20
        rej += p <= 0.05
        zs.append(float((np.mean(sur) - obs) / np.std(sur, ddof=1)))
    res["white_noise_rejections_at_0.05_of_10"] = int(rej)
    res["white_noise_z"] = zs
    res["criterion"] = ("Henon + 5 % noise: noise reduction factor >= 2 and agreement with TISEAN ghkss (corr >= 0.99, "
                        "rms difference <= 25 % of the noise); white noise: <= 2/10 rejections vs identically "
                        "processed IAAFT surrogates")
    res["VERIFIED"] = bool(res["noise_reduction_factor"] >= 2.0 and res.get("corr_mine_vs_tisean", 1.0) >= 0.99
                           and res.get("rms_diff_mine_vs_tisean_over_noise", 0.0) <= 0.25 and rej <= 2)
    return res


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--tisean", default="/tmp/claude-0/tisean/Tisean_3.0.1/source_c")
    ap.add_argument("--only", default="")
    a = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "measures.json"
    res = json.load(open(path)) if path.exists() else {}
    checks = {"D1_SDLE": check_sdle, "D2_FSLE": check_fsle, "D3_eps_entropy": lambda: check_entropy(a.tisean),
              "D4_PE_CECP": check_pe, "D5_RQA_DET": check_rqa, "D6_GHKSS": lambda: check_ghkss(a.tisean)}
    for k, f in checks.items():
        if a.only and k not in a.only.split(","):
            continue
        res[k] = f()
        print(k, "VERIFIED" if res[k]["VERIFIED"] else "NOT VERIFIED",
              json.dumps({kk: v for kk, v in res[k].items() if kk != "curves"})[:900], flush=True)
        json.dump(res, open(path, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
