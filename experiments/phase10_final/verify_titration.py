"""
Phase 10 Q2b: verification of the noise-titration implementation (titration10.py) against the published
behaviour in Poon & Barahona, PNAS 2001 (criteria V1-V6 fixed in METHODS.md 2.5 before this script was run).

    python -m experiments.phase10_final.verify_titration --workers 4

Output: results/verification/titration_verification.json and .md
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
from multiprocessing import Pool

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "results" / "verification"
N = 1000


def logistic(r, n=N, x0=0.4, burn=1000):
    x = x0
    out = np.empty(n)
    for i in range(burn + n):
        x = r * x * (1 - x)
        if i >= burn:
            out[i - burn] = x
    return out


def logistic_le(r, n=200000, x0=0.4):
    x = x0
    for _ in range(1000):
        x = r * x * (1 - x)
    s = 0.0
    for _ in range(n):
        x = r * x * (1 - x)
        s += np.log(max(abs(r * (1 - 2 * x)), 1e-300))
    return s / n


def henon(n=N, a=1.4, b=0.3, burn=1000):
    x, y = 0.1, 0.1
    out = np.empty(n)
    for i in range(burn + n):
        x, y = 1 - a * x * x + y, b * x
        if i >= burn:
            out[i - burn] = x
    return out


def lorenz(rho, n=N, ts=0.075, dt=0.0025, transient=100.0, sigma=10.0, beta=8.0 / 3.0):
    def f(s):
        x, y, z = s
        return np.array([sigma * (y - x), x * (rho - z) - y, x * y - beta * z])
    s = np.array([1.0, 1.0, 20.0])
    every = int(round(ts / dt))
    total = int(round(transient / dt)) + n * every
    out = []
    for i in range(1, total + 1):
        k1 = f(s); k2 = f(s + 0.5 * dt * k1); k3 = f(s + 0.5 * dt * k2); k4 = f(s + dt * k3)
        s = s + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        if i > int(round(transient / dt)) and (i - int(round(transient / dt))) % every == 0:
            out.append(s[0])
    return np.array(out[:n])


def mackey_glass(tau, n=N, ts=12.5, dt=0.05, transient=5000.0, seed=0):
    from experiments.phase8_cardiac.models import mackey_glass as MG
    rng = np.random.default_rng(seed)
    lag = int(round(tau / dt))
    every = int(round(ts / dt))
    steps = int(transient / dt) + n * every + 1
    x = MG.integrate(np.full(lag + 1, rng.uniform(0.5, 1.3)), steps, dt, lag, MG.BETA, MG.GAMMA, MG.NEXP, 0.0,
                     np.zeros(steps))
    x = x[lag + 1 + int(transient / dt):]
    return x[::every][:n]


def controls(kind, seed):
    rng = np.random.default_rng(10_000 + seed)
    k = np.arange(N)
    if kind == "white":
        return rng.standard_normal(N)
    if kind == "ar2":
        e = rng.standard_normal(N + 500)
        y = np.zeros(N + 500)
        for i in range(2, N + 500):
            y[i] = 1.5 * y[i - 1] - 0.75 * y[i - 2] + e[i]
        return y[500:]
    if kind == "limit_cycle":
        ph = rng.uniform(0, 2 * np.pi)
        w = 2 * np.pi / 23.7
        y = np.sin(w * k + ph) + 0.5 * np.sin(2 * w * k + 2 * ph) + 0.25 * np.sin(3 * w * k + 3 * ph)
        return y + 0.01 * y.std() * rng.standard_normal(N)
    if kind == "torus":
        p1, p2 = rng.uniform(0, 2 * np.pi, 2)
        y = (2 + np.sin(1.0 * k + p1)) * (2 + np.sin(np.sqrt(2) * k + p2))
        return y + 0.01 * y.std() * rng.standard_normal(N)
    if kind == "logistic_p3":
        y = logistic(3.83, x0=rng.uniform(0.2, 0.8))
        return y + 0.01 * y.std() * rng.standard_normal(N)
    raise ValueError(kind)


def task(t):
    from experiments.phase10_final import titration10 as T
    name, arg = t
    if name == "logistic":
        y, extra = logistic(arg), {"LE": logistic_le(arg)}
    elif name == "henon":
        y, extra = henon(), {}
    elif name == "lorenz":
        y, extra = lorenz(arg), {}
    elif name == "mg":
        y, extra = mackey_glass(arg), {}
    else:
        y, extra = controls(name, arg), {}
    import zlib
    rng = np.random.default_rng(np.random.SeedSequence([20261003, zlib.crc32(f"{name}|{arg}".encode())]))
    o = T.titrate(y, rng)
    return {"system": name, "arg": arg, "NL": o["NL"], "NL_i": o["NL_i"], "positive": o["positive"],
            "base": {k: v for k, v in o["base"].items()}, **extra}


def tasks():
    T = [("logistic", r) for r in (3.7, 3.575, 3.565)]
    T += [("logistic", round(3.5 + 0.005 * i, 3)) for i in range(101) if round(3.5 + 0.005 * i, 3) not in (3.7, 3.575, 3.565)]
    T += [(k, s) for k in ("white", "ar2", "limit_cycle", "torus", "logistic_p3") for s in range(20)]
    T += [("lorenz", 28.0), ("lorenz", 160.0), ("henon", 0)]
    T += [("mg", tau) for tau in (17.0, 23.0, 30.0, 14.0, 15.0, 16.0)]
    return T


def evaluate(R):
    from scipy import stats
    get = lambda s, a: next(r for r in R if r["system"] == s and r["arg"] == a)  # noqa: E731
    v = {}
    nl37, nl3575, nl3565 = get("logistic", 3.7)["NL"], get("logistic", 3.575)["NL"], get("logistic", 3.565)["NL"]
    v["V1"] = {"NL_3.7": nl37, "pass": abs(nl37 - 0.75) <= 0.15}
    v["V2"] = {"NL_3.575": nl3575, "pass": 0.03 <= nl3575 <= 0.20}
    v["V3"] = {"NL_3.565": nl3565, "pass": nl3565 <= 0.02}
    L = sorted([r for r in R if r["system"] == "logistic"], key=lambda r: r["arg"])
    ch = [r for r in L if r["LE"] > 0.02]
    pe = [r for r in L if r["LE"] < -0.02]
    fpos = np.mean([r["NL"] > 0 for r in ch]) if ch else float("nan")
    fzero = np.mean([r["NL"] == 0 for r in pe]) if pe else float("nan")
    rho = float(stats.spearmanr([r["NL"] for r in ch], [r["LE"] for r in ch])[0]) if len(ch) > 2 else float("nan")
    v["V4"] = {"n_chaotic": len(ch), "frac_NL_pos_chaotic": float(fpos), "n_periodic": len(pe),
               "frac_NL_zero_periodic": float(fzero), "spearman_NL_LE": rho,
               "pass": bool(fpos >= 0.9 and fzero >= 0.9 and rho >= 0.6),
               "scan": [(r["arg"], round(r["LE"], 4), round(r["NL"], 4)) for r in L]}
    c = {}
    for k in ("white", "ar2", "limit_cycle", "torus", "logistic_p3"):
        rr = [r for r in R if r["system"] == k]
        c[k] = {"n": len(rr), "positive": int(sum(r["positive"] for r in rr)),
                "NL_zero": int(sum(r["NL"] == 0 for r in rr))}
    v["V5"] = {"controls": c, "pass": bool(c["white"]["positive"] <= 2 and c["ar2"]["positive"] <= 2 and
                                          all(c[k]["NL_zero"] >= 18 for k in ("limit_cycle", "torus", "logistic_p3")))}
    mgc = [get("mg", t)["NL"] for t in (17.0, 23.0, 30.0)]
    mgp = [get("mg", t)["NL"] for t in (14.0, 15.0, 16.0)]
    v["V6"] = {"lorenz28": get("lorenz", 28.0)["NL"], "lorenz160": get("lorenz", 160.0)["NL"],
               "henon": get("henon", 0)["NL"], "mg_chaotic_17_23_30": mgc, "mg_periodic_14_15_16": mgp,
               "pass": bool(get("lorenz", 28.0)["NL"] > 0 and get("lorenz", 160.0)["NL"] == 0 and
                            get("henon", 0)["NL"] > 0 and sum(x > 0 for x in mgc) >= 2 and sum(x == 0 for x in mgp) >= 2)}
    v["VERIFIED"] = all(v[k]["pass"] for k in ("V1", "V2", "V3", "V4", "V5", "V6"))
    return v


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    with Pool(a.workers) as p:
        R = p.map(task, tasks(), chunksize=1)
    v = evaluate(R)
    json.dump({"records": R, "verification": v}, open(OUT / "titration_verification.json", "w"), indent=1,
              default=float)
    lines = ["# Titration verification (criteria: METHODS.md 2.5)", ""]
    for k in ("V1", "V2", "V3", "V4", "V5", "V6"):
        d = {kk: vv for kk, vv in v[k].items() if kk != "scan"}
        lines.append(f"- **{k}** {'PASS' if v[k]['pass'] else 'FAIL'}: {json.dumps(d, default=float)}")
    lines += ["", f"**VERIFIED: {v['VERIFIED']}**", "", "Logistic scan (r, LE, NL):", ""]
    lines += [f"| {r} | {le} | {nl} |" for r, le, nl in v["V4"]["scan"]]
    (OUT / "titration_verification.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[:12]))


if __name__ == "__main__":
    main()
