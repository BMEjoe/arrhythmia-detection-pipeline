"""
Final refinement, Part 1: independent verification of the key reported numbers.

    python -m experiments.final_refinement.verify_all          # writes VERIFICATION.md, verification.json

Every number is recomputed DIRECTLY from the committed raw result files (JSONL / CSV / JSON written by the
runners), with code written for this check. Nothing is imported from the project's analysis, decision or
pipeline modules: rates, Wilson / Clopper-Pearson bounds, cluster bootstraps, permutation tests, GEE
(independence and exchangeable working correlation, robust sandwich SE), Firth-free exclusion bounds,
leave-one-subject-out L2 logistic regression and AUC are all re-implemented here with numpy (scipy only for
the binomial CDF and the normal distribution).

Comparison rules (column "rule" in VERIFICATION.md):
  exact   counts, flags, labels and strings must be identical.
  round   the recomputed value, rounded to the precision printed in the report, equals the reported value
          (|recomputed - reported| <= half a unit of the last printed digit).
  mc      Monte Carlo quantities (cluster-bootstrap percentiles, permutation p-values) were recomputed with
          an independent bootstrap / permutation using a DIFFERENT seed; MATCH iff every bound is within the
          stated tolerance (absolute). The original analysis is seeded and therefore deterministic.
  text    a statement in a report checked against the recomputed numbers.

Each MISMATCH is handled under the error policy (experiments/final_refinement/ERRATA.md).
"""
from __future__ import annotations

import json
import math
import pathlib
import sys
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from scipy import stats

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = pathlib.Path(__file__).resolve().parent
E = ROOT / "experiments"
VSEED = 777_000_001          # verification RNG seed (deliberately different from every analysis seed)

CHECKS: list[dict] = []


# =================================================================================== comparison helpers
def _flat(v):
    if isinstance(v, (list, tuple)):
        out = []
        for x in v:
            out.extend(_flat(x))
        return out
    return [v]


def _fmt(v):
    if isinstance(v, (list, tuple)):
        return "(" + ", ".join(_fmt(x) for x in v) + ")"
    if isinstance(v, (bool, np.bool_)):
        return str(bool(v))
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        if v is None or not np.isfinite(v):
            return "nan"
        if v != 0 and (abs(v) < 1e-3 or abs(v) >= 1e5):
            return f"{v:.3g}"
        return f"{v:.6g}"
    return str(v)


def _decimals(x):
    """Number of printed decimals of a reported number given as a string or float."""
    s = x if isinstance(x, str) else repr(x)
    s = s.replace("−", "-")
    if "e" in s.lower():
        return None
    return len(s.split(".")[1]) if "." in s else 0


def check(cid, phase, quantity, reported, recomputed, source, where, rule="exact", tol=None, note="",
          scale=1.0):
    """reported: value(s) as printed in the report (numbers, or strings for numbers to keep the printed
    precision); recomputed: raw value(s); scale multiplies recomputed before comparison (e.g. 100 for %)."""
    rep = _flat(reported)
    rec = _flat(recomputed)
    ok = len(rep) == len(rec)
    if ok:
        for a, b in zip(rep, rec):
            if rule == "exact":
                ok &= (a == b) if not isinstance(a, str) else (str(a) == str(b))
            elif rule == "round":
                if b is None or (isinstance(b, float) and not np.isfinite(b)):
                    ok = False
                    continue
                d = _decimals(a)
                av = float(str(a).replace("−", "-"))
                bv = float(b) * scale
                if d is None:      # scientific notation: compare at 1 significant digit of the mantissa
                    ok &= abs(av - bv) <= 0.5 * 10 ** math.floor(math.log10(abs(av))) + 1e-300
                else:
                    ok &= abs(av - bv) <= 0.5 * 10 ** (-d) + 1e-9
            elif rule == "mc":
                av = float(str(a).replace("−", "-"))
                ok &= abs(av - float(b) * scale) <= tol + 1e-12
            elif rule == "text":
                ok &= bool(b) if isinstance(b, (bool, np.bool_)) else (a == b)
            else:
                raise ValueError(rule)
    rec_disp = recomputed
    if scale != 1.0:
        rec_disp = [x * scale if isinstance(x, (int, float, np.floating)) else x for x in _flat(recomputed)]
        if not isinstance(recomputed, (list, tuple)):
            rec_disp = rec_disp[0]
    CHECKS.append({"id": cid, "phase": phase, "quantity": quantity, "reported": _fmt(reported),
                   "recomputed": _fmt(rec_disp), "source": source, "where": where,
                   "rule": rule + (f" (±{tol:g})" if tol is not None else ""), "status": "MATCH" if ok else "MISMATCH",
                   "note": note})
    return ok


# =================================================================================== statistics (own code)
def wilson(k, n, z=1.959963984540054):
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def cp_upper_bisect(k, n, alpha=0.05):
    """One-sided (1 - alpha) Clopper-Pearson upper bound by bisection on the binomial CDF
    (P(X <= k; p) = alpha), independent of the beta-quantile formula used in the analysis."""
    if k >= n:
        return 1.0
    lo, hi = k / n, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if stats.binom.cdf(k, n, mid) > alpha:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def cluster_boot_means(vals, B, rng):
    vals = np.asarray(vals, float)
    n = len(vals)
    idx = rng.integers(0, n, size=(B, n))
    return vals[idx].mean(axis=1)


def pct(v, q):
    return float(np.percentile(np.asarray(v, float), q))


def logistic_irls(X, y, w=None, max_iter=100):
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    b = np.zeros(X.shape[1])
    for _ in range(max_iter):
        eta = X @ b
        mu = 1 / (1 + np.exp(-eta))
        W = mu * (1 - mu)
        H = X.T @ (X * W[:, None])
        g = X.T @ (y - mu)
        step = np.linalg.solve(H, g)
        b = b + step
        if np.max(np.abs(step)) < 1e-12:
            break
    return b


def gee_logit(X, y, groups, corr="independence", max_iter=200, tol=1e-10):
    """Own GEE (Liang & Zeger 1986), binomial / logit, robust sandwich covariance (no small-sample
    correction). Exchangeable dependence parameter: moment estimator with the scale and pair-count
    degrees-of-freedom corrections of statsmodels (ddof = number of coefficients)."""
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    groups = np.asarray(groups)
    p = X.shape[1]
    order = np.argsort(groups, kind="stable")
    X, y, groups = X[order], y[order], groups[order]
    _, starts = np.unique(groups, return_index=True)
    bounds = list(zip(starts, list(starts[1:]) + [len(y)]))
    b = logistic_irls(X, y)
    alpha = 0.0

    def vinv(n, a):
        # inverse of the exchangeable correlation matrix (1 - a) I + a 11'
        if a == 0.0 or corr == "independence":
            return np.eye(n)
        return (np.eye(n) - a / (1 + (n - 1) * a) * np.ones((n, n))) / (1 - a)

    for _ in range(max_iter):
        mu = 1 / (1 + np.exp(-(X @ b)))
        sd = np.sqrt(mu * (1 - mu))
        if corr == "exchangeable":
            r = (y - mu) / sd
            N = len(y)
            scale = np.sum(r * r) / (N - p)
            pair_sum, npairs = 0.0, 0.0
            for s, e in bounds:
                rr = r[s:e]
                pair_sum += (rr.sum() ** 2 - np.sum(rr * rr)) / 2
                npairs += 0.5 * (e - s) * (e - s - 1)
            alpha = (pair_sum / scale) / (npairs - p)
        A = np.zeros((p, p))
        U = np.zeros(p)
        for s, e in bounds:
            D = X[s:e] * (mu[s:e] * (1 - mu[s:e]))[:, None]       # d mu / d beta
            Ai = 1 / sd[s:e]
            Vi = Ai[:, None] * vinv(e - s, alpha) * Ai[None, :]      # V^-1 = A^-1/2 R^-1 A^-1/2
            A += D.T @ Vi @ D
            U += D.T @ Vi @ (y[s:e] - mu[s:e])
        step = np.linalg.solve(A, U)
        b = b + step
        if not np.all(np.isfinite(b)):
            return {"converged": False}
        if np.max(np.abs(step)) < tol:
            break
    mu = 1 / (1 + np.exp(-(X @ b)))
    sd = np.sqrt(mu * (1 - mu))
    A = np.zeros((p, p))
    Bm = np.zeros((p, p))
    for s, e in bounds:
        D = X[s:e] * (mu[s:e] * (1 - mu[s:e]))[:, None]
        Ai = 1 / sd[s:e]
        Vi = Ai[:, None] * vinv(e - s, alpha) * Ai[None, :]
        A += D.T @ Vi @ D
        u = D.T @ Vi @ (y[s:e] - mu[s:e])
        Bm += np.outer(u, u)
    Ainv = np.linalg.inv(A)
    cov = Ainv @ Bm @ Ainv
    se = np.sqrt(np.diag(cov))
    z = b / se
    pv = 2 * stats.norm.sf(np.abs(z))
    return {"converged": True, "beta": b, "se": se, "p": pv, "alpha": alpha,
            "OR": np.exp(b), "lo": np.exp(b - 1.96 * se), "hi": np.exp(b + 1.96 * se)}


def auc_mw(score, y):
    """Mann-Whitney AUC with mid-ranks for ties."""
    score = np.asarray(score, float)
    y = np.asarray(y, int)
    r = stats.rankdata(score)
    n1 = int(y.sum())
    n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def l2_logit_fit(X, y, C=1.0):
    """argmin 0.5 ||w||^2 + C * sum log-loss (intercept unpenalised), Newton iterations."""
    n, p = X.shape
    Xa = np.column_stack([np.ones(n), X])
    b = np.zeros(p + 1)
    P = np.eye(p + 1)
    P[0, 0] = 0.0
    for _ in range(200):
        mu = 1 / (1 + np.exp(-(Xa @ b)))
        g = C * Xa.T @ (mu - y) + P @ b
        H = C * Xa.T @ (Xa * (mu * (1 - mu))[:, None]) + P
        step = np.linalg.solve(H, g)
        b = b - step
        if np.max(np.abs(step)) < 1e-12:
            break
    return b


def loso_oof(df, feats, C=1.0):
    X = df[feats].to_numpy(float)
    y = df["label"].to_numpy(int)
    g = df["subject"].to_numpy()
    oof = np.full(len(df), np.nan)
    for s in np.unique(g):
        te = g == s
        tr = ~te
        if len(np.unique(y[tr])) < 2:
            continue
        mu, sdv = X[tr].mean(0), X[tr].std(0)          # StandardScaler (population SD) fit on the training fold
        sdv[sdv == 0] = 1.0
        b = l2_logit_fit((X[tr] - mu) / sdv, y[tr], C)
        oof[te] = 1 / (1 + np.exp(-(b[0] + ((X[te] - mu) / sdv) @ b[1:])))
    return oof


def jl(path):
    return [json.loads(line) for line in open(path)]


def rel(p):
    return str(pathlib.Path(p).relative_to(ROOT))


# =================================================================================== Phase 2E (quoted in Phase 10 S12)
def phase2e():
    src = E / "phase2e/results/core.jsonl"
    R = [json.loads(l)["row"] for l in open(src)]
    pos = Counter()
    lb = Counter()
    la = defaultdict(list)
    for r in R:
        s = r["system"]
        if r.get("lle_valid") and r.get("lle_per_beat") is not None and np.isfinite(r["lle_per_beat"]):
            pos[(s, "valid")] += 1
            pos[(s, "pos")] += r["lle_per_beat"] > 0
        if s in ("white_noise", "ar1"):
            k = r.get("significant_peak_count")
            lb[(s, "n")] += 1
            lb[(s, "k")] += bool(k is not None and np.isfinite(k) and k > 0)
            sp = r.get("source_peak_count")
            la[(s, r["window_length"] if "window_length" in r else r.get("n"))].append(
                bool(sp is not None and np.isfinite(sp) and sp > 0))
    where = "docs/PHASE10_FINAL_RESULTS.md S12 #1-2 (from PHASE2E 3.1-3.2)"
    check("2E-1", "2E", "positive production LLE, white noise (valid windows)", "194/195",
          f"{pos[('white_noise', 'pos')]}/{pos[('white_noise', 'valid')]}", rel(src), where)
    check("2E-2", "2E", "positive production LLE, AR(1)", "197/197",
          f"{pos[('ar1', 'pos')]}/{pos[('ar1', 'valid')]}", rel(src), where)
    check("2E-3", "2E", "Level-B false positives, white noise (128+256+512, 50 surrogates)", "14/200",
          f"{lb[('white_noise', 'k')]}/{lb[('white_noise', 'n')]}", rel(src), where)
    check("2E-4", "2E", "Level-B false positives, AR(1)", "15/200",
          f"{lb[('ar1', 'k')]}/{lb[('ar1', 'n')]}", rel(src), where)
    keys = [k for k in la if k[0] == "white_noise"]
    rates = [100 * np.mean(la[k]) for k in sorted(keys, key=str)]
    check("2E-5", "2E", "Level A fires on white-noise windows: range over lengths (%)", ("70", "98"),
          (min(rates), max(rates)), rel(src), where, rule="round",
          note="Level A = >= 1 source (So) peak; range over 128/256/512")


# =================================================================================== Phase 3
def phase3():
    d = E / "phase3_lle/results/test"
    where = "docs/PHASE3_LLE_VALIDATION.md 5.1-5.3; PHASE10 S12 #3"
    rec = {}
    for m in ("baseline", "c1_rosenstein_m2_iaaft", "c2_kantz_m3_sat_iaaft", "c3_eps_m2_iaaft",
              "c4_zero_one_iaaft"):
        rows = jl(d / f"{m}.jsonl")
        rec[m] = rows
    # detection recomputed from the stored p-value (p <= 0.05), independent of the stored flag
    def det(r):
        p = r["result"].get("p")
        return bool(p is not None and p <= 0.05)

    flag_mismatch = sum(det(r) != bool(r["result"]["detected"]) for m in rec for r in rec[m])
    check("3-0", "3", "stored 'detected' flag equals p <= 0.05 (all 12,000 test windows)", 0, flag_mismatch,
          rel(d / "*.jsonl"), where, note="0 = no window where the flag disagrees with its p-value")

    def cnt(m, sys_, snr, n):
        rows = [r for r in rec[m] if r["task"]["system"] == sys_ and r["task"]["snr_db"] == snr
                and r["task"]["window_length"] == n]
        return sum(det(r) for r in rows), len(rows)

    c1 = "c1_rosenstein_m2_iaaft"
    for i, (sys_, snr, rep) in enumerate([("white_noise", None, "4/100"), ("ar1", None, "3/100"),
                                          ("sinusoid", None, "0/100"), ("logistic_p4", None, "0/100"),
                                          ("logistic", None, "100/100"), ("henon", None, "100/100"),
                                          ("logistic", 10.0, "100/100"), ("henon", 10.0, "99/100")]):
        k, n = cnt(c1, sys_, snr, 256)
        check(f"3-1.{i}", "3", f"C1 detections, {sys_}{'' if snr is None else f' {snr:g} dB'}, 256", rep, f"{k}/{n}",
              rel(d / f"{c1}.jsonl"), where)
    k = sum(cnt(c1, s, 20.0, 256)[0] for s in ("logistic", "henon"))
    check("3-2", "3", "C1 pooled 20 dB detection (logistic + Henon), 256", "200/200", f"{k}/200",
          rel(d / f"{c1}.jsonl"), where)
    for i, (m, rep) in enumerate([("baseline", ("9/100", "5/100", "123/200")),
                                  ("c2_kantz_m3_sat_iaaft", ("9/100", "0/100", "200/200")),
                                  ("c3_eps_m2_iaaft", ("5/100", "3/100", "200/200")),
                                  ("c4_zero_one_iaaft", ("0/100", "2/100", "11/200"))]):
        wn, ar = cnt(m, "white_noise", None, 256)[0], cnt(m, "ar1", None, 256)[0]
        p20 = sum(cnt(m, s, 20.0, 256)[0] for s in ("logistic", "henon"))
        check(f"3-3.{i}", "3", f"{m}: WN FP, AR(1) FP, pooled 20 dB (256)", rep, (f"{wn}/100", f"{ar}/100", f"{p20}/200"),
              rel(d / f"{m}.jsonl"), where)
    lo_wn = wilson(4, 100)
    check("3-4", "3", "C1 white-noise FP Wilson 95 % (lower, upper)", ("0.016", "0.098"), lo_wn,
          rel(d / f"{c1}.jsonl"), where, rule="round")

    # LLE bias = median(LLE) - reference over finite estimates
    REF = {"logistic": math.log(2), "henon": 0.4192}

    def bias(m, sys_, snr, n):
        v = np.array([r["result"]["lle"] for r in rec[m] if r["task"]["system"] == sys_
                      and r["task"]["snr_db"] == snr and r["task"]["window_length"] == n
                      and r["result"].get("lle") is not None and np.isfinite(r["result"]["lle"])])
        return float(np.median(v) - REF[sys_]), float(np.percentile(v, 25)), float(np.percentile(v, 75))

    for i, (sys_, n, rep) in enumerate([("logistic", 256, ("-0.001", "0.684", "0.700")),
                                        ("henon", 256, ("-0.005", "0.400", "0.432")),
                                        ("logistic", 512, ("-0.000", "0.689", "0.697")),
                                        ("henon", 512, ("-0.001", "0.407", "0.430"))]):
        check(f"3-5.{i}", "3", f"C1 LLE bias (median - ref) and [Q1, Q3], clean {sys_}, {n}", rep, bias(c1, sys_, None, n),
              rel(d / f"{c1}.jsonl"), where, rule="round")
    for i, (sys_, snr, rep) in enumerate([("logistic", 20.0, "-0.213"), ("henon", 20.0, "-0.068"),
                                          ("logistic", 10.0, "-0.428"), ("henon", 10.0, "-0.183")]):
        check(f"3-6.{i}", "3", f"C1 LLE bias {sys_} {snr:g} dB, 256", rep, bias(c1, sys_, snr, 256)[0],
              rel(d / f"{c1}.jsonl"), where, rule="round")
    b0 = bias("baseline", "logistic", None, 256)[0], bias("baseline", "henon", None, 256)[0]
    check("3-7", "3", "baseline LLE bias clean logistic / Henon, 256", ("-0.567", "-0.302"), b0,
          rel(d / "baseline.jsonl"), where, rule="round")
    s = abs(bias(c1, "logistic", None, 256)[0]) + abs(bias(c1, "henon", None, 256)[0])
    check("3-8", "3", "C1 clean |bias| sum (tie-break)", "0.0057", s, rel(d / f"{c1}.jsonl"), where, rule="round")
    s3 = abs(bias("c3_eps_m2_iaaft", "logistic", None, 256)[0]) + abs(bias("c3_eps_m2_iaaft", "henon", None, 256)[0])
    s2 = abs(bias("c2_kantz_m3_sat_iaaft", "logistic", None, 256)[0]) + abs(bias("c2_kantz_m3_sat_iaaft", "henon", None, 256)[0])
    check("3-9", "3", "tie-break order C1 < C3 < C2 (|bias| sums 0.0057 < 0.0173 < 0.0236)",
          ("0.0173", "0.0236"), (s3, s2), rel(d), where, rule="round")


# =================================================================================== Phase 4
def phase4():
    src = E / "phase4_upo/results/test/test.jsonl"
    where = "docs/PHASE4_UPO_VALIDATION.md 4.1, 4.5; PHASE10 S12 #4"
    R = {}
    for r in jl(src):
        R[r["task"]["task_id"]] = r
    R = list(R.values())
    gates = {"baseline": ("mcao_M7_S50", None), "c1": ("mcao_M7_S50", "ext_lead"),
             "c2": ("m2_M15_S50", "src_median"), "c3": ("m2_M7_S50", "src_trim10")}

    def peaks(r, meth):
        run, stat = gates[meth]
        ro = r["runs"].get(run, {})
        out = []
        for p in ro.get("peaks", []) or []:
            J = p.get("J")
            if J is None or not (J < 0.05):          # Level B recomputed from J (production alpha 0.05)
                continue
            if stat is not None and not (p.get(stat) is not None and p[stat] >= 1.2):
                continue
            out.append(p)
        return out

    flag = 0
    for r in R:
        for ro in r["runs"].values():
            for p in ro.get("peaks", []) or []:
                flag += bool(p.get("J") is not None and (p["J"] < 0.05) != bool(p["significant"]))
    check("4-0", "4", "stored Level-B flag equals J < 0.05 (all peaks, all runs)", 0, flag, rel(src), where)

    def k(meth, sys_, snr, n, rule="det"):
        rows = [r for r in R if r["task"]["system"] == sys_ and r["task"]["snr_db"] == snr
                and r["task"]["window_length"] == n]
        if rule == "det":
            return sum(bool(peaks(r, meth)) for r in rows)
        lle = [bool(r["lle_chaos_test"]["detected"]) for r in rows]
        upo = [bool(peaks(r, meth)) for r in rows]
        if rule == "and":
            return sum(a and b for a, b in zip(lle, upo))
        if rule == "or":
            return sum(a or b for a, b in zip(lle, upo))
        return sum(lle)

    prim = {"baseline": (12, 7, 89, 67, 70, 40, 110), "c1": (0, 3, 0, 0, 65, 34, 99),
            "c2": (0, 0, 0, 1, 150, 143, 293), "c3": (4, 3, 0, 1, 149, 140, 289)}
    for meth, rep in prim.items():
        v = [k(meth, "white_noise", None, 256), k(meth, "ar1", None, 256), k(meth, "sinusoid", None, 256),
             k(meth, "two_tone", None, 256), k(meth, "henon", 30.0, 256), k(meth, "henon", 20.0, 256)]
        v.append(v[4] + v[5])
        check(f"4-1.{meth}", "4", f"{meth}: WN, AR1, sinusoid, two-tone, Hen 30 dB, Hen 20 dB, pooled (/150, 256)",
              rep, tuple(v), rel(src), where)
    # localization error, clean Henon, C2 (nearest gated peak to the analytic fixed point)
    a, b = 1.4, 0.3
    xstar = (-(1 - b) + math.sqrt((1 - b) ** 2 + 4 * a)) / (2 * a)
    errs = []
    for r in R:
        if r["task"]["system"] == "henon" and r["task"]["snr_db"] is None and r["task"]["window_length"] == 256:
            pk = peaks(r, "c2")
            if pk:
                errs.append(min(abs(p["loc"] - xstar) for p in pk))
    check("4-2", "4", "C2 median clean-Henon localization error (256)", "0.0102", float(np.median(errs)), rel(src),
          where, rule="round", note=f"analytic fixed point x* = {xstar:.4f}")
    # combined detector (AND / OR with lle_chaos_test)
    for n, rep in ((256, (0, 0, 0, 0, 48, 150, 143, 74)), (512, (0, 0, 0, 0, 89, 150, 148, 83))):
        v = (k("c2", "white_noise", None, n, "and"), k("c2", "ar1", None, n, "and"), k("c2", "sinusoid", None, n, "and"),
             k("c2", "two_tone", None, n, "and"), k("c2", "logistic", 10.0, n, "and"), k("c2", "henon", 30.0, n, "and"),
             k("c2", "henon", 20.0, n, "and"), k("c2", "henon", 10.0, n, "and"))
        check(f"4-3.{n}", "4", f"C2 AND lle_chaos_test: WN, AR1, sin, 2-tone, log10, Hen30, Hen20, Hen10 (/150, {n})",
              rep, v, rel(src), where)
    v = (k("c2", "white_noise", None, 256, "or"), k("c2", "ar1", None, 256, "or"), k("c2", "sinusoid", None, 256, "or"),
         k("c2", "two_tone", None, 256, "or"))
    check("4-4", "4", "C2 OR lle_chaos_test: WN, AR1, sin, 2-tone (/150, 256)", (11, 8, 0, 1), v, rel(src), where)
    v = (k("c2", "white_noise", None, 256, "lle"), k("c2", "ar1", None, 256, "lle"),
         k("c2", "white_noise", None, 512, "lle"), k("c2", "ar1", None, 512, "lle"))
    check("4-5", "4", "lle_chaos_test alone: WN, AR1 at 256; WN, AR1 at 512 (/150)", (11, 8, 4, 4), v, rel(src), where)


# =================================================================================== Phase 5
def phase5():
    src = E / "phase5_rr/results/test/test.jsonl"
    where = "docs/PHASE5_RR_STRESS_TEST.md 5.1-5.2"
    R = list({r["task"]["task_id"]: r for r in jl(src)}.values())
    by = defaultdict(list)
    for r in R:
        by[r["task"]["condition"]].append(r["m2"])

    def tab(c):
        rows = by[c]
        lle = [bool(x["lle"]) and not x.get("error") for x in rows]
        upo = [bool(x["upo"]) and not x.get("error") for x in rows]
        return (sum(a and b for a, b in zip(lle, upo)), sum(lle), sum(upo), sum(a or b for a, b in zip(lle, upo)),
                sum(x.get("error") is not None for x in rows), len(rows))

    nulls = {"N1_linear_rr": (1, 20, 4, 23, 0, 300), "N2_power_law": (1, 8, 4, 11, 0, 300),
             "N3_linear_rr_trend": (0, 24, 2, 26, 0, 300), "N4_linear_rr_step": (0, 14, 3, 17, 0, 300),
             "N5_linear_rr_warped": (1, 53, 6, 58, 0, 300), "N6_noisy_rsa": (0, 17, 1, 18, 0, 300)}
    for c, rep in nulls.items():
        check(f"5-1.{c[:2]}", "5", f"{c}: AND, LLE, UPO, OR, errors, N (m = 2)", rep, tab(c), rel(src), where)
    check("5-2", "5", "PASS: AND <= 21/300 on every N1-N6", True,
          all(tab(c)[0] <= 21 for c in nulls), rel(src), where, rule="text")
    s2 = {"S2_ectopic_2pct": (1, 194, 1, 194, 0, 200), "S2_ectopic_5pct": (2, 195, 2, 195, 0, 200),
          "S2_ectopic_10pct": (7, 197, 7, 197, 0, 200)}
    for c, rep in s2.items():
        check(f"5-3.{c}", "5", f"{c}: AND, LLE, UPO, OR, errors, N", rep, tab(c), rel(src), where)
    other = {"P3_henon_rr_trend": 41, "P3_logistic_rr_trend": 86, "G2_rossler_flow": 0, "G3_mackey_glass": 0,
             "P1_henon_rr": 198}
    for c, rep in other.items():
        check(f"5-4.{c}", "5", f"{c}: AND (m = 2)", rep, tab(c)[0], rel(src), where)
    check("5-5", "5", "total analysis errors (of 5,400)", 3, sum(tab(c)[4] for c in by), rel(src), where)
    g = [sum(bool(r[m]["lle"]) and bool(r[m]["upo"]) for r in R if r["task"]["condition"] == c)
         for c in ("G2_rossler_flow", "G3_mackey_glass") for m in ("m2", "m3", "m4")]
    check("5-6", "5", "K1-type AND on Rossler (G2) and Mackey-Glass (G3) at m = 2, 3, 4", (0, 0, 0, 0, 0, 0), tuple(g),
          rel(src), "docs/PHASE10_FINAL_RESULTS.md S12 #8")


# =================================================================================== Phase 6
def phase6():
    src = E / "phase6_robust/results/test/test.jsonl"
    where = "docs/PHASE6_ROBUSTNESS.md 3.1, 3.2, 4.4"
    R = list({r["task"]["task_id"]: r for r in jl(src)}.values())
    by = defaultdict(list)
    for r in R:
        by[r["task"]["condition"]].append(r["methods"])

    def t(c, m="baseline_k"):
        rows = [x[m] for x in by[c]]
        lle = [bool(x["lle"]) for x in rows]
        upo = [bool(x["upo"]) for x in rows]
        andv = [a and b for a, b in zip(lle, upo)]
        stored = sum(bool(x["and"]) for x in rows)
        return sum(andv), sum(lle), sum(upo), sum(a or b for a, b in zip(lle, upo)), len(rows), stored

    partb = {"N1_linear_rr": (0, 36, 0, 36, 400), "S2_ectopic_2pct": (1, 289, 1, 289, 300),
             "S2_ectopic_5pct": (5, 290, 5, 290, 300), "S2_ectopic_10pct": (7, 293, 7, 293, 300),
             "E1_bigeminy": (0, 37, 0, 37, 300), "E2_trigeminy": (0, 103, 0, 103, 300),
             "E3_couplets_5pct": (1, 292, 1, 292, 300), "E3_couplets_10pct": (3, 298, 3, 298, 300),
             "E4_runs": (0, 275, 0, 275, 300), "E5_atrial_5pct": (4, 269, 4, 269, 300),
             "E5_atrial_10pct": (11, 254, 12, 255, 300), "E6_ectopic10_trend": (10, 299, 10, 299, 300)}
    bad_and = 0
    for c, rep in partb.items():
        v = t(c)
        bad_and += v[0] != v[5]
        check(f"6-1.{c}", "6", f"BASELINE-K {c}: AND, LLE, UPO, OR, N", rep, v[:5], rel(src), where + " (3.1)")
    edited = {"S2_ectopic_5pct": (0, 23, 4), "S2_ectopic_10pct": (2, 46, 4), "E3_couplets_10pct": (1, 27, 24),
              "E2_trigeminy": (13, 291, 13), "E1_bigeminy": (0, 0, 0)}
    for c, rep in edited.items():
        v = t(c + "_edited")
        check(f"6-2.{c}", "6", f"BASELINE-K {c} edited: AND, LLE, UPO (/300)", rep, v[:3], rel(src), where + " (3.2)")
    worst = max(t(c)[0] / t(c)[4] for c in partb if c != "N1_linear_rr")
    check("6-3", "6", "max raw ectopy AND rate (%)", "3.7", worst, rel(src),
          where + "; PHASE10 S12 #6", rule="round", scale=100.0)
    # Part C (C4)
    meth = {"baseline_k": "BASELINE-K", "d1_linear_g05": "D1", "d2_linear_g07": "D2", "d3_smoothprior300_g07": "D3"}
    spec = {"baseline_k": (0, 0, 1, 1, 2, 0, 1, 5, 7), "d1_linear_g05": (0, 0, 0, 0, 2, 0, 1, 5, 7),
            "d2_linear_g07": (0, 0, 0, 0, 2, 0, 1, 5, 7), "d3_smoothprior300_g07": (0, 0, 0, 0, 2, 0, 1, 5, 7)}
    nulls = ["N1_linear_rr", "N2_power_law", "N3_linear_rr_trend", "N4_linear_rr_step", "N5_linear_rr_warped",
             "N6_noisy_rsa", "S2_ectopic_2pct", "S2_ectopic_5pct", "S2_ectopic_10pct"]
    for m, rep in spec.items():
        check(f"6-4.{m}", "6", f"{meth[m]} AND specificity N1-N6, S2 2/5/10 %", rep, tuple(t(c, m)[0] for c in nulls),
              rel(src), where + " (4.4)")
    powc = ["P1_henon_rr", "P2_logistic_rr", "P4_henon_rr_30dB", "P4_henon_rr_20dB", "P4_logistic_rr_30dB",
            "P4_logistic_rr_20dB", "G1_lorenz_maxima"]
    pw = {"baseline_k": (300, 294, 300, 296, 299, 284, 294), "d1_linear_g05": (300, 293, 300, 296, 298, 282, 279),
          "d2_linear_g07": (300, 294, 300, 296, 299, 284, 291), "d3_smoothprior300_g07": (300, 294, 300, 296, 299, 284, 292)}
    for m, rep in pw.items():
        check(f"6-5.{m}", "6", f"{meth[m]} AND power P1, P2, P4 H30/H20/L30/L20, G1 (/300)", rep,
              tuple(t(c, m)[0] for c in powc), rel(src), where + " (4.4)")
    partb_conds = ["E1_bigeminy", "E2_trigeminy", "E3_couplets_5pct", "E3_couplets_10pct", "E4_runs", "E5_atrial_5pct",
                   "E5_atrial_10pct", "E6_ectopic10_trend", "S2_ectopic_5pct_edited", "S2_ectopic_10pct_edited",
                   "E1_bigeminy_edited", "E2_trigeminy_edited", "E3_couplets_10pct_edited"]
    summ = {"baseline_k": (74, 121, 195, 45), "d1_linear_g05": (300, 281, 581, 44), "d2_linear_g07": (298, 273, 571, 45),
            "d3_smoothprior300_g07": (295, 244, 539, 45)}
    for m, rep in summ.items():
        h, lg = t("P3_henon_rr_trend", m)[0], t("P3_logistic_rr_trend", m)[0]
        check(f"6-6.{m}", "6", f"{meth[m]} pooled P3 AND (Hen, log, total /600) and pooled Part B AND (/3,900)", rep,
              (h, lg, h + lg, sum(t(c, m)[0] for c in partb_conds)), rel(src), where + " (4.4)")
    check("6-7", "6", "D1 fails (G1 more than 9 below BASELINE-K); D2, D3 pass the power rule", (True, True, True),
          (t("G1_lorenz_maxima", "baseline_k")[0] - t("G1_lorenz_maxima", "d1_linear_g05")[0] > 9,
           all(t(c, "baseline_k")[0] - t(c, "d2_linear_g07")[0] <= 9 for c in powc),
           all(t(c, "baseline_k")[0] - t(c, "d3_smoothprior300_g07")[0] <= 9 for c in powc)), rel(src), where, rule="exact")
    check("6-8", "6", "pooled P3 power BASELINE-K -> D2 (%)", ("32.5", "95.2"), (195 / 600, 571 / 600), rel(src),
          "PHASE6 4.4; PHASE10 S12 #7", rule="round", scale=100.0,
          note="fractions recomputed above (6-6); here rounded")
    check("6-9", "6", "LLE-alone N3, N4: BASELINE-K -> D2 (/400)", (23, 40, 17, 40),
          (t("N3_linear_rr_trend")[1], t("N3_linear_rr_trend", "d2_linear_g07")[1],
           t("N4_linear_rr_step")[1], t("N4_linear_rr_step", "d2_linear_g07")[1]), rel(src), where + " (4.4)")
    check("6-10", "6", "stored AND flag equals LLE AND UPO (BASELINE-K, Part B conditions)", 0, bad_and, rel(src), where)


# =================================================================================== Phase 7
def phase7():
    d = E / "phase7_mitbih/results"
    where = "docs/PHASE7_MITBIH_RESULTS.md"
    raw = pd.DataFrame(jl(d / "run/raw.jsonl"))
    A, N = raw[raw.label == 1], raw[raw.label == 0]
    check("7-1", "7", "Q1 raw AND: abnormal, normal; windows; subjects; subjects with abnormal windows",
          ("0/94", "0/211", 305, 43, 19),
          (f"{int(A.and_detected.sum())}/{len(A)}", f"{int(N.and_detected.sum())}/{len(N)}", len(raw),
           raw.subject.nunique(), A.subject.nunique()), rel(d / "run/raw.jsonl"), where + " 4")
    check("7-2", "7", "Wilson upper bound, abnormal / normal (%)", ("3.9", "1.8"),
          (wilson(0, 94)[1], wilson(0, 211)[1]), rel(d / "run/raw.jsonl"), where + " 4", rule="round", scale=100.0)
    comp = (f"{int(A.lle_detected.sum())}/94", f"{int(N.lle_detected.sum())}/211", f"{int(A.upo_detected.sum())}/94",
            f"{int(N.upo_detected.sum())}/211", f"{int((A.lle_detected | A.upo_detected).sum())}/94",
            f"{int((N.lle_detected | N.upo_detected).sum())}/211")
    check("7-3", "7", "components: LLE A, LLE N, UPO A, UPO N, OR A, OR N",
          ("47/94", "27/211", "3/94", "4/211", "50/94", "31/211"), comp, rel(d / "run/raw.jsonl"), where + " 7")
    check("7-4", "7", "UPO-gated windows that coincide with an LLE detection (raw)", 0,
          int((raw.lle_detected & raw.upo_detected).sum()), rel(d / "run/raw.jsonl"), where + " 7")
    ann = pd.DataFrame(jl(d / "run/annotation.jsonl"))
    aA, aN = ann[ann.label == 1], ann[ann.label == 0]
    check("7-5", "7", "sensitivity (annotation beat times): AND abnormal, normal", ("6/120", "0/248"),
          (f"{int(aA.and_detected.sum())}/{len(aA)}", f"{int(aN.and_detected.sum())}/{len(aN)}"),
          rel(d / "run/annotation.jsonl"), where + " 9")
    # Q3: refit the leave-one-subject-out L2 logistic models with own code
    HRV = ["sdnn_ms", "rmssd_ms", "pnn50"]
    models = {"M0": HRV, "M1": ["lle_z"], "M2": ["upo_score"], "M3": ["lle_z", "upo_score"],
              "M4": HRV + ["lle_z", "upo_score"]}
    cc = raw.dropna(subset=HRV + ["lle_z", "upo_score"]).reset_index(drop=True)
    aucs = {}
    for m, f in models.items():
        aucs[m] = auc_mw(loso_oof(cc, f), cc.label.values)
    check("7-6", "7", "Q3 out-of-fold AUC M0..M4 (own LOSO L2 logistic refit)",
          ("0.812", "0.685", "0.193", "0.675", "0.821"), tuple(aucs[m] for m in models), rel(d / "run/raw.jsonl"),
          where + " 6.1", rule="round", note=f"{len(cc)} complete-case windows")
    check("7-7", "7", "Q3 ΔAUC M3-M1, M3-M2, M4-M0", ("-0.010", "+0.482", "+0.009"),
          (aucs["M3"] - aucs["M1"], aucs["M3"] - aucs["M2"], aucs["M4"] - aucs["M0"]), rel(d / "run/raw.jsonl"),
          where + " 6.1", rule="round")
    oof = pd.read_csv(d / "analysis/raw_oof_predictions.csv")
    a2 = tuple(auc_mw(oof[f"oof_{m}"], oof.label) for m in models)
    check("7-8", "7", "Q3 AUC from the stored out-of-fold predictions (own AUC code)",
          ("0.812", "0.685", "0.193", "0.675", "0.821"), a2, rel(d / "analysis/raw_oof_predictions.csv"),
          where + " 6.1", rule="round")
    # spike-in table
    sp = pd.DataFrame(jl(d / "spike_in/spike_in.jsonl"))
    rep = {("henon", "a_clean"): ("210/211", "94/94"), ("henon", "c_ectopy"): ("205/211", "52/94"),
           ("henon", "d_ectopy_jitter"): ("208/211", "53/94"), ("logistic", "a_clean"): ("208/211", "94/94"),
           ("logistic", "c_ectopy"): ("204/211", "74/94"), ("logistic", "d_ectopy_jitter"): ("201/211", "79/94")}
    for (s, v), r in rep.items():
        g = sp[(sp.system == s) & (sp.variant == v)]
        gn, ga = g[g.label == 0], g[g.label == 1]
        check(f"7-9.{s}.{v}", "7", f"spike-in AND {s} {v}: normal, abnormal",
              r, (f"{int(gn.and_detected.sum())}/{len(gn)}", f"{int(ga.and_detected.sum())}/{len(ga)}"),
              rel(d / "spike_in/spike_in.jsonl"), where + " 12.2")
    check("7-10", "7", "spike-in: series, errors", (1830, 0), (len(sp), int(sp.error.notna().sum())),
          rel(d / "spike_in/spike_in.jsonl"), where + " 12.2")
    bins = [(0.10, 0.20), (0.20, 0.30), (0.30, 0.50), (0.50, 1.01)]
    reported = [(37, 26, 26, 35, 34), (28, 16, 18, 25, 27), (21, 8, 7, 9, 12), (8, 2, 2, 5, 6)]
    for (lo, hi), r in zip(bins, reported):
        ab = sp[(sp.label == 1)]

        def k(s, v):
            g = ab[(ab.system == s) & (ab.variant == v) & (ab.abnormal_fraction >= lo) & (ab.abnormal_fraction < hi)]
            return int(g.and_detected.sum()), len(g)

        nwin = k("henon", "c_ectopy")[1]
        check(f"7-11.{lo:.2f}", "7", f"spike-in abnormal windows, fraction [{lo:.2f}, {hi:.2f}): windows, Hen c/d, log c/d", r,
              (nwin, k("henon", "c_ectopy")[0], k("henon", "d_ectopy_jitter")[0], k("logistic", "c_ectopy")[0],
               k("logistic", "d_ectopy_jitter")[0]), rel(d / "spike_in/spike_in.jsonl"), where + " 12.2",
              note="bin edges: lower inclusive, upper exclusive (last bin > 0.50)")
    au = json.load(open(d / "spike_in/audit.json"))
    allok = all(all(v for k, v in w.items() if isinstance(v, bool)) for w in au["windows"])
    check("7-12", "7", "code audit: windows checked, discrepancies, every boolean check true", (10, 0, True),
          (au["n_checked"], len(au["discrepancies"]), allok), rel(d / "spike_in/audit.json"), where + " 12.1")


# =================================================================================== Phase 8
def phase8():
    d = E / "phase8_cardiac/results/test"
    where = "docs/PHASE8_CARDIAC_CHAOS.md 7.1-7.2"
    L = {"c1_frozen512": 512, "c2_nlp_iaaft512": 512, "c3_nlp_ep256": 256, "c4_titration256": 256,
         "baseline_frozen256": 256}
    rows = {}
    for n in (256, 512):
        for r in jl(d / f"test_{n}.jsonl"):
            rows[r["task_id"]] = r
    VI = ("vi_S2_5", "vi_E1", "vi_E3_10")
    out = {}
    for m, n in L.items():
        rr = [r for r in rows.values() if r["n"] == n and m in r["results"]]
        spec = defaultdict(lambda: [0, 0])
        for r in rr:
            if r["group"] in ("null", "noncha"):
                spec[r["regime"]][0] += bool(r["results"][m].get("detected"))
                spec[r["regime"]][1] += 1
        failed = sorted(c for c, (k, nn) in spec.items() if k > math.floor(0.07 * nn))
        vi = sum(bool(r["results"][m].get("detected")) for r in rr if r["group"] == "chaos_primary" and r["variant"] in VI)
        nvi = sum(1 for r in rr if r["group"] == "chaos_primary" and r["variant"] in VI)
        iv = sum(bool(r["results"][m].get("detected")) for r in rr if r["group"] == "chaos_primary" and r["variant"] == "iv_q")
        niv = sum(1 for r in rr if r["group"] == "chaos_primary" and r["variant"] == "iv_q")
        pn = sum(spec[c][0] for c in spec if c[0] in "ENS" and "_" in c)
        pr = sum(spec[c][0] for c in spec if ":" in c)
        out[m] = (len(failed) == 0, len(failed), f"{vi}/{nvi}", f"{iv}/{niv}", pn, pr, spec)
    rep = {"c1_frozen512": (True, 0, "66/990", "10/330", 8, 0), "c2_nlp_iaaft512": (False, 2, "220/990", "141/330", 2, 45),
           "c3_nlp_ep256": (False, 1, "148/990", "73/330", 1, 63), "c4_titration256": (False, 30, "966/990", "330/330", 1187, 1452),
           "baseline_frozen256": (True, 0, "52/990", "10/330", 7, 0)}
    for m, r in rep.items():
        check(f"8-1.{m}", "8", f"{m}: PASS, n failed, TEST chaotic (vi), (iv), pooled 17 nulls, pooled 18 non-chaotic",
              r, out[m][:6], rel(d / "test_*.jsonl"), where)
    passing = [m for m in ("c1_frozen512", "c2_nlp_iaaft512", "c3_nlp_ep256", "c4_titration256") if out[m][0]]
    check("8-2", "8", "winner", "c1_frozen512", passing[0] if len(passing) == 1 else str(passing), rel(d), where)
    s2, s3 = out["c2_nlp_iaaft512"][6], out["c3_nlp_ep256"][6]
    check("8-3", "8", "failing cells: C2 vdP (2, 5.6), (5.45, 5.6); C3 vdP (2, 5.6) (/100)", (24, 19, 62),
          (s2["coupled_vdp:rho=2.0,omega=5.6"][0], s2["coupled_vdp:rho=5.45,omega=5.6"][0],
           s3["coupled_vdp:rho=2.0,omega=5.6"][0]), rel(d), where)
    s4 = out["c4_titration256"][6]
    check("8-4", "8", "C4 titration on ectopy nulls: E1, E2, E3 5/10, E4, E5 5/10, E6, S2 5/10, N4, N5, S1 (/100)",
          (85, 100, 100, 100, 100, 89, 58, 100, 100, 100, 98, 100, 43),
          tuple(s4[c][0] for c in ("E1_bigeminy", "E2_trigeminy", "E3_couplets_5pct", "E3_couplets_10pct", "E4_runs",
                                   "E5_atrial_5pct", "E5_atrial_10pct", "E6_ectopic10_trend", "S2_ectopic_5pct",
                                   "S2_ectopic_10pct", "N4_linear_rr_step", "N5_linear_rr_warped", "S1_setar")),
          rel(d), where)
    fl = [r for r in rows.values() if r["group"] == "flows" and r["n"] == 512]
    check("8-5", "8", "C1 (K1) on Rossler and Mackey-Glass flows (/200)", 0,
          sum(bool(r["results"]["c1_frozen512"]["detected"]) for r in fl), rel(d / "test_512.jsonl"),
          "PHASE10 S12 #8")


# =================================================================================== Phase 9
def phase9():
    src = E / "phase9_waveform/results/test/test9.jsonl"
    where = "docs/PHASE9_WAVEFORM_NOISE_ROBUST.md 6"
    R = list({r["task_id"]: r for r in jl(src)}.values())
    Z, G, ZRR, GRR = 9.4521, 0.25, 12.1603, 0.85
    PASS = ("null_input", "null_real", "noncha_input", "noncha_real")
    cands = ("k1_frozen512", "k2_mnlp_ann", "k3_mnlp_growth_ann", "k4_mnlp_growth_rr", "k5_mnlp_fsle_ann")

    def det(r, c):
        res = (r.get("results") or {}).get(c, {})
        if c == "k3_mnlp_growth_ann":        # recomputed from the stored statistics and the frozen thresholds
            return bool(res.get("analysable") and res.get("zmax") is not None and res["zmax"] >= Z and res["G"] >= G)
        if c == "k4_mnlp_growth_rr":
            return bool(res.get("analysable") and res.get("zmax") is not None and res["zmax"] >= ZRR and res["G"] >= GRR)
        return bool(res.get("detected", False))

    flagdiff = sum(det(r, c) != bool(((r.get("results") or {}).get(c, {})).get("detected", False))
                   for r in R for c in ("k3_mnlp_growth_ann", "k4_mnlp_growth_rr"))
    check("9-0", "9", "stored k3/k4 flags equal (zmax >= Z and G >= G_min) recomputed", 0, flagdiff, rel(src), where)
    res = {}
    for c in cands:
        spec = defaultdict(lambda: [0, 0])       # preregistered conditions (PREREGISTRATION.md 3-4: 101)
        spec9 = defaultdict(lambda: [0, 0])      # decide9.py grouping (group|name)
        for r in R:
            if r["group"] in PASS:
                key9 = r["group"] + "|" + r["name"]
                # a KTz regime is (pacing period P, input pattern): ktz_labels.json; its name carries only P
                key = key9 + ("|" + r["ectopy"] if r.get("family") == "ktz" else "")
                spec[key][0] += det(r, c)
                spec[key][1] += 1
                spec9[key9][0] += det(r, c)
                spec9[key9][1] += 1
        failed = [k for k, (a, n) in spec.items() if a > math.floor(0.07 * n)]
        failed9 = [k for k, (a, n) in spec9.items() if a > math.floor(0.07 * n)]
        prim = [r for r in R if r["group"] == "chaos_primary"]
        res[c] = (len(failed) == 0, len(failed), sum(det(r, c) for r in prim), sum(v[0] for v in spec.values()),
                  len(spec), sum(v[1] for v in spec.values()), len(prim), spec, len(spec9), len(failed9),
                  len(failed9) == 0)
    rep = {"k1_frozen512": (False, 4, 74, 97), "k2_mnlp_ann": (False, 12, 304, 607),
           "k3_mnlp_growth_ann": (True, 0, 90, 0), "k4_mnlp_growth_rr": (True, 0, 2, 1), "k5_mnlp_fsle_ann": (True, 0, 0, 3)}
    for c, r in rep.items():
        check(f"9-1.{c}", "9", f"{c}: PASS, n failed (preregistered 101 conditions), primary chaotic (/1,320), "
              "PASS-condition detections (/10,100)", r, res[c][:4], rel(src), where,
              note=f"decide9.py grouping ({res[c][8]} conditions): PASS {res[c][10]}, {res[c][9]} failed")
    check("9-2", "9", "conditions, PASS windows, primary windows, total windows",
          (101, 10100, 1320, 12290), (res["k3_mnlp_growth_ann"][4], res["k3_mnlp_growth_ann"][5],
                                      res["k3_mnlp_growth_ann"][6], len(R)), rel(src), where)
    nerr = sum(1 for r in R if r.get("error"))
    check("9-3", "9", "window errors", 63, nerr, rel(src), where)
    prim = [r for r in R if r["group"] == "chaos_primary"]
    fam = defaultdict(lambda: [0, 0, 0])
    ect = Counter()
    for r in prim:
        fam[r["family"]][0] += det(r, "k3_mnlp_growth_ann")
        fam[r["family"]][1] += det(r, "k1_frozen512")
        fam[r["family"]][2] += 1
        if r["family"] == "coupled_vdp":
            ect[r["ectopy"]] += det(r, "k3_mnlp_growth_ann")
    check("9-4", "9", "primary by family k3/k1/N: coupled vdP, phase-reset, KTz",
          ((73, 41, 360), (17, 33, 630), (0, 0, 330)),
          (tuple(fam["coupled_vdp"]), tuple(fam["phase_reset"]), tuple(fam.get("ktz", fam.get("KTz", [None] * 3)))),
          rel(src), where, note=f"families present: {sorted(fam)}")
    check("9-5", "9", "k3 coupled vdP by ectopy S2_5, E3_10, E1", (46, 27, 0),
          (ect["S2_5"], ect["E3_10"], ect["E1"]), rel(src), where)
    k1v = Counter()
    k1n = Counter()
    for r in R:
        if r["group"] == "noncha_real" and r["family"] == "coupled_vdp":
            k1v[r["ectopy"]] += det(r, "k1_frozen512")
            k1n[r["ectopy"]] += 1
    check("9-6", "9", "k1 on non-chaotic vdP (realistic level) by ectopy: E1, E3_10, S2_5, none",
          ("43/150", "24/150", "11/150", "0/150"),
          tuple(f"{k1v[e]}/{k1n[e]}" for e in ("E1", "E3_10", "S2_5", "none")), rel(src), where)
    worst = max(res["k1_frozen512"][7][k][0] for k in res["k1_frozen512"][7]
                if res["k1_frozen512"][7][k][0] > math.floor(0.07 * res["k1_frozen512"][7][k][1]))
    check("9-7", "9", "k1 worst failing condition (/100)", 45, worst, rel(src), where + "; PHASE10 S12 #6")
    # Part F
    pf = E / "phase9_waveform/results/partf/partf_windows.jsonl"
    W = pd.DataFrame(jl(pf))
    W["det"] = W.analysable.astype(bool) & (W.zmax >= Z) & (W.G >= G)
    rows = []
    for db in ("nsrdb", "chfdb"):
        g = W[W.db == db]
        rows.append((g.subject.nunique(), len(g), int(g.det.sum()), int(g.analysable.sum())))
    check("9-8", "9", "Part F nsrdb, chfdb: subjects, windows, detected, analysable", ((18, 180, 1, 180), (15, 150, 0, 147)),
          tuple(rows), rel(pf), where + " 7 (Part F)")
    sr = W[W.db == "nsrdb"].groupby("subject").det.mean()
    check("9-9", "9", "Part F detected subject (nsrdb)", "16272", ",".join(sorted(W[W.det].subject.astype(str))),
          rel(pf), where + " 7")
    check("9-10", "9", "Part F NSR subject-level rate; W2 difference CHF - NSR", ("0.006", "-0.006"),
          (float(sr.mean()), float(W[W.db == "chfdb"].groupby("subject").det.mean().mean() - sr.mean())), rel(pf),
          where + " 7", rule="round")
    mx = E / "phase9_waveform/results/partf/mitbih_exploratory.jsonl"
    M = pd.DataFrame(jl(mx))
    M["det"] = M.analysable.astype(bool) & (M.zmax >= Z) & (M.G >= G)
    check("9-11", "9", "EXPLORATORY MIT-BIH: detected/windows, analysable, detected record",
          ("1/190", 141, "123"), (f"{int(M.det.sum())}/{len(M)}", int(M.analysable.sum()),
                                  ",".join(M[M.det].record.astype(str))), rel(mx), where + " 7")


# =================================================================================== Phase 10
P10 = E / "phase10_final/results"


def p10_real(phase="conf", exclude=()):
    R = list({r["id"]: r for r in jl(P10 / phase / "real.jsonl") if r.get("subject") not in exclude}.values())
    rows = []
    for r in R:
        if r.get("error"):
            continue
        clock = r["clock_hour"]
        if clock is None:
            clock = (9.9 + r["t_start"] / 3600.0) % 24.0
        rows.append({"subject": r["subject"], "group": r["group"], "j": r["j"], "burden": r["burden"],
                     "n_masked": r["n_masked"], "clock": clock, "sdnn": r["sdnn"],
                     "K1": bool(r["K1"]["lle_detected"] and r["K1"]["upo_detected"]),
                     "K1_stored": bool(r["K1"]["and_detected"]),
                     "LLE": bool(r["K1"]["lle_detected"]), "UPO": bool(r["K1"]["upo_detected"]),
                     "K3": bool(r["K3"].get("analysable") and r["K3"]["zmax"] >= 9.4521 and r["K3"]["G"] >= 0.25),
                     "K3_stored": bool(r["K3"]["detected"]), "K3_z": r["K3"].get("zmax"), "K3_G": r["K3"].get("G"),
                     "K3_an": bool(r["K3"].get("analysable")), "K3_nm": r["K3"].get("n_masked"),
                     "K4": bool(r["K4"]["detected"]), "K3RR": bool(r["K3RR"]["detected"]),
                     "TIT": bool(r["TIT_raw"]["NL"] > 0), "TIT_stored": bool(r["TIT_raw"]["positive"]),
                     "TITm": bool(r["TIT_masked"].get("positive", False)),
                     "TITm_ok": bool(r["TIT_masked"].get("analysable", False)),
                     "TITe": bool(r["TIT_edited"].get("positive", False)),
                     "TITe_ok": bool(r["TIT_edited"].get("eligible", False)),
                     "K1e_ok": r["K1_edited"].get("eligible", True),
                     "K1e": bool(r["K1_edited"].get("and_detected", False)),
                     "LLEe": bool(r["K1_edited"].get("lle_detected", False)),
                     "UPOe": bool(r["K1_edited"].get("upo_detected", False))})
    df = pd.DataFrame(rows)
    df["chf"] = (df.group == "CHF").astype(int)
    df["x"] = np.log2(1 + df.burden)
    df["night"] = (df.clock >= 0) & (df.clock < 5)
    return df


def p10_seg(phase="conf", exclude=()):
    R = list({r["id"]: r for r in jl(P10 / phase / "seg.jsonl") if r.get("subject") not in exclude}.values())
    rows = []
    for r in R:
        for s in r["segments"]:
            row = {"subject": r["subject"], "group": r["group"], "gap": s["gap"]}
            if not s["gap"]:
                clock = s.get("clock_hour")
                if clock is None:
                    clock = (9.9 + s["t_start"] / 3600.0) % 24.0
                row.update(burden=s["burden"], clock=clock, raw=bool(s["TIT_raw"]["NL"] > 0), raw_NL=s["TIT_raw"]["NL"],
                           masked=bool(s["TIT_masked"].get("positive", False)),
                           masked_ok=bool(s["TIT_masked"].get("analysable", False)),
                           edited=bool(s["TIT_edited"].get("positive", False)),
                           edited_ok=bool(s["TIT_edited"].get("eligible", False)),
                           wu=bool(s["TIT_wu"].get("positive", False)),
                           wu_ok=bool(s["TIT_wu"].get("eligible", True) and s["TIT_wu"].get("analysable", True)),
                           wu_NL=s["TIT_wu"].get("NL", 0.0))
            rows.append(row)
    return pd.DataFrame(rows)


def subj_rate(df, col):
    return df.groupby("subject")[col].mean()


def p10_q1(df):
    where = "docs/PHASE10_FINAL_RESULTS.md 5.1-5.2 (P1)"
    src = P10 / "conf/spike.jsonl"
    S = list({r["id"]: r for r in jl(src)}.values())
    nser = sum(len(r["designs"]) for r in S if "designs" in r)
    errs = sum(1 for r in S if r.get("error"))
    check("10-Q1-0", "10", "spike-in tasks, spiked series, errors", (1328, 10624, 0), (len(S), nser, errs), rel(src),
          "PHASE10 5 (Run)")
    cnt = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    lam = {}
    for r in S:
        key = (r["family"], r["param"])
        lam[key] = r["lam"]
        for dname, v in r["designs"].items():
            fkey = "rep" if dname == "replace" else float(v["f"])
            for det in ("K3", "K1"):
                hit = bool(v["K3"]["detected"]) if det == "K3" else bool(v["K1"]["and_detected"])
                cnt[(det, key)][fkey][0] += hit
                cnt[(det, key)][fkey][1] += 1
    FS = (0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9)
    # real-window bounds
    U = {}
    Ualt = {}
    rng = np.random.default_rng(VSEED)
    for det, kr, ur, cpr, bootr in (("K3", "8/988", "0.0146", "0.0146", "0.0141"),
                                    ("K1", "16/988", "0.0245", "0.0245", "0.0232")):
        k, n = int(df[det].sum()), len(df)
        ucp = cp_upper_bisect(k, n)
        rates = subj_rate(df, det).values
        ub = pct(cluster_boot_means(rates, 10000, rng), 95)
        U[det] = max(ucp, float(bootr))          # the reported (seeded) bootstrap term; CP recomputed exactly
        seeds = [pct(cluster_boot_means(rates, 10000, np.random.default_rng(VSEED + 100 + s)), 95) for s in range(50)]
        Ualt[det] = max(ucp, max(seeds))
        check(f"10-Q1-cp.{det}", "10", f"{det}: one-sided 95 % Clopper-Pearson bound on detected / windows", cpr, ucp,
              rel(P10 / "conf/real.jsonl"), where + "; tables.md", rule="round",
              note="own bisection on the binomial CDF")
        check(f"10-Q1-boot.{det}", "10", f"{det}: subject-cluster bootstrap 95th percentile of the mean subject rate",
              bootr, ub, rel(P10 / "conf/real.jsonl"), where + "; tables.md", rule="mc", tol=0.0015,
              note=f"independent bootstrap; over 50 further seeds the term takes values {min(seeds):.5f}-{max(seeds):.5f}")
        check(f"10-Q1-U.{det}.max", "10", f"{det}: U = max(bootstrap term, CP term) with the reported terms", ur,
              max(ucp, float(bootr)), rel(P10 / "conf/real.jsonl"), where, rule="round")
        check(f"10-Q1-U.{det}", "10", f"{det}: real detections / windows", kr, f"{k}/{n}",
              rel(P10 / "conf/real.jsonl"), where)
    kl = int(df.LLE.sum())
    ulle = max(cp_upper_bisect(kl, len(df)), pct(cluster_boot_means(subj_rate(df, "LLE").values, 10000, rng), 95))
    check("10-Q1-U.LLE", "10", "LLE alone: real detections", "401/988", f"{kl}/{len(df)}",
          rel(P10 / "conf/real.jsonl"), where)
    check("10-Q1-U2.LLE", "10", "LLE alone: U (bootstrap term decides; MC)", "0.45", ulle,
          rel(P10 / "conf/real.jsonl"), where, rule="mc", tol=0.01)

    # reported rows (detected / 83 at f; replacement; f_min 0.05; f_min 0.20)
    K3rep = {("henon", "a1.08"): ((0, 1, 0, 0, 0, 0, 0), 0, None, None),
             ("henon", "a1.14"): ((0, 0, 0, 2, 2, 9, 8), 4, None, 0.7),
             ("henon", "a1.22"): ((0, 0, 3, 6, 32, 58, 76), 33, 0.5, 0.5),
             ("henon", "a1.40"): ((0, 1, 4, 9, 29, 58, 76), 69, 0.5, 0.3),
             ("logistic", "r3.58"): ((0,) * 7, 0, None, None),
             ("logistic", "r3.65"): ((0, 0, 0, 0, 0, 1, 0), 0, None, None),
             ("logistic", "r3.88"): ((1, 1, 1, 5, 28, 59, 77), 62, 0.5, 0.5),
             ("logistic", "r4.00"): ((0, 0, 0, 1, 9, 25, 73), 76, 0.7, 0.5),
             ("coupled_vdp", "rho8_om3.3"): ((1, 2, 5, 13, 33, 48, 71), 31, 0.5, 0.3),
             ("coupled_vdp", "rho6_om3.3"): ((0, 0, 0, 3, 10, 13, 23), 45, None, 0.5),
             ("coupled_vdp", "rho2_om2.7"): ((0, 0, 0, 2, 6, 20, 63), 60, 0.9, 0.7),
             ("coupled_vdp", "rho6_om4.0"): ((0, 1, 2, 6, 13, 15, 27), 25, 0.9, 0.5)}
    K1rep = {("henon", "a1.08"): ((2, 2, 0, 0, 0, 0, 0), 0, None, None),
             ("henon", "a1.14"): ((0, 0, 4, 0, 7, 11, 24), 19, None, 0.7),
             ("henon", "a1.22"): ((0, 0, 2, 3, 8, 20, 38), 40, None, 0.7),
             ("henon", "a1.40"): ((0, 1, 0, 1, 6, 14, 42), 54, 0.9, 0.7),
             ("logistic", "r3.58"): ((0,) * 7, 0, None, None),
             ("logistic", "r3.88"): ((0, 0, 0, 0, 7, 20, 54), 38, 0.9, 0.7),
             ("logistic", "r4.00"): ((1, 0, 2, 6, 11, 30, 65), 53, 0.9, 0.5)}

    def fmin(det, key, th, UU=None):
        UU = U if UU is None else UU
        pe = [cnt[(det, key)][f][0] / cnt[(det, key)][f][1] for f in FS]
        pi = [min(1.0, UU[det] / p) if p > 0 else 1.0 for p in pe]
        ok = [i for i in range(len(FS)) if all(v < th for v in pi[i:])]
        return FS[ok[0]] if ok else None, pi

    an = json.load(open(P10 / "conf/analysis/analysis.json"))
    nmatch = 0
    ntot = 0
    for det, rep in (("K3", K3rep), ("K1", K1rep)):
        for key, (counts, rp, f5, f20) in rep.items():
            got = tuple(cnt[(det, key)][f][0] for f in FS)
            n83 = {cnt[(det, key)][f][1] for f in FS}
            check(f"10-Q1-{det}.{key[0]}.{key[1]}", "10",
                  f"{det} {key[0]} {key[1]}: detected/83 at f = 0.05..0.9, replacement, f_min(pi<0.05), f_min(pi<0.20)",
                  (counts, rp, str(f5), str(f20)),
                  (got, cnt[(det, key)]["rep"][0], str(fmin(det, key, 0.05)[0]), str(fmin(det, key, 0.20)[0])),
                  rel(src), where, note=f"base windows per f: {sorted(n83)}")
    # phase-reset (K3: <= 5 at any f, replacement 0-10, no bound) and K1 coupled vdP / phase-reset rows
    pr = [k for k in {k for (d, k) in cnt if d == "K3"} if k[0] == "phase_reset"]
    mx = max(cnt[("K3", k)][f][0] for k in pr for f in FS)
    rp = [cnt[("K3", k)]["rep"][0] for k in pr]
    nob = all(fmin("K3", k, th)[0] is None for k in pr for th in (0.05, 0.20))
    check("10-Q1-K3.phase_reset", "10", "K3 phase-reset (4 regimes): max detected at any f, replacement range, no bound",
          (5, (0, 10), True), (mx, (min(rp), max(rp)), nob), rel(src), where)
    vdp = [k for k in {k for (d, k) in cnt if d == "K1"} if k[0] == "coupled_vdp"]
    mx1 = max(cnt[("K1", k)][f][0] for k in vdp for f in FS)
    rp1 = [cnt[("K1", k)]["rep"][0] for k in vdp]
    nob1 = all(fmin("K1", k, th)[0] is None for k in vdp for th in (0.05, 0.20))
    check("10-Q1-K1.vdp", "10", "K1 coupled vdP (4 regimes): max detected at any f, replacement range, no bound",
          (3, (0, 7), True), (mx1, (min(rp1), max(rp1)), nob1), rel(src), where)
    pr1 = {k[1]: (cnt[("K1", k)][0.9][0], cnt[("K1", k)]["rep"][0], str(fmin("K1", k, 0.20)[0])) for k in pr}
    check("10-Q1-K1.pr", "10", "K1 phase-reset tau1.14 / tau1.16: detected at f = 0.9, replacement, f_min(pi<0.20)",
          ((13, 0, "0.9"), (8, 18, "None")), (pr1["tau1.14"], pr1["tau1.16"]), rel(src), where)
    lg365 = (max(cnt[("K1", ("logistic", "r3.65"))][f][0] for f in FS), cnt[("K1", ("logistic", "r3.65"))]["rep"][0])
    check("10-Q1-K1.r3.65", "10", "K1 logistic r3.65: max detected at any f (reported '<= 2'), replacement", (2, 1),
          lg365, rel(src), where)
    # pi_upper spot-checks against analysis.json (all K3 / K1 cells)
    for det in ("K3", "K1"):
        ex = an["Q1"]["detectors"][det]["exclusion_L0"]
        for key in [k for (d, k) in cnt if d == det]:
            mine = fmin(det, key, 0.05)[1]
            stored = ex[f"{key[0]}|{key[1]}"]["design"]["pi_upper"]
            ntot += 1
            nmatch += all(abs(a - b) < 5e-4 for a, b in zip(mine, stored))
    check("10-Q1-pi", "10", "pi_upper(f) at all 7 design f, all 32 K3/K1 (family, parameter) rows: agree with "
          "analysis.json (|diff| < 5e-4)", "32/32", f"{nmatch}/{ntot}", rel(P10 / "conf/analysis/analysis.json"), where,
          note="recomputed from the spike-in counts with U = max(CP, reported bootstrap term)")
    same = sum(fmin(d_, k_, th)[0] == fmin(d_, k_, th, Ualt)[0] for d_ in ("K3", "K1")
               for k_ in [k for (dd, k) in cnt if dd == d_] for th in (0.05, 0.20))
    check("10-Q1-robust", "10", "f_min (both thresholds, all 32 rows) unchanged when U uses the largest bootstrap term "
          "over 50 other seeds", "64/64", f"{same}/64", rel(src), where,
          note=f"U(K3) {U['K3']:.5f} -> {Ualt['K3']:.5f}; U(K1) {U['K1']:.5f} -> {Ualt['K1']:.5f}; robustness of the "
          "seeded analysis to its Monte Carlo term")
    spots = [("K3", ("henon", "a1.22"), 0.5), ("K3", ("henon", "a1.40"), 0.3), ("K3", ("logistic", "r4.00"), 0.7),
             ("K3", ("coupled_vdp", "rho8_om3.3"), 0.3), ("K3", ("coupled_vdp", "rho2_om2.7"), 0.9),
             ("K3", ("henon", "a1.14"), 0.9), ("K1", ("logistic", "r4.00"), 0.5), ("K1", ("henon", "a1.40"), 0.9),
             ("K1", ("henon", "a1.14"), 0.7), ("K3", ("coupled_vdp", "rho6_om3.3"), 0.9),
             ("K1", ("logistic", "r3.88"), 0.9), ("K3", ("logistic", "r3.88"), 0.5)]
    for det, key, f in spots:
        p = cnt[(det, key)][f][0] / cnt[(det, key)][f][1]
        st = an["Q1"]["detectors"][det]["exclusion_L0"][f"{key[0]}|{key[1]}"]["design"]
        stored = st["pi_upper"][st["f"].index(f)]
        check(f"10-Q1-spot.{det}.{key[1]}.{f:g}", "10", f"pi_upper spot-check {det} {key[0]} {key[1]} at f = {f:g} "
              "(= U / p_emp)", f"{stored:.4g}", min(1.0, U[det] / p), rel(src),
              where + "; analysis.json exclusion_L0 (f = 0.9 column printed in tables.md)", rule="round")


def p10_q2(sg, label="", exclude_note=""):
    where = "docs/PHASE10_FINAL_RESULTS.md 6 (P2)" + exclude_note
    src = rel(P10 / "conf/seg.jsonl")
    d = sg[~sg.gap].copy()
    rng = np.random.default_rng(VSEED + 1)
    out = {}

    def rates(col, ok=None):
        dd = d if ok is None else d[d[ok]]
        r = {g: subj_rate(dd[dd.group == g], col) for g in ("NSR", "CHF")}
        return r

    def boot_ci(v):
        b = cluster_boot_means(v, 10000, rng)
        return pct(b, 2.5), pct(b, 97.5)

    def diff_ci(a, b):
        da = cluster_boot_means(a, 10000, rng) - cluster_boot_means(b, 10000, rng)
        return pct(da, 2.5), pct(da, 97.5)

    def perm(a, b, B=10000):
        allv = np.concatenate([a, b])
        d0 = abs(a.mean() - b.mean())
        c = 0
        for _ in range(B):
            q = rng.permutation(allv)
            c += abs(q[:len(a)].mean() - q[len(a):].mean()) >= d0 - 1e-15
        return (1 + c) / (B + 1)

    arms = (("raw", None), ("masked", "masked_ok"), ("edited", "edited_ok"), ("wu", "wu_ok"))
    for col, ok in arms:
        r = rates(col, ok)
        nsr, chf = r["NSR"].values, r["CHF"].values
        out[col] = (nsr.mean(), chf.mean(), chf.mean() - nsr.mean(), boot_ci(nsr), boot_ci(chf), diff_ci(chf, nsr),
                    perm(chf, nsr), len(nsr), len(chf))
    return out, d


def p10_q2_checks(out, d, sens=False):
    where = "docs/PHASE10_FINAL_RESULTS.md 6 (P2)"
    src = rel(P10 / "conf/seg.jsonl")
    rep = {"raw": ("51.7", "74.2", "22.4", ("46.0", "57.7"), ("66.5", "81.4"), ("12.8", "31.7"), "0.0001"),
           "masked": ("25.7", "20.8", "-4.9", ("20.9", "30.9"), ("13.8", "28.5"), ("-13.7", "4.3"), "0.29"),
           "edited": ("28.9", "34.3", "5.4", ("23.9", "34.2"), ("26.6", "42.4"), ("-3.9", "14.7"), "0.24"),
           "wu": ("28.2", "31.1", "2.9", ("23.3", "33.4"), ("24.5", "38.1"), ("-5.3", "11.4"), "0.50")}
    for arm, r in rep.items():
        o = out[arm]
        check(f"10-P2-{arm}", "10", f"P2 {arm}: DR NSR, DR CHF, CHF - NSR (%; mean of subject DRs)", r[:3], o[:3], src,
              where, rule="round", scale=100.0, note=f"subjects NSR {o[7]}, CHF {o[8]}")
        check(f"10-P2-{arm}-ci", "10", f"P2 {arm}: subject-cluster bootstrap 95 % CIs (NSR, CHF, difference; %)",
              r[3:6], o[3:6], src, where, rule="mc", tol=1.0, scale=100.0,
              note="independent bootstrap, 10,000 resamples, different seed")
        tolp = 0.0005 if arm == "raw" else 0.02
        check(f"10-P2-{arm}-p", "10", f"P2 {arm}: two-sided permutation p (10,000)", r[6], o[6], src, where, rule="mc",
              tol=tolp)


def p10_q2_change(d):
    where = "docs/PHASE10_FINAL_RESULTS.md 6 (P2 change)"
    src = rel(P10 / "conf/seg.jsonl")
    rng = np.random.default_rng(VSEED + 2)

    def change(arm, ok, g):
        dd = d[d[ok] & (d.group == g)] if ok else d[d.group == g]
        a = subj_rate(dd, "raw")
        b = subj_rate(dd, arm)
        delta = (b - a).values
        bs = cluster_boot_means(delta, 10000, rng)
        both = int((dd.raw & dd[arm]).sum())
        return (len(dd), a.mean(), b.mean(), delta.mean(), (pct(bs, 2.5), pct(bs, 97.5)), 1 - both / dd.raw.sum())

    rep = {("masked", "masked_ok", "NSR"): (5671, "51.8", "25.7", "-26.1", ("-32.4", "-20.1"), "49"),
           ("masked", "masked_ok", "CHF"): (1865, "67.7", "20.8", "-46.8", ("-59.1", "-34.5"), "66"),
           ("edited", "edited_ok", "NSR"): (5910, None, None, "-23.0", ("-28.9", "-17.4"), "44"),
           ("edited", "edited_ok", "CHF"): (2763, None, None, "-36.7", ("-46.7", "-26.9"), "55"),
           ("wu", "wu_ok", "NSR"): (5946, None, None, "-23.6", ("-29.6", "-17.9"), "46"),
           ("wu", "wu_ok", "CHF"): (3017, None, None, "-43.1", ("-52.7", "-33.4"), "60")}
    for (arm, ok, g), r in rep.items():
        c = change(arm, ok, g)
        if r[1] is not None:
            check(f"10-P2c-{arm}.{g}", "10", f"change raw -> {arm}, {g}: segments; raw DR -> arm DR; mean change (%)",
                  (r[0], r[1], r[2], r[3]), (c[0], c[1] * 100, c[2] * 100, c[3] * 100), src, where, rule="round")
        else:
            check(f"10-P2c-{arm}.{g}", "10", f"change raw -> {arm}, {g}: segments; mean change (%)", (r[0], r[3]),
                  (c[0], c[3] * 100), src, where, rule="round")
        check(f"10-P2c-{arm}.{g}-ci", "10", f"change raw -> {arm}, {g}: 95 % CI (%)", r[4], c[4], src, where, rule="mc",
              tol=1.0, scale=100.0)
        check(f"10-P2c-{arm}.{g}-rm", "10", f"change raw -> {arm}, {g}: raw positives removed (%)", r[5], c[5] * 100,
              src, where, rule="round")
    # non-analysable counted negative
    for g, r in (("NSR", (5964, "-27.5", "53")), ("CHF", (3029, "-59.6", "81"))):
        dd = d[d.group == g].copy()
        dd["m2"] = dd.masked & dd.masked_ok
        a, b = subj_rate(dd, "raw"), subj_rate(dd, "m2")
        both = int((dd.raw & dd.m2).sum())
        check(f"10-P2c-masked0.{g}", "10", f"masked, non-analysable = negative, {g}: segments, change (%), removed (%)",
              r, (len(dd), (b - a).mean() * 100, (1 - both / dd.raw.sum()) * 100), src, where, rule="round")
    an = {g: d[d.group == g].masked_ok.mean() * 100 for g in ("NSR", "CHF")}
    el = {g: d[d.group == g].edited_ok.mean() * 100 for g in ("NSR", "CHF")}
    check("10-P2-af", "10", "masked analysable % NSR / CHF; edited eligible % NSR / CHF", ("95", "62", "99", "91"),
          (an["NSR"], an["CHF"], el["NSR"], el["CHF"]), src, where, rule="round")
    nl = tuple(d[(d.group == g) & d[c]][v].mean() for c, v in (("raw", "raw_NL"), ("wu", "wu_NL")) for g in ("NSR", "CHF"))
    check("10-P2-NL", "10", "mean NL among positive segments: raw NSR, raw CHF, Wu NSR, Wu CHF",
          ("0.37", "0.59", "0.19", "0.33"), nl, src, where, rule="round")
    is_night = (d.clock >= 0) & (d.clock < 5)
    nd = tuple(d[(d.group == g) & (is_night == night)].raw.mean() * 100
               for g in ("NSR", "CHF") for night in (True, False))
    check("10-P2-nd", "10", "raw DR night / day (pooled segments, approximate clock): NSR night, day, CHF night, day (%)",
          ("55.8", "50.4", "72.0", "75.1"), nd, src, where, rule="round")


def p10_q3(df):
    where = "docs/PHASE10_FINAL_RESULTS.md 7 (P3)"
    src = rel(P10 / "conf/real.jsonl")
    X = np.column_stack([np.ones(len(df)), df.x, df.chf])
    res = {}
    for col in ("TIT", "LLE", "UPO", "K1"):
        corr = "independence" if col == "TIT" else "exchangeable"
        res[col] = gee_logit(X, df[col].astype(float), df.subject.values, corr)
    rep = {"TIT": (455, ("2.31", "1.67", "3.19"), ("0.77", "0.46", "1.28")),
           "LLE": (401, ("1.51", "1.31", "1.75"), ("0.79", "0.49", "1.28")),
           "UPO": (31, ("1.17", "1.00", "1.37"), ("2.13", "1.00", "4.55")),
           "K1": (16, ("1.44", "1.09", "1.90"), ("2.33", "0.57", "9.41"))}
    for col, (npos, orb, orc) in rep.items():
        r = res[col]
        check(f"10-P3-{col}", "10", f"GEE {col}: positives; OR per doubling [95 % CI]; CHF OR adjusted [95 % CI] "
              f"({'independence (Amendment 2)' if col == 'TIT' else 'exchangeable'} working correlation)",
              (npos, orb, orc), (int(df[col].sum()), (r["OR"][1], r["lo"][1], r["hi"][1]), (r["OR"][2], r["lo"][2], r["hi"][2])),
              src, where, rule="round", note=f"own GEE; exchangeable alpha = {r['alpha']:.3f}" if col != "TIT" else "own GEE")
    check("10-P3-p", "10", "GEE p-values: TIT burden, TIT CHF, LLE burden, UPO burden, K1 burden",
          ("5e-7", "0.31", "2e-8", "0.053", "0.010"),
          (res["TIT"]["p"][1], res["TIT"]["p"][2], res["LLE"]["p"][1], res["UPO"]["p"][1], res["K1"]["p"][1]),
          src, where, rule="round")
    # Amendment 2 statement: exchangeable fit for TIT
    with np.errstate(all="ignore"):
        ex = gee_logit(X, df["TIT"].astype(float), df.subject.values, "exchangeable")
    check("10-P3-A2", "10", "Amendment 2: the exchangeable GEE for the titration P3 model is not estimable",
          "diverged", "converged" if ex["converged"] else "diverged", src, "PREREGISTRATION.md Amendment 2",
          note="own exchangeable implementation also diverges (dependence parameter estimate becomes non-finite)")
    # unadjusted
    Xu = np.column_stack([np.ones(len(df)), df.x])
    un = {}
    for col in ("TIT", "LLE", "UPO", "K1"):
        un[col] = gee_logit(Xu, df[col].astype(float), df.subject.values, "exchangeable")
    check("10-P3-un", "10", "unadjusted burden OR: TIT [CI], LLE, UPO, K1",
          (("3.85", "2.30", "6.46"), "1.48", "1.28", "1.56"),
          ((un["TIT"]["OR"][1], un["TIT"]["lo"][1], un["TIT"]["hi"][1]), un["LLE"]["OR"][1], un["UPO"]["OR"][1],
           un["K1"]["OR"][1]), src, where, rule="round", note="exchangeable working correlation (converged; as recorded "
          "in analysis.json)")
    Xc = np.column_stack([np.ones(len(df)), df.chf])
    gc = gee_logit(Xc, df["TIT"].astype(float), df.subject.values, "exchangeable")
    check("10-P3-3c", "10", "Q3c titration CHF OR unadjusted [95 % CI]", ("2.60", "1.64", "4.13"),
          (gc["OR"][1], gc["lo"][1], gc["hi"][1]), src, where + " 3c", rule="round")
    # burden table
    bins = [(0, 0), (1, 1), (2, 4), (5, 15), (16, 10 ** 9)]
    rep = [(562, 72, 123, 145, 15, 5, 7), (116, 55, 73, 63, 2, 0, 0), (116, 39, 93, 71, 3, 2, 0),
           (83, 26, 71, 42, 0, 0, 1), (111, 20, 95, 80, 11, 9, 0)]
    for (lo, hi), r in zip(bins, rep):
        g = df[(df.burden >= lo) & (df.burden <= hi)]
        check(f"10-P3-bt.{lo}", "10", f"burden {lo}-{hi if hi < 1e9 else 'inf'}: windows, subjects, TIT, LLE, UPO, K1, K3", r,
              (len(g), g.subject.nunique(), int(g.TIT.sum()), int(g.LLE.sum()), int(g.UPO.sum()), int(g.K1.sum()),
               int(g.K3.sum())), src, where)
    check("10-P3-K3", "10", "K3: positives, in burden-0 windows; not estimable (< 10 positives)", (8, 7, True),
          (int(df.K3.sum()), int(df[df.burden == 0].K3.sum()), int(df.K3.sum()) < 10), src, where)
    rng = np.random.default_rng(VSEED + 3)
    a = subj_rate(df[df.burden >= 1], "K3").values
    b = subj_rate(df[df.burden == 0], "K3").values
    bd = cluster_boot_means(a, 10000, rng) - cluster_boot_means(b, 10000, rng)
    check("10-P3-K3d", "10", "K3 rate difference burden >= 1 minus 0 (difference of mean subject rates; %) [CI]",
          ("-0.7", ("-1.6", "0.03")), (100 * (a.mean() - b.mean()), (100 * pct(bd, 2.5), 100 * pct(bd, 97.5))), src,
          where, rule="mc", tol=0.15)
    t0 = df[df.burden == 0].TIT.mean() * 100
    t1 = df[df.burden == 1].TIT.mean() * 100
    check("10-P3-22-63", "10", "titration positive rate, burden 0 vs 1 (%)", ("22", "63"), (t0, t1), src, where,
          rule="round")
    # 3b paired
    rng = np.random.default_rng(VSEED + 4)
    pairs = {"TIT_masked": ("TIT", "TITm", "TITm_ok", (281, 134, 62, "57.5", "-31.5")),
             "TIT_edited": ("TIT", "TITe", "TITe_ok", (504, 318, 114, "66.0", "-36.8")),
             "LLE_edited": ("LLE", "LLEe", "K1e_ok", (504, 263, 168, "47.1", "-19.5")),
             "K1_edited": ("K1", "K1e", "K1e_ok", (504, 7, 4, "71.4", "-0.5")),
             "UPO_edited": ("UPO", "UPOe", "K1e_ok", (504, 14, 18, "64.3", "0.6"))}
    for name, (a_, b_, ok, rep) in pairs.items():
        g = df[(df.n_masked >= 1) & df[ok].astype(bool)]
        both = int((g[a_] & g[b_]).sum())
        delta = (subj_rate(g, b_) - subj_rate(g, a_)).values
        check(f"10-P3b-{name}", "10", f"3b {name} (windows with >= 1 masked interval): windows, raw +, arm +, removed %, "
              "mean subject change %", rep,
              (len(g), int(g[a_].sum()), int(g[b_].sum()), 100 * (1 - both / g[a_].sum()), 100 * delta.mean()), src,
              where + " 3b", rule="round")
    g = df[(df.n_masked >= 1) & df.TITe_ok]
    g2 = df[(df.n_masked >= 1) & df.K1e_ok.astype(bool)]
    check("10-P3b-new", "10", "edited arm creates positives: new TIT-positive, new LLE-positive windows", (6, 29),
          (int((~g.TIT & g.TITe).sum()), int((~g2.LLE & g2.LLEe).sum())), src, where + " 3b")


def p10_q4(df):
    where = "docs/PHASE10_FINAL_RESULTS.md 8 (Q4)"
    src = P10 / "conf/robust.jsonl"
    rb = list({r["id"]: r for r in jl(src)}.values())
    m = df.set_index(["subject", "j"])
    rows = []
    for r in rb:
        if r.get("error"):
            continue
        w = m.loc[(r["subject"], r["j"])]
        row = {"group": r["group"], "K1": w.K1, "LLE": w.LLE, "UPO": w.UPO, "K3": w.K3,
               "K1_499": bool(r["K1_499"]["lle_detected"] and r["K1_499"]["upo_detected"]),
               "LLE_499": bool(r["K1_499"]["lle_detected"]), "UPO_499": bool(r["K1_499"]["upo_detected"]),
               "K3_499": bool(r["K3_499"]["detected"])}
        for L in (256, 1024):
            ok = r[f"L{L}"].get("available")
            row[f"ok{L}"] = bool(ok)
            if ok:
                row[f"K1_{L}"] = bool(r[f"L{L}"]["K1"]["and_detected"])
                row[f"LLE_{L}"] = bool(r[f"L{L}"]["K1"]["lle_detected"])
                row[f"UPO_{L}"] = bool(r[f"L{L}"]["K1"]["upo_detected"])
                row[f"K3_{L}"] = bool(r[f"L{L}"]["K3"]["detected"])
        rows.append(row)
    R = pd.DataFrame(rows)
    rep = {"K1": (9, 8, 8, 1, 0, "99.7", "0.94"), "LLE": (143, 143, 140, 3, 3, "98.2", "0.96"),
           "UPO": (12, 10, 10, 2, 0, "99.4", "0.91"), "K3": (1, 1, 1, 0, 0, "100", "1.00")}
    for d, r in rep.items():
        a, b = R[d].astype(bool).values, R[f"{d}_499"].values
        n11, n10, n01 = int((a & b).sum()), int((a & ~b).sum()), int((~a & b).sum())
        n00 = len(a) - n11 - n10 - n01
        po = (n11 + n00) / len(a)
        pe = ((n11 + n10) * (n11 + n01) + (n00 + n01) * (n00 + n10)) / len(a) ** 2
        kap = (po - pe) / (1 - pe)
        check(f"10-Q4a-{d}", "10", f"499 surrogates {d}: default +, 499 +, both, default only, 499 only, agreement %, kappa",
              r, (int(a.sum()), int(b.sum()), n11, n10, n01, 100 * po, kap), rel(src), where + " 8.1", rule="round")
    check("10-Q4a-n", "10", "Q4a windows", 331, len(R), rel(src), where + " 8.1")
    ok = R.ok256 & R.ok1024
    Rk = R[ok]
    rep = {"K1": ("0.3", "2.8", "2.8"), "LLE": ("28.4", "42.2", "53.1"), "UPO": ("1.9", "3.8", "4.7"),
           "K3": ("0.3", "0.3", "2.8")}
    check("10-Q4b-n", "10", "Q4b windows available at 256, 512 and 1024", 320, int(ok.sum()), rel(src), where + " 8.2")
    for d, r in rep.items():
        check(f"10-Q4b-{d}", "10", f"{d} rate at 256 / 512 / 1024 intervals (%)", r,
              (Rk[f"{d}_256"].mean() * 100, Rk[d].astype(bool).mean() * 100, Rk[f"{d}_1024"].mean() * 100), rel(src),
              where + " 8.2", rule="round")
    check("10-Q4b-K3n", "10", "K3 detections at 1024 / 512 / 256 (of 320)", (9, 1, 1),
          (int(Rk.K3_1024.sum()), int(Rk.K3.astype(bool).sum()), int(Rk.K3_256.sum())), rel(src), where + " 8.2")
    nd = {}
    for d in ("TIT", "LLE", "UPO", "K1", "K3", "K4"):
        nd[d] = tuple(f"{int(df[(df.group == g) & (df.night == n)][d].sum())}/{int(((df.group == g) & (df.night == n)).sum())}"
                      for g in ("NSR", "CHF") for n in (True, False))
    rep = {"TIT": ("66/152", "177/490", "47/82", "165/264"), "LLE": ("40/152", "190/490", "33/82", "138/264"),
           "UPO": ("1/152", "11/490", "3/82", "16/264"), "K1": ("1/152", "3/490", "2/82", "10/264"),
           "K3": ("1/152", "4/490", "0/82", "3/264"), "K4": ("0/152", "0/490", "0/82", "0/264")}
    for d, r in rep.items():
        check(f"10-Q4c-{d}", "10", f"{d} night/day counts: NSR night, NSR day, CHF night, CHF day", r, nd[d],
              rel(P10 / "conf/real.jsonl"), where + " 8.3")
    Xn = np.column_stack([np.ones(len(df)), df.x, df.chf, df.night.astype(float)])
    gn = gee_logit(Xn, df.TIT.astype(float), df.subject.values, "exchangeable")
    check("10-Q4c-gee", "10", "titration GEE night OR [95 % CI], p (adjusted for burden and group)",
          ("1.40", "0.97", "2.00", "0.069"), (gn["OR"][3], gn["lo"][3], gn["hi"][3], gn["p"][3]),
          rel(P10 / "conf/real.jsonl"), where + " 8.3", rule="round",
          note="exchangeable working correlation (converged; as recorded in analysis.json)")
    r3 = (int(df.K3.sum()), 100 * df[df.group == "NSR"].K3.mean(), 100 * df[df.group == "CHF"].K3.mean(),
          int(df.K4.sum()), int(df.K3RR.sum()), int((df.K3 != df.K4).sum()))
    check("10-Q4d", "10", "K3 detections, NSR %, CHF %; K4; K3RR; K3-K4 disagreements", (8, "0.78", "0.87", 0, 0, 8), r3,
          rel(P10 / "conf/real.jsonl"), where + " 8.4", rule="round")


def p10_k3_table(df):
    where = "docs/PHASE10_FINAL_RESULTS.md 8.5"
    src = rel(P10 / "conf/real.jsonl")
    R = {r["id"]: r for r in jl(P10 / "conf/real.jsonl")}
    rep = [("nsr2db:nsr029", 9, 0, 1, None, "0.40", "-+--"), ("nsr2db:nsr049", 3, 0, 1, "14.9", "0.53", "-+-+"),
           ("nsr2db:nsr049", 5, 0, 1, "11.1", "0.40", "-+--"), ("nsr2db:nsr049", 6, 0, 3, "13.8", "0.30", "-+--"),
           ("nsr2db:nsr052", 6, 0, 1, "11.1", "0.33", "-+--"), ("chf2db:chf211", 7, 0, 1, "10.0", "0.30", "++++"),
           ("chf2db:chf216", 2, 0, 0, "10.6", "0.40", "-+--"), ("chf2db:chf218", 12, 5, 11, "16.2", "0.27", "-+--")]
    det = df[df.K3]
    found = {f"{s}|{j}" for s, j in zip(det.subject, det.j)}
    order = [f"{s}|{j}" for s, j, *_ in rep]
    check("10-K3t-0", "10", "the 8 K3 detections (subject, window)", tuple(order),
          tuple(x for x in order if x in found) + tuple(sorted(found - set(order))), src, where)
    for s, j, bur, nm, z, G, flags in rep:
        db, sub = s.split(":")
        r = R[f"real|{db}|{sub}|{j}"]
        k1 = r["K1"]
        fl = "".join("+" if v else "-" for v in (k1["lle_detected"] and k1["upo_detected"], k1["lle_detected"],
                                                  k1["upo_detected"], r["TIT_raw"]["NL"] > 0))
        reported = (bur, nm, G, flags) if z is None else (bur, nm, z, G, flags)
        recomputed = (r["burden"], r["K3"]["n_masked"], r["K3"]["G"], fl) if z is None else \
            (r["burden"], r["K3"]["n_masked"], r["K3"]["zmax"], r["K3"]["G"], fl)
        ok = all((str(a) == str(b)) if isinstance(a, str) and not a.replace(".", "").replace("-", "").isdigit()
                 else abs(float(a) - float(b)) <= 0.5 * 10 ** (-_decimals(a)) + 1e-9
                 for a, b in zip(reported, recomputed))
        CHECKS.append({"id": f"10-K3t-{sub}.{j}", "phase": "10", "quantity": f"{s} window {j}: burden, masked, "
                       + ("" if z is None else "z_max, ") + "G, K1/LLE/UPO/TIT", "reported": _fmt(reported),
                       "recomputed": _fmt(recomputed), "source": src, "where": where, "rule": "round",
                       "status": "MATCH" if ok else "MISMATCH", "note": ""})
    r = R["real|nsr2db|nsr029|9"]
    check("10-K3t-nsr029", "10", "nsr029 w9: z_max (reported 2 x 10^11) and SDNN (ms)", ("2e11", "15"),
          (r["K3"]["zmax"], r["sdnn"] * 1000), src, where, rule="round")
    zz = df[df.K3_an & (df.K3_z >= 9.4521)]
    check("10-K3t-close", "10", "windows passing z >= 9.45; removed by the growth gate", (34, 26),
          (len(zz), int((zz.K3_G < 0.25).sum())), src, where)


def phase10():
    where = "docs/PHASE10_FINAL_RESULTS.md"
    df = p10_real()
    check("10-0", "10", "stored flags equal recomputed: K1 = LLE AND UPO; K3 = analysable & z >= 9.4521 & G >= 0.25; "
          "TIT = NL > 0 (disagreements)", (0, 0, 0),
          (int((df.K1 != df.K1_stored).sum()), int((df.K3 != df.K3_stored).sum()),
           int((df.TIT != df.TIT_stored).sum())), rel(P10 / "conf/real.jsonl"), where)
    per = df.groupby("subject").size()
    nfull = df[df.group == "NSR"].groupby("subject").size()
    cfull = df[df.group == "CHF"].groupby("subject").size()
    check("10-1", "10", "windows; subjects; NSR / CHF subjects with all 12 windows; subjects with 10-11",
          (988, 83, 49, 27, 6), (len(df), df.subject.nunique(), int((nfull == 12).sum()), int((cfull == 12).sum()),
                                 int(((per >= 10) & (per <= 11)).sum())), rel(P10 / "conf/real.jsonl"), where + " 3")
    p10_q1(df)
    sg = p10_seg()
    check("10-P2-n", "10", "12-min segments; analysable; NSR / CHF subjects", (9536, 8993, 54, 29),
          (len(sg), int((~sg.gap).sum()), sg[sg.group == "NSR"].subject.nunique(),
           sg[sg.group == "CHF"].subject.nunique()), rel(P10 / "conf/seg.jsonl"), where + " 6")
    out, d = p10_q2(sg)
    p10_q2_checks(out, d)
    p10_q2_change(d)
    gs = gee_logit(np.column_stack([np.ones(len(d)), np.log2(1 + d.burden), (d.group == "CHF").astype(float)]),
                   d.raw.astype(float), d.subject.values, "independence")
    check("10-P2-gee", "10", "segment GEE: burden OR [CI], CHF OR [CI] (independence fallback)",
          (("1.99", "1.51", "2.62"), ("0.87", "0.55", "1.36")),
          ((gs["OR"][1], gs["lo"][1], gs["hi"][1]), (gs["OR"][2], gs["lo"][2], gs["hi"][2])),
          rel(P10 / "conf/seg.jsonl"), where + " 6", rule="round")
    p10_q3(df)
    p10_q4(df)
    p10_k3_table(df)
    p10_q2e()
    p10_verification()
    p10_sensitivity()
    p10_dev()


def p10_q2e():
    where = "docs/PHASE10_FINAL_RESULTS.md 6.2 (Q2e)"
    src = P10 / "conf/synth.jsonl"
    R = list({r["id"]: r for r in jl(src)}.values())
    c = defaultdict(lambda: [0, 0, 0, 0, 0, 0, 0])
    for r in R:
        key = r["cond"] if r["kind"] == "null" else f"vdp{r['cond'][0]}_{r['cond'][1]}"
        v = c[(key, r["n"])]
        v[0] += 1
        v[1] += r["TIT_raw"]["positive"]
        v[2] += r.get("TIT_raw_q128", {}).get("positive", False)
        v[3] += r["TIT_masked"].get("positive", False)
        v[4] += r["TIT_masked"].get("analysable", False)
        v[5] += r["TIT_edited"].get("positive", False)
        v[6] += r["TIT_edited"].get("eligible", False)
    check("10-Q2e-n", "10", "synthetic windows; errors", (8200, 0), (len(R), sum(1 for r in R if r.get("error"))), rel(src),
          where)
    rep = {("N5_linear_rr_warped", 512): (100, 99, 100, 100), ("N5_linear_rr_warped", 800): (100, 98, None, None),
           ("E3_couplets_5pct", 512): (100, 100, 39, None), ("E3_couplets_10pct", 512): (100, 100, 23, None),
           ("E4_runs", 512): (96, 96, 11, 0), ("E4_runs", 800): (99, 99, 6, 0), ("N4_linear_rr_step", 512): (10, 0, 10, 10),
           ("S1_setar", 512): (0, 0, 0, 0), ("S1_setar", 800): (3, 3, 3, None), ("S2_ectopic_5pct", 512): (5, 5, 12, None),
           ("vdp26_none", 512): (100, 96, 100, 100), ("vdp26_none", 800): (100, 100, None, None),
           ("vdp26_E3_10", 512): (94, 93, 100, None), ("vdp27_E3_10", 512): (9, 8, 29, None),
           ("vdp27_E3_10", 800): (0, 0, 11, None), ("vdp29_E3_10", 512): (100, 100, 44, None),
           ("vdp29_E3_10", 800): (100, 100, 82, None), ("vdp27_S2_5", 512): (0, 0, 20, None),
           ("vdp30_S2_5", 512): (0, 0, 0, None)}
    for (cond, n), r in rep.items():
        v = c[(cond, n)]
        got = (v[1], v[2], v[5], v[3])
        check(f"10-Q2e-{cond}.{n}", "10", f"{cond} n = {n}: raw, raw 1/128 s, edited, masked positives (/100)",
              tuple(x for x in r if x is not None), tuple(g for g, x in zip(got, r) if x is not None), rel(src), where)
    lin = [k for k in c if k[0] in ("N1_linear_rr", "N2_power_law", "N3_linear_rr_trend", "N6_noisy_rsa")]
    check("10-Q2e-lin", "10", "N1, N2, N3, N6 (512 and 800): any positive in raw / q128 / edited / masked", 0,
          sum(c[k][1] + c[k][2] + c[k][3] + c[k][5] for k in lin), rel(src), where)
    other = [k for k in c if k[0] in ("vdp27_none", "vdp28_none", "vdp29_none", "vdp30_none", "vdp31_none")]
    check("10-Q2e-vdp", "10", "other 5 non-chaotic vdP regimes without ectopy: raw / q128 / edited positives (512 and 800)",
          0, sum(c[k][1] + c[k][2] + c[k][5] for k in other), rel(src), where,
          note="regime index -> (rho, omega) from phase8_cardiac/results/ground_truth/regimes.json: 26 = (5.45, 5.6), "
          "27 = (9.6, 2.1), 28 = (2, 5.6), 29 = (4, 5.6), 30 = (10, 3.3), 31 = (6, 5.6)")


def p10_verification():
    src = P10 / "verification/titration_verification.json"
    v = json.load(open(src))["verification"]
    scan = np.array(v["V4"]["scan"], float)
    ch = scan[scan[:, 1] > 0.02]          # predeclared margins (METHODS.md 2.5): chaotic LE > 0.02, periodic < -0.02
    pe = scan[scan[:, 1] < -0.02]
    rs = stats.spearmanr(ch[:, 2], ch[:, 1])[0]
    rho = np.corrcoef(stats.rankdata(ch[:, 2]), stats.rankdata(ch[:, 1]))[0, 1]
    a0, b0 = scan[scan[:, 1] > 0], scan[scan[:, 1] < 0]
    check("10-V4", "10", "titration V4: NL > 0 at chaotic r (LE > 0.02); NL = 0 at periodic r (LE < -0.02)",
          ("78/78", "18/18"), (f"{int((ch[:, 2] > 0).sum())}/{len(ch)}", f"{int((pe[:, 2] == 0).sum())}/{len(pe)}"),
          rel(src), "PHASE10 4.3; S12 #12; METHODS 2.5",
          note=f"without the ±0.02 margin: {int((a0[:, 2] > 0).sum())}/{len(a0)} and {int((b0[:, 2] == 0).sum())}/{len(b0)}")
    check("10-V4b", "10", "titration V4: Spearman(NL, LE) over chaotic r (own rank correlation)", "0.85", rho, rel(src),
          "PHASE10 4.3", rule="round", note=f"scipy spearmanr {rs:.4f}")
    check("10-V1-3", "10", "V1 NL(3.7), V2 NL(3.575), V3 NL(3.565) (%)", ("63.7", "8.7", "0.0"),
          (v["V1"]["NL_3.7"] * 100, v["V2"]["NL_3.575"] * 100, v["V3"]["NL_3.565"] * 100), rel(src), "PHASE10 4.3",
          rule="round")


def p10_sensitivity():
    ex = json.load(open(P10 / "conf/exclusions.json"))["originally_flagged_subjects"]
    where = "docs/PHASE10_FINAL_RESULTS.md 9.2"
    df = p10_real(exclude=set(ex))
    sg = p10_seg(exclude=set(ex))
    d = sg[~sg.gap]
    r = {g: subj_rate(d[d.group == g], "raw") for g in ("NSR", "CHF")}
    dm = d[d.masked_ok]
    m = {g: subj_rate(dm[dm.group == g], "masked") for g in ("NSR", "CHF")}
    X = np.column_stack([np.ones(len(df)), df.x, df.chf])
    g = gee_logit(X, df.TIT.astype(float), df.subject.values, "independence")
    check("10-S-1", "10", "without the 10 flagged subjects: subjects; P2 raw diff (%); masked diff (%); "
          "P3 TIT OR; CHF OR", (73, "23.3", "0.0", "2.21", "0.82"),
          (df.subject.nunique(), 100 * (r["CHF"].mean() - r["NSR"].mean()), 100 * (m["CHF"].mean() - m["NSR"].mean()),
           g["OR"][1], g["OR"][2]), rel(P10 / "conf/*.jsonl"), where, rule="round")


def p10_dev():
    where = "docs/PHASE10_FINAL_RESULTS.md 9.1 (development)"
    df = p10_real("dev")
    df = df[df.group.isin(["NSR", "CHF"])]
    sg = p10_seg("dev")
    d = sg[~sg.gap]
    r = {g: subj_rate(d[d.group == g], "raw") for g in ("NSR", "CHF")}
    dm = d[d.masked_ok]
    chg = (subj_rate(dm[dm.group == "CHF"], "masked") - subj_rate(dm[dm.group == "CHF"], "raw")).mean()
    m = {g: subj_rate(dm[dm.group == g], "masked") for g in ("NSR", "CHF")}
    X = np.column_stack([np.ones(len(df)), df.x, df.chf])
    g = gee_logit(X, df.TIT.astype(float), df.subject.values, "exchangeable")
    gl = gee_logit(X, df.LLE.astype(float), df.subject.values, "exchangeable")
    check("10-D-1", "10", "dev: raw DR NSR, CHF, diff (%); masked change CHF (%); masked diff (%)",
          ("47.5", "80.5", "33.0", "-58.3", "-23.6"),
          (100 * r["NSR"].mean(), 100 * r["CHF"].mean(), 100 * (r["CHF"].mean() - r["NSR"].mean()), 100 * chg,
           100 * (m["CHF"].mean() - m["NSR"].mean())), rel(P10 / "dev/seg.jsonl"), where, rule="round")
    check("10-D-2", "10", "dev: TIT OR per doubling, CHF OR adjusted, LLE OR per doubling",
          ("2.43", "0.78", "1.33"), (g["OR"][1], g["OR"][2], gl["OR"][1]), rel(P10 / "dev/real.jsonl"), where,
          rule="round")
    check("10-D-3", "10", "dev: K3 positives / K1 positives / windows (nsrdb + chfdb)", "0/3/396",
          f"{int(df.K3.sum())}/{int(df.K1.sum())}/{len(df)}", rel(P10 / "dev/real.jsonl"), where)


# =================================================================================== cross-report consistency
def cross():
    """Later reports quoting earlier phases' numbers: each quote vs the source table (already recomputed above
    where marked 'see')."""
    rows = [
        # (id, quoting report + place, quote, source report value, status-bool, note)
        ("X-1", "PHASE10 S12 #3", "White noise 4/100, AR(1) 3/100, 20 dB maps 200/200 (Phase 3 test seeds)",
         "PHASE3 5.1: 4/100, 3/100, 200/200 (see 3-1, 3-2)", True, ""),
        ("X-2", "PHASE10 S12 #4", "Sinusoid 0/150 and two-tone 1/150 (baseline 89 and 67/150); noisy Henon 293/300",
         "PHASE4 4.1: 0/150, 1/150; 89/150, 67/150; 293/300 (see 4-1)", True, ""),
        ("X-3", "PHASE10 S12 #2", "Level-B false positives 14/200 and 15/200 at nominal 0.059; Level A 70-98 %",
         "PHASE2E 3.1: 14/200, 15/200, nominal 0.0588; Level A 70-98 % (white noise) and 78-98 % (AR(1))", True,
         "the 70-98 % range is the white-noise range; AR(1) is 78-98 %"),
        ("X-4", "PHASE10 S12 #1", "positive in 194/195 white-noise and 197/197 AR(1) windows",
         "PHASE2E 3.2 Table H: 194/195, 197/197 (see 2E-1, 2E-2)", True, ""),
        ("X-5", "PHASE10 S12 #7", "Pooled power on trended maps rose from 32.5 % to 95.2 %",
         "PHASE6 4.4: 195/600 -> 571/600 (see 6-6, 6-8)", True, ""),
        ("X-6", "PHASE10 S12 #6", "At most 3.7 % false positives on synthetic ectopy patterns (Phase 6)",
         "PHASE6 3.3: at or below 3.7 % on every RAW ectopy pattern (11/300); edited trigeminy 13/300 = 4.3 %", True,
         "statement refers to K1 on raw RR, as in Phase 6 3.3"),
        ("X-7", "PHASE10 S12 #6", "Up to 45/100 on non-chaotic coupled vdP with ectopy (Phase 9 TEST)",
         "PHASE9 6: k1 failed 4 conditions, 8-45/100 (see 9-7)", True, ""),
        ("X-8", "PHASE10 S12 #9", "MIT-BIH, K1: 0/305; nsrdb / chfdb, K3: 1/180 and 0/150; nsr2db / chf2db, K3: 8/988",
         "PHASE7 4: 0/305 (see 7-1); PHASE9 7: 1/180, 0/150 (see 9-8); PHASE10 (see 10-Q1-U.K3)", True, ""),
        ("X-9", "PHASE10 S12 #10", "0/10,100 synthetic null and non-chaotic windows; 7 % pooled power on TEST chaos, "
         "0 % with bigeminy", "PHASE9 6: 0/10,100; 90/1,320 = 6.8 %; vdP E1 0 (see 9-1, 9-5)", True, ""),
        ("X-10", "PHASE10 S12 #12", "NL > 0 at 78/78 chaotic and NL = 0 at 18/18 periodic parameters",
         "PHASE10 4.3 V4 (see 10-V4)", True, ""),
        ("X-11", "PHASE10 S12 #8", "Rossler and Mackey-Glass were 0 % for K1 at every m (Phases 5, 8)",
         "PHASE5 G2/G3 AND 0/200 at m = 2, 3, 4 (see 5-6); PHASE8 C1 flows 0/200 (see 8-5)", True, ""),
        ("X-12", "PHASE10 analysis10.SYNTH_FP", "K3 0/10,100; K1 97/10,100 (Phase 9 TEST pooled PASS conditions)",
         "PHASE9 6: k3 0, k1 97 of 10,100 (see 9-1)", True, ""),
        ("X-13", "PHASE9 1", "C1 (frozen baseline) detected 74/1,320 ... failed specificity ... (up to 45/100)",
         "PHASE9 6 (see 9-1, 9-7)", True, ""),
        ("X-14", "PHASE7 4 and 5", "Phase 6 ectopy-only false-positive rates up to 3.7 %; synthetic nulls at most "
         "2/400", "PHASE6 3.1 (11/300 = 3.7 %); PHASE6 4.4 BASELINE-K N1-N6 max 2/400 (N5) (see 6-4)", True, ""),
        ("X-15", "PHASE7 12.3", "P1 Henon 99-100 %, P2 logistic 97-98 % (Phases 5-6)",
         "PHASE5: P1 198/200 = 99 %, P2 194/200 = 97 %; PHASE6: P1 300/300, P2 294/300 = 98 % (see 5-4, 6-5)", True, ""),
        ("X-16", "PHASE6 2", "Phase 5 N1-N6 AND counts identical: 1, 1, 0, 0, 1, 0", "PHASE5 5.1 (see 5-1)", True, ""),
        ("X-17", "PHASE10 S1 / S7", "22 % positive with no ectopic beat vs 63 % with one", "see 10-P3-22-63", True, ""),
        ("X-18", "PHASE9 6", "k1 on non-chaotic vdP with ectopy confirms the Phase 8 secondary warning (11/60)",
         "PHASE8 7.3 secondary table", None, "checked below from the Phase 8 TEST records"),
    ]
    for cid, where, quote, src, ok, note in rows:
        if ok is None:
            continue
        CHECKS.append({"id": cid, "phase": "X", "quantity": quote, "reported": "(quote)", "recomputed": src,
                       "source": "see linked checks", "where": where, "rule": "text", "status": "MATCH" if ok else
                       "MISMATCH", "note": note})
    # X-18: Phase 8 secondary: C1 on non-chaotic regimes with ectopy (noncha_ect), 512
    rows8 = [json.loads(l) for l in open(E / "phase8_cardiac/results/test/test_512.jsonl")]
    ne = [r for r in rows8 if r["group"] == "noncha_ect"]
    k = sum(bool(r["results"]["c1_frozen512"]["detected"]) for r in ne)
    vdp = [r for r in ne if "coupled_vdp" in r["regime"] and r["variant"] == "vi_E1"]
    kv = sum(bool(r["results"]["c1_frozen512"]["detected"]) for r in vdp)
    check("X-18", "X", "PHASE9 6 quotes the Phase 8 secondary warning '11/60' (C1 on non-chaotic vdP with bigeminy)",
          "11/60", f"{kv}/{len(vdp)}", rel(E / "phase8_cardiac/results/test/test_512.jsonl"), "PHASE9 6",
          note=f"all non-chaotic regimes with ectopy: {k}/{len(ne)}")
    # Phase 10 text statements checked against the recomputed numbers (see ERRATA.md for any MISMATCH)
    text_checks()


def data_checks():
    """Phase 10 Section 3 (data and overlap checks) against the committed check outputs, and two bounds quoted in
    Phases 7 and 9."""
    where = "docs/PHASE10_FINAL_RESULTS.md 3"
    a = json.load(open(P10 / "conf/data_checks.json"))
    b = json.load(open(P10 / "conf/data_checks_amended.json"))
    md = a["metadata"]["conf"]
    nsr = sum(1 for k in md if k.startswith("nsr2db"))
    chf = sum(1 for k in md if k.startswith("chf2db"))
    fs = sorted(set(a["resolution"].values()))
    starts = sum(1 for v in md.values() if v.get("base_time"))
    check("10-D3-1", "10", "records nsr2db / chf2db; annotation fs values; headers with a start time", (54, 29, "128.0", 0),
          (nsr, chf, ",".join(str(x) for x in fs), starts), rel(P10 / "conf/data_checks.json"), where)
    ov = a["rr_matching_summary"]["overlapping_pairs"]
    check("10-D3-2", "10", "original rule: flagged pairs, chf2db subjects involved; pairs compared", (13, 10, 2739),
          (len(ov), len({o["conf"] for o in ov}), a["rr_matching_summary"]["pairs"]), rel(P10 / "conf/data_checks.json"),
          where)
    mx = max(pp["n_consistent"] for pp in b["pairs"])
    val = [v for v in b["validation"] if v.get("same")]
    check("10-D3-3", "10", "amended rule: largest number of time-consistent probes; overlaps; validation duplicates "
          "detected 5/5", (1, 0, True), (mx, sum(bool(pp["overlap"]) for pp in b["pairs"]),
                                         all(v["n_consistent"] == 5 for v in val) and len(val) == 4),
          rel(P10 / "conf/data_checks_amended.json"), where, note=f"{len(val)} re-annotated duplicates in the validation")
    check("X-19", "X", "PHASE7 5: the Q2 rule needs a lower bound above 6.4 % (Phase 6 E5 atrial 10 %, 11/300, Wilson "
          "upper)", "6.4", wilson(11, 300)[1], "experiments/phase6_robust/results/test/test.jsonl", "PHASE7 5",
          rule="round", scale=100.0)
    check("X-20", "X", "PHASE9 7: synthetic false-positive upper bound for k3 (one-sided 95 % CP of 0/10,100)", "0.0003",
          cp_upper_bisect(0, 10100), "experiments/phase9_waveform/results/test/test9.jsonl", "PHASE9 7", rule="round")
    dev = P10 / "dev"
    check("10-D-4", "10", "development partial runs: spike-in tasks, robustness tasks", (96, 36),
          (sum(1 for _ in open(dev / "spike.jsonl")), sum(1 for _ in open(dev / "robust.jsonl"))), rel(dev),
          "PHASE10 9.1")
    df = p10_real("dev")
    mit = df[df.group == "MIT"]
    g = gee_logit(np.column_stack([np.ones(len(mit)), mit.x]), mit.TIT.astype(float), mit.subject.values, "exchangeable")
    check("10-D-5", "10", "development MIT-BIH exploratory titration OR per doubling [95 % CI]", ("1.23", "1.03", "1.46"),
          (g["OR"][1], g["lo"][1], g["hi"][1]), rel(dev / "real.jsonl"), "PHASE10 9.1", rule="round")


def text_checks():
    """Statements in docs/PHASE10_FINAL_RESULTS.md that summarise tables, checked against the tables'
    recomputed values (values from the checks above)."""
    src = "docs/PHASE10_FINAL_RESULTS.md"
    # S12 #13: 'Masking or editing ectopy removed 49-66 % of segment positives'
    rm = {c["id"]: c for c in CHECKS}
    vals = [float(rm[f"10-P2c-{a}.{g}-rm"]["recomputed"]) for a in ("masked", "edited") for g in ("NSR", "CHF")]
    check("T-1", "T", "S12 #13 'Masking or editing ectopy removed 49-66 % of segment positives': range over masked "
          "and edited arms (%)", ("49", "66"), (min(vals), max(vals)), "experiments/phase10_final/results/conf/seg.jsonl",
          src + " S12 #13", rule="round", note="edited arm removed 44.5 % (NSR) and 55.2 % (CHF)")
    # S12 #5 'fired in 94-99 % of linear-Gaussian windows with 2-10 % isolated ectopic beats (Phases 5-6)'
    p5 = [194 / 200, 195 / 200, 197 / 200]
    p6 = [289 / 300, 290 / 300, 293 / 300]
    check("T-2", "T", "S12 #5 'fired in 94-99 % of linear-Gaussian windows with 2-10 % isolated ectopic beats "
          "(Phases 5-6)': range of the LLE-alone rates (%)", ("94", "99"), (100 * min(p5 + p6), 100 * max(p5 + p6)),
          "experiments/phase5_rr, phase6_robust test.jsonl", src + " S12 #5", rule="round",
          note="94 % is the lowest Phase 6 Wilson lower bound, not a rate; rates are 96.3-98.5 %")
    # S1 'K3 never saw weak maps (lambda <= 0.26 per beat) ... so no bound is possible for them at any f'
    check("T-3", "T", "S1 'K3 never saw weak maps (λ <= 0.26 per beat) ... no bound at any f': maps with λ <= 0.26 "
          "and a K3 bound (f_min at π < 0.20)", "none", "Hénon a = 1.14 (λ 0.244): f_min(π < 0.20) = 0.7",
          "experiments/phase10_final/results/conf/spike.jsonl", src + " S1 (Q1)", rule="exact",
          note="Section 5 'Findings' states it correctly (weak = logistic λ <= 0.26 and Hénon λ = 0.14)")
    # S1 / S5: 'K1 reached only pi < 0.20 from f = 0.5-0.7, and only for the strong maps'
    check("T-4", "T", "S1 'K1 reached only π < 0.20 from f = 0.5-0.7, and only for the strong maps' (S5: 'π < 0.20 "
          "only from f = 0.5-0.7 for the strong maps')", "only π < 0.20; only strong maps",
          "π < 0.05 at f = 0.9 for Hénon a = 1.40, logistic r = 3.88, 4.00; π < 0.20 also for Hénon a = 1.14 "
          "(λ 0.244, f = 0.7) and phase-reset τ = 1.14 (f = 0.9)", "experiments/phase10_final/results/conf/spike.jsonl",
          src + " S1, S5 Findings", rule="exact",
          note="the K1 table in S5.2 itself is correct (see 10-Q1-K1.*)")
    # S1 / S8.1: '499 surrogates changed <= 3 decisions per detector'
    disc = {k: int(v["recomputed"].strip("()").split(", ")[3]) + int(v["recomputed"].strip("()").split(", ")[4])
            for k, v in ((d, rm[f"10-Q4a-{d}"]) for d in ("K1", "LLE", "UPO", "K3"))}
    check("T-5", "T", "S1 / S8.1 '499 surrogates changed <= 3 decisions per detector': largest number of windows whose "
          "decision changed (default-only + 499-only)", "3", max(disc.values()), "experiments/phase10_final/results/conf/"
          "robust.jsonl", src + " S1 (Q4), S8.1", rule="exact", note=f"per detector: {disc}; LLE 3 + 3")
    check("T-6", "T", "S5 Findings 'for the stronger Henon and logistic maps (λ >= 0.30) ... (f >= 0.5) ... at most 5 %': "
          "largest K3 f_min(π < 0.05) over those maps", "0.5",
          max(float(rm[f"10-Q1-K3.{f}.{p}"]["recomputed"].strip("()").split(", ")[-2]) for f, p in
              (("henon", "a1.22"), ("henon", "a1.40"), ("logistic", "r3.88"), ("logistic", "r4.00"))),
          "experiments/phase10_final/results/conf/spike.jsonl", src + " S5 Findings", rule="round",
          note="logistic r = 4.00 (λ 0.693): f_min(π < 0.05) = 0.7 (S1 states the range 0.5-0.7 correctly)")
    lle_rm = float(rm["10-P3b-LLE_edited"]["recomputed"].strip("()").split(", ")[3])
    check("T-7", "T", "S7 3b 'Masking or editing removes most titration and LLE positives': LLE raw positives removed "
          "by editing (%) > 50", True, lle_rm > 50, "experiments/phase10_final/results/conf/real.jsonl", src + " S7 (3b)",
          rule="exact", note=f"LLE: {lle_rm:.1f} % removed (titration 57.5 % masked, 66.0 % edited)")
    R = [json.loads(l) for l in open(P10 / "conf/robust.jsonl")]
    d1024 = {r["subject"] for r in R if r["L1024"].get("available") and r["L1024"]["K3"]["detected"]}
    d512 = {r["subject"] for r in (json.loads(l) for l in open(P10 / "conf/real.jsonl")) if r["K3"]["detected"]}
    check("T-8", "T", "S8.2 'Three subjects recur from the 512-interval detections': subjects with a 1,024-interval K3 "
          "detection that also had a 512-interval detection", 3, len(d1024 & d512),
          "experiments/phase10_final/results/conf/robust.jsonl", src + " S8.2", rule="exact",
          note=f"{sorted(d1024 & d512)} (the report lists these five in one sentence)")
    dev = json.load(open(P10 / "dev/analysis/analysis.json"))["Q2"]
    fr = [dev[k]["NSR"]["a_pos_removed_frac"] * 100 for k in ("P2_change_masked", "change_edited", "change_wu")]
    check("T-9", "T", "S6 'In development, nsrdb ... lost only 2-6 % of positives': range of NSR raw positives removed "
          "(masked, edited, Wu; %)", ("2", "6"), (min(fr), max(fr)), "experiments/phase10_final/results/dev/seg.jsonl",
          src + " S6", rule="round", note="from results/dev/analysis/analysis.json (development data)")


# =================================================================================== output
def write():
    df = pd.DataFrame(CHECKS)
    n_match = int((df.status == "MATCH").sum())
    n_mis = int((df.status == "MISMATCH").sum())
    n_info = int((df.status == "INFO").sum())
    lines = ["# Independent verification of the key results", "",
             "Generated by `python -m experiments.final_refinement.verify_all` from the committed raw result files "
             "only (no project analysis, decision or pipeline module is imported). Do not edit by hand.", "",
             f"**{len(df) - n_info} checks: {n_match} MATCH, {n_mis} MISMATCH** ({n_info} information-only row). "
             "Every MISMATCH is documented in `ERRATA.md`.", "",
             "Rules: `exact` = identical counts / flags / labels; `round` = the recomputed value rounded to the "
             "precision printed in the report equals the reported value; `mc (±t)` = Monte Carlo quantity "
             "(bootstrap percentile, permutation p) recomputed with an independent bootstrap / permutation and a "
             "different seed, MATCH iff within ±t (in the reported unit; the original analyses are seeded and "
             "deterministic); `text` = a statement in a report checked against recomputed numbers.", "",
             "Independent re-implementations used: Wilson interval; one-sided Clopper–Pearson bound by bisection on "
             "the binomial CDF; subject-cluster bootstrap and permutation test; GEE (binomial, logit; independence and "
             "exchangeable working correlation; robust sandwich SE) by Fisher scoring; L2 logistic regression "
             "(Newton) with per-fold standardisation for the Phase 7 leave-one-subject-out models; Mann–Whitney AUC; "
             "detector decisions re-derived from stored statistics (p ≤ 0.05; J < 0.05 and gate ≥ 1.2; "
             "z_max ≥ 9.4521 and G ≥ 0.25; NL > 0).", ""]
    for ph, title in (("2E", "Phase 2E (quoted in Phase 10)"), ("3", "Phase 3"), ("4", "Phase 4"), ("5", "Phase 5"),
                      ("6", "Phase 6"), ("7", "Phase 7"), ("8", "Phase 8"), ("9", "Phase 9"), ("10", "Phase 10"),
                      ("X", "Cross-report consistency (quotes of earlier phases)"),
                      ("T", "Report statements that summarise tables")):
        g = df[df.phase == ph]
        if g.empty:
            continue
        lines += [f"## {title}", "", "| id | quantity | reported | recomputed | rule | status | source | report | note |",
                  "|---|---|---|---|---|---|---|---|---|"]
        for _, r in g.iterrows():
            st = f"**{r['status']}**" if r["status"] != "MATCH" else r["status"]
            cells = [r["id"], r["quantity"], r["reported"], r["recomputed"], r["rule"], st, f"`{r['source']}`",
                     r["where"], r["note"]]
            lines.append("| " + " | ".join(str(c).replace("|", "/") for c in cells) + " |")
        lines.append("")
    (OUT / "VERIFICATION.md").write_text("\n".join(lines))
    json.dump({"n_checks": len(df) - n_info, "match": n_match, "mismatch": n_mis, "info": n_info,
               "checks": CHECKS}, open(OUT / "verification.json", "w"), indent=1, default=str)
    print(f"{len(df) - n_info} checks: {n_match} MATCH, {n_mis} MISMATCH, {n_info} INFO")
    for r in CHECKS:
        if r["status"] == "MISMATCH":
            print("MISMATCH", r["id"], r["quantity"], "| reported", r["reported"], "| recomputed", r["recomputed"])
    return n_mis


def main():
    for f in (phase2e, phase3, phase4, phase5, phase6, phase7, phase8, phase9, phase10, data_checks, cross):
        f()
        print(f"{f.__name__}: {len(CHECKS)} checks so far", file=sys.stderr)
    write()


if __name__ == "__main__":
    main()
