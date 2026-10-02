"""
Phase 9 Part E (development only): diagnosis statistics bank on DEVELOPMENT material.

    python -m experiments.phase9_waveform.dev_diag --workers 4

Material (seeds 9000-9499 only; no TEST family, no TEST seed):
  - the 17 PASS nulls at the input level (Phase 8 nulls.generate, quantized, 512) and at the
    realistic ECG level (Phase 9 null_window -> ECGSYN -> nstdb mix 12 dB (DEV signal 0) -> MIT-BIH
    quantization -> Pan-Tompkins beats -> Phase 7 V-offset jitter);
  - Mackey-Glass (DEV, 7 CHAOTIC, 3 NON-CHAOTIC): RR level (iv_q, vi_S2_5, vi_E1, vi_E3_10; Phase 8
    series.window) and ECG level (realistic, ectopy none / S2_5 / E1 / E3_10);
  - circle-map development controls (devmaps.py, NON-CHAOTIC: quasi-periodic and locked) at the ECG
    level, same four ectopy inputs;
  - Phase 5 positives P1, P2, G1, G2, G3 (native).
Per window (512 RR intervals) the data and 19 ectopy-preserving IAAFT surrogates (Phase 8
methods.ep_iaaft) get the statistics below; stored: data value, surrogate mean and SD.
  nlp_m{2..5}_h{1..5}  robust normalized NLP error (Phase 8 nl_pred_error, k 5, Theiler 5)
  sm_rho_m3_h{1..5}    Sugihara-May forecast correlation (Pearson r of prediction vs target), m 3
  sdle_m3_s{0..3}      SDLE (measures.sdle, m 3, 4 default shells, T 12): mean lambda over t = 4..10
  fsle_m3_l{..}        FSLE (measures.fsle, m 3, delta0 0.01 SD, r sqrt 2, 14 levels)
  epsent_m{2,3}_e{..}  eps-entropy h_m(eps) at eps = 0.25, 0.5, 1.0 SD (max norm)
  pe_H, pe_C           Bandt-Pompe d 4 (Rosso CECP)
  det                  RQA DET (m 3, recurrence rate 0.05, l_min 2, Theiler 1)
  gh_nlp_m{2,3}        NLP error (h 1) after GHKSS (m 5, q 2, k_min 20, 3 iterations) of the series
                       (surrogates denoised identically)
plus the Phase 8 IAAFT z of the NLP (C2 statistic, m 2..5, 39 surrogates) for reference.
Output: results/dev/diag.jsonl (resumable).
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import time
import zlib
from multiprocessing import Pool

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "results" / "dev"
N = 512
N_SUR = 19
FS = 360


def nlp_curve(x, m, hs=(1, 2, 3, 4, 5), k=5, theiler=5):
    """Robust NLP error and Sugihara-May correlation for several horizons (same neighbours)."""
    from experiments.phase8_cardiac import methods as MM
    x = np.asarray(x, float)
    E = MM.embed(x, m, 1)
    t0 = m - 1
    hmax = max(hs)
    nvec = len(E) - hmax
    A = E[:nvec]
    d = ((A[:, None, :] - A[None, :, :]) ** 2).sum(-1)
    idx = np.arange(nvec)
    d[np.abs(idx[:, None] - idx[None, :]) <= theiler] = np.inf
    nn = np.argpartition(d, k, axis=1)[:, :k]
    mad = max(np.median(np.abs(x - np.median(x))), 1e-12)
    err, rho = {}, {}
    for h in hs:
        tg = x[t0 + h: t0 + h + nvec]
        pr = tg[nn].mean(1)
        err[h] = float(np.median(np.abs(tg - pr)) / mad)
        rho[h] = float(np.corrcoef(tg, pr)[0, 1]) if np.std(pr) > 0 else 0.0
    return err, rho


def stats(x):
    from experiments.phase9_waveform import measures as M
    out = {}
    x = np.asarray(x, float)
    for m in (2, 3, 4, 5):
        e, r = nlp_curve(x, m)
        for h in e:
            out[f"nlp_m{m}_h{h}"] = e[h]
            if m == 3:
                out[f"sm_rho_m3_h{h}"] = r[h]
    try:
        s = M.sdle(x, m=3, L=1, T=12, rng=np.random.default_rng(0))
        for i, sh in enumerate(s["shells"]):
            out[f"sdle_m3_s{i}"] = float(np.mean(sh["lam"][3:9])) if "lam" in sh else float("nan")
    except Exception:                                   # noqa: BLE001
        pass
    dl, lam, cnt, _ = M.fsle(x, m=3, L=1, delta0=0.01, n_levels=14, maxt=50)
    for i in range(len(dl)):
        out[f"fsle_m3_l{i}"] = float(lam[i]) if cnt[i] >= 5 else float("nan")
    eps = np.array([0.25, 0.5, 1.0])
    _, C = M.correlation_sums(x, m_max=4, L=1, eps=eps, theiler=3)
    h = M.eps_entropy(C)
    for m in (2, 3):
        for j, e in enumerate(eps):
            out[f"epsent_m{m}_e{e}"] = float(h[m - 1, j])
    H, Cc = M.pe_cecp(x, d=4, tau=1)
    out["pe_H"], out["pe_C"] = float(H), float(Cc)
    out["det"] = float(M.rqa_det(x, m=3, tau=1, rr=0.05, lmin=2, theiler=1)[0])
    y = M.ghkss(x, m=5, q=2, k_min=20, iterations=3)
    for m in (2, 3):
        out[f"gh_nlp_m{m}"] = nlp_curve(y, m, hs=(1,))[0][1]
    return out


def series(task):
    from experiments.phase8_cardiac import nulls as NL8
    from experiments.phase8_cardiac import series as SR8
    from experiments.phase9_waveform import library as L
    kind, src, name, var, seed = task["kind"], task["source"], task["name"], task["variant"], task["seed"]
    if kind == "rr_null":
        return NL8.generate(name, seed, N), {}
    if kind == "rr_model":
        reg = {r["name"]: r for r in SR8.load_regimes()}[name]
        return SR8.window(reg["family"], reg["regime_idx"], reg["params"], N, var, seed), {}
    # ECG level, realistic: ectopy (models / devmaps) or built-in (nulls), mix12, detected, jitter
    ect = "none" if src == "null" else var
    w = L.build(src, name, ect, N, seed, "mix12", "DEV", "detected", True)
    if "error" in w:
        raise RuntimeError(w["error"])
    return np.asarray(w["rr"]), {"n_missed": w["n_missed"], "n_extra": w["n_extra"]}


def one(task):
    from experiments.phase8_cardiac import methods as MM
    t0 = time.time()
    out = dict(task)
    try:
        x, info = series(task)
        out.update(info)
        out["gen_s"] = time.time() - t0
        rng = np.random.default_rng(np.random.SeedSequence([20261005, task["seed"], zlib.crc32(task["tid"].encode())]))
        mask = MM.outlier_mask(x)
        out["n_flagged"] = int(mask.sum())
        d = stats(x)
        sur = [stats(MM.ep_iaaft(x, rng, mask)) for _ in range(N_SUR)]
        for k, v in d.items():
            sv = np.array([s.get(k, np.nan) for s in sur], float)
            out[k] = v
            out[k + "__sm"] = float(np.nanmean(sv)) if np.isfinite(sv).any() else float("nan")
            out[k + "__ss"] = float(np.nanstd(sv, ddof=1)) if np.isfinite(sv).sum() > 1 else float("nan")
            out[k + "__le"] = int(np.sum(sv <= v))
        c2 = MM._nlpred_scan(x)
        out["c2_zmax"] = float(max(c2[f"np_m{m}_z"] for m in (2, 3, 4, 5)))
        c3 = MM._nlpred_ep_scan(x)
        out["c3_zmax"] = float(max(c3[f"np_m{m}_z"] for m in (2, 3, 4, 5)))
    except Exception as exc:                            # noqa: BLE001
        out["error"] = f"{type(exc).__name__}: {exc}"
    out["t_total"] = time.time() - t0
    return out


def tasks(n_null=20, n_model=10):
    from experiments.phase8_cardiac import nulls as NL8
    from experiments.phase9_waveform import devmaps as DM
    from experiments.phase9_waveform import families as F
    T = []
    for c in NL8.PASS_NULLS:
        for s in range(9000, 9000 + n_null):
            T.append({"kind": "rr_null", "source": "null", "name": c, "variant": "input", "seed": s, "label": "NULL"})
            T.append({"kind": "ecg", "source": "null", "name": c, "variant": "real", "seed": s, "label": "NULL"})
    for c in ("P1_henon_rr", "P2_logistic_rr", "G1_lorenz_maxima", "G2_rossler_flow", "G3_mackey_glass"):
        for s in range(9000, 9000 + n_null):
            T.append({"kind": "rr_null", "source": "null", "name": c, "variant": "native", "seed": s, "label": "CHAOTIC"})
    for r in F.load_regimes():
        if r["family"] != "mackey_glass":
            continue
        for s in range(9000, 9000 + n_model):
            for v in ("iv_q", "vi_S2_5", "vi_E1", "vi_E3_10"):
                T.append({"kind": "rr_model", "source": "model", "name": r["name"], "variant": v, "seed": s,
                          "label": r["label"]})
            for e in ("none", "S2_5", "E1", "E3_10"):
                T.append({"kind": "ecg", "source": "model", "name": r["name"], "variant": e, "seed": s,
                          "label": r["label"]})
    for r in DM.REGIMES:
        for s in range(9000, 9000 + n_model):
            for e in ("none", "S2_5", "E1", "E3_10"):
                T.append({"kind": "ecg", "source": "devmap", "name": r["name"], "variant": e, "seed": s,
                          "label": "NON-CHAOTIC"})
    for t in T:
        t["tid"] = f"{t['kind']}|{t['name']}|{t['variant']}|{t['seed']}"
    return T


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "diag.jsonl"
    done = set()
    if path.exists():
        done = {json.loads(line)["tid"] for line in open(path)}
    T = [t for t in tasks() if t["tid"] not in done]
    if a.limit:
        T = T[:a.limit]
    print(len(T), "to run", flush=True)
    with open(path, "a") as fh, Pool(a.workers, maxtasksperchild=20) as pool:
        for i, r in enumerate(pool.imap_unordered(one, T, chunksize=1), 1):
            fh.write(json.dumps(r, default=float) + "\n")
            fh.flush()
            if i % 50 == 0:
                print(i, "done", flush=True)


if __name__ == "__main__":
    main()
