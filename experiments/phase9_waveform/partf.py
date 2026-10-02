"""
Phase 9 Part F (PREREGISTRATION_F.md): the adopted winner on nsrdb vs chfdb, then EXPLORATORY MIT-BIH.

    python -m experiments.phase9_waveform.partf --run --workers 4
    python -m experiments.phase9_waveform.partf --analyze
    python -m experiments.phase9_waveform.partf --mitbih --workers 4     # EXPLORATORY

Data: partf_data/{nsrdb,chfdb} (gitignored; headers + beat annotations, SHA-256 checked).
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
DATA = HERE / "partf_data"
OUT = HERE / "results" / "partf"
BEATS = set("NLRBAaJSVrFejnE/fQ?")
NORMAL = set("NLRBej")
ANN_EXT = {"nsrdb": "atr", "chfdb": "ecg"}
N = 512


def beat_series(db, rec):
    import wfdb
    a = wfdb.rdann(str(DATA / db / rec), ANN_EXT[db])
    fs = float(a.fs) if a.fs else float(wfdb.rdheader(str(DATA / db / rec)).fs)
    sym = np.array(a.symbol)
    smp = np.asarray(a.sample)
    keep = np.array([s in BEATS for s in sym])
    sym, smp = sym[keep], smp[keep]
    t = smp / fs
    lab = np.array(["N" if s in NORMAL else "X" for s in sym])
    rr = np.diff(t)
    bad = (rr < 0.25) | (rr > 2.5)
    lab[1:][bad] = "X"
    return t, lab


def windows(t, lab):
    out = []
    last_end = -1
    for j in range(1, 11):
        i0 = int(np.searchsorted(t, j * 3600.0))
        if i0 + N >= len(t) or i0 <= last_end:
            continue
        out.append((j, i0))
        last_end = i0 + N
    return out


def one(task):
    import final_pipeline as fp
    db, rec, j, i0 = task
    t, lab = beat_series(db, rec)
    rr = np.diff(t[i0:i0 + N + 1])
    lb = lab[i0:i0 + N + 1]
    r = fp.masked_growth_chaos_test(rr, lb, fp.CFG)
    return {"db": db, "subject": rec, "window": j, "start_s": float(t[i0]), "n_nonN": int(np.sum(lb != "N")),
            "rr_mean": float(rr.mean()), **{k: r.get(k) for k in ("detected", "analysable", "n_masked", "zmax", "G")}}


def run(workers):
    T = []
    for db in ("nsrdb", "chfdb"):
        for rec in (DATA / db / "RECORDS").read_text().split():
            t, lab = beat_series(db, rec)
            T += [(db, rec, j, i0) for j, i0 in windows(t, lab)]
    OUT.mkdir(parents=True, exist_ok=True)
    with Pool(workers) as pool:
        R = pool.map(one, T, chunksize=1)
    with open(OUT / "partf_windows.jsonl", "w") as fh:
        for r in R:
            fh.write(json.dumps(r, default=float) + "\n")
    print(len(R), "windows")


def _subject_rates(R, db):
    subs = sorted({r["subject"] for r in R if r["db"] == db})
    return np.array([np.mean([r["detected"] for r in R if r["subject"] == s]) for s in subs]), subs


def analyze():
    R = [json.loads(line) for line in open(OUT / "partf_windows.jsonl")]
    FP_UPPER = 1 - 0.05 ** (1 / 10100)
    res = {"synthetic_fp_upper95": FP_UPPER, "W1": {}}
    rates = {}
    for db, seed in (("nsrdb", 20261007), ("chfdb", 20261007)):
        p, subs = _subject_rates(R, db)
        rates[db] = p
        rng = np.random.default_rng(seed)
        bs = [np.mean(p[rng.integers(0, len(p), len(p))]) for _ in range(10000)]
        lo, hi = np.percentile(bs, [2.5, 97.5])
        w = [r for r in R if r["db"] == db]
        res["W1"][db] = {"n_subjects": len(p), "n_windows": len(w), "detected_windows": int(sum(r["detected"] for r in w)),
                         "analysable_windows": int(sum(bool(r["analysable"]) for r in w)),
                         "subject_rate": float(np.mean(p)), "ci95": [float(lo), float(hi)],
                         "exceeds_synthetic_fp": bool(lo > FP_UPPER),
                         "masked_fraction_median": float(np.median([r["n_masked"] / N for r in w])),
                         "per_subject": dict(zip(subs, map(float, p)))}
    a, b = rates["chfdb"], rates["nsrdb"]
    rng = np.random.default_rng(20261008)
    bs = [np.mean(a[rng.integers(0, len(a), len(a))]) - np.mean(b[rng.integers(0, len(b), len(b))]) for _ in range(10000)]
    lo, hi = np.percentile(bs, [2.5, 97.5])
    d0 = float(np.mean(a) - np.mean(b))
    rng = np.random.default_rng(20261009)
    allp = np.concatenate([a, b])
    perm = []
    for _ in range(10000):
        q = rng.permutation(allp)
        perm.append(np.mean(q[:len(a)]) - np.mean(q[len(a):]))
    pval = float((1 + np.sum(np.abs(perm) >= abs(d0) - 1e-15)) / 10001)
    res["W2"] = {"diff_chf_minus_nsr": d0, "ci95": [float(lo), float(hi)], "perm_p_two_sided": pval,
                 "supported": bool(lo > 0 or hi < 0)}
    zs = {db: [r["zmax"] for r in R if r["db"] == db and r["zmax"] is not None] for db in ("nsrdb", "chfdb")}
    gs = {db: [r["G"] for r in R if r["db"] == db and r["G"] is not None] for db in ("nsrdb", "chfdb")}
    res["secondary"] = {db: {"zmax_median": float(np.nanmedian(zs[db])) if zs[db] else None,
                             "zmax_ge_Z": int(np.sum(np.array(zs[db]) >= 9.4521)),
                             "G_median": float(np.nanmedian(gs[db])) if gs[db] else None} for db in zs}
    json.dump(res, open(OUT / "partf_results.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items()}, indent=1, default=float))


# --------------------------------------------------------------------------- EXPLORATORY MIT-BIH
def mitbih_one(task):
    import final_pipeline as fp
    rec, w0, rr, lab, cls = task
    r = fp.masked_growth_chaos_test(np.asarray(rr), np.asarray(lab), fp.CFG)
    return {"record": rec, "start_beat": w0, "class": cls,
            **{k: r.get(k) for k in ("detected", "analysable", "n_masked", "zmax", "G")}}


def mitbih(workers):
    import wfdb
    from experiments.phase7_mitbih import data as D7
    T = []
    for rec in D7.RECORDS:
        a = wfdb.rdann(str(pathlib.Path(D7.DATA_DIR) / rec), "atr")
        sym, smp = np.array(a.symbol), np.asarray(a.sample)
        keep = np.array([s in BEATS for s in sym])
        sym, smp = sym[keep], smp[keep]
        t = smp / float(a.fs or 360.0)
        lab = np.array(["N" if D7.effective_label(s) == 0 else "X" for s in sym])
        rr = np.diff(t)
        lab[1:][(rr < 0.25) | (rr > 2.5)] = "X"
        for w0 in range(0, len(t) - N - 1, N):
            ll = lab[w0:w0 + N + 1]
            T.append((rec, w0, rr[w0:w0 + N].tolist(), ll.tolist(), "abnormal" if np.mean(ll != "N") > 0.2 else "mostly_normal"))
    with Pool(workers) as pool:
        R = pool.map(mitbih_one, T, chunksize=1)
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "mitbih_exploratory.jsonl", "w") as fh:
        for r in R:
            fh.write(json.dumps(r, default=float) + "\n")
    print(len(R), "windows;", sum(r["detected"] for r in R), "detected;", sum(bool(r["analysable"]) for r in R), "analysable")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--analyze", action="store_true")
    ap.add_argument("--mitbih", action="store_true")
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args(argv)
    if a.run:
        run(a.workers)
    if a.analyze:
        analyze()
    if a.mitbih:
        mitbih(a.workers)


if __name__ == "__main__":
    main()
