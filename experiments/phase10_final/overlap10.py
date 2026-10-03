"""
Phase 10 data checks after the CONFIRMATORY download (PREREGISTRATION.md Section 1), run BEFORE any
detector on confirmatory data:

  1. annotation sampling rate of every nsr2db / chf2db record (resolution of the spike-ins);
  2. record metadata (header: fs, length, start time, comments) next to the development databases';
  3. subject overlap with nsrdb / chfdb: (a) documentation (recorded in METHODS.md / the report),
     (b) metadata (age / sex where the headers carry them), (c) RR-sequence matching: 5 probes of
     500 consecutive intervals from each confirmatory record (starting at 2, 6, 10, 14, 18 h, or the
     nearest available) are matched against every development record (all lags, FFT
     cross-correlation of the mean-removed series, then the median absolute interval difference at the
     best lag); a probe MATCHES iff that median is < 2 / 128 s; a pair OVERLAPS iff >= 2 probes match.
  Overlapping confirmatory subjects are excluded: results/conf/exclusions.json.

    python -m experiments.phase10_final.overlap10
"""
from __future__ import annotations

import json
import pathlib

import numpy as np

from experiments.phase10_final import data10 as D

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "results" / "conf"
PROBE_HOURS = (2, 6, 10, 14, 18)
PROBE_LEN = 500
TOL = 2.0 / 128.0


def rr_of(db, rec, root):
    t, lab, sym, fs = D.beat_series(db, rec, root)
    return t, np.diff(t), fs


def best_match(probe, rr):
    """Best lag of probe in rr by FFT cross-correlation of mean-removed series; returns (lag, median
    absolute difference at that lag, normalised correlation)."""
    n, m = len(rr), len(probe)
    if n < m:
        return None, np.inf, 0.0
    a = rr - rr.mean()
    b = probe - probe.mean()
    L = 1 << int(np.ceil(np.log2(n + m)))
    cc = np.fft.irfft(np.fft.rfft(a, L) * np.conj(np.fft.rfft(b, L)), L)[: n - m + 1]
    # local Pearson normalisation (b has zero mean, so sum (a_win - mean_win) b = sum a_win b)
    c1 = np.concatenate([[0], np.cumsum(a)])
    c2 = np.concatenate([[0], np.cumsum(a ** 2)])
    s1 = (c1[m:] - c1[:-m])[: n - m + 1]
    s2 = (c2[m:] - c2[:-m])[: n - m + 1]
    var_win = np.maximum(s2 - s1 ** 2 / m, 0.0)
    r = cc / (np.sqrt(var_win * np.sum(b ** 2)) + 1e-12)
    lag = int(np.argmax(r))
    mad = float(np.median(np.abs(rr[lag:lag + m] - probe)))
    return lag, mad, float(r[lag])


# ------------------------------------------------------------------ AMENDMENT 1 (2026-10-03; PREREGISTRATION.md)
OFFSET_TOL = 120.0        # s: a duplicate's probes match at one common time offset (annotation slack)
MIN_CONSISTENT = 3


def consistent_overlap(tc, rc, td, rd):
    """Amended rule: a pair overlaps iff >= MIN_CONSISTENT probes match (median |diff| < TOL at the
    best-correlation lag) AND those probes' matched positions share one time offset
    (t_dev(match) - t_conf(probe)) within OFFSET_TOL s of the cluster's median offset."""
    hits = []
    for h in PROBE_HOURS:
        i = int(np.searchsorted(tc, tc[0] + h * 3600.0))
        i = min(i, max(0, len(rc) - PROBE_LEN))
        lag, mad, r = best_match(rc[i:i + PROBE_LEN], rd)
        if lag is not None and mad < TOL:
            hits.append(float(td[lag] - tc[i]))
    best = 0
    for o in hits:
        best = max(best, sum(abs(o2 - o) <= OFFSET_TOL for o2 in hits))
    return best >= MIN_CONSISTENT, {"n_match": len(hits), "n_consistent": best, "offsets_s": hits}


def amended(validate_only=False):
    """Validation on development data (a known duplicate must be found, distinct subjects must not), then
    the amended rule on every confirmatory x development pair."""
    rng = np.random.default_rng(20261012)
    val = []
    for db, rec in (("chfdb", "chf09"), ("chfdb", "chf13"), ("nsrdb", "16265"), ("nsrdb", "19830")):
        t, lab, sym, fs = D.beat_series(db, rec, D.DEV_DIR)
        # duplicate: start 37 min later, re-annotated on a 1/128 s grid with +/- 1 sample jitter
        k0 = int(np.searchsorted(t, t[0] + 2220.0))
        tq = np.round(t[k0:] * 128.0 + rng.integers(-1, 2, len(t) - k0)) / 128.0
        tq = np.sort(tq) - tq[0]
        for db2, rec2 in (("chfdb", "chf09"), ("chfdb", "chf13"), ("nsrdb", "16265"), ("nsrdb", "19830")):
            t2, rr2, _ = rr_of(db2, rec2, D.DEV_DIR)
            ok, info = consistent_overlap(tq, np.diff(tq), t2, rr2)
            val.append({"duplicate_of": f"{db}:{rec}", "against": f"{db2}:{rec2}", "same": rec == rec2,
                        "overlap": ok, **info})
    v_ok = all(x["overlap"] == x["same"] for x in val)
    out = {"validation": val, "validation_pass": v_ok, "pairs": [], "excluded_subjects": []}
    if validate_only:
        return out
    dev = {}
    for db in ("nsrdb", "chfdb"):
        for rec in D.records(db, D.DEV_DIR):
            t, rr, _ = rr_of(db, rec, D.DEV_DIR)
            dev[f"{db}:{rec}"] = (t, rr)
    for db in D.CONF_DBS:
        for rec in D.records(db, D.CONF_DIR):
            tc, rc, _ = rr_of(db, rec, D.CONF_DIR)
            for did, (td, rd) in dev.items():
                ok, info = consistent_overlap(tc, rc, td, rd)
                out["pairs"].append({"conf": f"{db}:{rec}", "dev": did, "overlap": ok, **info})
                if ok and f"{db}:{rec}" not in out["excluded_subjects"]:
                    out["excluded_subjects"].append(f"{db}:{rec}")
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    report = {"resolution": {}, "metadata": {"conf": {}, "dev": {}}, "rr_matching": [], "excluded_subjects": []}
    dev = {}
    for db in ("nsrdb", "chfdb"):
        for rec in D.records(db, D.DEV_DIR):
            t, rr, fs = rr_of(db, rec, D.DEV_DIR)
            dev[f"{db}:{rec}"] = rr
            report["metadata"]["dev"][f"{db}:{rec}"] = D.header_info(db, rec, D.DEV_DIR) | {"ann_fs": fs, "n_beats": len(t),
                                                                                         "duration_h": float((t[-1] - t[0]) / 3600)}
    for db in D.CONF_DBS:
        for rec in D.records(db, D.CONF_DIR):
            t, rr, fs = rr_of(db, rec, D.CONF_DIR)
            sid = f"{db}:{rec}"
            report["resolution"][sid] = fs
            report["metadata"]["conf"][sid] = D.header_info(db, rec, D.CONF_DIR) | {"ann_fs": fs, "n_beats": len(t),
                                                                                    "duration_h": float((t[-1] - t[0]) / 3600)}
            probes = []
            for h in PROBE_HOURS:
                i = int(np.searchsorted(t, t[0] + h * 3600.0))
                i = min(i, max(0, len(rr) - PROBE_LEN))
                probes.append(rr[i:i + PROBE_LEN])
            for did, drr in dev.items():
                res = [best_match(p, drr) for p in probes]
                nmatch = int(sum(r[1] < TOL for r in res))
                report["rr_matching"].append({"conf": sid, "dev": did, "n_probe_match": nmatch,
                                              "min_median_abs_diff_s": float(min(r[1] for r in res)),
                                              "max_corr": float(max(r[2] for r in res))})
                if nmatch >= 2 and sid not in report["excluded_subjects"]:
                    report["excluded_subjects"].append(sid)
    fss = sorted(set(report["resolution"].values()))
    report["resolution_summary"] = {"annotation_fs_values": fss,
                                    "coarser_than_1_360": bool(min(fss) < 360.0),
                                    "spike_in_quantization": "1/fs of each record's annotation file (data10.beat_series)"}
    m = report["rr_matching"]
    report["rr_matching_summary"] = {"pairs": len(m), "overlapping_pairs": [x for x in m if x["n_probe_match"] >= 2],
                                     "smallest_median_abs_diff_s": float(min(x["min_median_abs_diff_s"] for x in m)),
                                     "largest_corr": float(max(x["max_corr"] for x in m))}
    json.dump({"excluded_subjects": report["excluded_subjects"]}, open(OUT / "exclusions.json", "w"), indent=1)
    json.dump(report, open(OUT / "data_checks.json", "w"), indent=1, default=str)
    print(json.dumps({"resolution": report["resolution_summary"], "excluded": report["excluded_subjects"],
                      "matching": {k: v for k, v in report["rr_matching_summary"].items() if k != "overlapping_pairs"},
                      "n_overlapping_pairs": len(report["rr_matching_summary"]["overlapping_pairs"])}, indent=1))


if __name__ == "__main__":
    import sys
    if "--amended" in sys.argv:
        res = amended()
        OUT.mkdir(parents=True, exist_ok=True)
        json.dump(res, open(OUT / "data_checks_amended.json", "w"), indent=1, default=str)
        prev = json.load(open(OUT / "data_checks.json"))["excluded_subjects"]
        json.dump({"excluded_subjects": res["excluded_subjects"], "rule": "amendment 1 (time-consistent matching)",
                   "originally_flagged_subjects": prev}, open(OUT / "exclusions.json", "w"), indent=1)
        print(json.dumps({"validation_pass": res["validation_pass"],
                          "validation": [(x["duplicate_of"], x["against"], x["overlap"], x["n_match"], x["n_consistent"])
                                         for x in res["validation"]],
                          "excluded": res["excluded_subjects"],
                          "max_consistent_any_pair": max(p["n_consistent"] for p in res["pairs"])}, indent=1))
    else:
        main()
