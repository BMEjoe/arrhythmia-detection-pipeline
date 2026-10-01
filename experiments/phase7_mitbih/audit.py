"""
Phase 7 EXPLORATORY code audit of the real-data path (chosen after the preregistered
results were seen).  For 5 normal and 5 abnormal raw (primary) windows, picked with a
fixed seed, everything is recomputed independently of the stored results and of
experiments/phase7_mitbih/data.py where possible:

  1. ECG read directly with wfdb.rdrecord (channel 0) == fp.load_mitbih_record_with_annotations
  2. fresh fp.detect_r_peaks == committed results/qc/detected_peaks.npz
  3. RR = diff(peaks / fs) computed by hand (the pipeline's arithmetic in
     fp.extract_rr_intervals); window = RR[start:start+256]; units (seconds), length,
     sum == (t1 - t0) / fs, SHA-256 == the stored rr_sha256 of the run.
     The algebraically identical diff(peaks) / fs differs in the last bit; it is run
     too, to show how far ulp-level input changes move the stored outputs (the
     lle_chaos_test surrogate RNG is seeded by a CRC of the data bytes).
  4. labels by an independent nearest-unused-annotation matcher (150 ms, beat codes) ==
     stored label and abnormal fraction
  5. detrend: least-squares trend ratio by np.linalg.lstsq == stored detrend_ratio;
     applied <=> ratio >= 0.7; analyze_segment's rr_dynamics == hand-detrended series
     (x - trend + mean(trend)) when applied, == x otherwise
  6. fresh detector.evaluate == stored decisions, p-value, statistic and scores

    python -m experiments.phase7_mitbih.audit
"""
from __future__ import annotations

import hashlib
import json
import os
import pathlib

import numpy as np
import wfdb

import final_pipeline as fp
from experiments.phase7_mitbih import analysis as A
from experiments.phase7_mitbih import data as D
from experiments.phase7_mitbih import detector as DT

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "results" / "spike_in"
BEATS = set("N L R B A a J S V r F e j n E / f Q ?".split())
NORMAL = set("N L R e j".split())
ABNORMAL = set("A a J S V E F f Q".split())


def own_labels(peaks, ann_s, ann_y, fs, tol_ms=150.0):
    keep = np.array([y in BEATS for y in ann_y])
    s, y = ann_s[keep], ann_y[keep]
    tol = int(round(tol_ms * fs / 1000.0))
    free = np.ones(len(s), bool)
    out = []
    for p in peaks:
        d = np.abs(s - p).astype(float)
        d[~free] = np.inf
        j = int(np.argmin(d)) if len(d) else -1
        if j < 0 or d[j] > tol:
            out.append(-1)
            continue
        free[j] = False
        out.append(0 if y[j] in NORMAL else 1 if y[j] in ABNORMAL else -1)
    return np.array(out)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    raw = A.load(HERE / "results" / "run", "raw")
    rng = np.random.default_rng(20261002)
    pick = []
    for lab in (0, 1):
        ids = sorted(raw[raw.label == lab].task_id)
        pick += list(rng.choice(ids, 5, replace=False))
    peaks_npz = np.load(HERE / "results" / "qc" / "detected_peaks.npz")
    rows, problems = [], []
    for tid in pick:
        st = raw[raw.task_id == tid].iloc[0]
        r, start = st.record, int(st.start_rr)
        rec = wfdb.rdrecord(os.path.join(D.DATA_DIR, r))
        ann = wfdb.rdann(os.path.join(D.DATA_DIR, r), "atr")
        sig_fp, fs, ann_s, ann_y = fp.load_mitbih_record_with_annotations(r, dl_dir=D.DATA_DIR)
        chk = {"task_id": tid, "record": r, "window": int(st.window), "label": int(st.label)}
        chk["ecg_equal"] = bool(np.array_equal(rec.p_signal[:, 0], sig_fp)) and float(rec.fs) == fs
        peaks = np.asarray(fp.detect_r_peaks(rec.p_signal[:, 0].astype(float), float(rec.fs), fp.CFG))
        chk["peaks_equal_committed"] = bool(np.array_equal(peaks, peaks_npz[r]))
        rr_all = np.diff(peaks / float(rec.fs))
        x = rr_all[start:start + 256]
        x_alt = (np.diff(peaks) / float(rec.fs))[start:start + 256]
        chk["n_intervals"] = int(len(x))
        chk["rr_min_s"], chk["rr_max_s"], chk["rr_mean_s"] = float(x.min()), float(x.max()), float(x.mean())
        chk["units_seconds"] = bool(0.2 < x.mean() < 2.0)
        chk["sum_equals_span"] = bool(np.isclose(x.sum(), (peaks[start + 256] - peaks[start]) / fs)
                                      and peaks[start] == st.t0 and peaks[start + 256] == st.t1)
        chk["sha256_equal_stored"] = hashlib.sha256(np.ascontiguousarray(x).tobytes()).hexdigest() == st.rr_sha256
        lab = own_labels(peaks, np.asarray(ann.sample), np.asarray(ann.symbol), fs)[start + 1:start + 257]
        chk["labels_all_usable"] = bool(np.all(lab >= 0))
        frac = float(np.mean(lab == 1))
        chk["abnormal_fraction_equal"] = bool(np.isclose(frac, st.abnormal_fraction))
        chk["label_equal"] = int(frac >= 0.10) == int(st.label)
        k = np.arange(256, dtype=float)
        Xd = np.column_stack([np.ones(256), k])
        beta = np.linalg.lstsq(Xd, x, rcond=None)[0]
        trend = Xd @ beta
        ratio = abs(beta[1]) * 255 / np.std(x - trend)
        chk["trend_ratio"] = float(ratio)
        chk["trend_ratio_equal_stored"] = bool(np.isclose(ratio, st.detrend_ratio, rtol=1e-9))
        chk["applied_equal_stored"] = bool((ratio >= 0.7) == bool(st.detrend_applied))
        out = fp.analyze_segment(x, DT.DETECTOR)
        expect = x - trend + trend.mean() if ratio >= 0.7 else x
        chk["rr_dynamics_matches_hand_detrend"] = bool(np.allclose(out["rr_dynamics"], expect, rtol=0, atol=1e-12))
        chk["rr_raw_is_input"] = bool(np.array_equal(out["rr_raw"], x))
        chk["config_detrend"] = out["rr_detrend"]["method"] == "linear" and out["rr_detrend"]["min_trend_sd"] == 0.7
        fresh = DT.evaluate(x)
        for key in ("and_detected", "lle_detected", "upo_detected", "lle_p", "lle_stat", "lle_z", "upo_score"):
            a, b = fresh[key], st[key]
            same = (a == b) if not isinstance(a, float) else bool(np.isclose(a, b, rtol=1e-12, atol=0))
            chk[f"fresh_{key}_equal"] = bool(same)
        chk["fresh"] = {k: fresh[k] for k in ("and_detected", "lle_detected", "upo_detected", "lle_p", "upo_score")}
        alt = DT.evaluate(x_alt)
        chk["ulp_variant"] = {"max_abs_rr_diff_s": float(np.max(np.abs(x - x_alt))),
                              "n_intervals_differing": int(np.sum(x != x_alt)),
                              **{k: alt[k] for k in ("and_detected", "lle_detected", "upo_detected", "lle_p",
                                                     "lle_stat", "upo_score")},
                              "stored_lle_p": float(st.lle_p), "stored_upo_score": float(st.upo_score),
                              "decisions_equal_stored": bool(alt["and_detected"] == st.and_detected
                                                             and alt["lle_detected"] == st.lle_detected
                                                             and alt["upo_detected"] == st.upo_detected)}
        bad = [k for k, v in chk.items() if isinstance(v, bool) and not v and k not in ("units_seconds",)]
        if not chk["units_seconds"]:
            bad.append("units_seconds")
        if bad:
            problems.append({"task_id": tid, "failed": bad})
        rows.append(chk)
        print(tid, "OK" if not bad else f"DISCREPANCY {bad}", flush=True)
    res = {"windows": rows, "discrepancies": problems, "n_checked": len(rows)}
    json.dump(res, open(OUT / "audit.json", "w"), indent=1, default=str)
    print(json.dumps({"n": len(rows), "discrepancies": problems}, indent=1))


if __name__ == "__main__":
    main()
