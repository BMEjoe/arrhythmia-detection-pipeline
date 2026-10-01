"""
Phase 7 Step 1: data-quality checks on MIT-BIH.  Uses only the ECG, the
Pan-Tompkins R-peak detector and the reference annotations.  No
analyze_segment and no chaos test is run here.

    python -m experiments.phase7_mitbih.qc

Writes results/qc/{records.csv, windows.csv, candidates.csv, symbols.csv,
summary.json, detected_peaks.npz} and QC.md is written by hand from them.
"""
from __future__ import annotations

import json
import os
from collections import Counter

import numpy as np
import pandas as pd
import wfdb

import final_pipeline as fp
from experiments.phase7_mitbih import data as D

HERE = os.path.dirname(__file__)
OUT = os.path.join(HERE, "results", "qc")
PEAK_CACHE = os.path.join(D.DATA_DIR, "_phase7_peaks")
W = fp.CLASSIFIER_WINDOW_RR
V_TYPE = frozenset({"V", "E", "F"})        # ventricular / ventricular fusion
S_TYPE = frozenset({"A", "a", "J", "S"})   # supraventricular
AF_RHYTHMS = ("(AFIB", "(AFL")


def rhythm_per_sample(ann):
    """Rhythm label in force at each annotation (from '+' aux notes)."""
    cur = "(N"
    rhythms = []
    for y, a in zip(ann.symbol, ann.aux_note):
        if y == "+" and a.strip("\x00").strip():
            cur = a.strip("\x00").strip()
        rhythms.append(cur)
    return np.asarray(rhythms)


def af_beats(record, beats):
    """True for each beat time lying in an AFIB/AFL rhythm episode."""
    ann = wfdb.rdann(os.path.join(D.DATA_DIR, record), "atr")
    s = np.asarray(ann.sample)
    rh = rhythm_per_sample(ann)
    plus = np.flatnonzero(np.asarray(ann.symbol) == "+")
    if not len(plus):
        return np.zeros(len(beats), bool)
    ps, pr = s[plus], rh[plus]
    k = np.searchsorted(ps, beats, side="right") - 1
    return np.array([(kk >= 0 and pr[kk].startswith(AF_RHYTHMS)) for kk in k])


def window_rows(series, extra_mask=None, missed_samples=None, sym_of_beat=None, af=None):
    rows = []
    labels = series["labels"]
    beats = series["beats"]
    for w in D.windows(series):
        s, e = w["start_rr"], w["end_rr"]
        syms = sym_of_beat[s + 1:e + 1] if sym_of_beat is not None else None
        n_v = int(np.isin(syms, list(V_TYPE)).sum()) if syms is not None else np.nan
        n_s = int(np.isin(syms, list(S_TYPE)).sum()) if syms is not None else np.nan
        n_q = int(np.sum(syms == "Q")) if syms is not None else np.nan
        hr = D.hrv(w["rr"])
        he = D.hrv(w["rr_edited"]) if w["rr_edited"] is not None else {k: np.nan for k in hr}
        n_missed = (int(np.sum((missed_samples > beats[s]) & (missed_samples < beats[e])))
                    if missed_samples is not None else 0)
        rows.append({
            "record": w["record"], "subject": w["subject"], "window": w["window"],
            "start_rr": s, "t0": w["t0"], "t1": w["t1"], "label": w["label"],
            "abnormal_fraction": w["abnormal_fraction"], "n_abnormal": w["n_abnormal"],
            "n_v_type": n_v, "n_s_type": n_s, "n_q": n_q,
            "beat_type": ("normal" if w["label"] == 0 else
                          ("ventricular" if n_v > n_s else "supraventricular")) if syms is not None else "",
            "nn_fraction": w["nn_fraction"], "edited_eligible": bool(w["nn_fraction"] >= 0.80),
            "n_missed_beats_in_span": n_missed,
            "start_beat_unusable": bool(labels[s] < 0),
            "af_fraction": float(np.mean(af[s + 1:e + 1])) if af is not None else np.nan,
            "rr_min": float(np.min(w["rr"])), "rr_max": float(np.max(w["rr"])),
            "rr_mean": float(np.mean(w["rr"])),
            "sdnn_ms": hr["sdnn_ms"], "rmssd_ms": hr["rmssd_ms"], "pnn50": hr["pnn50"],
            "sdnn_ms_edited": he["sdnn_ms"], "rmssd_ms_edited": he["rmssd_ms"], "pnn50_edited": he["pnn50"],
        })
    return rows


def _validate(peaks, ref_s, fs, tol):
    v = fp.validate_r_peaks_against_annotations(peaks, ref_s, fs, tol)
    return {k: v[k] for k in ("TP", "FP", "FN", "sensitivity", "PPV")}


def main():
    os.makedirs(OUT, exist_ok=True)
    rec_rows, win_rows, cand_rows, ann_win_rows = [], [], [], []
    sym_counts = Counter()
    all_peaks = {}
    for r in D.RECORDS:
        sig, fs, ann_s, ann_y = D.load(r)
        sym_counts.update(ann_y.tolist())
        beat = np.isin(ann_y, list(D.BEAT_SYMBOLS))
        ref_s, ref_y = ann_s[beat], ann_y[beat]
        peaks = D.detect(r, PEAK_CACHE)
        all_peaks[r] = peaks.astype(np.int32)
        rr, _, _ = fp.extract_rr_intervals(peaks, fs, fp.CFG.rr_min_seconds, fp.CFG.rr_max_seconds)
        af = af_beats(r, peaks)
        rec = {"record": r, "subject": D.SUBJECT[r], "fs": fs, "n_ref_beats": int(beat.sum()),
               "n_detected": int(len(peaks)), "n_rr": int(len(rr)),
               "paced": r in D.PACED_RECORDS, "af_record": bool(af.any()),
               "n_paced_beats": int(np.sum(ref_y == "/")), "n_Q_beats": int(np.sum(ref_y == "Q"))}
        per_tol = {}
        for tol, tag in ((D.TOL_PRIMARY_MS, "150"), (D.TOL_DEFAULT_MS, "75")):
            v = _validate(peaks, ref_s, fs, tol)
            v_all = fp.validate_r_peaks_against_annotations(peaks, ann_s, fs, tol)
            match = D.greedy_match(peaks, ref_s, fs, tol)
            extra = match < 0
            matched_ref = np.zeros(len(ref_s), bool)
            matched_ref[match[match >= 0]] = True
            missed_s = ref_s[~matched_ref]
            sym_of_beat = np.where(match >= 0, ref_y[np.maximum(match, 0)], "")
            lab = fp.match_detected_peaks_to_annotations(peaks, ref_s, ref_y, fs, tolerance_ms=tol)
            lab_all = fp.match_detected_peaks_to_annotations(peaks, ann_s, ann_y, fs, tolerance_ms=tol)
            assert np.array_equal(lab, np.array([-1 if j < 0 else D.effective_label(ref_y[j]) for j in match]))
            per_tol[tag] = (extra, missed_s, sym_of_beat, lab)
            for k, val in v.items():
                rec[f"{k}_{tag}"] = val
            rec[f"FN_vs_all_annotations_{tag}"] = int(v_all["FN"])
            rec[f"n_label_diff_beat_vs_all_ref_{tag}"] = int(np.sum(lab != lab_all))
            rec[f"n_windows_all_ref_{tag}"] = len(fp.make_rr_windows(rr, lab_all))
            # candidate windows (every start position, before label filtering)
            for start in range(0, len(rr) - W + 1, fp.CLASSIFIER_STEP_RR):
                end = start + W
                cand_rows.append({
                    "tol_ms": tag, "record": r, "subject": D.SUBJECT[r], "window": start // W,
                    "n_extra_detections": int(extra[start + 1:end + 1].sum()),
                    "n_missed_beats": int(np.sum((missed_s > peaks[start]) & (missed_s < peaks[end]))),
                    "n_excluded_symbol_beats": int(np.sum((~extra[start + 1:end + 1]) & (lab[start + 1:end + 1] < 0))),
                    "extra_at_start_beat": bool(extra[start]),
                    "retained": bool(np.all(lab[start + 1:end + 1] >= 0)),
                    "retained_all_annotations": bool(np.all(lab_all[start + 1:end + 1] >= 0))})
        # PRIMARY windows (150 ms) with the 75 ms label attached for the descriptive arm
        extra, missed_s, sym_of_beat, lab = per_tol["150"]
        series = {"record": r, "subject": D.SUBJECT[r], "fs": fs, "beats": peaks, "rr": rr, "labels": lab}
        rows = window_rows(series, extra, missed_s, sym_of_beat, af)
        lab75 = per_tol["75"][3]
        for row in rows:
            s0 = row["start_rr"]
            l75 = lab75[s0 + 1:s0 + W + 1]
            row["retained_tol75"] = bool(np.all(l75 >= 0))
            row["label_tol75"] = int(np.mean(l75 == 1) >= fp.CLASSIFIER_ABNORMAL_FRACTION) if row["retained_tol75"] else -1
        win_rows += rows
        rec["n_candidate_windows"] = len(range(0, len(rr) - W + 1, fp.CLASSIFIER_STEP_RR))
        rec["n_windows"] = len(rows)
        rec["n_windows_normal"] = sum(x["label"] == 0 for x in rows)
        rec["n_windows_abnormal"] = sum(x["label"] == 1 for x in rows)
        rec["n_windows_tol75"] = sum(x["retained_tol75"] for x in rows)
        # annotation-time series (sensitivity analysis)
        aser = D.record_series(r, source="annotation")
        bs, by, _ = D.annotation_beat_series(ann_s, ann_y)
        arows = window_rows(aser, None, None, by, af_beats(r, bs))
        ann_win_rows += arows
        rec["n_windows_annotation_times"] = len(arows)
        rec["n_windows_annotation_times_abnormal"] = sum(x["label"] == 1 for x in arows)
        rec_rows.append(rec)
        print(r, round(rec["sensitivity_150"], 4), round(rec["PPV_150"], 4), rec["n_windows"], flush=True)

    recs = pd.DataFrame(rec_rows)
    wins = pd.DataFrame(win_rows)
    cands = pd.DataFrame(cand_rows)
    awins = pd.DataFrame(ann_win_rows)
    recs.to_csv(os.path.join(OUT, "records.csv"), index=False)
    wins.to_csv(os.path.join(OUT, "windows.csv"), index=False)
    cands.to_csv(os.path.join(OUT, "candidates.csv"), index=False)
    awins.to_csv(os.path.join(OUT, "windows_annotation_times.csv"), index=False)
    np.savez_compressed(os.path.join(OUT, "detected_peaks.npz"), **all_peaks)
    pd.DataFrame([{"symbol": k, "count": n,
                   "is_beat": k in D.BEAT_SYMBOLS, "effective_label": D.effective_label(k)}
                  for k, n in sorted(sym_counts.items(), key=lambda kv: -kv[1])]
                 ).to_csv(os.path.join(OUT, "symbols.csv"), index=False)
    summary = {}
    for tag in ("75", "150"):
        tp, fpn, fn = recs[f"TP_{tag}"].sum(), recs[f"FP_{tag}"].sum(), recs[f"FN_{tag}"].sum()
        summary[f"tol_{tag}"] = {"TP": int(tp), "FP": int(fpn), "FN": int(fn),
                                 "sensitivity": float(tp / (tp + fn)), "PPV": float(tp / (tp + fpn))}
    summary.update({
        "n_candidate_windows": int(recs.n_candidate_windows.sum()),
        "n_windows": int(len(wins)), "n_normal": int((wins.label == 0).sum()),
        "n_abnormal": int((wins.label == 1).sum()),
        "n_subjects_with_windows": int(wins.subject.nunique()),
        "n_subjects_with_abnormal_windows": int(wins[wins.label == 1].subject.nunique()),
        "n_windows_tol75": int(wins.retained_tol75.sum()),
        "n_windows_annotation_times": int(len(awins)),
        "n_abnormal_annotation_times": int((awins.label == 1).sum())})
    json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
