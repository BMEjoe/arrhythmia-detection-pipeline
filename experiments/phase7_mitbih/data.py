"""
Phase 7 MIT-BIH data layer.  Reuses the existing pipeline helpers unchanged:
  fp.load_mitbih_record_with_annotations  (ECG channel 0, atr annotations)
  fp.detect_r_peaks                       (Phase 1 Pan-Tompkins, CFG)
  fp.extract_rr_intervals                 (RR = diff of peak times; no filtering)
  fp.match_detected_peaks_to_annotations  (greedy one-to-one, 75 ms)
  fp.make_rr_windows                      (256 / 256 / >= 10 % abnormal)
  fp.validate_r_peaks_against_annotations (sensitivity / PPV)
and the Phase 6 annotation edit experiments.phase6_robust.systems.edit_nn.

Nothing here runs analyze_segment or any chaos test.

Reference beats.  The atr files mix beat and non-beat annotations (rhythm
"+", noise "~", artifact "|", comment '"', flutter wave "!", VF onset/offset
"[" "]", non-conducted P "x").  The QRS (beat) annotation codes are the WFDB
beat set below (PhysioNet "MIT-BIH Arrhythmia Database Directory", table of
beat annotations; WFDB "isqrs").  wfdb-python 4.3.1's own `is_qrs` table is
misaligned with its `ann_labels` (it marks "[", "]", "x" as QRS and "f", "e"
as non-QRS), so it is NOT used.
"""
from __future__ import annotations

import os

import numpy as np

import final_pipeline as fp
from experiments.phase6_robust.systems import edit_nn

DATA_DIR = os.environ.get("PHASE7_DATA_DIR", "mitdb_data")
RECORDS = list(fp.MITDB_RECORDS)
PACED_RECORDS = ("102", "104", "107", "217")
BEAT_SYMBOLS = frozenset("N L R B A a J S V r F e j n E / f Q ?".split())
VF_SYMBOLS = frozenset({"!", "[", "]"})   # flutter waves and VF episode bounds (non-beat)

# MIT-BIH Arrhythmia Database directory: 48 records from 47 subjects; records
# 201 and 202 came from the same subject.  Every grouping uses SUBJECT.
# Label-transfer (matching) tolerance.  PRIMARY 150 ms = the ANSI/AAMI EC57
# beat-matching window; the pipeline default ANNOTATION_MATCH_TOLERANCE_MS = 75 ms
# leaves ~20 % of V beats unmatched because the refined Pan-Tompkins fiducial of a
# wide QRS lies ~90-110 ms after the annotation (QC.md Section 2).
TOL_PRIMARY_MS = 150.0
TOL_DEFAULT_MS = float(fp.ANNOTATION_MATCH_TOLERANCE_MS)

SUBJECT = {r: ("S201" if r == "202" else f"S{r}") for r in RECORDS}

# Effective label of each annotation symbol in the pipeline (fp._annotation_label:
# NORMAL is checked first, then ABNORMAL; EXCLUDED_BEAT_SYMBOLS is never consulted,
# so "Q" (in both ABNORMAL and EXCLUDED) is labelled ABNORMAL).


def effective_label(symbol):
    lab = fp._annotation_label(symbol)
    return -1 if lab is None else int(lab)


def load(record):
    sig, fs, ann_s, ann_y = fp.load_mitbih_record_with_annotations(record, dl_dir=DATA_DIR)
    return sig, fs, np.asarray(ann_s, dtype=int), np.asarray(ann_y, dtype=str)


def detect(record, cache_dir=None):
    """Pan-Tompkins R peaks (fp.detect_r_peaks on channel 0 with CFG), cached as .npy."""
    path = None if cache_dir is None else os.path.join(cache_dir, f"{record}.npy")
    if path is not None and os.path.exists(path):
        return np.load(path)
    sig, fs, _, _ = load(record)
    peaks = np.asarray(fp.detect_r_peaks(sig, fs, fp.CFG), dtype=np.int64)
    if path is not None:
        os.makedirs(cache_dir, exist_ok=True)
        np.save(path, peaks)
    return peaks


def greedy_match(detected, reference, fs, tolerance_ms=fp.ANNOTATION_MATCH_TOLERANCE_MS):
    """The exact matching rule of fp.match_detected_peaks_to_annotations (and of
    fp.validate_r_peaks_against_annotations), returning the matched reference
    index per detected peak (-1 = unmatched).  Used only for bookkeeping; the
    labels themselves always come from the fp function (tests check agreement)."""
    detected = np.asarray(detected, dtype=int)
    reference = np.asarray(reference, dtype=int)
    tol = int(round(float(tolerance_ms) * float(fs) / 1000.0))
    used = np.zeros(len(reference), dtype=bool)
    match = np.full(len(detected), -1, dtype=int)
    for i, p in enumerate(detected):
        cand = np.flatnonzero((np.abs(reference - p) <= tol) & (~used))
        if len(cand):
            j = cand[np.argmin(np.abs(reference[cand] - p))]
            used[j] = True
            match[i] = j
    return match


def nn_touched(labels):
    """Interval j (ending at beat j+1) is replaced by the Phase 6 edit iff it ends at
    an abnormal beat (labels[j+1] == 1) or follows one (labels[j] == 1).
    labels: per-beat effective labels (length n_rr + 1)."""
    labels = np.asarray(labels, dtype=int)
    return (labels[1:] == 1) | (labels[:-1] == 1)


def hrv(rr):
    """Time-domain HRV of the window exactly as passed to the detector (seconds in, ms out)."""
    rr = np.asarray(rr, dtype=float)
    d = np.diff(rr)
    return {"sdnn_ms": float(np.std(rr, ddof=1) * 1000.0),
            "rmssd_ms": float(np.sqrt(np.mean(d ** 2)) * 1000.0),
            "pnn50": float(np.mean(np.abs(d) > 0.050) * 100.0)}


def annotation_beat_series(ann_s, ann_y):
    """Annotation-time beat series for the sensitivity analysis: beats = WFDB beat
    annotations; per-beat label = effective_label; a beat whose preceding interval
    contains a VF/flutter annotation ("!", "[", "]") gets -1 (that interval spans
    an unannotated flutter episode)."""
    is_beat = np.isin(ann_y, list(BEAT_SYMBOLS))
    bs = ann_s[is_beat]
    by = ann_y[is_beat]
    labels = np.array([effective_label(y) for y in by], dtype=int)
    vf = ann_s[np.isin(ann_y, list(VF_SYMBOLS))]
    if len(vf) and len(bs) > 1:
        k = np.searchsorted(bs, vf, side="right")       # vf lies between beat k-1 and beat k
        k = k[(k >= 1) & (k < len(bs))]
        labels[k] = -1
    return bs, by, labels


def record_series(record, cache_dir=None, source="detected", tolerance_ms=TOL_PRIMARY_MS):
    """Per-record beat series used for windows.

    source="detected" (PRIMARY): Pan-Tompkins peaks; labels from
        fp.match_detected_peaks_to_annotations against the BEAT annotations only,
        at tolerance_ms (PRIMARY 150 ms; PREREGISTRATION.md B).
    source="annotation" (sensitivity): annotation beat times (annotation_beat_series).
    Returns dict(beats, fs, rr, labels, ...)."""
    sig, fs, ann_s, ann_y = load(record)
    is_beat = np.isin(ann_y, list(BEAT_SYMBOLS))
    if source == "detected":
        beats = detect(record, cache_dir)
        labels = fp.match_detected_peaks_to_annotations(beats, ann_s[is_beat], ann_y[is_beat], fs,
                                                        tolerance_ms=tolerance_ms)
    elif source == "annotation":
        beats, _, labels = annotation_beat_series(ann_s, ann_y)
    else:
        raise ValueError(source)
    rr, _, _ = fp.extract_rr_intervals(beats, fs, fp.CFG.rr_min_seconds, fp.CFG.rr_max_seconds)
    return {"record": record, "subject": SUBJECT[record], "fs": fs, "beats": np.asarray(beats),
            "rr": rr, "labels": np.asarray(labels, dtype=int), "ann_s": ann_s, "ann_y": ann_y,
            "is_beat": is_beat, "source": source}


def windows(series):
    """fp.make_rr_windows with the pipeline defaults, plus per-window bookkeeping:
    window index, beat-sample span, NN mask and fraction, edited series, HRV."""
    out = []
    labels = series["labels"]
    touched = nn_touched(labels)
    for w in fp.make_rr_windows(series["rr"], labels):
        s, e = w["start_rr"], w["end_rr"]
        mask = touched[s:e]
        nn_frac = float(1.0 - mask.mean())
        edited = None
        if (~mask).any():
            edited, _ = edit_nn(w["rr"], mask)
        out.append({**w, "window": s // fp.CLASSIFIER_STEP_RR, "record": series["record"],
                    "subject": series["subject"], "t0": int(series["beats"][s]),
                    "t1": int(series["beats"][e]), "nn_mask": mask, "nn_fraction": nn_frac,
                    "rr_edited": edited, "n_abnormal": int(np.sum(labels[s + 1:e + 1] == 1))})
    return out
