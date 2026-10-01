# Phase 7 Step 1: MIT-BIH data-quality checks

**No detector output was computed or viewed for this step.** Only the ECG, the
Pan-Tompkins R-peak detector (`fp.detect_r_peaks`, `CFG`, channel 0) and the
reference annotations were used. `analyze_segment`, `lle_chaos_test` and the UPO
analysis were not run on any MIT-BIH data.

- Code: `qc.py` (results), `qc_tables.py` (tables), `data.py` (thin layer over the
  reused `final_pipeline.py` helpers).
- Results: `results/qc/` (`records.csv`, `windows.csv`, `candidates.csv`,
  `windows_annotation_times.csv`, `symbols.csv`, `summary.json`, `tables.md`,
  `detected_peaks.npz`).
- Helper tests: `tests/test_phase7_mitbih_helpers.py` (43 tests, synthetic inputs).
- Data: 48 records of the MIT-BIH Arrhythmia Database v1.0.0. Every file matched
  PhysioNet's SHA-256 (`DATA_MANIFEST.json`).

## 1. R-peak quality (75 ms tolerance, as requested)

**Reference beats.** Only the 19 WFDB beat (QRS) codes `N L R B A a J S V r F e j n
E / f Q ?` are used as the reference. The non-beat annotations (`+ ~ | " ! [ ] x`) are
left out, as is standard for beat-by-beat evaluation.
- wfdb-python 4.3.1's own `is_qrs` table is misaligned with its label list. For
  example, it marks `[`, `]` and `x` as QRS and `e`, `f` as non-QRS. It was
  therefore not used.

**Counting.** `fp.validate_r_peaks_against_annotations` does the counting, with
greedy one-to-one matching.

| tolerance | TP | FP | FN | sensitivity | PPV |
|---|---:|---:|---:|---:|---:|
| **75 ms** (requested) | 107,157 | 2,931 | 2,337 | **0.9787** | **0.9734** |
| 150 ms (ANSI/AAMI EC57 match window) | 109,112 | 976 | 382 | 0.9965 | 0.9911 |
| 75 ms, 44 non-paced records | | | | 0.9821 | 0.9763 |
| 150 ms, 44 non-paced records | | | | 0.9964 | 0.9904 |

The per-record values are in `results/qc/tables.md`.

**Worst records at 75 ms.**

| record | sensitivity | PPV | cause |
|---|---:|---:|---|
| 208 | 0.681 | 0.683 | FP ≈ FN |
| 107 | 0.787 | 0.790 | paced; FP ≈ FN |
| 119 | 0.904 | 0.904 | FP ≈ FN |
| 207 | 0.973 | 0.798 | ventricular flutter episode with no beat annotations |
| 108 | 0.974 | 0.845 | noisy record, many extra detections |

### Most 75 ms "misses" are a timing offset, not missed beats

Records where FP ≈ FN at 75 ms (119, 124, 203, 208, 213, 221 and the paced records)
point to detections displaced from the annotation, not to missed beats.

- **Missed annotated beats by type (75 ms).**
  - V: 1,474 of 7,130, **20.7 %**;
  - paced (`/`): 490 of 7,028, 7.0 %;
  - N: 308 of 75,052, 0.4 %.
- **Where the unmatched detections are.** Among detections with no match at 75 ms:
  - 1,449 have a V beat as their nearest annotation. The distance is 89 / 103 /
    111 ms (10th / 50th / 90th percentile), and 98 % lie within (75, 150] ms;
  - 491 have a paced beat as their nearest annotation, with 99 % within
    (75, 150] ms.
- **Direction of the offset.** Among V beats matched at 150 ms, 19.9 % of the
  detected peaks lie more than 75 ms **after** the annotation; 0.2 % lie more than
  75 ms before it.

The refined Pan-Tompkins fiducial (±80 ms raw-ECG extremum, Phase 1 Section 4) lands
on a late deflection of a wide QRS complex. At 150 ms the totals rise to 0.9965 /
0.9911.

## 2. Missed and extra beats: how the reused code handles them

Code behaviour was verified by tests in `tests/test_phase7_mitbih_helpers.py`.

- **(a) Extra detections.** A detected peak with no annotation within tolerance gets
  label −1 from `match_detected_peaks_to_annotations`. `make_rr_windows` drops every
  window in which any of the 256 **ending** beats has label −1. Extra detections
  therefore never enter a window, except at a window's first beat (see below).
- **(b) Missed beats.** An annotated beat the detector missed leaves no trace in the
  labels. The two intervals around it merge into one long interval, labelled by the
  next detected beat, and **the window is kept**. The missed beat's own label is lost:
  a missed V beat is not counted as abnormal.
- **The window's starting beat is not checked.** Interval `rr[start]` begins at beat
  `start`, but `make_rr_windows` checks only beats `start+1 .. end`. An extra
  detection at the first beat therefore stays in the window.
- **Non-beat annotations can win a match.** `analyze_mitbih_record_windows` passes
  *all* annotations to the matcher. A detection that lies nearer a non-beat
  annotation (for example `+`) than its beat is then labelled −1 (test
  `test_non_beat_annotation_can_steal_a_match`).
  - On MIT-BIH this changes 0 beat labels at 75 ms and 5 at 150 ms, and **0
    windows** at either tolerance.
  - Phase 7 passes beat annotations only.

**Candidate windows** (every 256-interval position, step 256, before label filtering):

| | 75 ms | 150 ms (primary) |
|---|---:|---:|
| candidate windows | 403 | 403 |
| contain ≥ 1 extra detection | 174 | 69 |
| contain ≥ 1 missed annotated beat | 179 | 85 |
| contain ≥ 1 beat with an excluded label (paced `/`, unmatched to a usable symbol) | 32 | 32 |
| dropped for extra detections only / excluded symbols only / both | 147 / 5 / 27 | 66 / 29 / 3 |
| **retained** | **224** (187 N, 37 A) | **305** (211 N, 94 A) |
| retained windows containing ≥ 1 missed beat (merged interval) | 22 | **39** (18 N, 21 A; 111 missed beats) |
| retained windows with an extra detection at the starting beat | 2 | 1 |

At 75 ms, 81 of the 305 windows are lost, almost all through the late PVC
fiducial.
- 57 of them are abnormal windows, so 75 ms keeps only 37 of 94.
- 24 are normal windows.
- No window that is retained at both tolerances changes label (`windows.csv`,
  `label_tol75`).

The 75 ms default thus systematically removes windows with PVCs. That is a
selection bias against exactly the abnormal windows Q1 needs. **Decision (recorded
in PREREGISTRATION.md B, made before any detector run):**
- label transfer uses the **150 ms** EC57 window;
- the 75 ms arm is reported descriptively on the same windows. Its RR series are
  identical, so no extra detector runs are needed.

Missed beats inside retained windows are kept, as the reused code does. The
annotation-time sensitivity analysis measures whether such R-peak errors matter.

## 3. Labels

| symbol | count | effective label | note |
|---|---:|---|---|
| N | 75,052 | normal | |
| L, R | 8,075, 7,259 | normal | bundle-branch-block beats (AAMI class N) |
| e, j | 16, 229 | normal | atrial / nodal escape |
| V | 7,130 | abnormal | |
| A, a, J, S | 2,546, 150, 83, 2 | abnormal | supraventricular |
| E | 106 | abnormal | ventricular escape |
| F | 803 | abnormal | fusion of ventricular and normal |
| f | 982 | abnormal | fusion of paced and normal (paced records only) |
| **Q** | 33 | **abnormal** | see below |
| / | 7,028 | excluded (−1) | paced |
| B, r, n, ?, P | 0 | — | not present in mitdb |
| x | 193 | excluded (−1) | non-conducted P wave, not a beat |
| + ~ \| " ! [ ] | 1,291, 616, 132, 437, 472, 6, 6 | — | non-beat; not in the reference |

**"Q" is labelled abnormal.**
- `"Q"` is in both `ABNORMAL_BEAT_SYMBOLS` and `EXCLUDED_BEAT_SYMBOLS`.
- `_annotation_label` checks `NORMAL_BEAT_SYMBOLS`, then `ABNORMAL_BEAT_SYMBOLS`, and
  returns `None` otherwise. `EXCLUDED_BEAT_SYMBOLS` is **never consulted** anywhere in
  the code.
- So `"Q"` takes effect as abnormal (test `test_q_is_in_both_sets_and_abnormal_wins`).

**Impact of Q.** Only 2 retained windows contain a Q beat:
- 101 window 3: normal, 1 abnormal beat (the Q);
- 214 window 1: abnormal, 35 abnormal beats, 2 of them Q.

Excluding Q would drop the 101 window. Treating Q as normal would leave 214 window 1
abnormal (33/256 = 12.9 %). No label changes either way. **Decision:** keep the code's
effective label (abnormal).

**Paced records** (102, 104, 107, 217) contain 2,028, 1,380, 2,078 and 1,542 paced
beats. **All 32 of their candidate windows are dropped**, so they contribute no
window.

## 4. Subjects

- 48 records come from 47 subjects. Records 201 and 202 are one subject (`S201`).
  Every other record is its own subject (`data.SUBJECT`).
- **Subjects with retained primary windows: 43; with abnormal windows: 19.**
- The 4 subjects without windows are the paced records (102, 104, 107, 217).

## 5. Windows (256 RR intervals, step 256, abnormal iff ≥ 10 % abnormal beats)

| arm | windows | normal | abnormal | subjects | subjects with abnormal |
|---|---:|---:|---:|---:|---:|
| **PRIMARY**: detected peaks, 150 ms labels | **305** | **211** | **94** | 43 | 19 |
| detected peaks, 75 ms labels (descriptive) | 224 | 187 | 37 | | |
| annotation beat times (sensitivity analysis) | 368 | 248 | 120 | 43 | 19 |

**Abnormal windows.**
- Abnormal fraction: median 0.24, IQR 0.18–0.32, range 0.10–0.80.
- **Beat type:** 83 are mainly ventricular (V + E + F > A + a + J + S; 16 subjects)
  and 11 mainly supraventricular (4 subjects).

**Atrial fibrillation / flutter.**
- 34 windows overlap an AFIB or AFL rhythm episode: 21 labelled normal and 13
  abnormal.
- AF beats are annotated `N`, so the beat-based label calls AF windows "normal"
  unless they also contain enough ectopic beats.
- Records with AF/AFL annotations: 201, 202, 203, 210, 217, 219, 221, 222.

**RR extremes.** 19 windows contain an RR interval above 2 s, and 12 contain one
below 0.3 s. These are missed beats, pauses and noise; they are kept as the reused
code does.

The per-record counts are in `results/qc/tables.md`. **Abnormal windows per
subject:**
- S233 12, S213 10;
- S200, S214, S223, S221 8 each;
- S119 7, S201 6, S208 6, S232 5, S106 5, S209 3, S228 2;
- 1 each for S114, S203, S210, S205, S124 and S234.

The top six subjects hold 54 of the 94 abnormal windows. Subjects contribute 1–13
windows each (median 7).

## 6. Editing eligibility (NN fraction; Phase 6 Part B rule)

**Rule.** An interval is replaced if it ends at an abnormal beat **or** follows one
(`data.nn_touched`). It is replaced by linear interpolation using the Phase 6
function `experiments.phase6_robust.systems.edit_nn`, unchanged. The NN fraction is
the share of untouched intervals.

| NN fraction | [0, 0.2) | [0.2, 0.4) | [0.4, 0.6) | [0.6, 0.7) | [0.7, 0.8) | [0.8, 0.9) | [0.9, 0.95) | [0.95, 1] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| normal (211) | 0 | 0 | 0 | 0 | 0 | 26 | 34 | 151 |
| abnormal (94) | 11 | 12 | 30 | 16 | 18 | 7 | 0 | 0 |

**Windows with NN fraction ≥ 0.80:**
- **211 normal** (all);
- **7 abnormal**, from 7 subjects.

An abnormal window has at least 10 % abnormal beats, and each isolated abnormal beat
removes two intervals. So only windows near the 10 % threshold with paired
(couplet / run) ectopy reach 0.80. **The edited-NN secondary analysis therefore has
only 7 abnormal windows.** It can only be descriptive.

Annotation-time windows: 248 normal and 8 abnormal are eligible.

## 7. HRV baseline features (computed from the window passed to the detector)

SDNN (ddof = 1), RMSSD and pNN50 (|ΔRR| > 50 ms) are computed on the raw-RR window
(primary) and on the edited series (edited-NN arm).

| label | n | SDNN ms | RMSSD ms | pNN50 % | mean RR s |
|---|---:|---|---|---|---|
| normal | 211 | 52.7 [35.8, 88.6] | 64.6 [30.7, 107.3] | 12.5 [4.9, 35.9] | 0.772 [0.686, 0.922] |
| abnormal | 94 | 157.6 [109.7, 223.2] | 261.5 [148.2, 349.9] | 59.8 [46.3, 76.8] | 0.701 [0.592, 0.838] |
| normal, edited | 211 | 37.7 [25.8, 72.3] | 30.2 [20.9, 76.7] | 5.9 [1.2, 35.9] | |
| abnormal, edited (eligible) | 7 | 47.5 [21.8, 49.8] | 38.8 [17.3, 45.5] | 15.3 [0.2, 24.9] | |

Values are median [IQR]. On raw RR the ectopic beats inflate SDNN, RMSSD and pNN50,
as in Phase 5 S2.

## 8. Consequences for the preregistration

1. Label transfer at 150 ms (EC57), with beat annotations only. The 75 ms arm is
   descriptive.
2. The primary analysis has 305 windows (211 / 94) from 43 subjects, but only **19
   subjects with abnormal windows**. Cluster-bootstrap CIs will be wide.
3. The edited-NN arm has 7 abnormal windows and is descriptive only.
4. Missed beats are kept in primary windows (39 windows). The annotation-time
   sensitivity analysis and a descriptive split by missed-beat presence address them.
5. The paced records contribute nothing. Q stays abnormal (no label changes).
