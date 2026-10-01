# Phase 7: The Frozen Combined Chaos Detector on MIT-BIH

Research-development record in the style of the Phase 5 and Phase 6 reports. Every
statement is tied to a results file, table, preregistration or test under
`experiments/phase7_mitbih/`. Paths are relative to that directory unless they
start with `docs/`, `tests/` or `final_pipeline.py`.

This is the first analysis of real data. It is an **evaluation**: the detector,
its parameters and every analysis choice were fixed in `PREREGISTRATION.md`
before any detector output on MIT-BIH existed. Nothing was changed in response
to the results. Everything chosen afterwards is labelled **EXPLORATORY**.

## 1. Provenance and design

| Item | Value |
|---|---|
| Branch | `claude/amazing-cray-9jbx6p`, from `main` at `1b4b7a0` (Phases 2E–6) |
| Environment | Python 3.13.14, numpy 2.1.3, scipy 1.18.1, scikit-learn 1.9.1, wfdb 4.3.1; BLAS 1 thread, 4 workers; `results/run/manifest.jsonl` |
| Data | MIT-BIH Arrhythmia Database v1.0.0, all 48 records (`MITDB_RECORDS`) with `atr` annotations, from physionet.org. 144 files, every SHA-256 equal to PhysioNet's `SHA256SUMS.txt` (`DATA_MANIFEST.json`). Raw data not committed (`mitdb_data/` is gitignored) |
| Detector | `combined_chaos_detected(analyze_segment(rr, combined_chaos_config(replace(CFG, keep_upo_on_short_lle_embedding=True), detrend=True)))`, i.e. Phase 6 Section 5, unchanged. `final_pipeline.py` was not modified in Phase 7 |
| Windows | 256 RR intervals, step 256, abnormal iff ≥ 10 % abnormal beats (`fp.make_rr_windows`, unchanged) |
| Unit of analysis | subject (47 subjects; records 201 and 202 are one subject) for every grouping, split and bootstrap |
| Discipline | QC (`QC.md`, commit `2abd99b`) used no detector output. Code frozen in `573eefd`. `PREREGISTRATION.md` committed and pushed alone in `fb212da` before any MIT-BIH window was analysed. `run_phase7.py --phase test` reuses the Phase 3 guard |

**Summary.**
- On raw RR, the combined AND detector fired in **0 of 305** windows: 0/94
  abnormal and 0/211 normal. Q1 is not supported and Q2 is not met.
- With annotation beat times instead of detected peaks (sensitivity analysis), it
  fired in 6/120 abnormal and 0/248 normal windows. That is a difference of 5.0 %
  [1.9, 8.3], but the abnormal rate does not exceed the ectopy-only expectation.
- The LLE component alone separates the labels (50 % vs 13 %), as Phases 5–6
  predicted from ectopy alone. The UPO score carries no label information on raw RR.
- Standard HRV separates the labels better than either chaos score.

## 2. Data and quality control (`QC.md`, `results/qc/`; no detector output used)

### 2.1 R-peak quality (Pan-Tompkins, `fp.detect_r_peaks`, channel 0)

Reference: the 19 WFDB beat codes only.

| tolerance | TP | FP | FN | sensitivity | PPV |
|---|---:|---:|---:|---:|---:|
| 75 ms (requested) | 107,157 | 2,931 | 2,337 | 0.9787 | 0.9734 |
| 150 ms (ANSI/AAMI EC57 match window) | 109,112 | 976 | 382 | 0.9965 | 0.9911 |

**Most 75 ms errors are a late fiducial on wide QRS complexes, not missed beats:**
- 20.7 % of V beats and 7.0 % of paced beats are unmatched at 75 ms;
- in 98–99 % of these cases a detection lies 75–150 ms away, about 100 ms
  *after* the annotation;
- the ±80 ms raw-ECG refinement (Phase 1 Section 4) picks a late deflection of the
  wide complex.

The worst records at 75 ms:

| record | sensitivity | PPV |
|---|---:|---:|
| 208 | 0.68 | 0.68 |
| 107 (paced) | 0.79 | 0.79 |
| 119 | 0.90 | 0.90 |
| 207 | 0.97 | 0.80 (flutter episode) |
| 108 | 0.97 | 0.85 (noise) |

### 2.2 How the reused code handles R-peak errors (tests in `tests/test_phase7_mitbih_helpers.py`)

| case | what happens |
|---|---|
| **Extra detection** | label −1 drops every window in which it is an ending beat |
| **Missed beat** | leaves one merged interval, labelled by the next detected beat; the window is **kept** and the missed beat's label is lost |
| **Window's starting beat** | never label-checked |
| **`Q`** | in both `ABNORMAL_BEAT_SYMBOLS` and `EXCLUDED_BEAT_SYMBOLS`. `_annotation_label` never consults the excluded set, so `Q` is **abnormal** (2 windows contain Q; no label depends on it) |
| **Non-beat annotations** | can take a detection's match when all annotations are passed (as `analyze_mitbih_record_windows` does). This changes 0 windows here; Phase 7 passes beat annotations only |
| **Paced records** (102, 104, 107, 217) | all 32 candidate windows drop |

**Out of 403 candidate windows:**

| | 75 ms | 150 ms |
|---|---:|---:|
| windows with an extra detection | 174 | 69 |
| windows with a missed beat | 179 | 85 |
| retained windows | 224 (187 N / 37 A) | **305 (211 N / 94 A)** |

At 75 ms the late PVC fiducial alone would have removed 57 of 94 abnormal windows.
Label transfer therefore uses **150 ms** (decided in QC, before any detector run).

### 2.3 Windows, subjects, editing eligibility

| arm | windows | normal | abnormal | subjects | with abnormal |
|---|---:|---:|---:|---:|---:|
| **raw** (PRIMARY) | 305 | 211 | 94 | 43 | 19 |
| edited (NN fraction ≥ 0.80) | 218 | 211 | 7 | 38 | 7 |
| annotation beat times (SENSITIVITY) | 368 | 248 | 120 | 43 | 19 |

**Abnormal windows.**
- Abnormal fraction: median 0.24, range 0.10–0.80.
- 83 are mainly ventricular (16 subjects) and 11 mainly supraventricular (4
  subjects).
- Six subjects hold 54 of the 94 abnormal windows.
- 39 primary windows contain a missed beat (18 normal, 21 abnormal).

**Atrial fibrillation.** 34 windows overlap atrial fibrillation or flutter: 21
labelled normal and 13 abnormal. AF beats are annotated `N`.

**NN fraction.**
- Every normal window has an NN fraction ≥ 0.80.
- Only 7 abnormal windows do. The rest range from 0.02 to 0.80 (median 0.56).
- Each isolated ectopic beat removes two intervals, so the 10 % abnormal threshold
  and the 0.80 editing threshold are nearly incompatible.

**HRV on raw RR** (median [IQR]):

| label | SDNN ms | RMSSD ms | pNN50 % |
|---|---|---|---|
| normal | 52.7 [35.8, 88.6] | 64.6 [30.7, 107.3] | 12.5 [4.9, 35.9] |
| abnormal | 157.6 [109.7, 223.2] | 261.5 [148.2, 349.9] | 59.8 [46.3, 76.8] |

## 3. Preregistration summary (`PREREGISTRATION.md`, `fb212da`)

- **Detector:** as in Section 1.
- **Continuous scores**, fixed before the run:
  - **`lle_z`** = (Phase 3 statistic − surrogate mean) / surrogate SD, from
    `lle_chaos_test`'s 99 IAAFT surrogates;
  - **`upo_score`** = the largest per-peak rJ = deviation / W0 among period-1
    Level-A peaks passing the Phase 4 instability gate (robust median leading
    modulus ≥ 1.2), floored at 0.
- **Arms:** `raw` (primary); `edited` (NN fraction ≥ 0.80; Phase 6 `edit_nn`
  unchanged); `annotation` (sensitivity).
- **Q1:** Δ = AND rate (abnormal) − AND rate (normal). Supported iff the lower bound
  of its 95 % subject-cluster bootstrap CI (10,000 resamples) is > 0.
- **Q2:** "exceeds" iff the lower bound of the abnormal-rate CI is > 6.4 %.
- **Q3:** leave-one-subject-out L2 logistic regression, C = 1.0, scaler fit per
  fold. Models M0 HRV, M1 LLE, M2 UPO, M3 LLE + UPO, M4 all. Pooled out-of-fold
  AUC with cluster-bootstrap CIs and paired differences; DeLong secondary.
- **Amendments: none.** No code, window list or rule changed after the
  preregistration.
- **Run.** 891 windows in 630 s wall on 4 workers (preregistered estimate about 11
  min), in one pass. There were **0 errors** and 0 duplicate records. The UPO score
  reproduced the UPO decision in 891/891 windows (`upo_score_consistent`).

## 4. Q1 (PRIMARY): AND detection in abnormal vs normal windows, raw RR

Source: `results/analysis/results.md`, `results.json`.

| quantity | estimate | 95 % subject-cluster bootstrap CI | Wilson 95 % (ignores clustering) |
|---|---|---|---|
| AND rate, abnormal | **0/94 = 0.0 %** | [0.0, 0.0] | 0.0–3.9 % |
| AND rate, normal | **0/211 = 0.0 %** | [0.0, 0.0] | 0.0–1.8 % |
| **difference Δ** | **0.0 %** | **[0.0, 0.0]** | |
| sensitivity | 0.0 % | [0.0, 0.0] | |
| specificity | 100.0 % | [100.0, 100.0] | |
| PPV | undefined (no detection) | undefined in 10,000/10,000 resamples | |
| NPV | 69.2 % | [56.2, 81.8] | |
| balanced accuracy | 50.0 % | [50.0, 50.0] | |

Windows: 305 (211 N / 94 A), 43 subjects (19 with abnormal windows). TP 0, FN 94,
FP 0, TN 211.

**Q1 is not supported.** The detector made no detection at all, so the
preregistered cluster bootstrap is degenerate: every resample gives 0. **A
zero-width interval is not evidence of precision.** The informative bounds are the
Wilson upper limits:
- at most 3.9 % in abnormal windows;
- at most 1.8 % in normal windows.

The 0/211 on normal windows agrees with the synthetic specificity (Phases 5–6: at
most 2/400 on every null).

The descriptive 75 ms relabelled subset (224 windows, 37 abnormal) also has 0
detections. So do the 266 windows without a missed beat.

## 5. Q2: abnormal-window rate vs the ectopy-only expectation

- Rate: **0.0 %**, 95 % cluster-bootstrap CI [0.0, 0.0] (Wilson 0.0–3.9 %).
- The rule needs a lower bound above 6.4 %. **The rate does not exceed the Phase 6
  ectopy-only expectation.**
- It is in fact *below* the Phase 6 ectopy-only false-positive rates:
  - up to 3.7 % for isolated, couplet and atrial ectopy;
  - 0 % for bigeminy, trigeminy and runs.

The real abnormal windows produced fewer AND detections than synthetic linear RR
with added ectopy.

## 6. Q3: do UPO and the chaos measures add information? (raw RR)

### 6.1 Preregistered result

305 complete-case windows (0 excluded).

| model | features | out-of-fold AUC [95 % subject-cluster bootstrap] | DeLong 95 % CI (secondary) |
|---|---|---|---|
| M0 | SDNN, RMSSD, pNN50 | 0.812 [0.618, 0.949] | [0.746, 0.878] |
| M1 | `lle_z` | 0.685 [0.565, 0.797] | [0.618, 0.752] |
| M2 | `upo_score` | **0.193** [0.115, 0.290] | [0.125, 0.261] |
| M3 | `lle_z` + `upo_score` | 0.675 [0.553, 0.791] | [0.607, 0.743] |
| M4 | HRV + `lle_z` + `upo_score` | 0.821 [0.633, 0.951] | [0.758, 0.885] |

| difference | ΔAUC [95 % paired subject-cluster bootstrap] | DeLong z, p (secondary; windows not independent) |
|---|---|---|
| M3 − M1 (UPO beyond LLE) | −0.010 [−0.022, +0.002] | −2.29, 0.022 |
| M3 − M2 (LLE beyond UPO) | **+0.482 [+0.360, +0.605]** | 11.0, < 1e-27 |
| M4 − M0 (chaos beyond HRV) | +0.009 [−0.010, +0.028] | 1.37, 0.17 |

By the preregistered interpretation rules:

| claim | result |
|---|---|
| UPO adds beyond LLE | **no** (CI includes 0) |
| LLE adds beyond UPO | yes |
| Chaos measures add beyond HRV | **no** (CI includes 0) |

The secondary DeLong test calls M3 − M1 significantly *negative* (p = 0.022). It
ignores the clustering, and the subject-cluster CI includes 0.

### 6.2 Caveat: the pooled leave-one-subject-out AUC is biased here (EXPLORATORY diagnosis, Section 10)

An AUC of 0.19 for a single score is not an anti-signal. It is an artifact of
pooling out-of-fold predictions when labels cluster within subjects:
- many subjects are all-normal or nearly all-abnormal;
- the intercept fitted without a held-out subject is the training prevalence, which
  is strongly *anti-correlated* with that subject's own abnormal share (r = −0.92);
- an intercept-only model therefore scores a pooled AUC of **0.03**, and a pure
  noise feature a median of 0.15.

Weak features (M2) are pulled far below 0.5. Strong features (M0, M4) are pulled
down less. This bias was not anticipated in the preregistration.

**Fit-free AUCs** need no model and so have no fold effect (EXPLORATORY):

| score | AUC | 95 % subject-cluster bootstrap |
|---|---:|---|
| `lle_z` | 0.727 | [0.615, 0.825] |
| `upo_score` | 0.487 | [0.411, 0.563] |
| RMSSD | 0.865 | [0.729, 0.961] |
| SDNN | 0.821 | [0.661, 0.944] |
| pNN50 | 0.824 | [0.696, 0.925] |

**The qualitative conclusions do not change:**
- the UPO score carries no label information on raw RR;
- `lle_z` carries moderate information;
- HRV carries more than either chaos score.

## 7. Components alone (preregistered secondary, raw RR)

Detections / windows (Wilson 95 %), with the subject-cluster CI of the difference.

| decision | abnormal | normal | difference [cluster CI] | balanced accuracy |
|---|---|---|---|---|
| AND | 0/94 | 0/211 | 0.0 % | 0.50 |
| LLE alone | **47/94 (50.0 %; 40.1–59.9)** | 27/211 (12.8 %; 8.9–18.0) | +37.2 % [19.4, 54.8] | 0.69 [0.60, 0.77] |
| UPO alone | 3/94 (3.2 %; 1.1–9.0) | 4/211 (1.9 %; 0.7–4.8) | +1.3 % [−2.8, 6.5] | 0.51 [0.49, 0.53] |
| OR | 50/94 (53.2 %) | 31/211 (14.7 %) | +38.5 % [19.4, 57.0] | 0.69 [0.60, 0.79] |

- **The LLE component behaves as Phase 5 Section 8 predicted for ectopy.**
  - It fires in half the abnormal windows.
  - Its normal-window rate (12.8 %) is close to its synthetic null rates (7–18 %).
  - Its p-value sits at the minimum 0.01 in 24.5 % of abnormal vs 3.8 % of normal
    windows.
  - On synthetic ectopy it fired in 80–100 % of windows. Real abnormal windows
    trigger it less often, which fits Phase 6's finding that regular patterns
    (bigeminy, trigeminy) trigger it less.
- **The UPO gate passed in only 7 of 305 windows.** It never coincided with an LLE
  detection, so the AND fired nowhere.

## 8. Secondary results

### 8.1 Edited NN (NN fraction ≥ 0.80; descriptive, 7 abnormal windows)

- **Q1:** AND 0/7 abnormal and 0/211 normal; Δ = 0 (degenerate CI). Wilson upper
  bounds are 35 % and 1.8 %.
- **Components:** LLE 1/7 vs 17/211 (8.1 %); UPO 1/7 vs 5/211 (2.4 %).
- **The LLE component's normal-window rate falls from 12.8 % (raw) to 8.1 %
  (edited).** This is consistent with Phase 6: editing removes the ectopy that
  triggers it.
- **Q3 (preregistered):**

| model | AUC |
|---|---:|
| M0 | 0.230 |
| M1 | 0.114 |
| M2 | 0.579 |
| M3 | 0.376 |
| M4 | 0.438 |

  The CIs span about 0.1–0.8. With 7 abnormal windows from 7 subjects, and the
  same pooling artifact (intercept-only 0.04), these numbers carry no information.

### 8.2 Raw vs edited comparison

| window set | AND (abnormal) | AND (normal) | LLE alone (normal) | UPO alone (normal) |
|---|---|---|---|---|
| raw, all primary | 0/94 | 0/211 | 27/211 | 4/211 |
| edited, NN ≥ 0.80 | 0/7 | 0/211 | 17/211 | 5/211 |

### 8.3 By beat type and atrial fibrillation (raw RR, descriptive)

**Beat type:**
- mainly ventricular: AND 0/83, LLE 39/83 (47 %), UPO 3/83;
- mainly supraventricular: AND 0/11, LLE 8/11 (73 %), UPO 0/11.

**Records with AF / AFL annotations:**

| label | AND | LLE | UPO |
|---|---|---|---|
| normal | 0/25 | 2/25 | 2/25 |
| abnormal | 0/16 | 8/16 | 1/16 |

**Windows overlapping an AF / AFL episode:**

| label | AND | LLE | UPO |
|---|---|---|---|
| normal | 0/21 | 2/21 | 2/21 |
| abnormal | 0/13 | 7/13 | 1/13 |

No AF window was detected. The 21 AF-overlapping normal-labelled windows show the
LLE component at 9.5 %, no higher than in sinus windows (13.2 %). That fits AF's
irregularity being well described by the IAAFT (linear-Gaussian) null, but 21
windows from 5 subjects is too few to conclude anything.

### 8.4 Per record

The per-record table (`results/analysis/raw_per_record.csv`) shows no record with
an AND detection. The LLE component's abnormal-window detections are concentrated
in:

| record | LLE detections (abnormal windows) |
|---|---|
| 214 | 8/8 |
| 233 | 8/12 |
| 119 | 5/7 |
| 208 | 5/6 |
| 221 | 5/8 |
| 232 | 4/5 |

### 8.5 Missed beats and 75 ms labels (descriptive)

**Windows with a missed beat:**
- AND 0/39;
- LLE 3/18 normal and 7/21 abnormal, vs 24/193 and 40/73 without one.

Merged intervals do not inflate the LLE component; abnormal windows with a missed
beat were *less* often LLE-positive.

**75 ms labels:** AND 0/37 abnormal and 0/187 normal.

## 9. Sensitivity analysis: annotation beat times instead of detected peaks

368 windows (248 N / 120 A), 43 subjects (19 with abnormal windows), 0 errors.

| quantity | raw (primary) | annotation times |
|---|---|---|
| AND abnormal | 0/94 = 0.0 % | **6/120 = 5.0 %** [cluster 1.9, 8.3; Wilson 2.3–10.5] |
| AND normal | 0/211 = 0.0 % | 0/248 = 0.0 % [Wilson 0–1.5] |
| Δ (Q1) | 0.0 % [0, 0]: **not supported** | **+5.0 % [+1.9, +8.3]: supported by the Q1 rule** |
| Q2 (lower bound > 6.4 %) | no | **no** (lower bound 1.9 %) |
| LLE alone, A / N | 50.0 % / 12.8 % | 46.7 % / 13.7 % |
| UPO alone, A / N | 3.2 % / 1.9 % | **8.3 % / 1.2 %** |
| Q3 out-of-fold AUC M0 / M1 / M2 / M3 / M4 | 0.812 / 0.685 / 0.193 / 0.675 / 0.821 | 0.839 / 0.693 / 0.277 / 0.694 / 0.846 |
| ΔAUC M3 − M1 / M3 − M2 / M4 − M0 | −0.010 / +0.482 / +0.009 | +0.001 [−0.020, 0.022] / +0.417 [0.283, 0.537] / +0.007 [−0.005, 0.020] |

**R-peak errors matter for the UPO component, and therefore for the AND.**
- With exact annotation timing, the UPO gate fires in 10 of 120 abnormal windows
  instead of 3 of 94. The AND then fires in 6 abnormal windows from 6 subjects
  (EXPLORATORY E2):
  - 106, 119, 208 and 228 (mainly ventricular; abnormal fraction 0.20–0.41);
  - 232 (supraventricular, 73 % abnormal);
  - 233.
- Pairing windows by time (EXPLORATORY E3; 296 pairs at ≥ 90 % overlap):
  - the 4 paired AND detections are all **negative on raw RR**;
  - on raw RR, 3 of the 4 had LLE positive and UPO negative, with the UPO score
    0–2.0 on raw RR vs 2.2–3.4 with annotation timing;
  - the fourth (233 window 8) was UPO-positive and LLE-negative on raw RR.
- The late, variable PVC fiducial (Section 2.1) adds about 100 ms of jitter to
  exactly the premature and compensatory intervals. That jitter smears the
  So-transform structure that the UPO test needs.

**What the annotation-arm detections are (and are not).** Even there, the abnormal
AND rate (5.0 %) lies within the range Phase 6 measured for ectopy alone on
linear-Gaussian RR. Its lower bound (1.9 %) is far below the preregistered 6.4 %.
- Q1's rule is met in this arm: the rate is higher in abnormal windows.
- But by Q2 it is no more than ectopy alone could produce. These detections are
  therefore not evidence of deterministic chaos.
- Q1 in this arm is a sensitivity analysis, not the primary test.

## 10. Exploratory results (chosen after seeing the outputs; `exploratory.py`, `results/exploratory/`)

**E1. Pooled leave-one-subject-out AUC artifact** (Section 6.2).

| arm | intercept-only AUC | noise-feature AUC, median (range; 20 draws) | r(held-out share, training prevalence) |
|---|---:|---|---:|
| raw | 0.031 | 0.15 (0.03–0.41) | −0.92 |
| annotation | 0.027 | 0.11 | −0.97 |
| edited | 0.040 | 0.32 | −0.74 |

Fit-free AUCs, raw arm:

| score | AUC [95 % subject-cluster bootstrap] |
|---|---|
| `lle_z` | 0.727 [0.615, 0.825] |
| `lle_stat` | 0.720 [0.609, 0.814] |
| `upo_score` | 0.487 [0.411, 0.563] |
| `upo_source_rJ` | 0.591 [0.496, 0.679] |
| SDNN | 0.821 |
| RMSSD | 0.865 |
| pNN50 | 0.824 |

In the annotation arm: `lle_z` 0.726, `upo_score` 0.521.

**UPO score distribution.** `upo_score` is exactly 0 (no gate-passing peak) in
65 % of normal and 70 % of abnormal raw windows. Its median is 0 in both.

**E2. The 6 annotation-time detections** are listed in Section 9:
- all have LLE p ≤ 0.05 and one gated peak (modulus 1.26–5.95);
- none was detrended;
- none overlaps AF.

Of the 4 with a time-paired raw window:
- 3 were LLE-positive / UPO-negative on raw RR;
- 1 (233 window 8) was UPO-positive / LLE-negative.

**E3. Raw vs annotation pairs** (296; labels agree in 99.7 %):

| decision | both | annotation only | raw only |
|---|---:|---:|---:|
| LLE | 46 | 25 | 21 |
| UPO | 4 | 7 | 3 |
| AND | 0 | 4 | 0 |

Correlations between the arms: `lle_z` r = 0.71, `upo_score` r = 0.50. **The UPO
score is the component most sensitive to R-peak timing.**

## 11. Limitations

- **Few windows and subjects.**
  - 305 primary windows from 43 subjects, but only **19 subjects** contribute
    abnormal windows. Six of them contribute 54 of the 94.
  - With 0 detections, the preregistered cluster bootstrap is degenerate.
  - Only the Wilson bounds (≤ 3.9 % abnormal, ≤ 1.8 % normal) describe the
    uncertainty, and they ignore clustering.
  - The edited arm has 7 abnormal windows.
- **MIT-BIH was selected to contain arrhythmias.** Records 100–124 are a random
  sample and 200–234 were chosen for rare arrhythmias, from a 1975–1979 Boston
  population. The normal windows are mostly from arrhythmia patients, not healthy
  subjects. Rates do not transfer to screening populations.
- **The label is beat-based, not rhythm-based.**
  - "Abnormal" means ≥ 10 % abnormal beats.
  - AF beats are annotated `N`, so AF windows are mostly labelled normal.
  - Bundle-branch-block beats (L, R) count as normal.
  - Paced records are excluded entirely.
  - The question "is the arrhythmic rhythm chaotic?" is therefore only partly
    addressed.
- **R-peak detection errors.**
  - The Pan-Tompkins fiducial lands about 100 ms late on about 20 % of PVCs, which
    forced the 150 ms label tolerance.
  - 39 windows contain merged intervals from missed beats.
  - Section 9 shows these errors change the UPO component and the AND. **The
    primary null result is partly a consequence of R-peak timing error**, and
    annotation-timed RR gives a different (still sub-threshold) answer.
  - A better fiducial (for example a matched-filter or annotation-guided
    refinement) was not tested.
- **The edited-NN analysis is underpowered by construction.** The 10 % abnormal
  threshold and the 0.80 NN threshold leave almost no abnormal windows. Phase 6's
  limitation also still holds: power under editing was never measured.
- **Blindness to continuous-flow chaos (Phase 5 Section 8.6).** The LLE test missed
  sampled Rössler and Mackey–Glass flows at every m. A negative AND therefore means
  "no evidence of low-dimensional, map-like chaos detectable in 256 beats". It does
  not mean "no chaos". Heart-rate dynamics, if chaotic, may well be flow-like.
- **What a detection would and would not mean (Phase 5 Section 8).**
  - The IAAFT null is linear-Gaussian, so an LLE detection means "not a monotone
    transform of a linear Gaussian process". Ectopy, non-chaotic nonlinearity
    (SETAR 13 %) and static warping (18 %) all trigger it.
  - The UPO gate supplies specificity to unstable fixed points. Even so, AND
    detections at about 5 % in ectopic windows are within what ectopy alone
    produced in Phase 6.
  - No detection here can be read as evidence of deterministic chaos in the
    arrhythmic heart rhythm.
- **Nonstationarity costs power (Phases 5–6).** The D2 detrend was applied in 32 %
  of normal and 14 % of abnormal raw windows. A negative AND on a drifting window
  is weak evidence.
- **Q3 metric.** The preregistered pooled leave-one-subject-out AUC is biased
  downward when labels cluster within subjects (Section 6.2). The paired
  differences share the bias, so their direction and size must be read with that
  in mind. The exploratory fit-free AUCs are the cleaner description but were not
  preregistered.
- **Unvalidated fixed choices.** C = 1, untransformed HRV, the UPO score definition
  and m = 2 were fixed in advance without validation on RR data.
- **One machine, one run; no amendment.** The run was not interrupted.

## 12. Reproducibility

```bash
uv venv -p python3.13 /root/venv313 && uv pip install -p /root/venv313/bin/python -r requirements.txt wfdb
python - <<'PY'   # data (see DATA_MANIFEST.json; verify the SHA-256 values)
import wfdb, final_pipeline as fp
for r in fp.MITDB_RECORDS: wfdb.dl_database("mitdb", dl_dir="mitdb_data", records=[r], annotators=["atr"])
PY
python -m experiments.phase7_mitbih.qc && python -m experiments.phase7_mitbih.qc_tables   # Step 1
python -m experiments.phase7_mitbih.run_phase7 --phase synthetic                           # code check
python -m experiments.phase7_mitbih.run_phase7 --phase test --workers 4                    # guarded
python -m experiments.phase7_mitbih.analysis                                               # preregistered
python -m experiments.phase7_mitbih.exploratory                                            # EXPLORATORY
python -m pytest tests
```

Commits:
- `8f23a40` data manifest;
- `2abd99b` QC;
- `573eefd` code;
- `fb212da` preregistration;
- `e2e6964` run;
- `57d31bf` preregistered analyses;
- `39ba5d3` exploratory analyses.

`final_pipeline.py` is unchanged (SHA-256 prefix `27bf93808ab45487`).
