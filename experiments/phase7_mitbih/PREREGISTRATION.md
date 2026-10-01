# Phase 7 preregistration: the frozen combined chaos detector on MIT-BIH

This file is committed and pushed **on its own, before the detector is run on any
MIT-BIH window**. Up to this commit, nothing derived from `analyze_segment`,
`lle_chaos_test` or the UPO analysis has been computed on MIT-BIH data.

**What was seen before this file.** Only the Step 1 data-quality results
(`QC.md`, `results/qc/`): ECG, R peaks, annotations, window counts, NN fractions
and HRV summaries. The detector itself ran only on Phase 5 synthetic windows (code
check, `results/synthetic_check/`).

**What is frozen with this file.** After this commit, the following do not
change:
- `final_pipeline.py`;
- `data.py`, `detector.py`, `run_phase7.py`, `analysis.py`;
- the window lists in `results/qc/windows*.csv`.

Any unavoidable change goes in a dated amendment at the end of this file, saying
why and whether any detector output had been seen. `run_phase7.py --phase test`
reuses the Phase 3 guard: it refuses to start unless this file is committed,
unmodified, pushed, and lists the method below.

**Code at freeze** (git `573eefd`; SHA-256 prefixes):

| file | SHA-256 prefix |
|---|---|
| `final_pipeline.py` | `27bf93808ab45487` |
| `data.py` | `35fa25908551f64e` |
| `detector.py` | `0e38fea524b31652` |
| `run_phase7.py` | `8240a41da04d2cf9` |
| `analysis.py` | `8df93b7b0cb99133` |

## A. Detector

- method: `phase6_recommended_and`

The configuration is `docs/PHASE6_ROBUSTNESS.md` Section 5, unchanged:

```python
DETECTOR = combined_chaos_config(replace(CFG, keep_upo_on_short_lle_embedding=True), detrend=True)
decision = combined_chaos_detected(analyze_segment(rr, DETECTOR))
```

This is the Phase 3 `lle_chaos_test` (τ = 1, m = 2, 99 IAAFT surrogates,
p ≤ 0.05) **AND** the Phase 4 C2 instability-gated UPO test (m = 2, M = 15,
50 surrogates, median gate at δ = 0.2). It adds the Phase 6 D2 detrend: linear,
applied when the trend change is ≥ 0.7 residual SD. `analyze_segment` receives the
window exactly as built in Section B, with `use_corrected_rr_for_dynamics = False`
as tested. No parameter is changed for MIT-BIH.

**Recorded per window** (`detector.evaluate`; one JSON line per window):
- the decisions `and_detected`, `lle_detected` (`lle_chaos_test["detected"]`),
  `upo_detected` (`upo["instability_gate_detected"]`) and OR;
- `lle_stat` (`lle_chaos_test["lle"]`), `lle_p`, the surrogate mean and SD, and the
  number of undefined surrogates;
- UPO status, the significance status, Level-A / Level-B / gated peak counts,
  `upo_max_gated_modulus` (largest robust leading modulus among gated Level-B
  peaks; the "largest gated multiplier"), `upo_source_rJ`, and per-peak rJ, J,
  deviation, leading modulus and gate;
- the detrend ratio and whether it was applied, `lle_embedding_too_short`, runtime;
- the error message if `analyze_segment` raised.

**Errors** count as "not detected" in Q1 and Q2 and are reported separately, by
label.

### Continuous scores (for Q3), fixed in advance

**LLE score: `lle_z = (S − mean(S_sur)) / SD(S_sur)`.**
- S is the Phase 3 statistic `lle_chaos_test["lle"]`. S_sur are its finite IAAFT
  surrogate values (`lle_chaos_test["surrogate_lle"]`), and SD uses ddof = 1. It is
  NaN when S is undefined.
- **Why not p:** the p-value cannot rank windows. With 99 surrogates it has 100
  levels and saturates at its minimum 0.01 on most ectopic windows (Phase 5 S2:
  97–99 % of windows at 0.01).
- **Why not S itself:** S is not scale-free across windows. Its null level
  depends on the noise and spectrum of each window (Phase 3: white noise 0.12,
  AR(1) 0.20).
- The surrogate z-score is continuous, standardized against each window's own
  null, and monotone in the quantity the decision thresholds.

**UPO score: `upo_score`.** It is defined from these quantities:

| quantity | what it is | source |
|---|---|---|
| deviation | the peak's signed So-histogram excess over the surrogate mean | period-1 So surrogate test (`upo.periods[1].significance.per_peak`, aligned with `upo.periods[1].source_peak_candidates`) |
| W0 | the median surrogate maximum deviation | same test |
| rJ | deviation / W0 | the existing per-peak field |
| gate | robust (element-wise median) source-stability leading modulus ≥ 1.2 | `fp.robust_source_stability` on `upo.periods[1].jacobians`, exactly the computation of `fp._attach_instability_gate` |

`upo_score` is the largest rJ among period-1 Level-A peaks that pass the gate,
floored at 0. It is 0 when no peak passes the gate, or when the UPO analysis ended
in a failure status or without significance.

**Why this score:**
- The UPO decision is "a gate-passing period-1 peak has J < 0.05". J = fraction of
  surrogate maxima above the deviation, a decreasing function of it, and W0 is
  fixed within a window. So within a window, the decision is a threshold on this
  score: a detection implies a score above 1.
- It therefore carries both parts of the decision: significance, and instability,
  which supplies the specificity (Phases 4–6).
- The alternatives lose information:
  - **J itself** has only 51 levels with 50 surrogates;
  - **`upo_source_rJ`** ignores the gate, and the gate is what removed the
    periodic false positives in Phase 4;
  - **the largest gated multiplier** is undefined (0) in every window without a
    significant peak, which would tie almost all normal windows.
- `upo_max_gated_modulus` and `upo_source_rJ` are recorded and may be shown only
  as EXPLORATORY.
- `detector.py` records whether the recomputed gate reproduces the UPO decision
  (`upo_score_consistent`). Inconsistent windows are counted and reported.

## B. Inputs and handling decisions (from Step 1, `QC.md`)

All windows use the existing settings, through `fp.make_rr_windows`:
- 256 RR intervals, step 256;
- abnormal iff ≥ 10 % of the 256 ending beats are abnormal;
- a window is dropped if any ending beat has label −1.

RR = `fp.extract_rr_intervals` of the beat times in seconds (no RR correction).
Subject (Section B.6) is the unit of every grouping, split and bootstrap.

1. **PRIMARY, arm `raw` (305 windows: 211 normal, 94 abnormal; 43 subjects, 19 with
   abnormal windows).**
   - Peaks: Pan-Tompkins (`fp.detect_r_peaks`, `CFG`, ECG channel 0).
   - Labels: `fp.match_detected_peaks_to_annotations` against the **beat
     annotations only** (WFDB beat codes `N L R B A a J S V r F e j n E / f Q ?`),
     with a **150 ms** tolerance.
   - **Why 150 ms and not the 75 ms default:**
     - at 75 ms, 20.7 % of V beats and 7 % of paced beats are unmatched. In 98–99 %
       of these cases a detection lies 75–150 ms away, almost always after the
       annotation. The refined fiducial of a wide QRS lands about 100 ms late
       (QC.md 1);
     - every such PVC becomes an "extra" detection and drops its window. 75 ms
       keeps only 37 of the 94 abnormal windows, a selection bias against
       abnormal windows;
     - 150 ms is the ANSI/AAMI EC57 beat-matching window. At 150 ms sensitivity is
       0.9965 and PPV 0.9911.
   - **Why beat annotations only:** non-beat annotations (`+ ~ | " ! [ ] x`) are
     not beats. Passing them can let a non-beat annotation take a detection's
     match. On MIT-BIH this changes 0 windows at either tolerance (QC.md 2).

2. **75 ms labels (descriptive only).** The same detector runs as `raw`, restricted
   to the 224 windows retained at 75 ms and relabelled with 75 ms labels. Labels
   agree wherever both exist. No extra runs are needed: the RR series do not
   depend on the tolerance.

3. **SECONDARY, arm `edited` (218 windows: 211 normal, 7 abnormal).**
   - These are the primary windows with **NN fraction ≥ 0.80** (Phase 6's tested
     range).
   - Every interval that ends at an abnormal beat, and the interval that follows
     it, is replaced by the Phase 6 Part B edit:
     `experiments.phase6_robust.systems.edit_nn`, unchanged (linear interpolation
     by index; nearest NN value at the edges; no re-quantization).
   - Labels are unchanged.
   - HRV features are computed from the edited series.
   - Windows with NN fraction < 0.80 are not analysed in this arm.
   - **With 7 abnormal windows (7 subjects) this arm can only be descriptive.**

4. **SENSITIVITY, arm `annotation` (368 windows: 248 normal, 120 abnormal; 43
   subjects, 19 with abnormal windows).**
   - RR comes from the **annotation beat times** (beat codes above). Labels come
     from the annotation symbols.
   - A beat whose preceding interval contains a VF / flutter annotation (`! [ ]`)
     is set to −1, because that interval spans an unannotated flutter episode.
   - The windows are re-cut from the annotation series by the same
     `make_rr_windows`, so they are not paired one-to-one with `raw` windows.
   - Comparing `raw` with `annotation` shows whether R-peak errors matter: timing
     jitter, the late PVC fiducial, missed and extra beats.

5. **Missed and extra beats (all arms use the reused code's behaviour).**
   - An extra detection (no beat annotation within 150 ms) gets label −1 and drops
     every window in which it is an ending beat.
   - A missed beat leaves one merged interval, and the window is **kept**, as in
     `make_rr_windows`. This affects 39 primary windows: 18 normal, 21 abnormal,
     111 missed beats.
   - The starting beat of a window is not label-checked, as in `make_rr_windows`.
     This affects 1 primary window.
   - Reason: the primary arm is what the pipeline produces from the ECG. Their
     effect is measured by arm 4 and by a descriptive split of `raw` by
     missed-beat presence.

6. **Labels and subjects.**
   - **"Q" is ABNORMAL.** It is the effective label in the code: `_annotation_label`
     checks `ABNORMAL_BEAT_SYMBOLS` and never consults `EXCLUDED_BEAT_SYMBOLS`.
     Two windows contain Q, and no window label depends on it (QC.md 3).
   - **Paced beats (`/`) are excluded (−1).** All 32 candidate windows of the paced
     records 102, 104, 107 and 217 drop.
   - `x` (non-conducted P wave) is not a beat.
   - **Subject = record**, except records 201 and 202, which are the same subject
     (`S201`).

7. **No other exclusion.** RR extremes (19 windows with an interval > 2 s, 12 with
   one < 0.3 s) are kept.

## C. Analyses (`analysis.py`)

**Bootstrap.** Every confidence interval is a 95 % percentile **subject-cluster
bootstrap** with **10,000 resamples**:
- subjects are drawn with replacement, and all windows of a drawn subject enter
  with its multiplicity;
- fixed seeds (20261001 + offsets);
- a resample with no abnormal or no normal window is discarded, and the number of
  valid resamples is reported;
- a metric undefined in a resample (for example PPV with no detection) is dropped
  from that metric's percentiles, and the count is reported.

### Q1 (PRIMARY, arm `raw`)

- **Estimand.** Δ = (AND detections / abnormal windows) − (AND detections / normal
  windows), pooled over windows.
- **Q1 is supported iff the lower bound of the 95 % CI of Δ is > 0.** A CI entirely
  below 0 is reported as the opposite direction.
- **Also reported, each with its cluster-bootstrap CI:** sensitivity (= abnormal
  rate), specificity (= 1 − normal rate), PPV, NPV and balanced accuracy of the AND
  decision as a classifier of the window label. Wilson intervals, which ignore
  clustering, are shown for reference.

### Q2 (arm `raw`)

- The abnormal-window AND rate with its 95 % cluster-bootstrap CI (the same
  resamples as Q1).
- It is compared with the Phase 6 ectopy-only false-positive rates: highest 3.7 %
  (E5 atrial 10 %, 11/300), upper 95 % bound 6.4 %.
- **"Exceeds the ectopy-only expectation" iff the CI's lower bound is > 6.4 %.**

### Q3 (arm `raw`)

**Method:**
- leave-one-subject-out logistic regression (43 folds);
- L2 penalty with **C = 1.0** fixed (sklearn `LogisticRegression`, lbfgs,
  max_iter 10,000, no class weights);
- `StandardScaler` fit on each training fold only;
- features used untransformed.

**Models:**

| model | features |
|---|---|
| M0 | SDNN, RMSSD, pNN50 (HRV baseline) |
| M1 | `lle_z` |
| M2 | `upo_score` |
| M3 | `lle_z` + `upo_score` |
| M4 | SDNN, RMSSD, pNN50 + `lle_z` + `upo_score` |

**Features.** HRV is computed on the window passed to the detector: SDNN with
ddof = 1, RMSSD, and pNN50 = % |ΔRR| > 50 ms.

**Window set.** All models use the same complete-case windows: no error, all
features finite. Exclusions are reported by label.

**AUC estimates:**
- **Primary:** out-of-fold AUC of the pooled out-of-fold probabilities for each
  model, with a cluster-bootstrap CI. The fitted out-of-fold predictions are
  resampled; the models are not refitted.
- **Paired differences** M3 − M1, M3 − M2 and M4 − M0, with the same resamples for
  both models.
- **Secondary:** DeLong CIs and paired DeLong tests, noting that windows within a
  subject are not independent, so DeLong is anticonservative here.

**Interpretation rules:**
- "UPO adds information beyond the LLE" iff the lower bound of M3 − M1 is > 0.
- "The LLE adds beyond UPO" iff the lower bound of M3 − M2 is > 0.
- "The chaos measures add beyond HRV" iff the lower bound of M4 − M0 is > 0.

**Known before the run (QC.md 7, no detector output).** Ectopic beats inflate
RMSSD and pNN50 on raw RR: median RMSSD is 262 ms in abnormal vs 65 ms in normal
windows. M0 may therefore be near its ceiling, leaving little room for M4 − M0.

### SECONDARY and descriptive

Secondary and descriptive results support **no confirmatory claim**:
- Q1 and Q3 repeated on the `edited` arm. It has 7 abnormal windows, and its Q2 is
  not computed.
- Components (LLE alone, UPO alone, OR) as classifiers, by label, with CIs.
- A per-record table: windows by label, AND / LLE / UPO detections, errors, and
  median scores.
- Abnormal windows by beat type: **mainly ventricular** (V + E + F > A + a + J + S)
  vs **mainly supraventricular**.
- **Atrial fibrillation:**
  - windows of records with AFIB / AFL rhythm annotations (201, 202, 203, 210,
    217, 219, 221, 222) vs others, by label;
  - windows overlapping an AFIB / AFL episode vs not.
- `raw` by missed-beat presence. Q1 is also re-estimated without the 39
  missed-beat windows.
- The 75 ms relabelled subset (Section B.2).

### SENSITIVITY

Q1, Q2 and Q3 are rerun unchanged on the `annotation` arm and reported beside the
primary results. This involves no formal test of the difference.

Anything chosen after detector output exists is labelled **EXPLORATORY** and kept
separate.

## D. Runtime and outputs

**Runtime.**
- Tasks: 305 `raw` + 218 `edited` + 368 `annotation` = **891 windows**.
- Synthetic code check: 0.74 s wall per window on 4 workers (72 windows in 53 s),
  i.e. **about 11 minutes**, with 1 hour allowed.
- The runner is resumable; finished task ids are skipped.

**Command:**

```bash
OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /root/venv313/bin/python -m experiments.phase7_mitbih.run_phase7 --phase test --workers 4
/root/venv313/bin/python -m experiments.phase7_mitbih.analysis
```

**Outputs:**
- `results/run/raw.jsonl`, `results/run/edited.jsonl`,
  `results/run/annotation.jsonl` (one line per window), and
  `results/run/manifest.jsonl` (code hashes, git head, versions);
- `results/analysis/results.json` and `results/analysis/results.md` (all
  preregistered tables);
- `results/analysis/{arm}_oof_predictions.csv`;
- `results/analysis/{arm}_{by_label, abnormal_by_beat_type,
  by_label_and_af_record, by_label_and_af_overlap, by_label_and_missed_beat,
  per_record}.csv`;
- the report, `docs/PHASE7_MITBIH_RESULTS.md`.

## Amendments

(none)
