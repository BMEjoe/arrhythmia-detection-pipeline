**DRAFT — not yet the preregistration (thresholds pending); no TEST seed may be used while this line is present.**

# Phase 9 preregistration: annotation-masked chaos detectors on realistic ECG-derived RR

This file is committed and pushed **before any TEST seed (≥ 9500) is generated and before any
candidate is run on any TEST family** (phase_reset, coupled_vdp, av_node, KTz).

**What exists at that point** (development only, seeds 9000–9499):
- generator verification (ECGSYN, KTz);
- KTz ground-truth labels;
- realism statistics (realism seeds 980000–980004);
- surrogate, measure and feature verifications;
- the development diagnosis (`DIAGNOSIS.md`, `results/dev/diag*.{jsonl,md}`);
- the calibration (`calibrate9.py`, `results/dev/calib_raw.jsonl`, `results/dev/calibration9.json`).

**Prior knowledge** of the Phase 8 TEST families is disclosed in `HANDOFF.md`. In brief:
- the frozen detector is weak but specific;
- NLP vs IAAFT fails on quasi-periodic coupled vdP;
- ectopy harms non-chaotic vdP.

This motivated the DEV-only circle-map controls and the annotation mask.

**What is frozen after this file:** `candidates9.py` (thresholds included), `library.py`,
`families.py`, `ecgsyn.py`, `morphology.py`, `noise.py`, `regimes9.py`, `run_test9.py`,
`decide9.py`, the KTz labels, and `final_pipeline.py`. Any unavoidable change goes in a dated
amendment saying whether any TEST output had been seen. `run_test9.py` refuses to start unless
this file is committed, unmodified, pushed, and lists every candidate.

FROZEN_TABLE_PLACEHOLDER

## 1. Library and split (`HANDOFF.md`, PHASE 9 SPLIT)

- **TEST families:**
  - phase-resetting map: 7 CHAOTIC, 5 NON-CHAOTIC;
  - coupled modified vdP: 4 CHAOTIC, 6 NON-CHAOTIC;
  - AV node: 4 NON-CHAOTIC;
  - KTz paced-cell map → APD → T wave (morphology family): 11 CHAOTIC and 31 NON-CHAOTIC
    (P, input) regimes.
- **DEV family:** Mackey–Glass (7 CHAOTIC, 3 NON-CHAOTIC), used for development and calibration
  and included in the non-chaotic PASS conditions, as in Phase 8.
- **Labels:** Phase 8 `regimes.json`, `results/ground_truth/ktz_labels.json`.
- **Signal chain:** RR (+ ectopy) → ECGSYN with the predeclared ectopic morphology → 360 Hz →
  nstdb mix (bw + ma + em equal power) at 12 dB, nst definition, signal 1 for TEST → MIT-BIH
  quantization → Pan–Tompkins beats → Phase 7 V-offset jitter at V beats.
  - **Input level:** true beat times, no recording noise.
  - **Beat labels:** the generator's beat types of the matched beats (annotation quality, as in
    PhysioNet databases); 'X' for a detection matched to no beat.
- **The waveform arm and the T-wave / QT features are not candidates.** They were dropped on
  B3, B5 and C1 verification (`METHODS.md`).

## 2. Candidates (`candidates9.py`; 512-interval windows; every candidate has a null test)

- method: `k1_frozen512`
  - **K1** = the frozen Phase 6 detector (Phase 8 C1), baseline, on the window's RR.
  - Null control: LLE IAAFT test and the UPO surrogate test plus instability gate.
- method: `k2_mnlp_ann`
  - **K2** = annotation-masked nonlinear-prediction determinism test.
  - **Mask:** every ectopy-related interval (beat k or k+1 not 'N').
  - **Statistic:**
    - robust normalized local-average prediction error (k 5, Theiler 5, unit delay), h = 1;
    - computed only on delay vectors whose span to h = 5 contains no masked interval;
    - z_m = (mean − observed) / SD over 39 masked IAAFT surrogates (IAAFT of the interpolated
      series, statistic with the same mask), for m = 2, 3, 4, 5.
  - **Detected iff** analysable (≥ 100 valid vectors at every m) and zmax ≥ **Z_DET**.
  - **Null:** the unmasked intervals are a monotone static transform of a stationary linear
    Gaussian process; ectopy-related intervals are ignored.
- method: `k3_mnlp_growth_ann`
  - **K3** = K2 AND forecast-error growth.
  - G = E(h 5) − E(h 2) of the masked prediction error at m = 3 must satisfy G ≥ **G_MIN**.
- method: `k4_mnlp_growth_rr`
  - **K4** = K3 with the pipeline's RR outlier mask (`fp.correct_rr_intervals`) instead of
    annotations.
  - Thresholds **Z_DET_RR**, **G_MIN_RR**.
- method: `k5_mnlp_fsle_ann`
  - **K5** = K2 AND FSLE gate (Part D, D2).
  - FSLE at δ = 0.32 SD (m 3, δ0 0.01 SD, r √2, level 10) of the series with masked intervals
    linearly interpolated must satisfy FSLE ≥ **F_MIN**.

**Thresholds** (`calibrate9.py`; DEV seeds 9000–9219 only; rules in the calibrate9 docstring):
THRESHOLD_PLACEHOLDER

## 3. TEST seeds and windows (`run_test9.py`)

- TEST seeds start at 9500 and were never used before.
- Each window is generated once and evaluated by every candidate.
- All windows are 512 intervals. 12,290 windows in total.

| group | content | seeds | windows |
|---|---|---|---:|
| null_input | 17 PASS nulls, input level | 9500–9599 | 1,700 |
| null_real | 17 PASS nulls, realistic (nstdb + detected + jitter) | 9500–9599 | 1,700 |
| noncha_input | 18 Phase 8 NON-CHAOTIC regimes (15 TEST + 3 DEV), input level, no ectopy | 9500–9599 | 1,800 |
| noncha_real | 18 Phase 8 NON-CHAOTIC regimes, realistic, ectopy none / S2_5 / E1 / E3_10 by seed mod 4; 31 KTz NON-CHAOTIC regimes, realistic | 9500–9599 | 4,900 |
| chaos_primary | 11 Phase 8 TEST CHAOTIC × (S2_5, E1, E3_10), realistic; 11 KTz CHAOTIC, realistic | 9500–9529 | 1,320 |
| chaos_second | TEST CHAOTIC at input level (no ectopy / KTz input); DEV MG CHAOTIC × 3 ectopy, realistic (seeds 9500–9509) | 9500–9529 | 870 |

## 4. Primary rule (`decide9.py`)

**PASS** iff detections ≤ floor(0.07 N) = **7/100** in EACH condition of null_input,
null_real, noncha_input and noncha_real:
- 17 + 17 nulls;
- 18 + 49 non-chaotic regimes;
- 101 conditions in all.

An error or an unanalysable window counts as not detected.

**WINNER:** among passing candidates, the highest pooled detection count on chaos_primary
(1,320 windows, the most realistic variant).

**Tie-breaks:**
1. fewer pooled detections over all PASS conditions;
2. candidate index.

**If no candidate passes, that is the result.**

## 5. Secondary outcomes (report only)

- Detection by family, regime, ectopy pattern and level.
- K1 on the same windows.
- Analysable fraction (bigeminy).
- KTz chaos: expected undetectable from RR (its RR is the pacing input).
- DEV MG chaos at TEST seeds.
- Beat-detection errors.

## 6. Runtime and adoption

**Runtime.** Measured on DEV windows under 4-worker load: about 7–12 s per window for all
candidates (K1 ≈ 5–8 s), plus 0.1–3 s generation. 12,290 windows ≈ 1.2–1.5 × 10⁵ CPU-s ≈
8.5–10.5 h on 4 workers, within the ~16 h budget.

**Adoption.** If the winner is a new candidate (K2–K5), it is added to `final_pipeline.py` as
an opt-in `PipelineConfig` option (default off) with tests, and `groundrule_check.sh` must
pass. If K1 wins, nothing new is adopted, because K1 is already the pipeline's recommended
detector.

**Part F** (only if there is a winner) has its own preregistration, committed before any
nsrdb / chfdb download.

## Amendments

(none)
