# Phase 9: ECG waveform analysis, surrogates and noise-robust chaos measures

**Branch:** `claude/bold-dirac-fc69kl`. **Code:** `experiments/phase9_waveform/`. Living notes
in `HANDOFF.md`; sources and verification in `METHODS.md`; development diagnosis in
`DIAGNOSIS.md`; preregistrations in `PREREGISTRATION.md` (Part E) and
`PREREGISTRATION_F.md` (Part F). Part 0 (Phase 8 close-out) is
`docs/PHASE8_CARDIAC_CHAOS.md`.

## 1. Summary

**Question.** Can waveform analysis, multivariate beat features or noise-robust,
scale-resolved chaos measures break the sensitivity/specificity trade-off of the earlier
RR-only detectors?

**Answer.**
- **No for the waveform route.**
  - All three waveform surrogates are UNUSABLE: they reject a strictly periodic ECG
    5/5.
  - The published delineator failed its QT-database check.
  - T-wave features cannot be measured at the realistic noise level.
- **Partly for the beat route.** The preregistered winner `k3_mnlp_growth_ann` combines:
  - ectopy-related intervals masked using beat annotations;
  - a nonlinear-prediction determinism test against masked IAAFT surrogates;
  - a gate requiring the forecast error to keep growing with horizon.
- **Results of k3 on TEST** (new seeds ≥ 9500, nstdb noise, detected beats, V jitter):
  - specificity: **0 detections in 10,100 null and non-chaotic windows**, across 101
    conditions;
  - sensitivity: **90/1,320 pooled chaotic windows** at the most realistic variant
    (90/990 if the KTz morphology family, which is undetectable from RR, is set aside).
- **Comparison with the frozen baseline (Phase 8 C1).** C1 detected 74/1,320, but it
  failed specificity on non-chaotic coupled vdP with ectopy (up to 45/100).
- **The trade-off is improved, not broken.** The new detector is about as sensitive as C1
  and fully specific, but sensitivity is still low:
  - 7 % pooled;
  - 0 % with bigeminy;
  - 0 % for chaos visible only in the T wave.
- **On real long-term recordings (Part F),** the winner detected 1/180 nsrdb windows and
  0/150 chfdb windows. There is no excess over the synthetic false-positive rate and no
  group difference.

## 2. Sources and verification (`METHODS.md`)

| method | source used | verification | status |
|---|---|---|---|
| ECGSYN driven by external beat times | McSharry et al. 2003 + authors' `ecgsyn.c` (compiled) | vs C code r = 0.99999999994; Table I times; Figs 9–10 | VERIFIED |
| ectopic morphology (PVC / PAC) | LITFL PVC / PAC criteria, predeclared factors | calibrated widths and amplitudes | predeclared |
| KTz APD morphology family | Kinouchi–Tragtenberg z map, arXiv:2202.12406 | LE 0.00344 at P = 92 vs published ≈ 0.0035 | VERIFIED (TEST family) |
| modified LR1 EAD family | Tran et al. 2009 | chaos not reproduced | **DROPPED** |
| nstdb noise, nst SNR, MIT-BIH quantization | PhysioNet nstdb, nst docs | checksums; SNR definition | used |
| D1 SDLE | Gao et al. 2011 (PMC3264951) | Lorenz plateau 1.57–1.68 vs λ₁ 1.51; noise −γ ln ε | VERIFIED |
| D2 FSLE | Boffetta et al. 2002 Eq. 3.37; TISEAN estimator | Cencini map 0.888 vs 0.875; noise divergence | VERIFIED |
| D3 ε-entropy | Cencini et al. 2000; GP sums | logistic ln 2, Hénon 0.42, TISEAN d2 to 0.09 % | VERIFIED |
| D4 PE / CECP | Bandt–Pompe; Rosso; ordpy | exact d = 3 distributions; ordpy to 1e-15 | VERIFIED |
| D5 RQA DET (fixed recurrence rate) | Marwan 2011 | AR(3) 0.602 vs 0.6; Rössler +0.02–0.04 high | VERIFIED (deviation recorded) |
| D6 GHKSS noise reduction | TISEAN `ghkss.c`, Hegger et al. 1999 | 5.12 × noise reduction = TISEAN | VERIFIED |
| B2 PPS (noise-radius rule) | Small et al. 2001 via Luo et al. (open) | periodic + white noise rejected 10/10 | **NOT VERIFIED** |
| B2 cycle shuffle | Theiler 1995 via Luo et al. | periodic + white noise rejected 10/10 | **NOT VERIFIED** |
| B2 twin surrogates | Thiel et al. 2006/2008 | phase-synchronisation test 0/10 at ε 0.02, 10/10 at ε 0.045 | VERIFIED (null contains chaos: uninformative for chaos) |
| C1 delineator (NeuroKit2 DWT, Martínez 2004) | NK2 0.2.13; published QTDB accuracy (PMC3076264) | QRS onset −15.5 ± 37.2 ms vs 4.6 ± 7.7; T end −20.9 ± 34.9 vs −1.6 ± 18.1 | **NOT VERIFIED, dropped** |
| ecgpuwave | Laguna 1994 | no Fortran compiler | **DROPPED** |
| fixed-window RTp / Ta | predeclared windows | T apex +29.8 ± 79.6 ms at 12 dB (criterion ±13.9, SD ≤ 27.8) | **NOT VERIFIED at 12 dB** (passes post hoc at 24 dB, −1.5 ± 25.8) |
| C2 multivariate IAAFT | Schreiber & Schmitz 2000 Eqs. 18–20 | distributions exact; CCF error 0.05 vs 0.57; size 9/100 | VERIFIED (unused: no verified features) |
| ectopy masking | Peltola 2012 (PMC3358711): "both the ectopic beat and the following compensatory pause are edited" | equivalence test vs candidate code | used (k2, k3, k5) |
| forecast-error growth | Sugihara & May 1990 (primary not accessible; restated in PMC9760897) | development diagnosis | used (k3, k4) |

## 3. Synthetic library (Part A)

- **RR sources.** The RR series come from:
  - the 17 Phase 5–6 nulls, with their built-in ectopy, bit-identical to Phase 6;
  - every Phase 8 regime, with ectopy S2_5 / E1 / E3_10;
  - the KTz family, whose RR is the pacing input.
- **Signal chain.** ECGSYN at 360 Hz with predeclared ectopic morphology → nstdb
  bw + ma + em at 24 / 12 / 6 dB → MIT-BIH quantization → Pan–Tompkins beats → Phase 7
  V-offset jitter.
- **Realism.** Statistics and plots are in `results/realism/realism.md` and `plots/A_*`.
- **Known realism gaps:**
  - The coupled-vdP RR variability is far larger than human HRV (SDNN 300–450 ms), and its
    shortest intervals (≈ 0.2 s) are partly missed by beat detection.
  - ECGSYN PVCs enlarge S/T together with the compensatory interval.
  - Abrupt RR changes shift the baseline.
  - The dynamics are those of low-dimensional models, not human regulation.
- **The split** was recorded in `HANDOFF.md` before any measure was run:
  - Mackey–Glass is DEV;
  - phase-reset, coupled vdP, AV node and KTz are TEST;
  - KTz is the only morphology family, so development had none;
  - DEV seeds are 9000–9499 and TEST seeds start at 9500.
- **DEV-only control.** A sine circle-map family (quasi-periodic and locked, all
  non-chaotic) was added before any measure ran. It is the development stand-in for the
  quasi-periodic failure disclosed from Phase 8.

## 4. Surrogate false-positive validation (Part B3)

| surrogate | exact null (plain language) | FP on waveform nulls | verdict |
|---|---|---|---|
| PPS | a periodic orbit plus independent noise; no deterministic dependence between cycles | periodic ECG: 5/5 clean, 3/3 at 12 dB | **UNUSABLE** |
| cycle shuffle | cycles are independent draws of cycle shapes | periodic ECG: 5/5 clean, 3/3 at 12 dB | **UNUSABLE** |
| twin surrogates | another trajectory of the same system (chaotic or not) | periodic ECG: 5/5 clean, 0/3 at 12 dB | **UNUSABLE**; and its null contains chaos |
| masked IAAFT (RR) | unmasked intervals are a monotone transform of a stationary linear Gaussian process; ectopy-related intervals ignored | DEV nulls: thresholded at ≤ 2 % per condition + 0.5; TEST nulls 0/3,400 for k3 | used |
| multivariate IAAFT | stationary multivariate linear Gaussian process with the observed auto- and cross-spectra, monotone per channel | verified; size 9/100 | verified, not used |

Details for the waveform surrogates:
- Each PPS jump and each CS junction adds roughness that a prediction error sees. On a
  quantized periodic ECG, TS twins are not exact.
- In a smoke test, PPS and CS also rejected correlated HRV (N1).
- The run was interrupted after 8 of 65 windows by a container restart. It was not resumed,
  because every surrogate had already exceeded 7 % on one null.
- **B5 compute:** about 200 s per window, about 3.5 × the whole budget for the waveform arm.
  **The waveform arm was dropped.**

## 5. Development diagnosis (`DIAGNOSIS.md`)

1. **Determinism, not chaos.** Every Part D measure, and nonlinear prediction, separates
   deterministic from linear-stochastic signals. None separates chaotic from periodic or
   quasi-periodic dynamics. NLP z ≈ −9 to −14 for chaotic MG, non-chaotic MG,
   quasi-periodic and locked circle maps alike.
2. **Ectopy fakes chaos.** Ectopy in non-chaotic deterministic rhythms looks chaotic even
   after the pipeline's RR outlier mask (C3 zmax ≈ 20). Annotation masking removes the
   effect.
3. **The chaos signature that worked.** Forecast-error growth E(h5) − E(h2) is 0.37–0.40
   for chaos vs ≈ 0 for periodic and quasi-periodic dynamics. FSLE at 0.32 SD also separated
   them in medians.
4. **Measures that did not help.** SDLE at 512 beats and GHKSS at short N added nothing.

**DEV calibration** (`calibrate9.py`; 6,450 windows; rules predeclared):

| threshold | value |
|---|---|
| Z_DET | 9.4521 |
| G_MIN | 0.25 |
| Z_DET_RR | 12.1603 |
| G_MIN_RR | 0.85 |
| F_MIN | 1.6183 |

DEV results at these thresholds:

| candidate | DEV chaotic | DEV non-chaotic | DEV nulls |
|---|---|---|---|
| k2 | 586/800 | 943/2,250 | 3/3,400 |
| k3 | 411/800 | 1/2,250 | 0/3,400 |
| k4 | 16/800 | 1/2,250 | 0/3,400 |
| k5 | 3/800 | 2/2,250 | 0/3,400 |

## 6. Preregistered rule and TEST results (Part E)

- `PREREGISTRATION.md` was pushed (`31cd725`) before any TEST seed.
- **PASS:** ≤ 7/100 in each of 101 conditions:
  - 17 nulls at the input level and 17 at the realistic level;
  - 18 Phase 8 non-chaotic regimes at the input level, and the same 18 at the realistic
    level with ectopy cycling;
  - 31 KTz non-chaotic regimes.
- **WINNER:** the highest pooled count on chaos_primary: TEST chaotic regimes at the most
  realistic variant, 1,320 windows.
- **Run:** 12,290 windows, 30,052 s (8.3 h) on 4 workers.
- **Errors:** 63 window errors, all "fewer than n intervals" in fast coupled-vdP regimes,
  counted as not detected per the rule.

**[Erratum E1, final refinement]** `decide9.py` grouped the 31 preregistered KTz non-chaotic regimes by
pacing period only (81 conditions). Re-tabulated over the 101 preregistered conditions
(`experiments/final_refinement/errata/phase9_decision_101.md`), k2 fails 16 conditions instead of 12; every
PASS/FAIL and the winner are unchanged (`experiments/final_refinement/ERRATA.md`).

| candidate | PASS | failed conditions | primary chaotic (1,320) | detections in PASS conditions (10,100) |
|---|---|---|---|---|
| k1_frozen512 (baseline) | **no** | 4: non-chaotic coupled vdP with ectopy, 8–45/100 | 74 | 97 |
| k2_mnlp_ann | **no** | 16 [Erratum E1]: quasi-periodic / forced vdP (up to 100/100), MG τ = 15, 16 with ectopy, 4 KTz regimes without ectopy (8–12/100) | 304 | 607 |
| **k3_mnlp_growth_ann** | **yes** | none | **90** | **0** |
| k4_mnlp_growth_rr | yes | none | 2 | 1 |
| k5_mnlp_fsle_ann | yes | none | 0 | 3 |

**WINNER: k3_mnlp_growth_ann.**

Primary detections by family (k3 / k1):

| family | k3 | k1 | notes |
|---|---|---|---|
| coupled vdP | 73/360 | 41/360 | k3 by ectopy: S2_5 46, E3_10 27, E1 0 |
| phase-reset | 17/630 | 33/630 | |
| KTz | 0/330 | 0/330 | expected: chaos in the T wave only, RR is the pacing input |

Regimes:
- k3 detects all four chaotic vdP regimes, at 12–29 of 90 each.
- Most phase-reset regimes are near 0 for both detectors. The best is τ = 0.60: k3 14,
  k1 24.

Secondary results (report only):
- **Input level, no ectopy:** k3 34/120 vdP and 13/210 phase-reset; k2 118/120 and 164/210.
- **DEV MG chaos at TEST seeds:** k3 81/140 with S2_5 / E3_10 and 0/70 with bigeminy.
- **k1 on non-chaotic vdP with ectopy:** E1 43/150, E3_10 24/150, S2_5 11/150, none 0/150.
  This confirms the Phase 8 secondary warning (11/60) on new seeds.

**Adoption.**
- Form: `final_pipeline.masked_growth_chaos_test(rr, beat_labels, config)`, a standalone
  opt-in function with parameters in `PipelineConfig.masked_growth_*`.
- `analyze_segment` and every default are unchanged. A label argument on `analyze_segment`
  would break the label-free contract test, so it was not added.
- Tests: `tests/test_phase9_masked_growth.py`, 7 tests, including bit-equivalence with the
  preregistered candidate code.
- `groundrule_check.sh` passed:
  - pytest 462 passed + 1 pre-existing failure (default) and 463 passed (AVX-512 off);
  - replicability 44/44;
  - `replicability.json` identical.

## 7. Part F: nsrdb vs chfdb (`PREREGISTRATION_F.md`, pushed `89d4ba9` before download)

**Data and windows.**
- Annotation beat series; up to 10 windows of 512 intervals per subject, starting at hours
  1–10.
- One amendment before any detector run: the adopted form is the standalone function.

| group | subjects | windows | detected | analysable | subject-level rate (95 % cluster bootstrap CI) | exceeds synthetic FP (upper bound 0.0003)? |
|---|---|---|---|---|---|---|
| NSR (nsrdb) | 18 | 180 | 1 (subject 16272) | 180 | 0.006 (0 – 0.017) | no |
| CHF (chfdb) | 15 | 150 | 0 | 147 | 0.000 (0 – 0) | no |

- **W2 (CHF − NSR):** −0.006, 95 % CI (−0.017, 0), permutation p = 1.0. **Not supported.**
- **Secondary:** median zmax 2.6 (NSR) and 2.0 (CHF). zmax ≥ Z in 1 and 4 windows; the
  CHF ones then fail the growth gate.

**EXPLORATORY MIT-BIH** (48 records, 512-interval windows of the annotation series):
- 1/190 windows detected (record 123, a mostly normal record);
- 141/190 analysable;
- 41 of the 53 windows with > 20 % non-normal beats are not analysable.

## 8. Limitations and what I am unsure about

- **Synthetic realism.**
  - TEST chaos comes from low-dimensional models whose RR variability is unlike human HRV
    (coupled vdP SDNN 300–450 ms).
  - ECGSYN morphology is a sum of Gaussians.
  - The ectopy patterns, noise and jitter are realistic in kind, not in every detail.
  - A detector tuned on this library may behave differently on human data.
- **Phase 8 TEST disclosure.** The Phase 8 TEST families had been evaluated before (seeds
  5000–5099). Their failure modes (quasi-periodic vdP, ectopy in non-chaotic vdP) motivated
  two Phase 9 choices: the DEV circle-map controls and annotation masking. Phase 9 uses new
  seeds and never used a TEST family for development. Still, the design was not blind to
  the kind of TEST failure.
- **High-dimensional chaos vs noise.** At 512 beats, no measure here can tell
  high-dimensional chaos from correlated stochastic dynamics.
  - The determinism gate needs strong low-dimensional predictability.
  - The growth gate needs that predictability to decay within 5 beats.
  - Weak or high-dimensional chaos (MG τ = 16.5–17; most phase-reset regimes) goes
    undetected.
  - Real-data detections (Part F) are therefore not proof of chaos. Nonstationarity,
    sleep/wake transitions and nonlinear stochastic regulation are not among the controlled
    nulls.
- **Annotation dependence.** k3 needs beat labels of annotation quality:
  - On synthetic data these were generator types matched to detected beats.
  - With imperfect automatic labels, ectopy would leak into the test. k4, the RR-rule
    version, shows how much power is lost without labels.
  - Bigeminy and other ectopy-dense windows are not analysable.
- **KTz.** Morphology chaos was undetectable by every candidate. No verified T-wave
  measurement exists at 12 dB, so the morphology question remains open; at 24 dB a feature
  route may be feasible (exploratory accuracy only).
- **Interrupted B3 run.** It was stopped after 8 windows on a decisive outcome. The PPS / CS
  verdicts also rest on their own failed Rössler verification.
- **Uncertainties:**
  - **Sugihara–May citation.** The paper itself was not readable here; the open
    restatement was used.
  - **Peltola PMC number.** It was corrected during writing: PMC3358711.
  - **Development diagnosis bank.** It covered 379 of 1,820 windows.
  - **D5 Rössler deviation.** It is recorded rather than resolved.
  - **Part F window counts.** The Part F rates rest on very few detections (one window), so
    the CIs are wide in relative terms.
