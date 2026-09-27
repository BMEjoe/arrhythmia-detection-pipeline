# Phase 3: LLE as a Statistically Meaningful, Noise-Robust Chaos Detector

Research-development record in the style of the Phase 1 and Phase 2E reports.
Every statement is tied to a results file, table, preregistration or test under
`experiments/phase3_lle/` (paths below are relative to it unless they start with
`docs/`, `tests/` or `final_pipeline.py`). Synthetic data only; no MIT-BIH or
real ECG data.

## 1. Provenance and design

| Item | Value |
|---|---|
| Branch | `claude/amazing-allen-p4ssx3` |
| Environment | Python 3.13.12, numpy 2.1.3, scipy 1.18.1 (`requirements.txt`), BLAS 1 thread; manifests in `results/{dev,test}/manifest.jsonl` |
| Systems | white noise, AR(1) φ = 0.8 (nulls); sinusoid period 7.3, logistic r = 3.5 (periodic); logistic r = 4, Hénon (chaotic); noisy logistic / Hénon at 30, 20, 10 dB. Phase 2E generators (`experiments/phase2e/systems.py`) |
| Windows | 256 (primary; matches `CLASSIFIER_WINDOW_RR`) and 512 |
| Development seeds | Phase 2E core / noise seeds: nulls 0–99 (256) and 0–49 (512), all other conditions 0–29 |
| Test seeds | 1000–1099 for every condition (100 per condition and length) |
| Seed discipline | Test seeds were first run after `PREREGISTRATION.md` was committed and pushed (commit `25ec207`); `run_phase3.py` refuses `--phase test` otherwise. The decision rule was coded (`decide.py`, commit `7c4a496`) before any test result was read |
| References (bias only) | logistic ln 2 = 0.693, Hénon 0.4192 (natural log per iteration) |

Rates are detections / **all** windows (a failed window counts as "not
detected"), with 95 % Wilson intervals. For null and periodic systems the rate is
the false-positive rate; for chaotic systems it is the power. Bias = median LLE −
reference, over windows with a finite estimate.

## 2. B1: baseline on development seeds

The baseline is the existing LLE path, unchanged:
- `fp.takens_embed`: TDMI τ, Cao m at τ;
- `fp.rosenstein_lle`: production early fit, Theiler = mean period;
- `fp.run_surrogate_analysis`: 199 AAFT surrogates at fixed τ and m;
- detection ⇔ `lle_p_upper` ≤ 0.05.

(`estimators.py::baseline`; results in `results/dev/baseline.jsonl`,
`results/dev/tables/B1_baseline.md`.)

| 256 samples (dev) | rate | median LLE | bias |
|---|---|---|---|
| white noise | 7/100 = 0.070 [0.034, 0.137] | 0.075 | |
| AR(1) | 4/100 = 0.040 [0.016, 0.098] | 0.077 | |
| sinusoid | 0/30 | 0.011 | |
| logistic r = 3.5 | 0/30 (LLE undefined in 30/30) | NA | |
| logistic | 14/30 = 0.467 | 0.122 | −0.571 |
| Hénon | 12/30 = 0.400 | 0.118 | −0.301 |
| logistic 30 / 20 / 10 dB | 16/30, 19/30, 8/30 | 0.117, 0.117, 0.051 | −0.577, −0.576, −0.642 |
| Hénon 30 / 20 / 10 dB | 15/30, 14/30, 16/30 | 0.114, 0.118, 0.091 | −0.305, −0.301, −0.328 |

At 512 samples the baseline detects 30/30 (logistic) and 27/30 (Hénon), and 10/30
and 22/30 at 10 dB. Its null rates are 3/50 and 1/50. The LLE values reproduce
Phase 2E Table H exactly (for example logistic 256: 0.122 [0.109, 0.151]).

So the existing test holds its size but has low power at 256 samples. Its LLE
is 3–6× too small, and white noise (0.075) and chaos (0.12) overlap
(`docs/PHASE2E_SYNTHETIC_VALIDATION.md` 3.2).

## 3. B2: diagnosis on development seeds

Source: `results/dev/diagnosis.jsonl`, `results/dev/tables/B2_diagnosis.md`,
`plots/B2_divergence_curves_N{256,512}.png`; first 30 development seeds per
condition. Every curve is the frozen `fp.rosenstein_lle` mean log-divergence
curve, computed on four embeddings:
- V0 production;
- V1 τ = 1 with the production m;
- V2 τ = 1 with Cao at lag 1;
- V3 τ = 1 with the true map dimension (2 for non-map systems).

| Suspected cause | Test | Measured effect (256, median) |
|---|---|---|
| **Delay** | V0 → V1: only τ changes (TDMI τ* = 8 logistic, 11 Hénon → 1) | bias logistic −0.571 → −0.400; Hénon −0.301 → −0.141 |
| **Embedding dimension** | V1 → V2 → V3 at τ = 1 (m* ≈ 10 → Cao lag-1 m = 4 → true m) | logistic −0.400 → −0.022 → −0.000; Hénon −0.141 → −0.046 → −0.020 |
| Dimension under noise | Cao lag-1 m at 30 / 20 / 10 dB | logistic 8 / 9 / 9 and Hénon 7 / 8 / 8 (clean 4). V2 bias at 20 dB: −0.440 (logistic), −0.192 (Hénon) |
| **Fit region** | step k95 at which the curve reaches 95 % of its rise, vs the fixed fit k = 0..5 | production V0: k95 = 2 in 30/30 chaotic windows, so the fit spans the plateau. V3 clean: k95 = 9 (logistic), 12 (Hénon), so the fit is in the linear region |
| **Noise floor** | increment y[1] − y[0] at τ = 1 (V3) | clean logistic 0.69, Hénon 0.29. Logistic 30 / 20 / 10 dB: 2.98 / 3.88 / 4.69. Hénon: 0.67 / 1.12 / 1.68. White noise 1.97, AR(1) 1.72 |
| Noise-floor consequence | slope over k = 0..5 vs k = 1..6 (V3) | white noise 0.424 → 0.088; AR(1) 0.464 → 0.158; logistic 20 dB 0.898 (+0.205) → 0.308 (−0.385); Hénon 20 dB 0.480 → 0.323; clean logistic 0.693 → 0.696 (unchanged) |

Conclusions:
1. The production bias comes from the TDMI delay (about 30–50 % of it) and from
   the large Cao dimension at that delay (the rest). Together they saturate the
   curve inside the fit window.
2. With τ = 1 and a small m the Rosenstein curve is unbiased for clean maps.
   `results/posthoc/lle.json` in Phase 2E confirms this up to N = 2048.
3. With τ = 1, noise and stochastic data produce a large first-step jump. A fit
   that includes k = 0 turns this jump into a large spurious "LLE" (white noise
   0.42, noisy logistic 0.90). A fit starting at k = 1 removes most of it.
4. A positive LLE is still not specific. Even starting at k = 1, white noise
   (0.09) and AR(1) (0.16) give positive slopes, so a null-hypothesis test is
   required.

## 4. B3: candidates and preregistration

Tuning used development seeds only (`tune.py`, `results/dev/tuning.jsonl`,
39 surrogates, 256 samples). The grid:
- m ∈ {2, 3};
- neighbour rule ∈ {1-NN, ε-exclusion, 5-NN mean};
- surrogate ∈ {AAFT, IAAFT};
- fit ∈ {0..5, 1..5, saturation 50 %, saturation 30 %};
- plus the 0–1 test.

Findings:
- a fit from k = 1 at τ = 1 gave development power 1.00 at 30, 20 and 10 dB for
  most LLE combinations;
- IAAFT kept AR(1) nearer nominal than AAFT;
- the 0–1 test had power ≤ 0.07.

Four candidates were frozen in `PREREGISTRATION.md` (commit `25ec207`). All use
τ = 1, Theiler = mean period, 99 IAAFT surrogates from a window-specific RNG,
and a one-sided p = (1 + #{surrogate ≥ observed or undefined}) / 100, with
detection ⇔ p ≤ 0.05 (exact size 0.05):

| id | method name | statistic |
|---|---|---|
| C1 | `c1_rosenstein_m2_iaaft` | m = 2, Rosenstein 1-NN curve, slope of y[1..5] |
| C2 | `c2_kantz_m3_sat_iaaft` | m = 3, log of the mean distance to 5 nearest neighbours; slope of y[1..k_end], k_end at 30 % of the rise to the plateau (3..12) |
| C3 | `c3_eps_m2_iaaft` | m = 2, 1-NN beyond ε = median NN distance (noise-scale exclusion), slope of y[1..5] |
| C4 | `c4_zero_one_iaaft` | 0–1 test K (independent comparison; no LLE) |

**Primary rule** (test seeds, 256 samples):
- **PASS:** the Wilson lower bound of the false-positive rate is ≤ 0.05 on
  white noise and on AR(1);
- **WINNER:** among passing candidates, the highest pooled detection on
  logistic and Hénon at 20 dB;
- **tie-break:** the lowest |bias| summed over clean logistic and Hénon.

The baseline is reported but not eligible. Development results of the frozen
methods are in `PREREGISTRATION.md` Section 7. On development seeds all three
LLE candidates detected every chaotic window, so the test-seed winner was
expected to be decided by the tie-break. This was stated in the preregistration
before the test run.

## 5. B4: test-seed results (all methods, as preregistered)

(to be completed from `results/test/tables/summary.md` and `results/test/tables/primary.md`)
