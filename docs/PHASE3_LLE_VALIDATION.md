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

Sources: `results/test/*.jsonl` (12,000 windows: 5 methods × 12 conditions ×
2 lengths × 100 seeds), `results/test/tables/summary.{md,csv}` and
`results/test/tables/primary.md` (`decide.py`). The runner completed all 12,000
tasks. Failures are only undefined LLEs on logistic r = 3.5 (200 per LLE
method; the curve has no finite values because of exact duplicate vectors),
counted as "not detected".

### 5.1 Primary decision (256 samples)

| method | white noise FP (Wilson) | AR(1) FP (Wilson) | PASS | 20 dB pooled detection | clean \|bias\| sum (tie-break) |
|---|---|---|---|---|---|
| baseline (not eligible) | 9/100 = 0.090 [0.048, 0.162] | 5/100 = 0.050 [0.022, 0.112] | (yes) | 123/200 = 0.615 | 0.870 |
| **C1** `c1_rosenstein_m2_iaaft` | 4/100 = 0.040 [0.016, 0.098] | 3/100 = 0.030 [0.010, 0.085] | yes | 200/200 = 1.000 | **0.006** |
| C2 `c2_kantz_m3_sat_iaaft` | 9/100 = 0.090 [0.048, 0.162] | 0/100 = 0.000 [0.000, 0.037] | yes | 200/200 = 1.000 | 0.024 |
| C3 `c3_eps_m2_iaaft` | 5/100 = 0.050 [0.022, 0.112] | 3/100 = 0.030 [0.010, 0.085] | yes | 200/200 = 1.000 | 0.017 |
| C4 `c4_zero_one_iaaft` | 0/100 = 0.000 [0.000, 0.037] | 2/100 = 0.020 [0.006, 0.070] | yes | 11/200 = 0.055 | NA (no LLE) |

**All four candidates pass. C1, C2 and C3 tie at 200/200 on the 20 dB
criterion, so the preregistered tie-break decides:**
- C1 |bias| sum 0.0057;
- C3 0.0173;
- C2 0.0236.

**Winner: C1.** Against the baseline, C1 raises pooled 20 dB detection from
123/200 (0.615) to 200/200. It also lowers the clean |bias| sum from 0.870 to
0.006. Its margin over the runner-up C3 is only on the tie-break (0.012 in
summed |bias|); C1 and C3 do not differ on detection at 20 dB.

C2 detects 9/100 white-noise windows, the largest count that still passes
(k = 10 would fail).

### 5.2 Detection and false-positive rates (secondary outcomes), detections / 100

| method | N | WN | AR(1) | sinusoid | logistic r 3.5 | logistic | Hénon | log 30 dB | Hén 30 dB | log 20 dB | Hén 20 dB | log 10 dB | Hén 10 dB |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| baseline | 256 | 9 | 5 | 0 | 0 | 58 | 49 | 53 | 50 | 75 | 48 | 28 | 47 |
| C1 | 256 | 4 | 3 | 0 | 0 | 100 | 100 | 100 | 100 | 100 | 100 | 100 | 99 |
| C2 | 256 | 9 | 0 | 0 | 0 | 100 | 100 | 100 | 99 | 100 | 100 | 100 | 100 |
| C3 | 256 | 5 | 3 | 0 | 0 | 100 | 100 | 100 | 100 | 100 | 100 | 100 | 99 |
| C4 | 256 | 0 | 2 | 0 | 0 | 3 | 6 | 3 | 6 | 4 | 7 | 5 | 4 |
| baseline | 512 | 7 | 8 | 0 | 0 | 98 | 89 | 98 | 92 | 96 | 88 | 54 | 72 |
| C1 | 512 | 6 | 6 | 0 | 0 | 100 | 100 | 100 | 100 | 100 | 100 | 100 | 100 |
| C2 | 512 | 3 | 3 | 0 | 0 | 100 | 100 | 100 | 99 | 100 | 100 | 100 | 100 |
| C3 | 512 | 2 | 4 | 0 | 0 | 100 | 100 | 100 | 100 | 100 | 100 | 100 | 100 |
| C4 | 512 | 5 | 7 | 0 | 0 | 3 | 4 | 4 | 5 | 4 | 8 | 4 | 5 |

- Pooled logistic + Hénon detection at 256 (30 / 20 / 10 dB):
  - baseline 103 / 123 / 75;
  - C1 200 / 200 / 199;
  - C2 199 / 200 / 200;
  - C3 200 / 200 / 199;
  - C4 9 / 11 / 9.
- Periodic rejection: no method flags the sinusoid or logistic r = 3.5 in any of
  its 200 windows (0/100 at each length).
- At 512 every null rate is within its Wilson interval of 0.05. The highest is
  the baseline AR(1) at 8/100 [0.041, 0.150].

### 5.3 LLE bias (median − reference; [Q1, Q3] of the estimate)

| method | N | logistic | Hénon | log 30 dB | Hén 30 dB | log 20 dB | Hén 20 dB | log 10 dB | Hén 10 dB |
|---|---:|---|---|---|---|---|---|---|---|
| baseline | 256 | −0.567 [0.107, 0.141] | −0.302 [0.098, 0.134] | −0.575 | −0.295 | −0.587 | −0.297 | −0.640 | −0.323 |
| **C1** | 256 | **−0.001** [0.684, 0.700] | **−0.005** [0.400, 0.432] | −0.096 | −0.054 | −0.213 | −0.068 | −0.428 | −0.183 |
| C2 | 256 | +0.007 [0.693, 0.713] | +0.017 [0.424, 0.449] | −0.042 | −0.011 | −0.095 | −0.023 | −0.244 | −0.043 |
| C3 | 256 | −0.011 [0.674, 0.691] | −0.007 [0.396, 0.428] | −0.109 | −0.058 | −0.221 | −0.079 | −0.439 | −0.190 |
| baseline | 512 | −0.539 | −0.279 | −0.546 | −0.273 | −0.556 | −0.273 | −0.632 | −0.289 |
| **C1** | 512 | **−0.000** [0.689, 0.697] | **−0.001** [0.407, 0.430] | −0.091 | −0.048 | −0.193 | −0.059 | −0.419 | −0.171 |
| C2 | 512 | +0.008 | +0.014 | −0.052 | −0.021 | −0.076 | −0.001 | −0.215 | −0.014 |
| C3 | 512 | −0.002 | −0.005 | −0.101 | −0.059 | −0.204 | −0.067 | −0.426 | −0.179 |

- On clean data the bias drops from −0.567 / −0.302 (baseline, 256) to
  −0.001 / −0.005 (C1). The clean-bias target is met.
- Under noise C1 still underestimates: −0.213 (logistic) and −0.068 (Hénon) at
  20 dB, and −0.428 / −0.183 at 10 dB. C2 (5-NN Kantz-style, saturation-aware
  fit) has the smallest noisy bias at every SNR, for example −0.095 / −0.023 at
  20 dB. It lost only the clean-bias tie-break.

### 5.4 Runtime (median s per window, 4 workers sharing 4 cores)

| method | 256 | 512 |
|---|---:|---:|
| baseline | 6.9 | 18.2 |
| C1 | 1.0 | 3.2 |
| C2 | 1.6 | 4.4 |
| C3 | 1.2 | 3.8 |
| C4 | 3.8 | 7.0 |

The whole test run took 14,854 s wall time, as reported by the runner's progress
output. Per-task runtimes are in `results/test/*.jsonl` (`runtime_s`).

## 6. B5: adoption of the winner

C1 is in `final_pipeline.py` as an opt-in option, **default off**:
- `PipelineConfig.lle_chaos_test = False`, with parameter fields fixed at the
  preregistered values: m = 2, fit steps 1..5, 99 surrogates, α = 0.05,
  seed entropy 20260927;
- new functions `lle_chaos_test`, `iaaft_surrogate` and `_nn_divergence_curve`;
- when enabled, `analyze_segment` adds `result["lle_chaos_test"]` (LLE estimate,
  p-value, surrogate values, decision). The production LLE, its AAFT test and
  every default output are unchanged.

Checks:
- `tests/test_lle_chaos_test.py` (11 tests): the pipeline function equals
  `estimators.c1_rosenstein_m2_iaaft` bitwise on seven signal types, and the
  defaults are off.
- Recomputing 240 stored test-seed C1 results (10 per condition and length)
  through the pipeline gave 0 mismatches.
- Ground rules (`experiments/phase3_lle/groundrule_check.sh`):
  - pytest 367 passed plus the pre-existing stability failure (368/368 with
    numpy AVX-512 dispatch disabled; Phase 2E report 5.3);
  - same-machine replicability 44/44 bitwise identical;
  - `run_phase2e --replicate` unchanged (13/44 on this machine, identical file;
    see `experiments/phase3_lle/HANDOFF.md` for why this machine gives 13/44
    with the unmodified pipeline).

## 7. Limitations

- **Synthetic maps are not RR intervals.** The chaotic systems are one-step maps
  sampled once per iteration, which is why τ = 1 is correct for them. RR
  intervals are event series with physiological variability, non-stationarity,
  artefact corrections and 1/360 s quantization (Phase 2E report 6). The
  winning choices (τ = 1, m = 2, fit steps 1..5) were justified and tested only
  on these maps. None of the rates here transfer to MIT-BIH without separate
  validation.
- **Quantization.** Logistic r = 3.5 has exact duplicate vectors, so every
  candidate returns an undefined LLE (counted as "not detected"). RR data
  quantized at 1/360 s will contain many exact ties (Phase 2E 3.7). The
  1e-15 floor then removes neighbour pairs, and the effect on C1 was not
  measured.
- **The LLE estimate is not specific without the test.** C1's estimate on the
  nulls is positive (test, 256: white noise 0.122 [0.109, 0.134], AR(1) 0.204
  [0.187, 0.221]; sinusoid 0.019). Only the surrogate decision separates them
  from chaos. The number must not be read as "λ > 0 ⇒ chaos".
- **The surrogate null is linear-Gaussian (IAAFT).** Rejection means "not a
  monotone transform of a linear Gaussian process". It is not proof of chaos:
  nonlinear stochastic or non-stationary signals could also be rejected, and
  none were tested. Periodic signals were not rejected here only because their
  surrogates are themselves nearly periodic.
- **Noise bias remains.** C1 underestimates by 0.21 (logistic) and 0.07 (Hénon)
  at 20 dB. C2 halves this but was not the preregistered winner. A
  noise-robust LLE estimate would need a separate preregistered comparison.
- **Power saturated**, so the winner was decided by the tie-break. The test
  distinguishes the candidates only on clean bias (a 0.012 difference between C1
  and C3). This design was fixed in advance, but the ranking among C1–C3 is
  weakly informed.
- **Multiplicity.** Null rates were checked in 10 method × null cells at 256 and
  10 at 512 without adjustment. The primary rule uses 2 cells per method, as
  preregistered.
- **Development/test independence.** Test seeds come from the same generators
  and parameters as development seeds, so the split protects against tuning to
  noise, not against a mismatch between models and real data.
- **Machine dependence.** Results were produced on one machine. Phase 2E
  replicability across machines is 13/44 bitwise (Phase 2E report 5.3); the
  C1 statistic uses `argsort` and `polyfit`, whose ulp-level behaviour could
  differ elsewhere.

## 8. Reproducibility

```bash
uv venv -p python3.13 /root/venv313 && uv pip install -p /root/venv313/bin/python -r requirements.txt
python -m experiments.phase3_lle.diagnose && python -m experiments.phase3_lle.diagnose --summary   # B2
python -m experiments.phase3_lle.tune && python -m experiments.phase3_lle.tune --summary           # B3 tuning
python -m experiments.phase3_lle.run_phase3 --phase dev  --methods all                              # B1, B3 dev
python -m experiments.phase3_lle.run_phase3 --phase test --methods all                              # B4 (guarded)
python -m experiments.phase3_lle.analysis --phase test && python -m experiments.phase3_lle.decide   # tables, rule
python -m pytest tests                                                                              # 368 tests
```

Commits:
- `655233a` scaffold and B2;
- `fd3c751` B1;
- `25ec207` preregistration;
- `7c4a496` decision rule;
- `5cd69d5` B4 results;
- `f117480` B5.
