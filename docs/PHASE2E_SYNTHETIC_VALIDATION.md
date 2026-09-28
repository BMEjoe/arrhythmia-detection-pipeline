# Phase 2E: Synthetic End-to-End Validation of the Frozen UPO Pipeline

Research-development record. Factual and technical; every statement below is
tied to a table, results file, audit, or test committed in `88835d5`. Paths are
relative to `experiments/phase2e/` unless they start with `docs/`, `tests/` or
`final_pipeline.py`.

No MIT-BIH or real ECG data were used. `final_pipeline.py`, `config.py` and every
predeclared grid are unchanged by this report.

## 1. Baseline and provenance

| Item | Value |
|---|---|
| Pipeline commit analysed | `dd43fd0` (Reconstruct and validate source-faithful UPO pipeline) |
| Results / tables / audits commit | `88835d5` (Complete Phase 2E analysis, diagnostics, and validation outputs) |
| `final_pipeline.py` SHA-256 | `fc7a7db5…a014a0`, identical in both runs in `results/manifest.json` and in the working tree at `88835d5` |
| `config.py`, `systems.py`, `metrics.py`, `run_phase2e.py` SHA-256 | identical between `results/manifest.json` and `88835d5` |
| Tasks | 3,643 of 3,643 (`results/*.jsonl`: core 769, noise 300, tau_step 183, surrogate_sens 320, algo_seed 100, oracle_dim 60, fp_sens 360, henon_loc 1,452, henon_loc_oracleJ 72, stability 21, flag_independence 6) |
| Run environment | Python 3.13.15, numpy 2.1.3, BLAS 1 thread (`results/manifest.json`) |
| Total task runtime | 31,379 s (`results/audits.json`, `total_runtime_s`) |

The design (systems, seeds, arms) is in `README.md` and `config.py`. In brief:
`core` runs 7 systems × {128, 256, 512} samples with the production `one_sample`
map and 50 surrogates. The null systems use 50/100/50 seeds per window length, the
deterministic systems 30, and `constant` 3. The other arms are paired with the
core seeds at 256 samples unless stated otherwise.

Terminology (from `final_pipeline.py` and Table I):
**Level A** = source peaks (So et al. detection),
**Level B** = Level-A peaks with surrogate `J < 0.05`,
**Level C** = PROJECT EXTENSION verification gates (residual ≤ 0.05, R² ≥ 0.9,
support ≥ 2), never part of source detection.

## 2. Verification of the HANDOFF.md preliminary findings

`HANDOFF.md` was written at 1,988 of 3,643 tasks. Each finding is checked below
against the final data.

| # | Handoff finding | Verdict | Evidence |
|---|---|---|---|
| 1 | Hénon ~0.01 localization error is finite resolution, not a bug | **Held**, with one wording correction and a wider full-grid range | See 2.1 |
| 2 | Off-attractor fixed point: standard Hénon 0/3, skewed Hénon 3/3 at the production point | **Held** | See 2.2 |
| 3 | Period-2 Level-A peak count depends on resolution; true 2-cycle stays nearest | **Held**; the error range is wider over the full grid | See 2.3 |
| 4 | Cao at lag 1 overestimates dimension (Hénon 4–6, logistic 4) | **Held** for 256–512 samples; the 1024–2048 part **cannot be verified** from committed files | See 2.4 |
| 5 | `embedding_not_saturated` conflates "undefined" and "not saturated" | **Held**; count updated from "about 14 so far" to 43 core windows | See 2.5 |
| 6 | Sinusoid LLE at the floating-point floor, one 128-sample run +0.24; white noise LLE positive | **Changed** for the 128-sample sinusoid; **held** for 256/512 and for white noise (range updated) | See 2.6 |
| 7 | Nominal size 3/51 at n = 50; W_i mean includes surrogate i | **Held** | See 2.7 and 3.1 |
| 8 | Shared random streams across windows; tested by `algo_seed` | **Held** (design fact); no evidence the null rate depends on the fixed seed | See 2.8 and 3.1 |

### 2.1 Finding 1: Hénon localization error

Source: `results/tables/table_F_henon_localization.{md,csv}`. N = 1024, d = 2,
3 realizations per cell, production tube 10 %, range `pad10`, estimated Jacobian, M = 2.

| bins | bin width (std / skewed) | standard Hénon \|err\| (err / bin) | skewed Hénon \|err\| (err / bin) |
|---:|---|---|---|
| 15 | 0.2044 / 0.3323 | 0.0118 (0.058) | 0.0059 (0.018) |
| 20 | 0.1533 / 0.2492 | 0.0152 (0.099) | 0.0378 (0.152) |
| 30 | 0.1022 / 0.1662 | 0.0097 (0.094) | 0.0141 (0.085) |
| 45 | 0.0681 / 0.1108 | 0.0132 (0.193) | 0.0296 (0.267) |
| 60 | 0.0511 / 0.0831 | 0.0009 (0.018) | 0.0124 (0.150) |
| 90 | 0.0341 / 0.0554 | 0.0034 (0.100) | 0.0133 (0.240) |

- On this bin sweep the error is 0.0009–0.038 and at most 0.27 bin widths. It is
  not monotone in bin count. This matches the handoff.
- **Over the full 121-setting grid** (all tube percentiles and range policies) the
  error reaches 0.084 (0.37 bin widths) for standard Hénon and 0.144
  (0.50 bin widths) for skewed Hénon. The largest errors all occur at the 40 %
  tube percentile.
- **Analytic Jacobian (`henon_loc_oracleJ`).** The handoff said the error was
  "identical". It is close but **not identical**. On the bin sweep the
  estimated-minus-analytic difference in |err| is at most 0.0022 (standard) and
  0.00065 (skewed). The signed error has the same sign in all 12 paired cells.
  The Jacobian estimate is therefore not the source of the error.
- **Histogram-free KDE mode** of the same tube points: |err| 3.7e-5 to 2.7e-4
  (standard) and 6.0e-4 to 1.0e-3 (skewed) over the whole grid.
- Conclusion (unchanged): the offset comes from the pipeline's location estimator
  (mean of tube points in the peak cell) under histogram quantization.

The same mechanism is visible in the core arm (Cao dimension, N = 256/512);
see 3.3.

### 2.2 Finding 2: off-attractor fixed points

Source: Table F, column `fp2 (off-attractor) detected <=2 bins`.

- Standard Hénon (−1.131): 0/3 at the production point (30 bins, tube 10 %,
  `pad10`) with M = 2 and with M = 7. Over the grid: 0/3 in 117 of 121 settings.
  1/3 in 4 settings, all at tube 2.5 %. 0/3 in all 6 analytic-Jacobian settings.
- Skewed Hénon (−1.903): 3/3 at the production point with M = 2 and M = 7. By tube
  percentile, 3/3 in 24/24 settings at 2.5 %, 24/24 at 5 %, 22/25 at 10 %,
  13/24 at 20 % and 0/24 at 40 %.

### 2.3 Finding 3: period-2 resolution dependence

Source: Table F, period-2 columns (estimated Jacobian, tube 10 %, `pad10`, M = 2).

- Mean number of minimal-period-2 Level-A peaks, from 15 to 90 bins:
  standard 4.0 → 8.0 → 15.3 → 33.0 → 44.3 → 55.7;
  skewed 4.0 → 7.3 → 14.7 → 41.3 → 58.0 → 78.3.
- Error from the true 2-cycle to the nearest period-2 peak on this sweep:
  standard 0.0005–0.0078, skewed 0.016–0.033. Over the full grid: standard
  0.0004–0.025, skewed 0.011–0.073.
- The strongest period-2 peak lies within one bin of the true 2-cycle in 3/3
  realizations in all 121 standard-Hénon settings. For skewed Hénon: 73 settings
  3/3, 42 settings 2/3, 1 setting 1/3, 5 settings 0/3.

### 2.4 Finding 4: Cao dimension at lag 1

Sources: Table A (`cao_dim (UPO)`), `results/core.jsonl` (`embedding_dimension`),
`results/tables/table_oracle_dimension.md`.

- Hénon (true dimension 2): at 256 samples d = 4 in 25/30 windows (range 4–7); at
  512, d = 4 in 26/30 (range 4–8).
- Logistic (true dimension 1): d = 4 in 30/30 at both 256 and 512.
- At 128 samples the estimates spread further (logistic 4–13, Hénon 4–8 among
  windows that completed).
- The handoff's "256–2048 samples" range relied on the post-hoc LLE diagnostic
  (`posthoc_diagnostics.py`, `lle_diagnostics`, N up to 2048). Its output file
  `results/posthoc/lle.json` **is not in the repository** (only `tdmi.json` and
  `levelc.json` are). The 1024–2048 part therefore cannot be verified here.
- Oracle dimension (diagnostic only, never substituted): fixed d = 1 (logistic) and
  d = 2 (Hénon) gave Level-B detection 30/30 and 30/30, against 30/30 and 29/30
  with production Cao (`table_oracle_dimension.md`).

### 2.5 Finding 5: `embedding_not_saturated` status semantics

Sources: Table G, Table A, `results/tables/table_fp_sensitivity.md`,
`final_pipeline.py::cao_method`.

- `cao_method` returns the sentinel `max_dim + 1` (= 13) when E1 never reaches a
  two-step plateau. This happens both when E1 is genuinely not saturated and
  when E1 is NaN. E1 is NaN when every nearest-neighbour distance is ≤ 1e-15 (exact
  duplicate vectors) or when too few vectors exist at large m.
  `run_upo_analysis` maps both cases to `embedding_not_saturated`.
- Logistic r = 3.5 (4-cycle): 30/30 at 128, 256 and 512 (Table A, Table G).
- RR-quantized (1/360 s) sinusoid: 30/30 on both quantization routes
  (`table_fp_sensitivity.md`). Logistic: 2/30. Hénon: 0/30.
- **Updated count, white noise and AR(1), core arm:** white noise 10/50 (128),
  11/100 (256), 1/50 (512), 22/200 in total. AR(1) 8/50, 12/100, 1/50, 21/200 in
  total. Together 43 distinct windows (handoff: "about 14 so far"). Every one has
  Cao dimension 13 (Table A). The stored rows do not record which branch
  produced the sentinel. The other null-system arms (algo_seed, surrogate_sens,
  tau_step) reuse core data seeds, so they add no distinct windows for the
  `one_sample` map.

### 2.6 Finding 6: Rosenstein LLE on the sinusoid and white noise

Sources: Table H, `results/core.jsonl` (`lle_per_beat`, `lle_r2`).

| Sinusoid window | LLE median (range) | R² median | R² ≥ 0.9 |
|---:|---|---:|---:|
| 128 | 0.188 (0.098–0.508) | 0.964 | 26/30 |
| 256 | 0.0109 (0.0041–0.0415) | 0.036 | 0/30 |
| 512 | 0.0097 (0.0040–0.0416) | 0.017 | 0/30 |

- At 256/512 the handoff description (slope ≈ 0.008–0.01, R² ≈ 0.03) holds.
- **At 128 it did not hold.** The handoff described a single +0.24 outlier. In
  fact all 30/30 windows are positive, the median is 0.188, and 26/30 pass the
  R² ≥ 0.9 fit-quality threshold. No chaotic window at any length passes that
  threshold (Table H, `R2>=0.9 rate`: 0/N for logistic and Hénon).
- The handoff's statement that nearest-neighbour distances are about 1e-14 is not
  recorded in any committed results file. It is not verified here.
- White noise: 194 of 195 valid estimates are positive (44/45, 100/100, 50/50).
  Medians are 0.059 / 0.075 / 0.091; the range is −0.0009 to 0.146 (handoff:
  0.01–0.08). AR(1): 197/197 positive. See 3.2.

### 2.7 Finding 7: significance arithmetic

Sources: `final_pipeline.py::so_surrogate_significance`, `so_J`;
`results/audits.json` (`significance_audit`).

- `J(W') = mean(W_i > W')`. A peak is Level B when its J < 0.05.
- A peak's deviation is at most W, so its J is at least J(W). A window with at
  least one Level-B peak therefore always has J(W) < 0.05. The window-level
  size of J(W) < 0.05 bounds the window-level Level-B false-positive rate.
- Under exact exchangeability, J(W) < 0.05 ⇔ #{W_i > W} ≤ ⌈0.05 n⌉ − 1. The
  nominal size is ⌈0.05 n⌉ / (n + 1): **1/21 = 0.0476** (n = 20),
  **3/51 = 0.0588** (n = 50), **5/101 = 0.0495** (n = 100).
- The audit confirms the self-inclusion asymmetry. W_i uses a surrogate mean that
  contains surrogate i; the observed W does not. The scale ratio is
  √((n+1)/(n−1)) = 1.0202 at n = 50 (anti-conservative, source-faithful).
- Empirical rates are in 3.1.

### 2.8 Finding 8: shared random streams

- `final_pipeline.py` seeds the So R draws and surrogate phases from
  `config.random_seed` (lines 1599, 1698, 2083, 2124). All windows therefore share
  those streams.
- `algo_seed` (random_seed = CFG + 1000 + seed, the same 50 data seeds as the
  paired core windows): white noise 1/50 vs 3/50 paired core; AR(1) 5/50 vs 5/50
  (`table_surrogate_count.md`). At this sample size there is no evidence that
  the null rate depends on the fixed seed.

## 3. Results

### 3.1 Empirical false-positive rate on the null systems

Definition: fraction of **all** windows (failures count as "no peak") with at
least one Level-B peak. 95 % Wilson intervals are from the tables. Tail
probabilities are exact binomial P(X ≥ k) under the nominal size. They are
computed from the table counts and are **not** adjusted for multiple comparisons.

**50 surrogates (production)**, `table_B_white_noise_null.md`,
`table_C_ar1_null.md`, nominal 0.0588:

| System | 128 | 256 | 512 | Pooled (independent data) | P(X ≥ k) |
|---|---|---|---|---|---:|
| White noise | 2/50 = 0.040 [0.011, 0.135] | 7/100 = 0.070 [0.034, 0.137] | 5/50 = 0.100 [0.043, 0.214] | 14/200 = 0.070 [0.042, 0.114] | 0.29 |
| AR(1), φ = 0.8 | 2/50 = 0.040 [0.011, 0.135] | 11/100 = 0.110 [0.063, 0.186] | 2/50 = 0.040 [0.011, 0.135] | 15/200 = 0.075 [0.046, 0.120] | 0.20 |

Among windows where significance was actually assessed, the rates are higher.
For example, white noise at 256 is 7/89 (`assessed_N`, Table B).
AR(1) at 256 is the only cell whose one-sided tail probability is below 0.05
(P = 0.033) and whose Wilson interval (0.063–0.186) excludes the nominal 0.0588.
This is not adjusted for the number of cells compared. The pooled rates are
consistent with the nominal size.

**Surrogate count** (256 samples, paired data seeds 0–49),
`table_surrogate_count.md` (also Tables B/C after the A1 fix):

| System | n = 20 (nominal 0.0476) | n = 50 (nominal 0.0588) | n = 100 (nominal 0.0495) |
|---|---|---|---|
| White noise | 1/50 = 0.020 [0.004, 0.105] | 3/50 = 0.060 [0.021, 0.162] | 2/50 = 0.040 [0.011, 0.135] |
| AR(1) | 3/50 = 0.060 [0.021, 0.162] | 5/50 = 0.100 [0.043, 0.214] | 4/50 = 0.080 [0.032, 0.188] |
| Algorithm seed varied (n = 50) | WN 1/50 = 0.020 | AR 5/50 = 0.100 | |

Every Wilson interval in this table contains its nominal size. Pooling white noise and AR(1)
gives 4/100 (n = 20), 8/100 (n = 50), 6/100 (n = 100).

**Constant:** 9/9 windows return `constant_data` and are
`not_assessed_detector_failure` (Table A, Table G). There are 0 Level-B peaks,
but the surrogate test never runs, so these windows give no information about the
false-positive rate.

**Level A on nulls:** at least one Level-A peak in 70–98 % of white-noise windows
and 78–98 % of AR(1) windows (Table A). Level A alone does not separate noise
from structure. Level C gated **0** null-system peaks at every window length
(Table I).

### 3.2 Largest Lyapunov exponent by system (Table H)

Frozen Rosenstein estimator; LLE in natural log per sample (= per map iteration).
The theoretical references are not used by the pipeline.

| System | Reference | 128 | 256 | 512 | Fraction > 0 (valid) | R² ≥ 0.9 |
|---|---|---|---|---|---|---|
| Logistic r = 4 | ln 2 = 0.693 | 0.091 (20/30 valid) | 0.122 | 0.141 | 80/80 | 0/80 |
| Hénon | ≈ 0.42 | 0.115 (6/30 valid) | 0.118 | 0.150 | 66/66 | 0/66 |
| Sinusoid | ≤ 0 | 0.188 | 0.011 | 0.010 | 90/90 | 26/90 (all at 128) |
| White noise | undefined | 0.059 | 0.075 | 0.091 | 194/195 | 0/195 |
| AR(1) | undefined | 0.048 | 0.077 | 0.119 | 197/197 | 0/197 |
| Logistic r = 3.5 | < 0 | NaN 30/30 | NaN 30/30 | NaN 30/30 | — | — |
| Constant | — | NaN 3/3 | NaN 3/3 | NaN 3/3 | — | — |

(Medians. The NaN reason is "insufficient divergence points" for logistic r = 3.5
and constant, and "analyze_segment raised" for the 128-sample failures.)

- **Positive LLE on white noise and AR(1) means a positive LLE alone is not
  evidence of chaos.** The chaotic and stochastic distributions overlap. At 256
  the IQRs are 0.109–0.151 (logistic), 0.103–0.137 (Hénon), 0.058–0.098 (white
  noise) and 0.050–0.098 (AR(1)).
- The estimates for the chaotic maps are far below their references (logistic
  0.12–0.14 vs 0.693; Hénon 0.12–0.15 vs 0.42). The LLE path embeds with the TDMI
  delay: median τ 6–15 samples for the chaotic maps (Table H, `LLE tau`), against
  a one-step map. `results/posthoc/tdmi.json` shows TDMI selects τ = 10 / 11
  (logistic, N = 512 / 4096) and 18 / 21 (Hénon). The post-hoc LLE diagnostic,
  regenerated in Phase 3 (Section 5.3), shows that τ = 1 with the true dimension
  removes most of the bias (logistic 0.69, Hénon 0.40–0.41).
- Floating-point sensitivity: a 1-ulp perturbation (`nextafter`) leaves Cao, UPO
  status and source peak count identical in 30/30. The LLE is identical in only
  17/30 (Hénon), 10/30 (logistic) and 30/30 (sinusoid)
  (`table_fp_sensitivity.md`).

### 3.3 Chaotic systems (logistic r = 4, Hénon)

Source: Table A (`one_sample`, 50 surrogates), `table_peak_locations.md`, Table I.

| System | Window | Status ok | ≥ 1 Level-A | ≥ 1 Level-B | r_J median | Level C gated | Verified unstable |
|---|---:|---|---|---|---:|---:|---:|
| Logistic | 128 | 18/30 | 18/30 | 17/30 = 0.567 | 4.99 | 0 | 0 |
| Logistic | 256 | 30/30 | 30/30 | 30/30 = 1.000 [0.886, 1.000] | 9.41 | 0 | 0 |
| Logistic | 512 | 30/30 | 30/30 | 30/30 = 1.000 | 9.54 | 0 | 0 |
| Hénon | 128 | 6/30 | 6/30 | 4/30 = 0.133 | 3.46 | 0 | 0 |
| Hénon | 256 | 30/30 | 30/30 | 29/30 = 0.967 [0.833, 0.994] | 6.07 | 0 | 0 |
| Hénon | 512 | 30/30 | 30/30 | 30/30 = 1.000 | 7.85 | 0 | 0 |

- **128-sample failures.** 10/30 logistic and 24/30 Hénon windows are
  `analysis_error`, with the message "Insufficient embedded points after
  reconstruction" (`results/core.jsonl`). This comes from the guard in
  `analyze_segment` that requires at least 50 LLE-embedded points
  (`final_pipeline.py` line 2351). It fires before the UPO result is returned.
  The Phase 2E diagnostic field `upo_direct_status` shows `run_upo_analysis`
  alone would have returned `ok` for 24/24 Hénon and 10/10 logistic (Table G).
  The same guard fails 5 white-noise and 3 AR(1) windows at 128, and 2
  RR-quantized Hénon windows in `fp_sens`.
- **Location, logistic.** Level-B peaks lie within 0.05 of a fixed point (0.75 or
  0) in 54/54 (256) and 58/58 (512) cases. The median location is split between
  about 0.001 and about 0.757 because both fixed points are detected. The table
  title lists only 0.75, but `analysis.py::peak_location_summary` uses
  {0.75, 0}.
- **Location, Hénon (Cao d = 4).** Level-B peaks lie within 0.05 of 0.631 in only
  12/29 (256) and 12/30 (512) cases. The per-peak detail in `results/core.jsonl`
  at 512 samples shows why. The bin width is 0.102 and a bin edge sits at about
  0.607, which is 0.024 (0.23 bin widths) below the fixed point. The peak falls
  in the correct cell in 12 windows (error +0.010 to +0.016) and in the adjacent
  lower cell in 17 windows (error −0.069 to −0.076, about 0.7 bin widths). One
  window with d = 8 is two cells away (−0.183). This is the same bin-quantization
  mechanism as finding 1.
- **Level C.** No chaotic peak passes the Level-C gates at production dimension.
  The residual gate fails for every peak (Table I, `C failed gates`). The post-hoc
  diagnostic `results/posthoc/levelc.json` (3 seeds, N = 512) shows the following.
  For Hénon at Cao m = 4 the pipeline residual is 0.115–0.689. At forced m = 2
  it is 0.027–0.050 and all gates pass in 3/3. For logistic at Cao m = 4 the
  residual is 6.9–11.9, and at forced m = 1 it is 0.063–0.073, still above 0.05.
  With oracle dimension, `verified_unstable` peaks appear at logistic x = 0
  (30 windows), logistic x ≈ 0.75 (1) and Hénon x ≈ 0.63 (7)
  (`results/oracle_dim.jsonl`).
- **Surrogate count** does not change chaotic detection. Logistic is 30/30 and
  Hénon 29/30 at n = 20, 50 and 100 (`table_surrogate_count.md`).

### 3.4 Periodic systems (sinusoid, logistic r = 3.5)

- **Sinusoid** (period 7.3 samples, exactly 73-periodic; the linear map has a
  neutral fixed point at the centre). The status is `ok` in 89/90 windows
  (1 `no_valid_transforms` at 128). At least one Level-B peak appears in
  11/30 = 0.367 (128), 21/30 = 0.700 (256) and 19/30 = 0.633 (512) windows
  (Table A). These rates are far above the null rates. Level-B peaks lie within
  0.05 of the centre in 0/15, 0/42 and 0/38 cases. The Level-A location IQR is
  [−0.201, 0.202] (`table_peak_locations.md`). Level C gates 0, 1 and 3 peaks, and
  verified_unstable is 0 at every length (Table I). **A deterministic periodic
  signal therefore produces surrogate-significant "UPO" peaks, not at its fixed
  point**, so Level B does not separate chaotic from periodic structure.
- **Logistic r = 3.5** (stable 4-cycle, 4 distinct values): 90/90 windows are
  `embedding_not_saturated` with Cao dimension 13, and there are 0 Level-A peaks
  (Table A). The LLE is NaN in 90/90 (Table H). The pipeline reports a detector
  failure, not "no UPO", for this input (see 5.1).

### 3.5 Noise robustness (Table D; 256 samples, seeds 0–29 paired with clean)

| SNR (dB) | Hénon ≥ 1 Level-B | Logistic ≥ 1 Level-B | Hénon Cao d median | Logistic Cao d median | Hénon r_J | Logistic r_J |
|---|---|---|---:|---:|---:|---:|
| clean | 29/30 = 0.967 | 30/30 = 1.000 | 4 | 4 | 6.07 | 9.41 |
| 30 | 11/30 = 0.367 | 29/30 = 0.967 | 7 | 8 | 2.10 | 3.93 |
| 20 | 8/30 = 0.267 | 27/30 = 0.900 | 8 | 9 | 1.53 | 3.44 |
| 10 | 2/30 = 0.067 | 22/30 = 0.733 | 8 | 9 | 1.26 | 2.69 |
| 5 | 3/30 = 0.100 | 13/30 = 0.433 | 8 | 8 | 0.87 | 1.69 |
| 0 | 1/30 = 0.033 | 3/30 = 0.100 | 9 | 8 | 0.97 | 1.08 |

- Hénon detection drops from 0.967 to 0.367 at 30 dB. From 10 dB down it is
  within the null range (Section 3.1). Logistic stays ≥ 0.90 down to 20 dB and
  reaches 0.10 at 0 dB.
- Median source coverage rises from 0.107 / 0.103 (clean) to 1.0 by 10 dB
  (Hénon) and 20 dB (logistic). This matches the nulls (Table B/C coverage 1.0).
  Coverage near 1 therefore indicates loss of structure, not better coverage.
- `embedding_not_saturated` appears in 0–4 of 30 noisy windows per SNR level
  (Table G). There are no `analysis_error` windows in the noise arm.
- The median LLE falls from 0.118 / 0.122 (clean) to 0.072 / 0.060 at 0 dB
  (Table H), into the white-noise range.

### 3.6 Hénon localization

Summarized in 2.1–2.3 (Table F; plots `10_henon_fixed_point_localization.png`,
`11_henon_period2_localization.png`). Additional points:

- Period-1 significance at the three one-factor sweeps through the production
  point: the on-attractor fixed point is Level B in 3/3 realizations in every
  assessed setting, for both systems (Table F, `fp1 significant`).
- The on-attractor fixed point is the strongest period-1 peak (rank 1) in 3/3 in
  all 121 settings for both systems.
- Stability audit: the hybrid source-stability estimate classifies the
  on-attractor peak as unstable in 12/12 detector runs. Disabling the stability
  diagnostic leaves the peaks identical in 12/12 detector runs and 9/9 pipeline
  runs, with no differing fields (`results/audits.json`, `stability_detector`,
  `stability_pipeline`).

### 3.7 Map mode (Table E) and floating-point / RR quantization

- `tau_step` (PROJECT EXTENSION) vs `one_sample` (SOURCE), 256 samples.
  Logistic Level-B drops from 30/30 to 0/30 and Hénon from 29/30 to 0/30.
  The null rates are 1/30 (white noise) and 1/30 (AR(1)), against 3/30 each for
  `one_sample`. The sinusoid drops from 21/30 to 12/30. The two modes give the same
  Level-B count in 0/30 (logistic) and 1/30 (Hénon) paired windows.
- RR quantization at 1/360 s (`table_fp_sensitivity.md`; rr = 0.8 + 0.05 z). It
  leaves a median of 37.5 / 51 / 61 distinct values in 256 samples (sinusoid /
  logistic / Hénon). Compared with the unquantized signal, embedding dimension is
  identical in 0/30, 5/30 and 9/30 windows, and UPO status in 0/30, 26/30 and 29/30.
- The two quantization routes differ by at most 2.2e-16 in a median of 76–83
  samples per window (`results/audits.json`, `rr_route_bit_differences`). Even so,
  embedding dimension is identical between the routes in only 10/30 (logistic) and
  15/30 (Hénon) windows. Source peak count is identical in 14/30 and 24/30, and the
  LLE in 3/30 and 4/30. On quantized data the results depend on the arithmetic
  route at the level of 1 ulp. Unquantized data are not affected: `nextafter`
  leaves Cao, status and peaks identical in 30/30.

## 4. Replicability and audits

All from `results/audits.json` and `results/replicability.json` unless stated.

| Audit | Result |
|---|---|
| Replicability rerun (Step 25; first 4 tasks of each of 11 experiments) | 44/44 row-identical and 44/44 detail-identical |
| Significance recomputation from stored surrogate W_i (1,378 windows, 2,911 peaks) | 0 mismatches in W, W0, r_J, J(W) and per-peak J; hand check passed; signed (not absolute) deviation confirmed |
| J(W0) | 0.50 in 1,351 windows; 0.48 / 0.49 / 0.45 in 27 (ties in W_i) |
| Status consistency violations | 0 |
| Harness errors | 0 |
| `analyze_segment` exceptions | 44, all "Insufficient embedded points after reconstruction" (42 core at 128 samples, 2 fp_sens) |
| Runtime warnings | 196 across 2,032 window rows |
| Flag independence (`compute_surrogates` False/True; 6 windows) | all UPO fields, W, W0, J(W), r_J and surrogate W_i identical in 6/6 |
| Stability independence | peaks identical without stability in 21/21 |
| Source integrity | `final_pipeline.py` and the four Phase 2E runtime sources hash-identical to `results/manifest.json` |

**Test suite.** Under Python 3.11.15 and, after the Phase 3 A2 follow-up, under
Python 3.13.12 (numpy 2.1.3, scipy 1.18.1, pandas 3.0.6; `requirements.txt`):
**334 passed, 1 failed**. The failing test is
`tests/test_upo_stability.py::test_hybrid_period1_stability_matches_analytic_henon_multipliers[1-prl_norm]`.
For the skewed-Hénon off-attractor peak (−1.903) the leading Lyapunov number is
17.24 against an analytic 2.09. The `pre_tensor` variant of the same case passes
(1.81). `HANDOFF.md` records 335/335 at `dd43fd0` under Python 3.13.15.
Diagnosis (A2, Section 5.3): the failure is a floating-point-path effect
amplified by one ill-conditioned local Jacobian, not a Python or scipy version
effect.

## 5. Findings about the analysis outputs

### 5.1 `embedding_not_saturated` is not a structure verdict

See 2.5. The label covers exact-duplicate (undefined E1), too-few-vectors and
genuinely non-saturated cases, and all of them yield `not_assessed_detector_failure`.
A stable periodic orbit (logistic r = 3.5) is therefore reported as a detector
failure in 90/90 windows.

**Phase 3 (A5) diagnostic.** The status is unchanged. The opt-in option
`PipelineConfig.upo_report_cao_diagnostics` (default False) adds these fields to
the `run_upo_analysis` output: `cao_e1_undefined`,
`cao_e1_undefined_duplicate_vectors`, `cao_e1_undefined_too_few_vectors`,
`cao_not_saturated`, and a per-dimension `cao_diagnostics` dict
(`cao_method(..., return_diagnostics=True)`; `tests/test_cao_e1_diagnostics.py`).
All 846 Phase 2E core windows of the six non-constant systems were recomputed
with Cao at lag 1 on the Phase 3 machine:

- all 90 logistic r = 3.5 `embedding_not_saturated` windows have E1 undefined
  because of duplicate vectors;
- all 22 white-noise, 22 AR(1) and 1 logistic `embedding_not_saturated` windows
  have E1 defined, so they are genuinely not saturated. The 22 AR(1) windows
  include the one hidden behind `analysis_error` at 128 samples.

### 5.2 Tables B and C: surrogate_sens row pooled two surrogate counts (fixed)

`analysis.py::null_table` grouped `surrogate_sens` only by `window_length`. That
file holds both the 20- and the 100-surrogate runs (50 white-noise and 50 AR(1)
windows each, seeds 0–49; `results/surrogate_sens.jsonl`). The original Tables B
and C (commit `88835d5`) therefore showed one pooled row, labelled 20 surrogates,
N = 100, nominal 0.0476: white noise 3/100 and AR(1) 7/100.

**Fix (Phase 3, item A1).** `null_table` now groups by `(window_length,
surrogates)`. After rerunning `python -m experiments.phase2e.analysis`, the only
content changes are in Tables B and C. Each pooled row is replaced by two rows:

| Table | n = 20 (nominal 0.0476) | n = 100 (nominal 0.0495) |
|---|---|---|
| B, white noise | 1/50 = 0.020 [0.004, 0.105]; assessed 44 | 2/50 = 0.040 [0.011, 0.135]; assessed 44 |
| C, AR(1) | 3/50 = 0.060 [0.021, 0.162]; assessed 45 | 4/50 = 0.080 [0.032, 0.188]; assessed 45 |

These match `table_surrogate_count.md` and Section 3.1. The core and algo_seed
rows are unchanged. The rerun also rewrote Table E (`.md` column-alignment marker
only), two CSVs (Table H and `table_coverage_vs_count`, last-digit float
differences at the 1e-16 relative level) and the PNG plots. These are
library-version effects, not content changes, so those files were restored to
their `88835d5` versions. `results/audits.json` was regenerated byte-identical.
No pipeline result is affected.

### 5.3 Items not verifiable from committed files

- **`results/posthoc/lle.json` (regenerated in Phase 3, A3).** It was absent at
  `88835d5`. It is written last by `posthoc_diagnostics.py::main`. It was
  regenerated by calling `posthoc_diagnostics.lle_diagnostics()` unchanged; the
  committed `tdmi.json` and `levelc.json` were not rewritten. The data come from
  the deterministic Phase 2E generators (`systems.generate`, same seeds). The
  run took 198 s on the Phase 3 machine, which is not the Phase 2E run machine
  (see the replicability item below). Results, as the median over seeds 0–4 of
  the Rosenstein LLE per sample:

  | System / N | production (TDMI τ*, Cao m* at τ*) | τ = 1, Cao lag-1 m | τ = 1, true m | τ*, true m | reference |
  |---|---:|---:|---:|---:|---:|
  | logistic 256 / 512 / 1024 / 2048 | 0.108 / 0.172 / 0.157 / 0.184 | 0.667 / 0.682 / 0.692 / 0.691 | 0.690 / 0.688 / 0.694 / 0.691 | 0.690 / 0.688 / 0.694 / 0.691 | 0.693 |
  | Hénon 256 / 512 / 1024 / 2048 | 0.132 / 0.167 / 0.174 / 0.186 | 0.379 / 0.382 / 0.385 / 0.387 | 0.411 / 0.404 / 0.407 / 0.404 | 0.385 / 0.474 / 0.515 / 0.553 | 0.419 |
  | white noise 256 / 512 / 1024 / 2048 | 0.060 / 0.094 / 0.105 / 0.077 | 0.110 / 0.122 / 0.118 / 0.147 | 0.422 / 0.474 / 0.527 / 0.577 (m = 2) | 0.319 / 0.401 / 0.451 / 0.458 | undefined |

  - The production LLE path underestimates by a factor of 3.8–6.4 (logistic)
    and 2.3–3.2 (Hénon). With τ = 1 and the true dimension the bias nearly
    vanishes (logistic within 0.005 of ln 2; Hénon within 0.015 of 0.419).
    The TDMI delay is the dominant cause for these one-step maps. The
    production dimension is 7–13 at N = 512 (`production_vs_tau1_N512_seeds0_9`).
  - At N = 2048 with τ = 1 and the true m, the production fit gives 0.695
    (logistic) and 0.404 (Hénon) for every Theiler window in {0, 1, 5, 20}
    (`theiler_fit_range_N2048_seed0`). For Hénon the first-K slope rises from
    0.30 (K = 2) to 0.41 (K = 10).
  - **White noise is not separated by bias correction alone.** With τ = 1 and m = 2
    its "LLE" is 0.42–0.58, above the Hénon value.
  - `lle.json` stores LLE medians, not the lag-1 Cao dimensions themselves. The
    handoff's Cao claim at N = 1024–2048 therefore remains unverified. At
    256–512 it is verified from `core.jsonl` (Section 2.4).
- **Stability-test failure (diagnosed in Phase 3, A2).** The test source is
  byte-identical to the one compiled when the test passed (compiled bytecode
  in `tests/__pycache__` at `dd43fd0`, same source size and constants). In this
  container the failing value (17.241) is identical under Python 3.11 and 3.13,
  scipy 1.14.1–1.18.1 and every forced OpenBLAS core type. Its cause:
  - the off-attractor peak averages 108 member Jacobians
    (`source_period1_stability`, arithmetic mean);
  - one member is an ill-conditioned M = 2 neighbour fit with Frobenius norm
    4,078; the next largest is 22.9;
  - that single outlier moves the leading modulus from about 2.0 (element-wise
    median 2.011, trimmed mean 1.977) to 17.24.

  Disabling numpy's AVX-512 dispatch (`NPY_DISABLE_CPU_FEATURES`) changes the
  ulp-level arithmetic. The member set then has 106 points and the largest member
  norm is 93.4 instead of 4,078. The mean gives 1.917 and the full suite passes
  335/335. The test outcome therefore depends on the
  CPU's floating-point path. The run machine for `dd43fd0` is not recorded, so
  which path it used is unknown. Neither the test nor the pipeline was changed.
- **Bitwise replicability across machines.** In this container,
  `run_phase2e --replicate` gives 13/44 bitwise-identical tasks against the stored
  records, with the unmodified pipeline and at every numpy/glibc SIMD level
  tried. The differing fields are at the 1e-15 relative level (LLE, R², extension
  moduli, orbit coordinates). They are amplified only in ill-conditioned
  quantities: hybrid source-stability Lyapunov numbers (up to 98 % relative),
  Level-C residuals (up to 0.8 %), and a few histogram counts and peak
  locations (up to 0.9 %). Reruns within one machine are bitwise identical
  (44/44). The Phase 2E replicability claim (Section 4) holds for the original
  machine only.

## 6. Limitations

- **`embedding_not_saturated` label conflation.** One status covers "E1 undefined"
  (duplicate vectors or too few vectors) and "E1 did not plateau". Downstream it
  becomes `not_assessed_detector_failure` with NaN features. It must not be read
  as either "no structure" or "high-dimensional".
- **RR quantization at 1/360 s (MIT-BIH).** MIT-BIH R-peak times are sampled at
  360 Hz, so RR intervals are integer multiples of 1/360 s. In the synthetic
  RR-style test (Section 3.7) this quantization:
  - turned every sinusoid window into `embedding_not_saturated`;
  - changed the Cao dimension in 25/30 logistic and 21/30 Hénon windows;
  - made embedding dimension, peak count and LLE depend on 1-ulp differences in
    how the quantized values are computed.

  On MIT-BIH, (a) a share of `embedding_not_saturated` segments is expected from
  duplicate delay vectors alone, (b) Cao, peak and LLE outputs should be treated
  as sensitive to arithmetic route, and (c) the quantization step (2.8 ms)
  relative to RR variability sets how many distinct values a segment has. None
  of this was measured on real RR data.
- **Null false-positive rate.** Empirical rates are 0.02–0.11 per cell and 0.070 /
  0.075 pooled at n = 50. They are consistent with the nominal sizes within
  binomial uncertainty in every cell except AR(1) at 256 (11/100). The cells have
  only 50–100 windows. The shared random streams
  (Section 2.8) make windows not fully independent. Rates among assessed windows
  are higher than among all windows.
- **Level B is not specific to chaos.** A periodic sinusoid gives Level-B peaks in
  37–70 % of windows (3.4), and Level A fires on 70–98 % of noise windows (3.1).
- **LLE.** It is positive for white noise, AR(1) and the periodic sinusoid, and
  well below the reference values for the chaotic maps. At 128 samples the sinusoid
  gives large positive LLEs with good fits. The LLE also changes under a 1-ulp
  perturbation in 13/30 (Hénon) and 20/30 (logistic) windows.
- **Short windows.** At 128 samples the LLE-embedding guard in `analyze_segment`
  (at least 50 points) fails 24/30 Hénon and 10/30 logistic windows. It discards
  a UPO result that `run_upo_analysis` would have returned. This is recorded
  behaviour of the production wrapper, reproduced by `metrics.py::error_fields`,
  not a Phase 2E harness error. Phase 3 (A4) added the opt-in option
  `PipelineConfig.keep_upo_on_short_lle_embedding` (default False = unchanged).
  With it on, all 42 core 128-sample `analysis_error` windows return the UPO
  status that `run_upo_analysis` gives alone (41 `ok`, 1
  `embedding_not_saturated`). Only the LLE is marked failed
  (`tests/test_lle_short_embedding.py`).
- **Cao at lag 1 overestimates dimension** (Hénon 4, logistic 4). Level C
  therefore gates no chaotic peak at production dimension (3.3).
- **Localization.** The bin-mean location estimator has errors up to 0.27 bin
  widths on the production sweep and 0.50 over the full grid (Table F). In the
  core arm about 0.7 bin widths occurs when the fixed point sits near a bin edge
  (3.3).
- **Coverage** rises to about 1 for noise and noisy chaos (3.5). It is not a
  monotone measure of structure.
- **Synthetic only.** No statement here transfers to MIT-BIH without a separate
  validation.
- **Not re-examined in Phase 2E:** the Pan-Tompkins detector, RR correction, TDA
  and classifier logic.

## 7. Reproducibility

- Raw results: `results/*.jsonl` / `*.csv`; tables: `results/tables/`; audits:
  `results/audits.json`, `results/replicability.json`; post-hoc diagnostics:
  `results/posthoc/`; figures: `plots/`.
- Reproduction commands: `experiments/phase2e/README.md`. The runner and
  analysis were not rerun for this report. All numbers are read from the files
  committed in `88835d5`. The binomial tail probabilities in 3.1 and the per-peak
  location breakdown in 3.3 were computed directly from those files.
- Later corrections (Phase 3, Part A) are marked in the text: Section 5.2
  (Tables B/C fix, analysis rerun); Sections 4 and 5.3 (test failure and
  cross-machine replicability, A2); Sections 3.2 and 5.3 (`lle.json`
  regenerated, A3).
- `experiments/phase2e/HANDOFF.md` is superseded by this report. Its progress
  section and resume instructions no longer apply.
