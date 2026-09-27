# Phase 3 preregistration: LLE-based chaos test

Written and committed **before any test-seed run**. After this file is pushed,
the candidates, their parameters and the decision rule are frozen. Any change
goes in a dated amendment at the end of this file that says why. The runner
(`run_phase3.py --phase test`) refuses to start unless this file is committed,
unmodified and pushed, and it only runs the methods listed below with
`method:` lines.

## 1. Data

- **Systems** (Phase 2E generators and seed scheme,
  `experiments/phase2e/systems.py::generate`):
  - `white_noise`, `ar1` (φ = 0.8): nulls;
  - `sinusoid` (period 7.3), `logistic_p4` (r = 3.5): periodic;
  - `logistic` (r = 4), `henon` (a = 1.4, b = 0.3): chaotic;
  - `logistic` and `henon` with white Gaussian observation noise at SNR 30, 20
    and 10 dB.
- **Window lengths:** 256 (primary) and 512.
- **Test seeds:** 1000–1099 for every condition, i.e. 100 windows per
  (condition, length) and 2,400 windows per method. Data seed =
  `SeedSequence([20260926, system_code, window_length, seed])`; noise from the
  Phase 2E child stream per (seed, SNR).
- **Development seeds** (Phase 2E core/noise seeds) were used for all tuning; see
  Section 5.

## 2. Methods (all parameters fixed)

Common to C1–C4:
- τ = 1 delay embedding (`fp._embed_backward`);
- Theiler window = `fp._mean_period_beats(x)` (as in production);
- **99 IAAFT surrogates** (Schreiber & Schmitz 1996; at most 200 iterations,
  stop when the rank order is unchanged; `estimators.iaaft_surrogate`) from a
  window-specific RNG `SeedSequence([20260927, crc32(x bytes), 2])`;
- test statistic computed identically on the data and on every surrogate;
- one-sided p = (1 + #{surrogate statistic ≥ observed, or undefined}) / 100;
- **α = 0.05**; detected ⇔ p ≤ 0.05, i.e. at most 4 of the 99 surrogates reach
  the observed value. Under exchangeability the size is exactly 5/100 = 0.05;
- a window whose observed statistic is undefined (NaN) is "not detected";
- the LLE estimate reported is the observed statistic (natural log per sample).

The divergence curve y[k], k = 0..39, is the mean log distance between each
point and its neighbour(s) k steps later. It uses the same Theiler exclusion,
1e-15 floor and count rule (≥ max(10, n/4) finite terms) as
`fp.rosenstein_lle`. With one neighbour it equals `fp.rosenstein_lle`'s curve to
1e-12 (`tests/test_phase3_estimators.py`).

- method: `c1_rosenstein_m2_iaaft`
  - C1: m = 2; single nearest neighbour (Rosenstein); statistic = least-squares
    slope of y[1..5].
- method: `c2_kantz_m3_sat_iaaft`
  - C2: m = 3; y[k] = mean over points of log(mean distance to its 5 nearest
    admissible neighbours, k steps later) (Kantz-style);
  - statistic = slope of y[1..k_end]. k_end = first k in 1..12 with
    y[k] ≥ y[1] + 0.3·(plateau − y[1]), plateau = mean of y[25..39], clipped
    to [3, 12]. If plateau ≤ y[1], k_end = 3.
- method: `c3_eps_m2_iaaft`
  - C3: m = 2; single nearest neighbour among points farther than ε, where
    ε = the median over points of the nearest-neighbour distance (noise-scale
    exclusion); statistic = slope of y[1..5].
- method: `c4_zero_one_iaaft`
  - C4 (independent comparison, not LLE-based): 0–1 test for chaos,
    modified-correlation K (Gottwald & Melbourne 2009); 100 frequencies
    c ~ U(π/5, 4π/5) from `default_rng(0)`, n ≤ N/10; the same one-sided
    surrogate rule on K. It produces no LLE estimate.
- method: `baseline`
  - The existing pipeline, evaluated for reference and **not eligible to win**:
    TDMI τ, Cao m at τ, `fp.rosenstein_lle` (production early fit),
    `fp.run_surrogate_analysis` (199 AAFT, rng = random_seed + 7919);
  - detected ⇔ `lle_p_upper` ≤ 0.05 (size 10/200 = 0.05).

## 3. Decision rule (PRIMARY, fixed in advance)

All on **test seeds at 256 samples**. Rate = detected / all 100 windows; failed
windows count as not detected.

1. **PASS.** A candidate (C1–C4) passes if the false-positive rate on
   `white_noise` and on `ar1` each has a 95 % Wilson interval that contains
   α = 0.05 or lies entirely below it (equivalently, the lower Wilson bound is
   ≤ 0.05; with n = 100 this means at most 9 detections).
2. **WINNER.** Among passing candidates, the winner has the highest detection
   rate on noisy Hénon plus noisy logistic at 20 dB, pooled (k out of 200
   windows).
3. **Tie-break** (equal pooled count): the lowest absolute LLE bias on clean
   logistic plus clean Hénon, defined as
   |median LLE − ln 2| + |median LLE − 0.4192|, over windows with a finite
   estimate. A candidate without an LLE estimate (C4) ranks last in the
   tie-break. If still tied, prefer the lower pooled false-positive count on
   `white_noise` + `ar1`, then the lower method index (C1 < C2 < C3).
4. If no candidate passes, there is no winner, and this is reported plainly.
5. **Adoption (B5).** The winner, if any, is added to `final_pipeline.py` as an
   opt-in `PipelineConfig` option, default off, which reproduces the estimator
   in `estimators.py` exactly (checked by a unit test). No default changes.
6. **Frozen code.** The methods are the functions in `estimators.py` at the
   commit that adds this file. The runner records their SHA-256 in
   `results/test/manifest.jsonl`. `estimators.py` is not edited between that
   commit and the end of B4.

## 4. Secondary outcomes (reported for every method, including the baseline)

- Rejection (false-positive) rate on `sinusoid` and `logistic_p4`.
- Detection rate at 30 and 10 dB (each system and pooled).
- Clean detection rate and LLE bias (median − reference; also the IQR) for
  logistic and Hénon; bias at 30 / 20 / 10 dB.
- All of the above at 512 samples.
- Runtime per window.
- Failure counts and statuses.

## 5. Development evidence behind the choices (DEVELOPMENT seeds only)

- **B2 diagnosis** (`results/dev/tables/B2_diagnosis.md`). The production LLE
  bias comes from:
  - the TDMI delay: logistic −0.571 → −0.400 when τ* → 1 at the same m;
  - the Cao dimension at that delay: → −0.022 with the lag-1 Cao m, ≈ 0 with the
    true m;
  - the resulting early saturation (k95 = 2 in 30/30 production curves).

  At τ = 1 noise adds a large jump at k = 0 → 1 (white noise 2.0; logistic
  20 dB 3.9 vs 0.69 clean). This is why every LLE candidate starts its fit at
  k = 1. Cao inflates m to 7–9 under noise, so m is fixed (2 or 3), not
  estimated.
- **B3 tuning** (`results/dev/tuning.jsonl`, 39 surrogates, 256 samples). Grid:
  m {2, 3} × neighbour rule {1-NN, ε-exclusion, 5-NN} × surrogate {AAFT, IAAFT}
  × fit {0..5, 1..5, saturation 50 %, 30 %}, plus the 0–1 test.
  - IAAFT kept AR(1) closer to nominal than AAFT (for example m2/5-NN/1..5:
    AR(1) 0.21 with AAFT vs 0.05 with IAAFT).
  - The 0–1 test had power ≤ 0.07 on every chaotic condition. It is kept only
    as the requested independent comparison.
  - The candidates were chosen for low development null rates, power, and low
    clean / noisy bias.
- **Development results of the frozen candidates** (99 surrogates, all dev seeds,
  256 and 512): see Section 7, added before commit.

## 6. Runtime estimate (single process, 256 / 512 samples, measured on the Phase 3 machine)

| method | s per window | CPU-hours for 2,400 test windows |
|---|---|---|
| baseline | 5–9 / 16–19 | about 8.2 |
| C1 | 1.1 / 3.2 | about 1.4 |
| C2 | 1.8 / 4.5 | about 2.1 |
| C3 | 1.5 / 4.1 | about 1.9 |
| C4 | 4.1 / 7.0 | about 3.7 |

The total is about 17 CPU-hours: 4.5–6 hours wall time on 4 workers, since
development runs showed contention (baseline median 14 s / 32 s per window).
99 surrogates give an exact size of 0.05 with a p-value resolution of 0.01. 19
surrogates also give exact size 0.05, but only p ∈ {0.05, 0.10, …}: detection
then needs every surrogate below the observed value, so power is lower. 199
would double the candidate cost.

## 7. Development results of the frozen methods

From `results/dev/tables/summary.md` (development seeds only; 99 surrogates; nulls 100 windows at 256 and 50 at 512, other conditions 30).

**256 samples, detections / windows** (development seeds)

| method | WN | ar1 | sinusoid | log r3.5 | logistic | henon | logistic@30dB | henon@30dB | logistic@20dB | henon@20dB | logistic@10dB | henon@10dB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| baseline | 7/100 | 4/100 | 0/30 | 0/30 | 14/30 | 12/30 | 16/30 | 15/30 | 19/30 | 14/30 | 8/30 | 16/30 |
| c1_rosenstein_m2_iaaft | 1/100 | 5/100 | 0/30 | 0/30 | 30/30 | 30/30 | 30/30 | 30/30 | 30/30 | 30/30 | 30/30 | 30/30 |
| c2_kantz_m3_sat_iaaft | 6/100 | 0/100 | 0/30 | 0/30 | 30/30 | 30/30 | 30/30 | 30/30 | 30/30 | 30/30 | 30/30 | 30/30 |
| c3_eps_m2_iaaft | 6/100 | 5/100 | 0/30 | 0/30 | 30/30 | 30/30 | 30/30 | 30/30 | 30/30 | 30/30 | 30/30 | 30/30 |
| c4_zero_one_iaaft | 3/100 | 1/100 | 0/30 | 0/30 | 0/30 | 1/30 | 0/30 | 1/30 | 0/30 | 1/30 | 0/30 | 2/30 |

**256 samples, median LLE bias** (median − reference)

| method | logistic | henon | logistic@30dB | henon@30dB | logistic@20dB | henon@20dB | logistic@10dB | henon@10dB |
|---|---|---|---|---|---|---|---|---|
| baseline | -0.571 | -0.301 | -0.577 | -0.305 | -0.576 | -0.301 | -0.642 | -0.328 |
| c1_rosenstein_m2_iaaft | -0.002 | -0.002 | -0.100 | -0.054 | -0.203 | -0.065 | -0.424 | -0.176 |
| c2_kantz_m3_sat_iaaft | +0.015 | +0.018 | -0.040 | -0.011 | -0.095 | -0.028 | -0.245 | -0.045 |
| c3_eps_m2_iaaft | -0.007 | -0.004 | -0.124 | -0.062 | -0.224 | -0.078 | -0.433 | -0.185 |

**512 samples, detections / windows** (development seeds)

| method | WN | ar1 | sinusoid | log r3.5 | logistic | henon | logistic@30dB | henon@30dB | logistic@20dB | henon@20dB | logistic@10dB | henon@10dB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| baseline | 3/50 | 1/50 | 0/30 | 0/30 | 30/30 | 27/30 | 30/30 | 27/30 | 30/30 | 27/30 | 10/30 | 22/30 |
| c1_rosenstein_m2_iaaft | 4/50 | 3/50 | 0/30 | 0/30 | 30/30 | 30/30 | 30/30 | 30/30 | 30/30 | 30/30 | 30/30 | 30/30 |
| c2_kantz_m3_sat_iaaft | 1/50 | 0/50 | 0/30 | 0/30 | 30/30 | 30/30 | 30/30 | 30/30 | 30/30 | 30/30 | 30/30 | 30/30 |
| c3_eps_m2_iaaft | 3/50 | 2/50 | 0/30 | 0/30 | 30/30 | 30/30 | 30/30 | 30/30 | 30/30 | 30/30 | 30/30 | 30/30 |
| c4_zero_one_iaaft | 2/50 | 4/50 | 0/30 | 0/30 | 6/30 | 4/30 | 4/30 | 3/30 | 4/30 | 5/30 | 1/30 | 2/30 |

**512 samples, median LLE bias** (median − reference)

| method | logistic | henon | logistic@30dB | henon@30dB | logistic@20dB | henon@20dB | logistic@10dB | henon@10dB |
|---|---|---|---|---|---|---|---|---|
| baseline | -0.552 | -0.269 | -0.542 | -0.277 | -0.548 | -0.271 | -0.628 | -0.293 |
| c1_rosenstein_m2_iaaft | -0.003 | -0.000 | -0.090 | -0.052 | -0.193 | -0.057 | -0.419 | -0.169 |
| c2_kantz_m3_sat_iaaft | +0.007 | +0.012 | -0.055 | -0.030 | -0.074 | -0.000 | -0.215 | -0.008 |
| c3_eps_m2_iaaft | -0.005 | -0.004 | -0.099 | -0.062 | -0.204 | -0.064 | -0.427 | -0.175 |

On development seeds every LLE candidate detected all chaotic windows, so the primary winner on test seeds is expected to be decided by the tie-break. This was known when the rule was fixed; the rule is the one specified in the task and was not altered.

## 8. Amendments

(none)

