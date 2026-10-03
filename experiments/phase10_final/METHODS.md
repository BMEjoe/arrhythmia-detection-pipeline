# Phase 10 methods and sources

Every method is implemented from a primary published source (or, where the primary is not
accessible, an open restatement, stated as such). Every implementation choice and deviation is
recorded here. Sections are filled in as each method is built.

## 1. Detectors (frozen; no retuning)
- K1 `k1_frozen512`: `experiments.phase7_mitbih.detector.evaluate` (Phase 6 Section 5 combined
  detector, `keep_upo_on_short_lle_embedding = True`, D2 detrend). Components LLE alone, UPO alone.
- K3: `final_pipeline.masked_growth_chaos_test(rr, labels, CFG)`.

## 2. Noise titration (Question 2), `titration10.py`

### 2.1 Sources (accessed 2026-10-03; copies or extracted text in `sources/` where licences allow)
| key | source | access | used for |
|---|---|---|---|
| [PB01] | Poon CS, Barahona M. Titration of chaos with added noise. PNAS 98:7107-7112 (2001), PMC34630 | full text read through the PMC page (PMC XML and PDF blocked: publisher / CAPTCHA) | titration algorithm; Fig. 1 legend (1,000-point series, unit-variance white noise xi, y + alpha xi, bisection for the largest alpha still nonlinear at the 1 % level, NL = alpha / sigma_y, 10 noise realisations averaged); Fig. 2 values (logistic r = 3.7: NL ~ 75 %, r = 3.575: ~ 9 %, r = 3.565 (periodic): ~ 0 %); Figs. 3-5 qualitative results (Mackey-Glass, Lorenz, NL = 0 in periodic windows) |
| [Wu09] | Wu GQ et al. Chaotic signatures of heart rate variability and its power spectrum in health, aging and heart failure. PLoS ONE 4:e4323 (2009) | open (XML, PDF, equation images `sources/wu2009_e001.png`, `_e002.png`) | Eq. 1 (VWK model, term order), Eq. 2 C(r) = log eps(r) + r/N with r = number of LEADING terms; F-test at 1 %; NL in 5-10 repetitions; HRV data, preprocessing and 12-min segments (Section 2.6) |
| [Wy06] | Wysocki M et al. Chaotic dynamics of resting ventilatory flow in humans assessed through noise titration. arXiv:nlin/0606032 (Respir Physiol Neurobiol 153:54, 2006) | open | best linear model = d = 1, kappa minimising the Akaike criterion; best nonlinear d > 1; F-test and Whitney-Mann at 1 %; routine defaults kappa = 6, d = 3; Fig. 3: C(r) curves of linear and nonlinear models plotted on the same r axis |
| [P10] | Poon CS, Li C, Wu GQ. A unified theory of chaos ... arXiv:1004.1427 | open | same restatement of Eqs. 1-2 as [Wu09] |
| [BP96] | Barahona M, Poon CS. Detection of nonlinear dynamics in short, noisy time series. Nature 381:215 (1996) | **not accessible** (nature.com login redirect); only the abstract | primary definition of the indicator; its details are taken from [PB01], [Wu09], [Wy06], [P10] |
| [PM97] | Poon CS, Merrill CK. Decrease of cardiac chaos in congestive heart failure. Nature 389:492 (1997) | **not accessible**; PubMed abstract only | the CHF claim; PhysioNet's chfdb page lists it as a user of chfdb |

### 2.2 Implementation (as published)
- Model (Eq. 1 [Wu09]): y_n = a_0 + sum_{j<=kappa} a_j y_{n-j} + all products of the lags up to degree d, ordered
  by degree and lexicographically in the lag indices; M = (kappa + d)! / (kappa! d!) terms.
- Coefficients: least squares on every leading set of r terms, by modified Gram-Schmidt in the term order
  (Korenberg's recursive orthogonal estimation; numerically dependent terms add nothing).
- eps(r): normalised RMS one-step prediction error, eps^2 = RSS(r) / sum (y_n - mean)^2.
- C(r) = log eps(r) + r / N (Eq. 2), N = number of fitted prediction equations.
- Best linear model: minimum of C over the linear family (d = 1). Best nonlinear model: minimum of C over the
  leading truncations that contain at least one nonlinear term, for kappa = 1..6 and d up to 3 ([Wy06] defaults).
- Nonlinear iff C_nl < C_lin and the F-test at the 1 % level rejects linearity.
- Titration ([PB01] Fig. 1 legend): xi unit-variance white Gaussian (standardised sample), bisection on alpha in
  units of sigma_y (bracket [0, 2], doubled while still nonlinear, cap 16; 10 halvings, resolution 0.2 % of
  sigma_y); NL_i = largest alpha still nonlinear; NL = mean of 10 realisations; NL = 0 if y itself is not
  nonlinear. Positive iff NL > 0.
- Series standardised before the fit (the span of every leading term set is invariant to affine changes of y,
  because translation of a degree-k monomial produces only lower-degree terms, all of which precede it).
- All models are fitted on the same target rows n = K .. N-1, K = the longest lag of any candidate model.

### 2.3 Implementation choices where the accessible sources are ambiguous (each disclosed)
1. **eps is the RMS (square root of the normalised residual variance)**, not the variance: only then is
   Eq. 2 Akaike's criterion (AIC / 2N + const), the name used by [PB01], [Wy06], [Wu09]. Alternative recorded.
2. **F-test = one-sided variance-ratio F-test of the two models' residuals**, df (N - r_lin, N - r_nl):
   [PB01], [Wy06], [P10] pair the F-test with a Whitney-Mann (Mann-Whitney) test, which compares two residual
   samples, so the F-test is read as its parametric counterpart. The nested partial F-test is recorded as an
   alternative.
3. **Linear memory.** [Wy06]'s default kappa = 6 is applied to the nonlinear family. For the linear family the
   memory is allowed up to M - 1 = 83 lags (as many terms as the largest nonlinear model). Reason: [PB01]
   states that periodic and quasi-periodic signals give NL ~ 0 "because a nonlinear limit cycle can be described
   both by nonlinear models and by linear models of enough memory" (memory 2M for M harmonics, 3N for an
   N-torus), and [Wy06] Fig. 3 plots the linear and nonlinear C(r) curves over the same r range. With linear
   kappa <= 6 the periodic example of [PB01] Fig. 2 is NOT reproduced (Section 2.4).
4. Masked variant (Question 2d): a row is used only if its whole span (target and every lag up to K) is free of
   ectopy-related intervals (the K3 rule); >= 300 usable rows required, else not analysable.

### 2.4 Readings tried on the three numeric values of [PB01] Fig. 2 (logistic map, 1,000 points after 1,000
transient iterations from x0 = 0.4; 5 noise realisations, rng seed 1). These values were USED to choose
between readings 2-3 (they are not an independent check).

| reading | r = 3.7 (pub. ~0.75) | r = 3.575 (pub. ~0.09) | r = 3.565 periodic (pub. ~0) |
|---|---|---|---|
| linear kappa <= 6, F variance ratio | 0.666 | 0.631 | **0.638** (fails) |
| linear kappa <= 6, nested F | 1.364 | 1.412 | 1.418 (fails) |
| linear kappa <= 10 / 16 / 20 / 30, F variance ratio | 0.646 / 0.643 / 0.642 / 0.641 | 0.375 / 0.170 / 0.156 / 0.112 | 0.348 / 0.104 / 0.092 / 0.036 |
| linear kappa <= 50, F variance ratio | 0.641 | 0.092 | 0.000 |
| **linear kappa <= 83, F variance ratio (PRIMARY)** | **0.637** | **0.088** | **0.000** |
| linear kappa <= 83, nested F | 1.298 | 0.161 | 0.009 |
| linear kappa <= 83, Mann-Whitney on squared residuals | 0.695 | 0.095 | 0.000 |

### 2.5 Independent verification (criteria fixed here BEFORE running `verify_titration.py`)
Primary reading only; 10 noise realisations; 1,000-point series; seeds fixed in the script.
- V1-V3 (Fig. 2, not independent): |NL(3.7) - 0.75| <= 0.15; 0.03 <= NL(3.575) <= 0.20; NL(3.565) <= 0.02.
- V4 (Fig. 1, logistic bifurcation scan r = 3.50, 3.505, ..., 4.00; LE from the equation, 2e5 iterations):
  NL > 0 in >= 90 % of r with LE > 0.02; NL = 0 in >= 90 % of r with LE < -0.02; Spearman(NL, LE) >= 0.6 over
  r with LE > 0.02.
- V5 (non-chaotic controls, 20 series each): white Gaussian noise and Gaussian AR(2) (y_n = 1.5 y_{n-1} -
  0.75 y_{n-2} + e_n): positive in <= 2/20 each; periodic signals (3-harmonic limit cycle + 1 % white noise;
  logistic r = 3.83 period 3) and the quasi-periodic 2-torus y_k = prod_i [2 + sin(w_i k)], w = (1, sqrt 2),
  + 1 % noise: NL = 0 in >= 18/20 each.
- V6 (flows, [PB01] Figs. 3c, 5): Lorenz (sigma 10, b 8/3) x sampled every 0.075: r = 28 (chaotic) NL > 0;
  r = 160 (limit cycle) NL = 0. Mackey-Glass (beta 0.2, gamma 0.1, n 10; Phase 8 verified integrator and
  ground-truth labels) sampled every 12.5: chaotic tau = 17, 23, 30: NL > 0 in >= 2 of 3; periodic tau = 14, 15,
  16: NL = 0 in >= 2 of 3. Henon (1.4, 0.3) x: NL > 0.
- **VERIFIED** iff V1-V6 all hold. Otherwise titration claims are not made (Question 2 reports the failure).

### 2.6 Verification result (`verify_titration.py`, `results/verification/titration_verification.{json,md}`; 59 s)
**VERIFIED (V1-V6 all pass).**
- V1-V3 (used to choose the reading): NL = 0.637 (r = 3.7; published ~0.75), 0.087 (r = 3.575; ~0.09),
  0.000 (r = 3.565 periodic; ~0).
- V4 logistic scan: NL > 0 at 78/78 r with LE > 0.02; NL = 0 at 18/18 r with LE < -0.02; Spearman(NL, LE) 0.85.
- V5 controls (20 each): white noise 0/20 positive, AR(2) 0/20, limit cycle NL = 0 in 20/20, 2-torus 20/20,
  logistic period 3 20/20.
- V6: Lorenz r = 28 NL 0.61, r = 160 NL 0; Henon NL 0.77; Mackey-Glass chaotic tau 17 / 23 / 30: NL 0.33 / 0.68 /
  0.66; periodic tau 14 / 15 / 16: NL 0 / 0 / 0.
- Remaining difference: r = 3.7 gives 0.64 vs the published ~0.75 (a value read from a figure). The
  Mann-Whitney reading gives 0.70. Titration claims are made with the primary reading.

### 2.7 HRV studies, data and preprocessing (Question 2c)
- [Wu09] (the "transient chaos" interpretation tied to ectopic beats): young group n = 13 (32 +/- 8 y) and CHF
  n = 14 (NYHA III-IV, no beta-blockers) "from the PhysioNet database"; elderly n = 16 from the authors' laboratory.
  These match nsrdb (18 subjects, 20-50 y) and chfdb (15 subjects, NYHA 3-4), i.e. this project's DEVELOPMENT
  data. Subjects "selected on the basis of stability of the mean heart rate and limited number of ectopic beats
  and undetected beats" (not reproducible). RR: PhysioNet software with "elimination (without interpolation) of
  premature or missing beats and other ectopic beats", then manual removal of residual premature beats with
  abnormal QRS (not reproducible without the ECG; nsr2db / chf2db have no ECG). Segments: 24-h data divided into
  120 segments of 12 min (~800 beats). Outcomes: nonlinear detection rate DR (% of segments with nonlinearity
  detected, per 3-h window) and NL (mean over detected segments only). Kappa and d are not reported; the
  routine defaults of [Wy06] (kappa 6, d 3) are used. Fig. 2C: a CHF segment with NL 134.2 % fell to 0 % after
  manual removal of RR "spikes" (undetected ectopic beats, compensatory pauses); Fig. 2D: 12 segments in 7 CHF
  subjects.
- [PM97] ("decrease of cardiac chaos in CHF"): healthy subjects and severe CHF; PhysioNet lists it as a user of
  chfdb. Window length and preprocessing not accessible; [Wu09] states its 12-min / ~800-beat segmentation is
  "as with previous studies [PM97]".
- Neither study used nsr2db or chf2db (the CONFIRMATORY databases): chf2db is NYHA I-III (Columbia), nsr2db
  healthy (Washington University / Columbia), both 128 Hz annotations with manual review.

## 3. Data handling (`data10.py`, `harness10.py`)
- Downloads: PhysioNet `files/<db>/1.0.0/` RECORDS, headers and the beat-annotation file only; every file's
  SHA-256 equals PhysioNet's SHA256SUMS.txt (`MANIFEST.json` per database; data gitignored).
  Development: mitdb (atr, 48 records), nsrdb (atr, 18), chfdb (ecg, 15).
- Beats and labels: the Phase 9 Part F rule (wfdb beat symbols; 'N' iff N L R B e j; an interval < 0.25 s or
  > 2.5 s marks its ending beat non-normal). Ectopy-related interval = the K3 rule = the Phase 7 editing rule.
- Development finding that shaped the window rule: nsrdb contains long unannotated stretches (intervals up to
  207 s; 407 intervals > 3 s in record 16272) and chfdb spurious short intervals (< 0.25 s, up to 70 in chf03).
  These are data gaps / annotation artefacts, so windows are moved past intervals outside the pipeline's hard
  limits [0.25, 3.0] s, and 12-min segments drop such intervals and require >= 90 % coverage.
- Editing: Phase 6 `edit_nn` (linear interpolation by beat index) with the K3 mask; eligible iff NN
  fraction >= 0.80 (Phases 6-7).
- Clock time: header start time + elapsed time. Development start times (33 long-term records) range
  08:00-14:35, median 09:54 (used only as the predeclared fallback when a header has no start time).

## 4. Question 1 spike-ins
- Lyapunov exponents of the maps from the equations (`/tmp` scan reproduced in the preregistration table):
  logistic mean log|r (1 - 2x)| over 4e5 iterations after 2,000; Henon (b = 0.3) Benettin / QR on the
  Jacobian [[-2 a x, 1], [b, 0]] over 4e5 iterations. Grid chosen where the exponent stays > 0 at +/- 0.001
  and +/- 0.002 of the parameter (avoids periodic windows): logistic r = 3.58 (0.105), 3.65 (0.255),
  3.88 (0.464), 4.00 (0.693); Henon a = 1.08 (0.136), 1.14 (0.244), 1.22 (0.303), 1.40 (0.419).
- Phase 8 models (`experiments.phase8_cardiac.series.model_rr`, verified in Phase 8) with their ground-truth
  exponents (Phase 8 `regimes.json`): coupled vdP chaotic regimes (all four) and four of the seven chaotic
  phase-resetting regimes spanning their exponent range.
- Additive design: x' = x + s c, s^2 = f / (1 - f) var_NN (independent-sum definition of the chaos fraction);
  beat times rounded to the record's annotation resolution. Replacement design: Phase 7 `spike_in`
  (`groups`, `insert_ectopy`) reused unchanged.
- Curves: Firth (1993) penalised logistic regression (Jeffreys prior; modified score
  X'(y - p + h (1/2 - p)), h = hat-matrix diagonal), finite under complete separation; subject-cluster
  bootstrap bands.
- Bound: U = max(cluster-bootstrap 95th percentile, one-sided Clopper-Pearson 95 %); pi_upper = U / p.

## 5. Question 3 model
- GEE (Liang & Zeger 1986) as implemented in statsmodels 0.15.0: binomial, logit, exchangeable working
  correlation, robust sandwich covariance; predictor log2(1 + burden) and CHF indicator. Not fitted with
  < 10 positives or negatives (sparse-outcome rule).
