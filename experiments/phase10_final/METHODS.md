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
