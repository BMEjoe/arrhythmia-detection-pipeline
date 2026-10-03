# Conclusions checklist (Phase 10 Section 12)

A record for the author: each conclusion of `docs/PHASE10_FINAL_RESULTS.md` Section 12, with its precise scope,
the results that support it, and the results that limit it. Wording of the conclusions is as in the report
(after the final-refinement corrections marked `[Erratum En]`). Check ids refer to
`experiments/final_refinement/VERIFICATION.md`; tables to `experiments/final_refinement/tables/`.

**Status flags.** ✔ = supported within the stated scope; ⚠ = supported, but a reported result limits it in a way
the conclusion's wording does not state; ✖ = a reported result contradicts part of the wording (author's decision
required; see `ERRATA.md`).

Detector names: **LLE alone** = Phase 3 `lle_chaos_test` (τ = 1, m = 2 Rosenstein slope vs 99 IAAFT surrogates);
**UPO alone** = Phase 4 C2 instability-gated UPO test; **K1** = LLE AND UPO with the Phase 6 D2 detrend
(`k1_frozen512` at 512 intervals; tested at 256 in Phases 5-7); **K3** = annotation-masked nonlinear-prediction
determinism test with forecast-error growth gate (`final_pipeline.masked_growth_chaos_test`); **titration** = the
verified Poon-Barahona noise titration of Phase 10 (the Phase 8 C4 version was approximate).

---

### 1. A positive largest Lyapunov exponent from the production Rosenstein estimator is not evidence of chaos. ✔
- **Scope.** The production pipeline LLE (TDMI delay, Cao dimension, early fit); synthetic series of 128-512
  samples; no ectopy.
- **Supports.** Positive in 194/195 white-noise and 197/197 AR(1) windows (2E-1, 2E-2); underestimates known
  exponents (Phase 2E Table H; Phase 3 baseline bias −0.567 logistic, −0.302 Hénon at 256, check 3-7).
- **Limits.** Concerns the production estimator only. The Phase 3 C1 slope is nearly unbiased on clean maps
  (−0.001 / −0.005, check 3-5) but biased under noise (−0.213 logistic at 20 dB, check 3-6). The production LLE
  is never a detection rule in later phases.

### 2. The So et al. surrogate test holds its nominal size on linear nulls, but Level-A peaks alone do not separate noise from structure. ✔
- **Scope.** Production UPO detector (period 1, 50 surrogates); white noise and AR(1), 128-512 samples.
- **Supports.** Level-B false positives 14/200 and 15/200 at nominal 0.059 (2E-3, 2E-4); Level A in 70-98 %
  of white-noise windows (2E-5; 78-98 % for AR(1)).
- **Limits.** One cell (AR(1), 256: 11/100, P = 0.033 unadjusted) has a Wilson interval excluding the nominal
  size. The test is not specific on periodic signals: baseline sinusoid 89/150 and two-tone 67/150 (check 4-1),
  the problem Conclusion 4 addresses.

### 3. The τ = 1 Rosenstein slope against IAAFT surrogates (Phase 3 C1) is a size-controlled test of "not a monotone transform of a linear Gaussian process", with full power on noisy logistic and Hénon maps. ✔
- **Scope.** Phase 3 test seeds, 256 samples; nulls white noise and AR(1); chaotic logistic and Hénon maps
  (clean and 30/20/10 dB).
- **Supports.** White noise 4/100, AR(1) 3/100, 20 dB maps 200/200 (3-1, 3-2); also 99-100/100 at 10 dB.
- **Limits.** The null is linear-Gaussian, not "non-chaotic": the test fires on static nonlinearity (Phase 5 N5
  53/300), SETAR (26/200) and ectopy (Conclusion 5). It misses sampled flows (Phase 5 G2/G3 0/200, p = 1.00).

### 4. The instability-gated UPO test (Phase 4 C2) removes the periodic and quasi-periodic false positives of the source UPO detector. ✔
- **Scope.** Phase 4 test seeds, 256 samples, discrete maps, period 1.
- **Supports.** Sinusoid 0/150 and two-tone 1/150 (baseline 89 and 67/150); noisy Hénon 293/300 (check 4-1).
- **Limits.** At 512 samples C2 still flags two-tone 2/150. Power falls under drift (Phase 5 P3 trend: 41/200
  and 86/200) until the D2 detrend (Conclusion 7) and at 10 dB (logistic 48/150). In Phase 9 the UPO component
  of K1 fired on non-chaotic coupled vdP with ectopy (Conclusion 6).

### 5. The LLE test alone cannot be interpreted on RR intervals that contain ectopic beats. ✔
- **Scope.** LLE alone at 256 intervals on synthetic linear-Gaussian RR (Phases 5-6) and at 512 on real windows
  (Phase 10).
- **Supports.** Fired in 96-99 % of windows with 2-10 % isolated ectopic beats ([Erratum E5]; checks 5-3,
  6-1.S2*); 80-100 % for couplets, runs and atrial-type ectopy (6-1.*); in MIT-BIH, 50 % of abnormal vs 12.8 % of
  normal windows (7-3); real-data dose-response OR 1.51 [1.31, 1.75] per doubling of 1 + burden (10-P3-LLE).
- **Limits.** Strictly regular patterns trigger it less (bigeminy 37/300, trigeminy 103/300; 6-1). Editing
  reduces firing to 8-15 % when ≥ 80 % of intervals are normal but raises it on edited trigeminy (291/300; 6-2).

### 6. The AND of the LLE and UPO tests (K1) is protected from isolated ectopy only partially, and not from ectopy in non-chaotic deterministic rhythms. ✔
- **Scope.** K1 configuration at 256 (Phase 6, raw synthetic ectopy) and K1 at 512 (Phases 9-10).
- **Supports.** At most 3.7 % false positives on raw synthetic ectopy patterns (11/300 atrial 10 %; 6-3); up to
  45/100 on non-chaotic coupled vdP with ectopy (9-7; failed 4 Phase 9 conditions); real positives track burden,
  OR 1.44 [1.09, 1.90] (10-P3-K1).
- **Limits.** "Partially" is quantified only for the tested patterns; edited trigeminy gave 13/300 (4.3 %;
  6-2). The real-data dose-response rests on 16 positives, 9 of them in windows with ≥ 16 ectopic beats
  (10-P3-bt.16).

### 7. Linear detrending gated at trend/residual ratio 0.7 (D2) recovers power lost to drift without adding false positives. ✔
- **Scope.** Phase 6 test seeds, 256 intervals, linear trends of 5-10 % of the mean RR; the AND detector.
- **Supports.** Pooled power on trended maps 32.5 % → 95.2 % (195 → 571 of 600; 6-6, 6-8); AND false positives
  unchanged on N1-N6 and S2 (6-4).
- **Limits.** The LLE component alone fires more often on trend and step nulls after detrending (N3 23 → 40,
  N4 17 → 40 of 400; 6-9); the AND is unaffected. G1 (Lorenz maxima) power −3/300. Only linear drift was tested.

### 8. None of the project's detectors sees continuous-flow chaos sampled as RR intervals at realistic lengths. ✖ (see ERRATA.md E14; not corrected, author's decision)
- **Scope as evidenced.** Flows sampled at fixed time steps: Phase 5 G2 (Rössler, about 6 samples per orbit)
  and G3 (Mackey-Glass τ = 17, step 6.0), 256 intervals; Phase 8 the same flows at 512.
- **Supports.** K1-type AND 0/200 on G2 and G3 at m = 2, 3, 4 (5-6); K1 0/200 on Phase 8 flows (8-5); LLE alone
  p = 1.00 on every G2/G3 window (Phase 5 report 5.3).
- **Contradicting results (reported in Phases 8 and 9).** Mackey-Glass chaos represented as intervals between
  successive maxima (the Phase 8 model's beat definition) **was** detected: Phase 8 C1 (= K1) 26/70 at the input
  variant and 23/70, 18/70 with isolated ectopy or couplets (0/70 with bigeminy); Phase 9 TEST seeds, realistic
  level, DEV regimes: K3 81/140 and K1 42/140 with isolated ectopy or couplets, 0/70 with bigeminy; by regime, K3
  20/20 at τ = 18 and 20, none at τ = 16.5-17 (VERIFICATION.md X-21). The statement holds for fixed-step sampling
  and for weak / near-threshold Mackey-Glass chaos (τ ≤ 17), not for flows in general.
- **Also limits.** "K3 needs strong low-dimensional predictability that decays within 5 beats (Phase 9)" is a
  design property of the growth gate, not a measured result.

### 9. No preregistered real-data analysis found evidence of deterministic chaos beyond what the detectors' known false-positive sources explain. ⚠
- **Scope.** K1 on MIT-BIH (Phase 7, 305 windows of 256), K3 on nsrdb/chfdb (Phase 9 Part F, 330 windows of 512)
  and on nsr2db/chf2db (Phase 10, 988 windows of 512).
- **Supports.** MIT-BIH K1 0/305 (7-1); nsrdb / chfdb K3 1/180 and 0/150, no excess over the synthetic
  false-positive bound under the preregistered Part F rule (9-8, 9-10); nsr2db/chf2db K3 8/988 (10-Q1-U.K3).
- **Limits.** K3's 8/988 (0.8 %) is above its synthetic false-positive rate (0/10,100). The report attributes the
  detections to sources outside the synthetic nulls: one degenerate z statistic at 1/128 s and low variability
  (nsr029 window 9); three windows from one subject; sensitivity to window length (9/320 at 1,024 intervals vs
  1/320 at 512). Nonstationarity, sleep/wake transitions and nonlinear stochastic regulation were not among the
  tested nulls (Phase 9 report 8), so for 7 of the 8 detections the "known false-positive source" is a
  hypothesis, not a demonstrated cause. Phase 7's annotation-time arm gave K1 6/120 abnormal windows (within the
  ectopy-only expectation; not the primary arm).

### 10. The annotation-masked, growth-gated nonlinear-prediction test (K3) is the most specific detector built in the project, at a large cost in power. ✔
- **Scope.** Phase 9 TEST (512 intervals; 101 conditions; true beat types as labels) and Phase 10 spike-ins.
- **Supports.** 0/10,100 synthetic null and non-chaotic windows (9-1); 90/1,320 = 6.8 % pooled power on TEST
  chaos, 0 % with bigeminy (9-1, 9-5); additive-chaos exclusion in real windows only from chaos fractions of about
  0.5 (π < 0.05; 0.3 for π < 0.20) for strongly chaotic maps (10-Q1-K3.*).
- **Limits.** K3 depends on annotation-quality labels (Conclusion 11); bigeminy and other ectopy-dense windows
  are not analysable; it was calibrated at 512 intervals only and its real-data rate rose at 1,024 (0.3 % → 2.8 %;
  10-Q4b-K3); at 1/128 s and low variability its z statistic can degenerate (10-K3t-nsr029).

### 11. K3's real-data behaviour depends on annotation-quality beat labels. ✔
- **Scope.** Phase 10 confirmatory windows (manually reviewed nsr2db/chf2db annotations).
- **Supports.** K4 and K3RR (RR-rule labels) 0/988 vs K3 8/988 (10-Q4d).
- **Limits.** With RR-rule labels power is also far lower on synthetic chaos (Phase 9 k4 2/1,320), so the
  disappearance of detections does not by itself show they were label artefacts.

### 12. A faithful re-implementation of Poon-Barahona noise titration reproduces the published behaviour on standard test systems. ✔
- **Scope.** Logistic bifurcation scan, non-chaotic controls, Lorenz, Mackey-Glass and Hénon (Phase 10 Q2b).
- **Supports.** NL > 0 at 78/78 chaotic and NL = 0 at 18/18 periodic r (predeclared ±0.02 margin; 10-V4);
  Spearman(NL, LE) 0.85 (10-V4b); V1-V3 within criteria (10-V1-3).
- **Limits.** The primary method papers (Barahona & Poon 1996; Poon & Merrill 1997) were not accessible; three
  implementation readings were chosen from restatements, two of them partly by matching the PNAS 2001 Fig. 2
  values (V1-V3 are therefore not independent checks). V1 gives 63.7 % against ≈ 75 % read from the figure.

### 13. Noise titration on heart-rate data mostly detects ectopy, not chaos. ✔
- **Scope.** Verified titration; Phase 10 confirmatory windows (512) and 12-min segments; Q2e synthetic checks.
- **Supports.** Window positive rate 22 % with no ectopic beat vs 63 % with one (10-P3-22-63); OR 2.31
  [1.67, 3.19] per doubling (10-P3-TIT); masking removed 49-66 % and editing 44-55 % of segment positives
  ([Erratum E4]; 10-P2c-*); couplets on linear RR 100/100 and runs 96-99/100 (10-Q2e-*).
- **Limits.** The masked arm is not analysable in ectopy-dense segments (62 % of CHF segments analysable;
  10-P2-af), so the masked comparison favours CHF segments with little ectopy. Unannotated ectopic beats would
  remain in every masked arm. About a quarter to a third of segments stay positive after ectopy removal
  (Conclusion 15).

### 14. The "decreased cardiac chaos in CHF" claim is not supported by titration on the confirmatory databases. ✔
- **Scope.** nsr2db (54 healthy, 58-76 y) vs chf2db (29 CHF, NYHA I-III); 12-min segments; window-level GEE.
- **Supports.** Raw RR: CHF more positives, +22.4 % [12.8, 31.7] (10-P2-raw); after masking −4.9 % [−13.7,
  +4.3], editing +5.4 % [−3.9, +14.7], Wu-style +2.9 % [−5.3, +11.4] (10-P2-*); CHF OR adjusted for burden 0.77
  [0.46, 1.28] (10-P3-TIT); unchanged without the 10 subjects flagged by the original overlap rule (10-S-1).
- **Limits.** Different subjects and preprocessing from the original studies: chf2db is NYHA I-III (Poon &
  Merrill: severe CHF), the healthy group is older than Wu's young group, manual removal of premature beats
  could not be reproduced on RR-only data, and Poon & Merrill's preprocessing was not accessible. The adjusted
  (independence) and unadjusted (exchangeable) titration ORs use different GEE working correlations (ERRATA.md
  O3).

### 15. About a quarter of ectopy-free healthy and CHF segments remain titration-positive, but this cannot be read as chaos. ✔
- **Scope.** Segments after masking, editing or Wu-style elimination of ectopy-related intervals.
- **Supports.** DR after ectopy removal 21-34 % in both groups (masked 25.7 % / 20.8 %, edited 28.9 % / 34.3 %,
  Wu 28.2 % / 31.1 %; 10-P2-*); the verified indicator is positive on a static monotone transform of linear noise
  (N5 100/100) and on the non-chaotic forced vdP (5.45, 5.6) rhythm (100/100) (10-Q2e-*).
- **Limits.** "Ectopy-free" means after the ectopy-handling arms, not segments without any ectopic beat; the
  remaining positives were not tested against window-level surrogates (listed as future work).

---

## Other scope notes for writing

- **Exclusion bounds (P1)** refer only to an independent additive chaotic component at the tested families and
  exponents, with real ectopy kept and 1/128 s quantization; they are 95 % upper bounds, not prevalence estimates
  (F4, T4). No bound exists for logistic λ ≤ 0.26, Hénon λ = 0.14 and any phase-resetting regime (K3).
- **Phase 7 MIT-BIH** used the K1 configuration at 256 intervals with detected R peaks; its spike-in analysis is
  exploratory (F6).
- **Night/day** results are approximate: no confirmatory header has a start time (09:54 assumed).
