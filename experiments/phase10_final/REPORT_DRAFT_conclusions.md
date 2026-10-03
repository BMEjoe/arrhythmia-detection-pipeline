## 12. Conclusions supported by Phases 2E–10

One sentence per conclusion, each with its evidence. "Detected" always means the stated statistical
test rejected its null hypothesis; that is not proof of chaos.

1. **A positive largest Lyapunov exponent from the production Rosenstein estimator is not evidence of
   chaos.** It was positive in 194/195 white-noise and 197/197 AR(1) windows and underestimated
   known exponents (Phase 2E, Table H).
2. **The So et al. surrogate test for unstable periodic orbits holds its nominal size on linear
   nulls, but Level-A peaks alone do not separate noise from structure.** Level-B false positives
   were 14/200 and 15/200 at nominal 0.059; Level A fired in 70–98 % of null windows (Phase 2E, 3.1).
3. **The τ = 1 Rosenstein slope tested against IAAFT surrogates (Phase 3 C1) is a size-controlled
   test of "not a monotone transform of a linear Gaussian process", with full power on noisy
   logistic and Hénon maps.** White noise 4/100, AR(1) 3/100, 20 dB maps 200/200 (Phase 3 test
   seeds).
4. **The instability-gated UPO test (Phase 4 C2) removes the periodic and quasi-periodic false
   positives of the source UPO detector.** Sinusoid 0/150 and two-tone 1/150 (baseline 89 and
   67/150); noisy Hénon 293/300 (Phase 4 test seeds).
5. **The LLE test alone cannot be interpreted on RR intervals that contain ectopic beats.** It fired in
   94–99 % of linear-Gaussian windows with 2–10 % isolated ectopic beats (Phases 5–6). On real data
   its positives track ectopy burden: OR 1.51 [1.31, 1.75] per doubling (Phase 10, Q3).
6. **The AND of the LLE and UPO tests (K1) is protected from isolated ectopy only partially, and not from
   ectopy in non-chaotic deterministic rhythms.**
   - At most 3.7 % false positives on synthetic ectopy patterns (Phase 6).
   - Up to 45/100 on non-chaotic coupled vdP with ectopy (Phase 9 TEST).
   - Its 16 real-data positives also track burden: OR 1.44 [1.09, 1.90] (Phase 10).
7. **Linear detrending gated at trend/residual ratio 0.7 (D2) recovers power lost to drift without
   adding false positives.** Pooled power on trended maps rose from 32.5 % to 95.2 % (Phase 6).
8. **None of the project's detectors sees continuous-flow chaos sampled as RR intervals at realistic
   lengths.** Rössler and Mackey–Glass were 0 % for K1 at every m (Phases 5, 8). K3 needs strong
   low-dimensional predictability that decays within 5 beats (Phase 9).
9. **No preregistered real-data analysis found evidence of deterministic chaos beyond what the
   detectors' known false-positive sources explain.**
   - MIT-BIH, K1: 0/305 windows (Phase 7).
   - nsrdb / chfdb, K3: 1/180 and 0/150 (Phase 9).
   - nsr2db / chf2db, K3: 8/988, one of them degenerate (Phase 10).
10. **The annotation-masked, growth-gated nonlinear-prediction test (K3) is the most specific detector
    built in the project, at a large cost in power.**
    - 0/10,100 synthetic null and non-chaotic windows.
    - 7 % pooled power on TEST chaos, 0 % with bigeminy (Phase 9).
    - Additive-chaos detection in real windows only at chaos fractions ≥ ~0.5 for strongly chaotic
      maps (Phase 10, Q1).
11. **K3's real-data behaviour depends on annotation-quality beat labels.** With RR-rule masks (K4,
    K3RR) it detected 0/988 confirmatory windows, against 8/988 with annotations (Phase 10, Q4d).
12. **A faithful re-implementation of Poon–Barahona noise titration reproduces the published
    behaviour on standard test systems.** The logistic scan gave NL > 0 at 78/78 chaotic and NL = 0
    at 18/18 periodic parameters; Lorenz, Mackey–Glass and Hénon behaved as published (Phase 10,
    Q2b). The approximate Phase 8 version did not: it was positive on SETAR and step nulls.
13. **Noise titration on heart-rate data mostly detects ectopy, not chaos.**
    - Window positive rate rose from 22 % with no ectopic beat to 63 % with one (Phase 10, Q3).
    - Masking or editing ectopy removed 49–66 % of segment positives (Phase 10, P2).
    - Isolated couplets and runs on linear RR were positive in every development window (Phase 10,
      Q2e).
14. **The "decreased cardiac chaos in CHF" claim is not supported by titration on the confirmatory
    databases.**
    - On raw RR, CHF had more positives than healthy subjects: +22 % [13, 32].
    - After masking, editing or Wu-style elimination of ectopy, there was no group difference:
      −4.9 % [−13.7, +4.3] masked.
    - The window-level group effect vanished after adjusting for burden: OR 0.77 [0.46, 1.28]
      (Phase 10, P2 and Q3c).
15. **About a quarter of ectopy-free healthy and CHF segments remain titration-positive, but this
    cannot be read as chaos.** The verified indicator is also positive on a static monotone
    transform of linear noise (N5) and on a non-chaotic forced vdP rhythm (Phase 10, Q2e).

## 11. Limitations

- **Titration sources.** The primary method papers (Barahona & Poon 1996; Poon & Merrill 1997) were
  not accessible.
  - Three implementation readings (ε as RMS, the variance-ratio F-test, linear memory up to 83 lags)
    were chosen from restatements. The last two were chosen partly by matching the three PNAS 2001
    Fig. 2 values, so those values are not an independent check.
  - The independent checks (bifurcation scan, controls, flows) all passed.
  - The κ_lin ≤ 83 reading makes the K3-style masked arm unanalysable in ectopy-dense segments.
  - A different reading (for example linear κ ≤ 6) would give many more positives, including on
    periodic signals.
- **Study replication is approximate.**
  - Wu et al. selected subjects with stable heart rate and few ectopic beats, and removed residual
    premature beats by hand from the ECG. Neither step can be reproduced on RR-only databases.
  - Their κ and d were not reported; the routine defaults (6, 3) were used.
  - Poon & Merrill's preprocessing was not accessible.
- **The confirmatory databases differ from the original studies' databases.** chf2db is NYHA I–III
  (Wu: III–IV); nsr2db subjects are older (58–76 y) than Wu's young group.
- **The additive spike-in is a model of hidden chaos, not of a chaotic heart.**
  - The exclusion bounds apply only to an independent low-dimensional chaotic component added to real
    variability, at the tested families and exponents.
  - Phase-reset and weak maps were rarely detected, so no bound is possible for them.
- **Real-rate bound.** U is driven by the Clopper–Pearson term when real detections are rare. A
  subject-level (rather than window-level) prevalence would need the more conservative bound reported
  as secondary.
- **No start times.** All confirmatory night/day labels assume a 09:54 recording start.
- **K3 at 1/128 s.** K3 was validated at 1/360 s. At 1/128 s and low variability its z statistic can
  degenerate (one of the 8 detections).
- **Annotation quality.** chfdb annotations (development) were not manually corrected; nsr2db / chf2db
  were reviewed. Unannotated ectopic beats would leak into every "masked" or "NN" arm.
- **Overlap rule.** The preregistered rule produced false matches and was amended after the download
  (before any detector output). The sensitivity analysis without the 10 flagged subjects is reported.
- **Model fallback.** The exchangeable GEE diverged for titration (Amendment 2). The independence
  working correlation was used; it gives valid robust SEs but is less efficient.
- **One machine, one run.** Detector and titration Monte Carlo variability (surrogate draws, noise
  realisations) is not repeated except for the 499-surrogate subset (Q4a).

## 13. Future work (not pursued, by design)

- Obtain Barahona & Poon (1996) and run the original Poon-lab titration routine on the same segments to
  settle the three implementation readings.
- Repeat the dose-response with ECG-based ectopy review (databases with waveforms, e.g. the MIT-BIH Long
  Term or the BIDMC CHF waveforms) to catch unannotated ectopic beats.
- A masked titration variant whose linear memory adapts to the clean-span length, so that ectopy-dense
  segments are analysable.
- Detection limits for flow-like and high-dimensional chaos, which no current detector sees.
- Spike-in designs that model chaos *in* the sinus node (e.g. a chaotic modulation of a physiological
  integrate-and-fire model) rather than an additive component.
- Window-level surrogate-null calibration of titration (titration against IAAFT surrogates of each
  window) as an explicit specificity control.
