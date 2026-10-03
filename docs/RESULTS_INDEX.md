# Results index

Key numbers likely to be used in a paper, each with its value, uncertainty, sample size, the exact file it comes
from and the script that produces it. Every count, rate and point estimate here was recomputed independently from
the raw result files (`experiments/final_refinement/VERIFICATION.md`, check id in the last column) and matches the
reports after the final-refinement errata. Intervals marked † are quoted from the phase's own (seeded) analysis
output and were not recomputed independently; all other intervals were (Wilson and Wald exactly, bootstrap
intervals within Monte Carlo tolerance). Paths are relative to the repository root; `P10 = experiments/phase10_final/results`.

Intervals: **W** = 95 % Wilson score interval; **CB** = 95 % subject-cluster bootstrap percentile interval
(10,000 resamples, seeded); **Wald** = 95 % Wald interval from the GEE robust (sandwich) standard error; **CP** =
one-sided 95 % Clopper-Pearson upper bound.

## Phase 2E: frozen pipeline on synthetic systems

| quantity | value | uncertainty | n | source file | script | check |
|---|---|---|---|---|---|---|
| production LLE > 0, white noise | 194/195 valid windows | | 195 | `experiments/phase2e/results/core.jsonl` | `experiments/phase2e/analysis.py` (table H) | 2E-1 |
| production LLE > 0, AR(1) | 197/197 | | 197 | same | same | 2E-2 |
| Level-B false positives, white noise (nominal 0.0588) | 14/200 = 7.0 % | W 4.2-11.4 % | 200 | same (`table_B_white_noise_null`) | same | 2E-3 |
| Level-B false positives, AR(1) | 15/200 = 7.5 % | W 4.6-12.0 % | 200 | same (`table_C_ar1_null`) | same | 2E-4 |

## Phase 3: LLE surrogate test (test seeds 1000-1099, 256 samples)

| quantity | value | uncertainty | n | source file | script | check |
|---|---|---|---|---|---|---|
| C1 (`lle_chaos_test`) false positives, white noise / AR(1) | 4/100 / 3/100 | W 1.6-9.8 % / 1.0-8.5 % | 100 each | `experiments/phase3_lle/results/test/c1_rosenstein_m2_iaaft.jsonl` | `phase3_lle/analysis.py`, `decide.py` | 3-1, 3-4 |
| C1 detection, logistic + Hénon at 20 dB | 200/200 | W ≥ 98.1 % | 200 | same | same | 3-2 |
| C1 clean LLE bias (median − reference), logistic / Hénon | −0.001 / −0.005 | IQR of estimate [0.684, 0.700] / [0.400, 0.432] | 100 each | same | same | 3-5 |
| C1 LLE bias at 20 dB, logistic / Hénon | −0.213 / −0.068 | | 100 each | same | same | 3-6 |
| production baseline: false positives WN / AR(1); 20 dB detection | 9/100, 5/100; 123/200 | | | `results/test/baseline.jsonl` | same | 3-3 |

## Phase 4: instability-gated UPO test (test seeds 2000-2149, 256 samples)

| quantity | value | uncertainty | n | source file | script | check |
|---|---|---|---|---|---|---|
| C2 false positives: white noise, AR(1), sinusoid, two-tone | 0, 0, 0, 1 /150 | W upper 2.5 %; two-tone 0.1-3.7 % | 150 each | `experiments/phase4_upo/results/test/test.jsonl` | `phase4_upo/analysis.py`, `decide.py` | 4-1 |
| baseline UPO false positives: sinusoid, two-tone | 89/150, 67/150 | W 51.3-66.9 %, 36.9-52.7 % | 150 each | same | same | 4-1 |
| C2 detection, Hénon 30 / 20 dB (pooled) | 150, 143 /150 (293/300) | W 97.5-100 %, 90.7-97.7 % | 300 | same | same | 4-1 |
| C2 AND `lle_chaos_test`: false positives on the 4 controls | 0 on each | | 150 each | same | same | 4-3 |

## Phase 5: combined detector on RR-like series (test seeds 3000-3299, 256 intervals, m = 2)

| quantity | value | uncertainty | n | source file | script | check |
|---|---|---|---|---|---|---|
| AND on nulls N1-N6 (PASS limit 21/300) | 1, 1, 0, 0, 1, 0 /300 | W upper 1.9 % (1/300), 1.3 % (0/300) | 300 each | `experiments/phase5_rr/results/test/test.jsonl` | `phase5_rr/analysis.py`, `decide.py` | 5-1 |
| LLE alone on isolated ectopy 2 / 5 / 10 % | 194, 195, 197 /200 | | 200 each | same | same | 5-3 |
| AND on isolated ectopy 10 % | 7/200 = 3.5 % | W 1.7-7.0 % | 200 | same | same | 5-3 |
| LLE alone on static warp N5 | 53/300 = 17.7 % | | 300 | same | same | 5-1 |
| AND on maps with trend (P3 Hénon / logistic) | 41/200, 86/200 | | 200 each | same | same | 5-4 |
| AND on Rössler / Mackey-Glass sampled flows, m = 2-4 | 0/200 each | | 200 each | same | same | 5-6 |

## Phase 6: robustness (test seeds 4000-4399, 256 intervals)

| quantity | value | uncertainty | n | source file | script | check |
|---|---|---|---|---|---|---|
| LLE alone on raw ectopy patterns (isolated, couplets, runs, atrial) | 80-100 % (254-299/300) | | 300 each | `experiments/phase6_robust/results/test/test.jsonl` | `phase6_robust/analysis.py` | 6-1 |
| AND (K1 configuration) worst raw ectopy pattern (atrial 10 %) | 11/300 = 3.7 % | W 2-6 % | 300 | same | same | 6-1, 6-3 |
| pooled AND power on trended maps, BASELINE-K → D2 | 195/600 (32.5 %) → 571/600 (95.2 %) | | 600 | same | `decide.py` | 6-6, 6-8 |
| pooled AND on 13 Part B ectopy conditions, D2 | 45/3,900 | | 3,900 | same | same | 6-6 |
| edited trigeminy: AND, LLE alone | 13/300 (4.3 %, W 2.5-7.2 %), 291/300 | | 300 | same | `analysis.py` | 6-2 |

## Phase 7: MIT-BIH (preregistered, 256 intervals, detected R peaks)

| quantity | value | uncertainty | n | source file | script | check |
|---|---|---|---|---|---|---|
| AND, abnormal / normal windows | 0/94 / 0/211 | W upper 3.9 % / 1.8 % | 305 windows, 43 subjects | `experiments/phase7_mitbih/results/run/raw.jsonl` | `phase7_mitbih/analysis.py` | 7-1, 7-2 |
| LLE alone, abnormal vs normal | 47/94 (50.0 %) vs 27/211 (12.8 %); difference +37.2 % | W 40.1-59.9 % / 8.9-18.0 %; CB [19.4, 54.8]† | 305 | same | same | 7-3 |
| out-of-fold AUC M0 (HRV) / M1 (lle_z) / M2 (upo_score) | 0.812 / 0.685 / 0.193 | CB [0.618, 0.949] / [0.565, 0.797] / [0.115, 0.290]† | 305 | same | same | 7-6 |
| annotation-time sensitivity arm: AND abnormal / normal | 6/120 (5.0 %) / 0/248 | CB [1.9, 8.3]† | 368 | `results/run/annotation.jsonl` | same | 7-5 |
| spike-in (EXPLORATORY), Hénon clean: normal / abnormal | 210/211 / 94/94 | | | `results/spike_in/spike_in.jsonl` | `spike_in.py --summary` | 7-9 |
| spike-in with real ectopy, abnormal: Hénon / logistic | 52/94 (55.3 %) / 74/94 (78.7 %) | W 45.3-65.0 % / 69.4-85.8 % | 94 | same | same | 7-9 |
| code audit | 10/10 windows reproduced bit for bit, 0 discrepancies | | 10 | `results/spike_in/audit.json` | `audit.py` | 7-12 |

## Phase 8: cardiac-model RR (TEST seeds 5000-5099)

| quantity | value | uncertainty | n | source file | script | check |
|---|---|---|---|---|---|---|
| C1 (= K1, 512): PASS; TEST chaotic detected at (vi) / (iv) | yes; 66/990 (6.7 %) / 10/330 (3.0 %) | | 35 conditions × 100 | `experiments/phase8_cardiac/results/test/test_512.jsonl` | `phase8_cardiac/decide.py` | 8-1 |
| approximate titration (C4, 256): failed conditions | 30 of 35 (e.g. SETAR 43/100, step 98/100, every ectopy pattern 58-100/100) | | | `test_256.jsonl` | same | 8-1, 8-4 |
| C1 on Rössler / Mackey-Glass flows | 0/200 | | 200 | `test_512.jsonl` | same | 8-5 |

## Phase 9: noise-robust measures (TEST seeds ≥ 9500, 512 intervals)

| quantity | value | uncertainty | n | source file | script | check |
|---|---|---|---|---|---|---|
| K3 (k3_mnlp_growth_ann): detections in the 101 PASS conditions | 0/10,100 | CP upper 0.03 % | 10,100 | `experiments/phase9_waveform/results/test/test9.jsonl` | `phase9_waveform/decide9.py`; `final_refinement/errata_phase9_decision.py` | 9-1, X-20 |
| K3 pooled primary chaotic detection | 90/1,320 = 6.8 % (coupled vdP 73/360, phase-reset 17/630, KTz 0/330) | | 1,320 | same | same | 9-1, 9-4 |
| K1 (k1_frozen512): failed conditions; worst | 4; 45/100 (non-chaotic vdP (9.6, 2.1) with ectopy) | | | same | same | 9-1, 9-7 |
| k2 failed conditions (101 preregistered) | 16 (erratum E1; decide9.py tabulated 12) | | | same | `errata_phase9_decision.py` | 9-1 |
| Part F: nsrdb / chfdb K3 detected windows | 1/180 / 0/150 (147 analysable) | subject rate 0.006 (CB 0-0.017†) | 18 + 15 subjects | `results/partf/partf_windows.jsonl` | `partf.py --analyze` | 9-8, 9-10 |

## Phase 10: confirmatory (nsr2db 54 healthy, chf2db 29 CHF)

| quantity | value | uncertainty | n | source file | script | check |
|---|---|---|---|---|---|---|
| windows / segments | 988 windows of 512 intervals; 9,536 12-min segments (8,993 analysable) | | 83 subjects | `P10/conf/real.jsonl`, `seg.jsonl` | `phase10_final/run10.py` | 10-1, 10-P2-n |
| titration verification V4 | NL > 0 at 78/78 chaotic r, NL = 0 at 18/18 periodic r; Spearman 0.85 | | 101 r values | `P10/verification/titration_verification.json` | `verify_titration.py` | 10-V4 |
| **P1** real detections, K3 / K1 | 8/988 / 16/988 | U = 0.0146 / 0.0245 (max of CB 95th percentile and CP) | 988 | `P10/conf/real.jsonl` | `analysis10.py` (`analysis.json` Q1) | 10-Q1-U.*, 10-Q1-cp.* |
| **P1** K3 smallest f with π < 0.05 / π < 0.20 (Hénon a = 1.22, 1.40; logistic r = 3.88, 4.00; vdP (8, 3.3)) | 0.5/0.5; 0.5/0.3; 0.5/0.5; 0.7/0.5; 0.5/0.3 | π_upper is itself a 95 % upper bound | 83 base windows per point | `P10/conf/spike.jsonl`, `analysis.json` | `analysis10.py`; table `tables/T4b_exclusion_fmin` | 10-Q1-K3.* |
| **P2** titration DR, healthy vs CHF (raw) | 51.7 % vs 74.2 %; difference +22.4 % | CB [46.0, 57.7], [66.5, 81.4]; [12.8, 31.7]; permutation p = 0.0001 | 54 / 29 subjects | `P10/conf/seg.jsonl` | `analysis10.py` (Q2) | 10-P2-raw |
| **P2** change after masking, healthy / CHF | −26.1 % / −46.8 % (raw positives removed 49 % / 66 %) | CB [−32.4, −20.1] / [−59.1, −34.5] | 5,671 / 1,865 segments | same | same | 10-P2c-masked.* |
| **P2** CHF − healthy after masking / editing / Wu-style | −4.9 % / +5.4 % / +2.9 % | CB [−13.7, 4.3] / [−3.9, 14.7] / [−5.3, 11.4]; p 0.29 / 0.24 / 0.50 | | same | same | 10-P2-* |
| **P3** titration OR per doubling of (1 + burden) | 2.31 | Wald [1.67, 3.19], p = 5 × 10⁻⁷ | 988 windows, 455 positive | `P10/conf/real.jsonl` | `analysis10.py` (Q3; independence GEE, Amendment 2) | 10-P3-TIT |
| **P3** titration CHF OR adjusted for burden (unadjusted) | 0.77 (2.60) | Wald [0.46, 1.28], p = 0.31 ([1.64, 4.13]) | 988 | same | same | 10-P3-TIT, 10-P3-3c |
| **P3** K3 | 8 positives, 7 in ectopy-free windows; not estimable (< 10) | rate difference burden ≥ 1 − 0: −0.7 %, CB [−1.6, +0.03] | 988 | same | same | 10-P3-K3, 10-P3-K3d |
| titration positive rate, 0 vs 1 ectopic beat | 22 % vs 63 % | | 562 / 116 windows | same | same | 10-P3-22-63 |
| LLE alone / UPO alone / K1 OR per doubling | 1.51 / 1.17 / 1.44 | Wald [1.31, 1.75] / [1.00, 1.37] / [1.09, 1.90] | 988 (401 / 31 / 16 positive) | same | same | 10-P3-* |
| Q2e: verified titration on static warp N5, couplets, runs, forced vdP (5.45, 5.6) | 100/100, 100/100, 96/100, 100/100 (512) | | 100 each | `P10/conf/synth.jsonl` | `analysis10.py` (Q2e) | 10-Q2e-* |
| Q4a 499 surrogates: changed window decisions | LLE 6 (3 + 3), K1 1, UPO 2, K3 0; κ 0.91-1.00 | | 331 windows | `P10/conf/robust.jsonl` | `analysis10.py` (Q4) | 10-Q4a-* |
| Q4b K3 rate at 256 / 512 / 1,024 intervals | 0.3 % / 0.3 % / 2.8 % | | 320 windows | same | same | 10-Q4b-K3 |
| Q4d K4 / K3RR (RR-rule labels) | 0/988 / 0/988 | | 988 | `P10/conf/real.jsonl` | same | 10-Q4d |
| sensitivity without the 10 originally flagged subjects: P2 raw difference; P3 titration OR; CHF OR | +23.3 %; 2.21; 0.82 | CB [11.1, 34.6]†; Wald [1.55, 3.15]†; [0.41, 1.65]† | 73 subjects | `P10/conf/*.jsonl` | `analysis10.py --sensitivity-flagged` | 10-S-1 |

The full tables behind these numbers are in `experiments/final_refinement/tables/` (CSV and Markdown) and in each
phase's `results/.../tables*`.
