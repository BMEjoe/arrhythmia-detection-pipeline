# Phase 6 Part C primary decision (test seeds, file test)

PASS: AND <= floor(0.07 N) on EACH of N1-N6 and each S2 rate, AND power on each of P1, P2, P4 (4) and G1 at most floor(0.03 N) below BASELINE-K on the same windows. WINNER: most pooled AND on P3; tie-break fewer pooled AND on the Part B ectopy conditions, then candidate index.

## Specificity (AND detections; limit)

| method | N1_linear_rr | N2_power_law | N3_linear_rr_trend | N4_linear_rr_step | N5_linear_rr_warped | N6_noisy_rsa | S2_ectopic_2pct | S2_ectopic_5pct | S2_ectopic_10pct |
|---|---|---|---|---|---|---|---|---|---|
| baseline_k | 0/400 (≤28) | 0/400 (≤28) | 1/400 (≤28) | 1/400 (≤28) | 2/400 (≤28) | 0/400 (≤28) | 1/300 (≤21) | 5/300 (≤21) | 7/300 (≤21) |
| d1_linear_g05 | 0/400 (≤28) | 0/400 (≤28) | 0/400 (≤28) | 0/400 (≤28) | 2/400 (≤28) | 0/400 (≤28) | 1/300 (≤21) | 5/300 (≤21) | 7/300 (≤21) |
| d2_linear_g07 | 0/400 (≤28) | 0/400 (≤28) | 0/400 (≤28) | 0/400 (≤28) | 2/400 (≤28) | 0/400 (≤28) | 1/300 (≤21) | 5/300 (≤21) | 7/300 (≤21) |
| d3_smoothprior300_g07 | 0/400 (≤28) | 0/400 (≤28) | 0/400 (≤28) | 0/400 (≤28) | 2/400 (≤28) | 0/400 (≤28) | 1/300 (≤21) | 5/300 (≤21) | 7/300 (≤21) |

## Power on untrended positives (AND; BASELINE-K; max loss)

| method | P1_henon_rr | P2_logistic_rr | P4_henon_rr_30dB | P4_henon_rr_20dB | P4_logistic_rr_30dB | P4_logistic_rr_20dB | G1_lorenz_maxima |
|---|---|---|---|---|---|---|---|
| baseline_k | 300 (base 300, ≥291) | 294 (base 294, ≥285) | 300 (base 300, ≥291) | 296 (base 296, ≥287) | 299 (base 299, ≥290) | 284 (base 284, ≥275) | 294 (base 294, ≥285) |
| d1_linear_g05 | 300 (base 300, ≥291) | 293 (base 294, ≥285) | 300 (base 300, ≥291) | 296 (base 296, ≥287) | 298 (base 299, ≥290) | 282 (base 284, ≥275) | 279 (base 294, ≥285) **FAIL** |
| d2_linear_g07 | 300 (base 300, ≥291) | 294 (base 294, ≥285) | 300 (base 300, ≥291) | 296 (base 296, ≥287) | 299 (base 299, ≥290) | 284 (base 284, ≥275) | 291 (base 294, ≥285) |
| d3_smoothprior300_g07 | 300 (base 300, ≥291) | 294 (base 294, ≥285) | 300 (base 300, ≥291) | 296 (base 296, ≥287) | 299 (base 299, ≥290) | 284 (base 284, ≥275) | 292 (base 294, ≥285) |

## Summary

| method | PASS | pooled P3 AND | pooled Part B ectopy AND |
|---|---|---|---|
| baseline_k | yes (not eligible) | 195/600 | 45/3900 |
| d1_linear_g05 | no | 581/600 | 44/3900 |
| d2_linear_g07 | yes | 571/600 | 45/3900 |
| d3_smoothprior300_g07 | yes | 539/600 | 45/3900 |

Ranking of passing candidates: d2_linear_g07, d3_smoothprior300_g07

**WINNER: d2_linear_g07**
