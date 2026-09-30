# Phase 6 Part A: BASELINE-K (keep_upo_on_short_lle_embedding = True) vs Phase 5

## The three Phase 5 windows that raised

| window | Phase 5 | BASELINE-K error | LLE embedding too short | LLE | UPO | AND |
|---|---|---|---|---|---|---|
| P1_henon_rr 3008 | raised (ValueError); AND counted 0 | none | True | True | True | True |
| P1_henon_rr 3154 | raised (ValueError); AND counted 0 | none | True | True | True | True |
| P4_henon_rr_30dB 3014 | raised (ValueError); AND counted 0 | none | True | True | True | True |

## N1-N6, Phase 5 test seeds 3000-3299

| condition | AND Phase 5 | AND BASELINE-K | LLE P5 / K | UPO P5 / K | windows with any difference | errors (K) |
|---|---|---|---|---|---|---|
| N1_linear_rr | 1/300 | 1/300 | 20 / 20 | 4 / 4 | 0 | 0 |
| N2_power_law | 1/300 | 1/300 | 8 / 8 | 4 / 4 | 0 | 0 |
| N3_linear_rr_trend | 0/300 | 0/300 | 24 / 24 | 2 / 2 | 0 | 0 |
| N4_linear_rr_step | 0/300 | 0/300 | 14 / 14 | 3 / 3 | 0 | 0 |
| N5_linear_rr_warped | 1/300 | 1/300 | 53 / 53 | 6 / 6 | 0 | 0 |
| N6_noisy_rsa | 0/300 | 0/300 | 17 / 17 | 1 / 1 | 0 | 0 |

Windows with any decision difference over N1-N6: **0** / 1800.
