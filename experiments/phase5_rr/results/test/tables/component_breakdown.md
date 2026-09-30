# Phase 5 test records: component breakdown at m = 2 (descriptive)

| condition | N | >= 1 Level-B peak | UPO gate passed | LLE p quartiles | LLE p <= 0.05 | AND | production TDMI tau / Cao m (median) |
|---|---|---|---|---|---|---|---|
| N1_linear_rr | 300 | 29 | 4 | 0.18 / 0.43 / 0.70 | 20 | 1 | 3 / 8 |
| N2_power_law | 300 | 24 | 4 | 0.33 / 0.55 / 0.79 | 8 | 1 | 6 / 9 |
| N3_linear_rr_trend | 300 | 29 | 2 | 0.19 / 0.38 / 0.67 | 24 | 0 | 4 / 8 |
| N4_linear_rr_step | 300 | 26 | 3 | 0.26 / 0.53 / 0.78 | 14 | 0 | 4 / 8 |
| N5_linear_rr_warped | 300 | 31 | 6 | 0.09 / 0.27 / 0.52 | 53 | 1 | 3 / 9 |
| N6_noisy_rsa | 300 | 16 | 1 | 0.27 / 0.54 / 0.80 | 17 | 0 | 5 / 9 |
| S1_setar | 200 | 15 | 2 | 0.12 / 0.34 / 0.59 | 26 | 1 | 4 / 9 |
| S2_ectopic_2pct | 200 | 19 | 1 | 0.01 / 0.01 / 0.01 | 194 | 1 | 3 / 10 |
| S2_ectopic_5pct | 200 | 21 | 2 | 0.01 / 0.01 / 0.01 | 195 | 2 | 3 / 11 |
| S2_ectopic_10pct | 200 | 18 | 7 | 0.01 / 0.01 / 0.01 | 197 | 7 | 3 / 10 |
| S3_linear_rr_unq | 200 | 19 | 2 | 0.23 / 0.41 / 0.68 | 15 | 0 | 3 / 8 |
| S3_linear_rr_trend_unq | 200 | 15 | 2 | 0.17 / 0.44 / 0.71 | 25 | 0 | 4 / 9 |
| S3_linear_rr_step_unq | 200 | 19 | 3 | 0.29 / 0.56 / 0.77 | 8 | 0 | 4 / 8 |
| P1_henon_rr | 200 | 200 | 200 | 0.01 / 0.01 / 0.01 | 200 | 200 | 11 / 10 |
| P2_logistic_rr | 200 | 194 | 194 | 0.01 / 0.01 / 0.01 | 200 | 194 | 8 / 10 |
| P3_henon_rr_trend | 200 | 44 | 41 | 0.01 / 0.01 / 0.01 | 200 | 41 | 7 / 9 |
| P3_logistic_rr_trend | 200 | 86 | 86 | 0.01 / 0.01 / 0.01 | 200 | 86 | 5 / 9 |
| P4_henon_rr_30dB | 200 | 200 | 200 | 0.01 / 0.01 / 0.01 | 200 | 200 | 11 / 10 |
| P4_henon_rr_20dB | 200 | 199 | 198 | 0.01 / 0.01 / 0.01 | 200 | 198 | 10 / 10 |
| P4_logistic_rr_30dB | 200 | 200 | 200 | 0.01 / 0.01 / 0.01 | 200 | 200 | 7 / 10 |
| P4_logistic_rr_20dB | 200 | 190 | 188 | 0.01 / 0.01 / 0.01 | 200 | 188 | 6 / 10 |
| G1_lorenz_maxima | 200 | 199 | 199 | 0.01 / 0.01 / 0.01 | 200 | 199 | 8 / 9 |
| G2_rossler_flow | 200 | 200 | 0 | 1.00 / 1.00 / 1.00 | 0 | 0 | 5 / 7 |
| G3_mackey_glass | 200 | 106 | 62 | 1.00 / 1.00 / 1.00 | 0 | 0 | 6 / 7 |

Counts include the 3 windows where analyze_segment raised (component decisions from the direct calls); the preregistered AND counts those windows as not detected.
