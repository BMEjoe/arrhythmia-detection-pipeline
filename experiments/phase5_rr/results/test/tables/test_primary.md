# Phase 5 primary decision (test seeds, file test)

PASS iff the AND detector (m = 2, analyze_segment) fires on <= floor(0.07 N) windows on EACH of N1-N6 (quantized, 256 RR intervals). Errors count as not detected.

| condition | AND detections | limit | within limit | analyze_segment errors | AND if errors were not counted |
|---|---|---|---|---|---|
| N1_linear_rr | 1/300 | 21 | yes | 0 | 1 |
| N2_power_law | 1/300 | 21 | yes | 0 | 1 |
| N3_linear_rr_trend | 0/300 | 21 | yes | 0 | 0 |
| N4_linear_rr_step | 0/300 | 21 | yes | 0 | 0 |
| N5_linear_rr_warped | 1/300 | 21 | yes | 0 | 1 |
| N6_noisy_rsa | 0/300 | 21 | yes | 0 | 0 |

**PASS**

Power (AND, reported, not pass/fail):

| condition | AND | LLE alone | UPO alone | OR |
|---|---|---|---|---|
| P1_henon_rr | 198/200 | 198 | 198 | 198 |
| P2_logistic_rr | 194/200 | 200 | 194 | 200 |
| P3_henon_rr_trend | 41/200 | 200 | 41 | 200 |
| P3_logistic_rr_trend | 86/200 | 200 | 86 | 200 |
| P4_henon_rr_30dB | 199/200 | 199 | 199 | 199 |
| P4_henon_rr_20dB | 198/200 | 200 | 198 | 200 |
| P4_logistic_rr_30dB | 200/200 | 200 | 200 | 200 |
| P4_logistic_rr_20dB | 188/200 | 200 | 188 | 200 |
| G1_lorenz_maxima | 199/200 | 200 | 199 | 200 |
| G2_rossler_flow | 0/200 | 0 | 0 | 0 |
| G3_mackey_glass | 0/200 | 0 | 62 | 62 |
