# Phase 5 primary decision (dev seeds, file runtime)

PASS iff the AND detector (m = 2, analyze_segment) fires on <= floor(0.07 N) windows on EACH of N1-N6 (quantized, 256 RR intervals). Errors count as not detected.

| condition | AND detections | limit | within limit | analyze_segment errors | AND if errors were not counted |
|---|---|---|---|---|---|
| N1_linear_rr | 0/10 | 0 | yes | 0 | 0 |
| N2_power_law | 0/10 | 0 | yes | 0 | 0 |
| N3_linear_rr_trend | 0/10 | 0 | yes | 0 | 0 |
| N4_linear_rr_step | 0/10 | 0 | yes | 0 | 0 |
| N5_linear_rr_warped | 0/10 | 0 | yes | 0 | 0 |
| N6_noisy_rsa | 0/10 | 0 | yes | 0 | 0 |

**PASS**

Power (AND, reported, not pass/fail):

| condition | AND | LLE alone | UPO alone | OR |
|---|---|---|---|---|
| P1_henon_rr | 10/10 | 10 | 10 | 10 |
| P2_logistic_rr | 10/10 | 10 | 10 | 10 |
| P3_henon_rr_trend | 2/10 | 10 | 2 | 10 |
| P3_logistic_rr_trend | 5/10 | 10 | 5 | 10 |
| P4_henon_rr_30dB | 10/10 | 10 | 10 | 10 |
| P4_henon_rr_20dB | 10/10 | 10 | 10 | 10 |
| P4_logistic_rr_30dB | 10/10 | 10 | 10 | 10 |
| P4_logistic_rr_20dB | 10/10 | 10 | 10 | 10 |
| G1_lorenz_maxima | 10/10 | 10 | 10 | 10 |
| G2_rossler_flow | 0/10 | 0 | 0 | 0 |
| G3_mackey_glass | 0/10 | 0 | 2 | 2 |
