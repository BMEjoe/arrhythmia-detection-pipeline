# Phase 4 primary decision (test seeds, 256 samples, file test)

PASS: detections <= floor(0.07 N) on EACH of white_noise, ar1, sinusoid, two_tone. WINNER: most pooled detections on henon@30dB + henon@20dB; tie-break: median localization error on clean henon, then fewer pooled false positives, then candidate index. The baseline is not eligible.

| method | WN | AR(1) | sinusoid | two-tone | PASS | Hénon 30+20 dB | clean Hénon loc. error |
|---|---|---|---|---|---|---|---|
| baseline | 12/150 (≤10) | 7/150 (≤10) | 89/150 (≤10) | 67/150 (≤10) | no (not eligible) | 110/300 | 0.0679 |
| c1_prod_extgate | 0/150 (≤10) | 3/150 (≤10) | 0/150 (≤10) | 0/150 (≤10) | yes | 99/300 | 0.0679 |
| c2_m2M15_mediangate | 0/150 (≤10) | 0/150 (≤10) | 0/150 (≤10) | 1/150 (≤10) | yes | 293/300 | 0.0102 |
| c3_m2M7_trimgate | 4/150 (≤10) | 3/150 (≤10) | 0/150 (≤10) | 1/150 (≤10) | yes | 289/300 | 0.0100 |

Ranking of passing candidates: c2_m2M15_mediangate, c3_m2M7_trimgate, c1_prod_extgate

**WINNER: c2_m2M15_mediangate**
