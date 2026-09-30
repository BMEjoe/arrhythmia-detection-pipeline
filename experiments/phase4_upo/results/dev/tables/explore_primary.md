# Phase 4 primary decision (dev seeds, 256 samples, file explore)

PASS: detections <= floor(0.07 N) on EACH of white_noise, ar1, sinusoid, two_tone. WINNER: most pooled detections on henon@30dB + henon@20dB; tie-break: median localization error on clean henon, then fewer pooled false positives, then candidate index. The baseline is not eligible.

| method | WN | AR(1) | sinusoid | two-tone | PASS | Hénon 30+20 dB | clean Hénon loc. error |
|---|---|---|---|---|---|---|---|
| baseline | 7/100 (≤7) | 11/100 (≤7) | 21/30 (≤2) | 15/30 (≤2) | no (not eligible) | 19/60 | 0.0679 |
| c1_prod_extgate | 0/100 (≤7) | 2/100 (≤7) | 0/30 (≤2) | 0/30 (≤2) | yes | 18/60 | 0.0679 |
| c2_m2M15_mediangate | 0/100 (≤7) | 1/100 (≤7) | 0/30 (≤2) | 0/30 (≤2) | yes | 60/60 | 0.0094 |
| c3_m2M7_trimgate | 0/100 (≤7) | 3/100 (≤7) | 0/30 (≤2) | 0/30 (≤2) | yes | 57/60 | 0.0097 |

Ranking of passing candidates: c2_m2M15_mediangate, c3_m2M7_trimgate, c1_prod_extgate

**WINNER: c2_m2M15_mediangate**
