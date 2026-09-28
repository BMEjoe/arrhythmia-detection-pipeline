# Primary decision rule (test seeds, 256 samples)

PASS: lower 95 % Wilson bound of the false-positive rate <= 0.05 on white_noise AND ar1. WINNER: highest pooled detection on logistic@20dB + henon@20dB; tie-break: lowest |bias logistic| + |bias henon| (clean), then fewer null false positives, then method index. The baseline is reported but not eligible.

| method | WN FP (Wilson lower) | AR1 FP (Wilson lower) | PASS | 20 dB pooled | clean |bias| sum |
|---|---|---|---|---|---|
| baseline | 9/100 (0.048) | 5/100 (0.022) | yes (not eligible) | 123/200 | 0.8695 |
| c1_rosenstein_m2_iaaft | 4/100 (0.016) | 3/100 (0.010) | yes | 200/200 | 0.0057 |
| c2_kantz_m3_sat_iaaft | 9/100 (0.048) | 0/100 (0.000) | yes | 200/200 | 0.0236 |
| c3_eps_m2_iaaft | 5/100 (0.022) | 3/100 (0.010) | yes | 200/200 | 0.0173 |
| c4_zero_one_iaaft | 0/100 (0.000) | 2/100 (0.006) | yes | 11/200 | NA (no LLE) |

Ranking of passing candidates: c1_rosenstein_m2_iaaft, c3_eps_m2_iaaft, c2_kantz_m3_sat_iaaft, c4_zero_one_iaaft

**WINNER: c1_rosenstein_m2_iaaft**
