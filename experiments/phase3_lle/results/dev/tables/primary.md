# Primary decision rule (dev seeds, 256 samples)

PASS: lower 95 % Wilson bound of the false-positive rate <= 0.05 on white_noise AND ar1. WINNER: highest pooled detection on logistic@20dB + henon@20dB; tie-break: lowest |bias logistic| + |bias henon| (clean), then fewer null false positives, then method index. The baseline is reported but not eligible.

| method | WN FP (Wilson lower) | AR1 FP (Wilson lower) | PASS | 20 dB pooled | clean |bias| sum |
|---|---|---|---|---|---|
| baseline | 7/100 (0.034) | 4/100 (0.016) | yes (not eligible) | 33/60 | 0.8718 |
| c1_rosenstein_m2_iaaft | 1/100 (0.002) | 5/100 (0.022) | yes | 60/60 | 0.0047 |
| c2_kantz_m3_sat_iaaft | 6/100 (0.028) | 0/100 (0.000) | yes | 60/60 | 0.0331 |
| c3_eps_m2_iaaft | 6/100 (0.028) | 5/100 (0.022) | yes | 60/60 | 0.0110 |
| c4_zero_one_iaaft | 3/100 (0.010) | 1/100 (0.002) | yes | 1/60 | NA (no LLE) |

Ranking of passing candidates: c1_rosenstein_m2_iaaft, c3_eps_m2_iaaft, c2_kantz_m3_sat_iaaft, c4_zero_one_iaaft

**WINNER: c1_rosenstein_m2_iaaft**
