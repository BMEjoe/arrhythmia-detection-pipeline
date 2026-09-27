### Step 23 -- source_peak_coverage vs source_peak_count (core, non-failed windows)

| system           |   N |   spearman_rho |           p |   coverage_SD_within_count (median over counts) | coverage==0 when count==0   |
|:-----------------|----:|---------------:|------------:|------------------------------------------------:|:----------------------------|
| ALL (non-failed) | 583 |     0.319201   | 2.82848e-15 |                                      0.0978399  | 1/1                         |
| ar1              | 176 |     0.370496   | 4.16071e-07 |                                      0.00643977 | 0/0                         |
| henon            |  66 |     0.123077   | 0.324856    |                                      0.174033   | 0/0                         |
| logistic         |  79 |    -0.00956271 | 0.933341    |                                      0.192299   | 1/1                         |
| sinusoid         |  89 |     0.0573917  | 0.59319     |                                      0.00451447 | 0/0                         |
| white_noise      | 173 |    -0.0443463  | 0.562359    |                                      0          | 0/0                         |
