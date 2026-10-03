# Phase 10 tables (conf, analysis)



Windows 988, subjects 83.



## Q1 detection limits



## Q2 titration (12-min segments)



## Q2e synthetic titration checks



## Q3 ectopy dose-response

| outcome | windows | positive | OR per doubling of 1 + burden [95 % CI] | p | CHF OR (adjusted for burden) |
|---|---|---|---|---|---|
| **titration (P3)** | 988 | 455 | 2.31 [1.67, 3.19] | 4.6e-07 | 0.77 [0.46, 1.28], p = 0.31 |
| **K3 (P3)** | 988 | 8 | not estimable (fewer than 10 positives or negatives (or no burden variation)) | | |
| LLE alone | 988 | 401 | 1.51 [1.31, 1.75] | 1.6e-08 | 0.79 [0.49, 1.28], p = 0.34 |
| UPO alone | 988 | 31 | 1.17 [1.00, 1.37] | 0.053 | 2.13 [1.00, 4.55], p = 0.05 |
| K1 | 988 | 16 | 1.44 [1.09, 1.90] | 0.0096 | 2.33 [0.57, 9.41], p = 0.24 |

| burden (beats) | windows | subjects | titration | LLE | UPO | K1 | K3 |
|---|---|---|---|---|---|---|---|
| 0 | 562 | 72 | 123 (22 %) | 145 (26 %) | 15 (3 %) | 5 (1 %) | 7 (1 %) |
| 1 | 116 | 55 | 73 (63 %) | 63 (54 %) | 2 (2 %) | 0 (0 %) | 0 (0 %) |
| 2-4 | 116 | 39 | 93 (80 %) | 71 (61 %) | 3 (3 %) | 2 (2 %) | 0 (0 %) |
| 5-15 | 83 | 26 | 71 (86 %) | 42 (51 %) | 0 (0 %) | 0 (0 %) | 1 (1 %) |
| >=16 | 111 | 20 | 95 (86 %) | 80 (72 %) | 11 (10 %) | 9 (8 %) | 0 (0 %) |

| paired comparison (3b) | subset | windows | raw positive | other arm positive | both | raw positives removed | mean subject-level change [95 % CI] |
|---|---|---|---|---|---|---|---|
| TIT_masked | all | 720 | 231 | 159 | 154 | 33.3 % | -13.9 % [-18.7, -9.3] |
| TIT_masked | n_masked_ge1 | 281 | 134 | 62 | 57 | 57.5 % | -31.5 % [-40.8, -22.4] |
| TIT_edited | all | 943 | 415 | 211 | 205 | 50.6 % | -22.4 % [-28.0, -17.1] |
| TIT_edited | n_masked_ge1 | 504 | 318 | 114 | 108 | 66.0 % | -36.8 % [-44.5, -29.3] |
| K1_edited | all | 943 | 10 | 7 | 5 | 50.0 % | -0.4 % [-1.0, 0.2] |
| K1_edited | n_masked_ge1 | 504 | 7 | 4 | 2 | 71.4 % | -0.5 % [-1.4, 0.2] |
| LLE_edited | all | 943 | 365 | 270 | 241 | 34.0 % | -11.4 % [-15.8, -7.2] |
| LLE_edited | n_masked_ge1 | 504 | 263 | 168 | 139 | 47.1 % | -19.5 % [-25.8, -13.3] |
| UPO_edited | all | 943 | 25 | 29 | 16 | 36.0 % | 0.6 % [-0.6, 1.9] |
| UPO_edited | n_masked_ge1 | 504 | 14 | 18 | 5 | 64.3 % | 0.6 % [-0.9, 2.3] |

Q3c: titration CHF OR unadjusted 2.60 [1.64, 4.13], p = 4.9e-05.

## Q4 robustness


| detector | NSR night | NSR day | CHF night | CHF day |
|---|---|---|---|---|
| TIT | 66/152 | 177/490 | 47/82 | 165/264 |
| LLE | 40/152 | 190/490 | 33/82 | 138/264 |
| UPO | 1/152 | 11/490 | 3/82 | 16/264 |
| K1 | 1/152 | 3/490 | 2/82 | 10/264 |
| K3 | 1/152 | 4/490 | 0/82 | 3/264 |
| K4 | 0/152 | 0/490 | 0/82 | 0/264 |

Titration GEE night effect (adjusted for burden and group): OR 1.40 [0.97, 2.00], p = 0.069.

Q4d rates: {"K3": {"NSR": 0.00778816199376947, "CHF": 0.008670520231213872}, "K4": {"NSR": 0.0, "CHF": 0.0}, "K3RR": {"NSR": 0.0, "CHF": 0.0}}; K3 vs K4 agreement {"both": 0, "a_only": 8, "b_only": 0, "neither": 980}; K3 vs K3RR {"both": 0, "a_only": 8, "b_only": 0, "neither": 980}.
