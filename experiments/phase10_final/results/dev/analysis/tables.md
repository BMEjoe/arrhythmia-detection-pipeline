# Phase 10 tables (dev, analysis)



Windows 485, subjects 77.



## Q1 detection limits

#### K3: real detections 0/396 windows (33 subjects); U = 0.0075 (cluster bootstrap 0.0000, Clopper–Pearson 0.0075)

| family | parameter | λ | detected / base windows at f = 0.1, 0.3, 0.5, 0.9 | replacement | π_upper at f = 0.9 | smallest f, π < 0.05 | smallest f, π < 0.20 | (secondary: Wilson-lower p) π < 0.20 | (descriptive: Firth fit) π < 0.20 |
|---|---|---|---|---|---|---|---|---|---|
| Hénon | a1.08 | 0.136 | 0/6 · 0/6 · 0/6 · 0/6 | 0/6 | 1.000 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | 0.48 |
| Hénon | a1.14 | 0.244 | 0/6 · 0/6 · 0/6 · 0/6 | 0/6 | 1.000 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | 0.48 |
| Hénon | a1.22 | 0.303 | 0/6 · 1/6 · 3/6 · 6/6 | 2/6 | 0.008 | 0.3 | 0.3 | 0.5 | 0.09 |
| Hénon | a1.40 | 0.419 | 0/6 · 0/6 · 3/6 · 6/6 | 6/6 | 0.008 | 0.5 | 0.5 | 0.5 | 0.21 |
| logistic | r3.58 | 0.105 | 0/6 · 0/6 · 0/6 · 0/6 | 0/6 | 1.000 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | 0.48 |
| logistic | r3.65 | 0.255 | 0/6 · 0/6 · 0/6 · 0/6 | 0/6 | 1.000 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | 0.48 |
| logistic | r3.88 | 0.464 | 0/6 · 0/6 · 3/6 · 6/6 | 5/6 | 0.008 | 0.5 | 0.5 | 0.5 | 0.21 |
| logistic | r4.00 | 0.693 | 0/6 · 0/6 · 0/6 · 6/6 | 6/6 | 0.008 | 0.9 | 0.9 | 0.9 | 0.45 |
| coupled vdP | rho8_om3.3 | 0.0751 | 0/6 · 1/6 · 3/6 · 6/6 | 2/6 | 0.008 | 0.3 | 0.3 | 0.5 | 0.09 |
| coupled vdP | rho6_om3.3 | 0.0788 | 0/6 · 0/6 · 0/6 · 1/6 | 4/6 | 0.045 | 0.9 | 0.9 | none ≤ 0.9 | 0.4 |
| coupled vdP | rho2_om2.7 | 0.085 | 0/6 · 0/6 · 0/6 · 6/6 | 6/6 | 0.008 | 0.9 | 0.9 | 0.9 | 0.45 |
| coupled vdP | rho6_om4.0 | 0.101 | 0/6 · 1/6 · 2/6 · 2/6 | 4/6 | 0.023 | 0.3 | 0.3 | 0.5 | 0.05 |
| phase-reset | tau1.14 | 0.0354 | 0/6 · 0/6 · 0/6 · 0/6 | 0/6 | 1.000 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | 0.48 |
| phase-reset | tau0.58 | 0.081 | 0/6 · 0/6 · 1/6 · 0/6 | 0/6 | 1.000 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | 0.05 |
| phase-reset | tau1.20 | 0.138 | 0/6 · 0/6 · 0/6 · 0/6 | 0/6 | 1.000 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | 0.48 |
| phase-reset | tau1.16 | 0.185 | 0/6 · 0/6 · 0/6 · 0/6 | 0/6 | 1.000 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | 0.48 |

#### K1: real detections 3/396 windows (33 subjects); U = 0.0195 (cluster bootstrap 0.0152, Clopper–Pearson 0.0195)

| family | parameter | λ | detected / base windows at f = 0.1, 0.3, 0.5, 0.9 | replacement | π_upper at f = 0.9 | smallest f, π < 0.05 | smallest f, π < 0.20 | (secondary: Wilson-lower p) π < 0.20 | (descriptive: Firth fit) π < 0.20 |
|---|---|---|---|---|---|---|---|---|---|
| Hénon | a1.08 | 0.136 | 0/6 · 0/6 · 0/6 · 0/6 | 0/6 | 1.000 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |
| Hénon | a1.14 | 0.244 | 0/6 · 0/6 · 0/6 · 2/6 | 2/6 | 0.058 | none ≤ 0.9 | 0.9 | none ≤ 0.9 | 0.6 |
| Hénon | a1.22 | 0.303 | 0/6 · 0/6 · 2/6 · 3/6 | 4/6 | 0.039 | 0.9 | 0.5 | 0.9 | 0.3 |
| Hénon | a1.40 | 0.419 | 0/6 · 0/6 · 0/6 · 2/6 | 5/6 | 0.058 | none ≤ 0.9 | 0.9 | none ≤ 0.9 | 0.6 |
| logistic | r3.58 | 0.105 | 0/6 · 0/6 · 0/6 · 0/6 | 0/6 | 1.000 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |
| logistic | r3.65 | 0.255 | 0/6 · 0/6 · 0/6 · 0/6 | 0/6 | 1.000 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |
| logistic | r3.88 | 0.464 | 0/6 · 0/6 · 0/6 · 4/6 | 5/6 | 0.029 | 0.9 | 0.9 | 0.9 | 0.55 |
| logistic | r4.00 | 0.693 | 0/6 · 0/6 · 0/6 · 4/6 | 5/6 | 0.029 | 0.9 | 0.9 | 0.9 | 0.55 |
| coupled vdP | rho8_om3.3 | 0.0751 | 0/6 · 0/6 · 0/6 · 0/6 | 0/6 | 1.000 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |
| coupled vdP | rho6_om3.3 | 0.0788 | 0/6 · 0/6 · 0/6 · 0/6 | 0/6 | 1.000 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |
| coupled vdP | rho2_om2.7 | 0.085 | 0/6 · 0/6 · 0/6 · 1/6 | 1/6 | 0.117 | none ≤ 0.9 | 0.9 | none ≤ 0.9 | 0.67 |
| coupled vdP | rho6_om4.0 | 0.101 | 0/6 · 0/6 · 0/6 · 0/6 | 0/6 | 1.000 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |
| phase-reset | tau1.14 | 0.0354 | 0/6 · 0/6 · 0/6 · 1/6 | 0/6 | 0.117 | none ≤ 0.9 | 0.9 | none ≤ 0.9 | 0.67 |
| phase-reset | tau0.58 | 0.081 | 0/6 · 0/6 · 0/6 · 0/6 | 0/6 | 1.000 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |
| phase-reset | tau1.20 | 0.138 | 0/6 · 0/6 · 1/6 · 1/6 | 0/6 | 0.117 | none ≤ 0.9 | 0.5 | none ≤ 0.9 | 0.47 |
| phase-reset | tau1.16 | 0.185 | 1/6 · 0/6 · 0/6 · 1/6 | 2/6 | 0.117 | none ≤ 0.9 | 0.9 | none ≤ 0.9 | 0.19 |

#### LLE: real detections 126/396 windows (33 subjects); U = 0.3889 (cluster bootstrap 0.3889, Clopper–Pearson 0.3589)

| family | parameter | λ | detected / base windows at f = 0.1, 0.3, 0.5, 0.9 | replacement | π_upper at f = 0.9 | smallest f, π < 0.05 | smallest f, π < 0.20 | (secondary: Wilson-lower p) π < 0.20 | (descriptive: Firth fit) π < 0.20 |
|---|---|---|---|---|---|---|---|---|---|
| Hénon | a1.08 | 0.136 | 1/6 · 1/6 · 2/6 · 4/6 | 1/6 | 0.583 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |
| Hénon | a1.14 | 0.244 | 1/6 · 4/6 · 5/6 · 6/6 | 4/6 | 0.389 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |
| Hénon | a1.22 | 0.303 | 1/6 · 5/6 · 5/6 · 6/6 | 6/6 | 0.389 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |
| Hénon | a1.40 | 0.419 | 1/6 · 4/6 · 5/6 · 6/6 | 6/6 | 0.389 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |
| logistic | r3.58 | 0.105 | 3/6 · 1/6 · 2/6 · 1/6 | 1/6 | 1.000 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |
| logistic | r3.65 | 0.255 | 2/6 · 3/6 · 2/6 · 6/6 | 0/6 | 0.389 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |
| logistic | r3.88 | 0.464 | 2/6 · 4/6 · 6/6 · 6/6 | 6/6 | 0.389 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |
| logistic | r4.00 | 0.693 | 0/6 · 5/6 · 5/6 · 6/6 | 6/6 | 0.389 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |
| coupled vdP | rho8_om3.3 | 0.0751 | 2/6 · 6/6 · 6/6 · 6/6 | 1/6 | 0.389 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |
| coupled vdP | rho6_om3.3 | 0.0788 | 1/6 · 5/6 · 5/6 · 6/6 | 6/6 | 0.389 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |
| coupled vdP | rho2_om2.7 | 0.085 | 1/6 · 4/6 · 5/6 · 6/6 | 6/6 | 0.389 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |
| coupled vdP | rho6_om4.0 | 0.101 | 4/6 · 5/6 · 5/6 · 6/6 | 3/6 | 0.389 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |
| phase-reset | tau1.14 | 0.0354 | 0/6 · 5/6 · 6/6 · 6/6 | 5/6 | 0.389 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |
| phase-reset | tau0.58 | 0.081 | 2/6 · 5/6 · 5/6 · 6/6 | 1/6 | 0.389 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |
| phase-reset | tau1.20 | 0.138 | 2/6 · 1/6 · 3/6 · 1/6 | 0/6 | 1.000 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |
| phase-reset | tau1.16 | 0.185 | 2/6 · 3/6 · 6/6 · 6/6 | 3/6 | 0.389 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |

#### UPO: real detections 8/396 windows (33 subjects); U = 0.0362 (cluster bootstrap 0.0328, Clopper–Pearson 0.0362)

| family | parameter | λ | detected / base windows at f = 0.1, 0.3, 0.5, 0.9 | replacement | π_upper at f = 0.9 | smallest f, π < 0.05 | smallest f, π < 0.20 | (secondary: Wilson-lower p) π < 0.20 | (descriptive: Firth fit) π < 0.20 |
|---|---|---|---|---|---|---|---|---|---|
| Hénon | a1.08 | 0.136 | 0/6 · 0/6 · 0/6 · 0/6 | 0/6 | 1.000 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |
| Hénon | a1.14 | 0.244 | 0/6 · 0/6 · 0/6 · 2/6 | 2/6 | 0.108 | none ≤ 0.9 | 0.9 | none ≤ 0.9 | 0.74 |
| Hénon | a1.22 | 0.303 | 0/6 · 0/6 · 2/6 · 3/6 | 4/6 | 0.072 | none ≤ 0.9 | 0.5 | 0.9 | 0.48 |
| Hénon | a1.40 | 0.419 | 0/6 · 0/6 · 0/6 · 2/6 | 5/6 | 0.108 | none ≤ 0.9 | 0.9 | none ≤ 0.9 | 0.74 |
| logistic | r3.58 | 0.105 | 0/6 · 0/6 · 0/6 · 0/6 | 0/6 | 1.000 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |
| logistic | r3.65 | 0.255 | 0/6 · 0/6 · 0/6 · 0/6 | 0/6 | 1.000 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |
| logistic | r3.88 | 0.464 | 0/6 · 0/6 · 0/6 · 4/6 | 5/6 | 0.054 | none ≤ 0.9 | 0.9 | 0.9 | 0.64 |
| logistic | r4.00 | 0.693 | 0/6 · 0/6 · 0/6 · 4/6 | 5/6 | 0.054 | none ≤ 0.9 | 0.9 | 0.9 | 0.64 |
| coupled vdP | rho8_om3.3 | 0.0751 | 0/6 · 0/6 · 0/6 · 0/6 | 0/6 | 1.000 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |
| coupled vdP | rho6_om3.3 | 0.0788 | 0/6 · 0/6 · 0/6 · 0/6 | 0/6 | 1.000 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |
| coupled vdP | rho2_om2.7 | 0.085 | 0/6 · 0/6 · 0/6 · 1/6 | 1/6 | 0.217 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | 0.87 |
| coupled vdP | rho6_om4.0 | 0.101 | 0/6 · 0/6 · 0/6 · 0/6 | 0/6 | 1.000 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |
| phase-reset | tau1.14 | 0.0354 | 0/6 · 0/6 · 0/6 · 1/6 | 0/6 | 0.217 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | 0.87 |
| phase-reset | tau0.58 | 0.081 | 0/6 · 0/6 · 0/6 · 0/6 | 0/6 | 1.000 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |
| phase-reset | tau1.20 | 0.138 | 0/6 · 0/6 · 1/6 · 2/6 | 1/6 | 0.108 | none ≤ 0.9 | 0.9 | none ≤ 0.9 | 0.63 |
| phase-reset | tau1.16 | 0.185 | 2/6 · 0/6 · 0/6 · 1/6 | 3/6 | 0.217 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 | none ≤ 0.9 |


## Q2 titration (12-min segments)

Segments: 3395 (31 not analysable: coverage < 90 % or < 400 intervals).

| arm | NSR mean subject DR [95 % CI] | CHF mean subject DR [95 % CI] | CHF − NSR [95 % CI] | permutation p |
|---|---|---|---|---|
| **raw (P2)** | 47.5 % [40.8, 55.1] (n = 18) | 80.5 % [65.4, 93.4] (n = 15) | 33.0 % [16.6, 47.3] | 0.0002 |
| masked (analysable) | 45.2 % [38.2, 53.0] (n = 18) | 21.6 % [11.6, 34.4] (n = 15) | -23.6 % [-36.1, -9.0] | 0.0015 |
| edited (eligible) | 46.0 % [38.8, 53.9] (n = 18) | 37.8 % [24.3, 52.4] (n = 15) | -8.2 % [-23.4, 8.4] | 0.3142 |
| Wu NN preprocessing | 46.0 % [38.9, 53.9] (n = 18) | 35.9 % [22.9, 49.6] (n = 15) | -10.1 % [-24.8, 5.5] | 0.1945 |

| change vs raw | group | segments | raw DR | arm DR | mean change [95 % CI] | raw positives removed |
|---|---|---|---|---|---|---|
| **masked (P2)** | NSR | 1889 | 47.6 % | 45.2 % | -2.4 % [-3.5, -1.4] | 6.1 % |
| **masked (P2)** | CHF | 954 | 79.9 % | 21.6 % | -58.3 % [-72.2, -43.2] | 74.0 % |
| masked, non-analysable = negative | NSR | 1891 | 47.5 % | 45.1 % | -2.4 % [-3.5, -1.4] | 6.1 % |
| masked, non-analysable = negative | CHF | 1473 | 80.5 % | 11.9 % | -68.6 % [-82.1, -53.7] | 85.2 % |
| edited | NSR | 1891 | 47.5 % | 46.0 % | -1.6 % [-2.6, -0.7] | 3.6 % |
| edited | CHF | 1406 | 80.5 % | 37.8 % | -42.7 % [-56.6, -28.7] | 54.4 % |
| Wu NN | NSR | 1891 | 47.5 % | 46.0 % | -1.5 % [-2.5, -0.6] | 3.5 % |
| Wu NN | CHF | 1462 | 80.5 % | 35.9 % | -44.7 % [-58.0, -31.4] | 56.1 % |

Masked arm analysable fraction: {"NSR": 0.999, "CHF": 0.648}; edited eligible fraction: {"NSR": 1.0, "CHF": 0.955}.
Mean NL among positive segments, raw: {"NSR": 0.172, "CHF": 0.706}; Wu: {"NSR": 0.168, "CHF": 0.36}.
Raw DR night / day (clock time approximate): {"NSR": {"night": 0.7008928571428571, "day": 0.41094941094941095}, "CHF": {"night": 0.8189910979228486, "day": 0.8001760563380281}}.
Segment GEE (raw positive ~ log2(1 + burden) + CHF): burden OR 3.14 [1.69, 5.83], p = 0.00029; CHF OR 2.91 [0.10, 88.75], p = 0.54.

## Q2e synthetic titration checks

| condition | n | windows | raw | raw at 1/128 s | masked (analysable) | edited (eligible) |
|---|---|---|---|---|---|---|
| E1_bigeminy | 512 | 5 | 0 | 0 | 0 (0) | 0 (0) |
| E1_bigeminy | 800 | 5 | 0 | 0 | 0 (0) | 0 (0) |
| E2_trigeminy | 512 | 5 | 0 | 0 | 0 (0) | 0 (0) |
| E2_trigeminy | 800 | 5 | 0 | 0 | 0 (0) | 0 (0) |
| E3_couplets_10pct | 512 | 5 | 5 | 0 | 0 (0) | 2 (5) |
| E3_couplets_10pct | 800 | 5 | 5 | 0 | 0 (0) | 2 (5) |
| E3_couplets_5pct | 512 | 5 | 5 | 0 | 0 (0) | 2 (5) |
| E3_couplets_5pct | 800 | 5 | 5 | 0 | 0 (0) | 3 (5) |
| E4_runs | 512 | 5 | 5 | 0 | 0 (2) | 0 (5) |
| E4_runs | 800 | 5 | 5 | 0 | 0 (5) | 0 (5) |
| E5_atrial_10pct | 512 | 5 | 0 | 0 | 0 (0) | 1 (5) |
| E5_atrial_10pct | 800 | 5 | 0 | 0 | 0 (0) | 1 (5) |
| E5_atrial_5pct | 512 | 5 | 0 | 0 | 0 (0) | 0 (5) |
| E5_atrial_5pct | 800 | 5 | 0 | 0 | 0 (0) | 0 (5) |
| E6_ectopic10_trend | 512 | 5 | 0 | 0 | 0 (0) | 0 (5) |
| E6_ectopic10_trend | 800 | 5 | 0 | 0 | 0 (0) | 0 (5) |
| N1_linear_rr | 512 | 5 | 0 | 0 | 0 (5) | 0 (5) |
| N1_linear_rr | 800 | 5 | 0 | 0 | 0 (5) | 0 (5) |
| N2_power_law | 512 | 5 | 0 | 0 | 0 (5) | 0 (5) |
| N2_power_law | 800 | 5 | 0 | 0 | 0 (5) | 0 (5) |
| N3_linear_rr_trend | 512 | 5 | 0 | 0 | 0 (5) | 0 (5) |
| N3_linear_rr_trend | 800 | 5 | 0 | 0 | 0 (5) | 0 (5) |
| N4_linear_rr_step | 512 | 5 | 1 | 0 | 1 (5) | 1 (5) |
| N4_linear_rr_step | 800 | 5 | 0 | 0 | 0 (5) | 0 (5) |
| N5_linear_rr_warped | 512 | 5 | 5 | 0 | 5 (5) | 5 (5) |
| N5_linear_rr_warped | 800 | 5 | 5 | 0 | 5 (5) | 5 (5) |
| N6_noisy_rsa | 512 | 5 | 0 | 0 | 0 (5) | 0 (5) |
| N6_noisy_rsa | 800 | 5 | 0 | 0 | 0 (5) | 0 (5) |
| S1_setar | 512 | 5 | 0 | 0 | 0 (5) | 0 (5) |
| S1_setar | 800 | 5 | 0 | 0 | 0 (5) | 0 (5) |
| S2_ectopic_10pct | 512 | 5 | 0 | 0 | 0 (0) | 1 (5) |
| S2_ectopic_10pct | 800 | 5 | 0 | 0 | 0 (0) | 0 (5) |
| S2_ectopic_5pct | 512 | 5 | 0 | 0 | 0 (0) | 2 (5) |
| S2_ectopic_5pct | 800 | 5 | 0 | 0 | 0 (0) | 2 (5) |
| vdp26_E1 | 512 | 5 | 1 | 0 | 0 (0) | 0 (0) |
| vdp26_E1 | 800 | 5 | 2 | 0 | 0 (0) | 0 (0) |
| vdp26_E3_10 | 512 | 5 | 5 | 0 | 0 (0) | 5 (5) |
| vdp26_E3_10 | 800 | 5 | 5 | 0 | 0 (0) | 5 (5) |
| vdp26_S2_5 | 512 | 5 | 0 | 0 | 0 (0) | 5 (5) |
| vdp26_S2_5 | 800 | 5 | 0 | 0 | 0 (0) | 5 (5) |
| vdp26_none | 512 | 5 | 5 | 0 | 5 (5) | 5 (5) |
| vdp26_none | 800 | 5 | 5 | 0 | 5 (5) | 5 (5) |
| vdp27_E1 | 512 | 5 | 0 | 0 | 0 (0) | 0 (0) |
| vdp27_E1 | 800 | 5 | 0 | 0 | 0 (0) | 0 (0) |
| vdp27_E3_10 | 512 | 5 | 0 | 0 | 0 (0) | 2 (5) |
| vdp27_E3_10 | 800 | 5 | 0 | 0 | 0 (0) | 1 (5) |
| vdp27_S2_5 | 512 | 5 | 0 | 0 | 0 (0) | 1 (5) |
| vdp27_S2_5 | 800 | 5 | 0 | 0 | 0 (0) | 0 (5) |
| vdp27_none | 512 | 5 | 0 | 0 | 0 (5) | 0 (5) |
| vdp27_none | 800 | 5 | 0 | 0 | 0 (5) | 0 (5) |
| vdp28_E1 | 512 | 5 | 0 | 0 | 0 (0) | 0 (0) |
| vdp28_E1 | 800 | 5 | 0 | 0 | 0 (0) | 0 (0) |
| vdp28_E3_10 | 512 | 5 | 5 | 0 | 0 (0) | 5 (5) |
| vdp28_E3_10 | 800 | 5 | 5 | 0 | 0 (0) | 5 (5) |
| vdp28_S2_5 | 512 | 5 | 0 | 0 | 0 (0) | 2 (5) |
| vdp28_S2_5 | 800 | 5 | 0 | 0 | 0 (0) | 4 (5) |
| vdp28_none | 512 | 5 | 0 | 0 | 0 (5) | 0 (5) |
| vdp28_none | 800 | 5 | 0 | 0 | 0 (5) | 0 (5) |
| vdp29_E1 | 512 | 5 | 0 | 0 | 0 (0) | 0 (0) |
| vdp29_E1 | 800 | 5 | 0 | 0 | 0 (0) | 0 (0) |
| vdp29_E3_10 | 512 | 5 | 5 | 0 | 0 (0) | 3 (5) |
| vdp29_E3_10 | 800 | 5 | 5 | 0 | 0 (0) | 3 (5) |
| vdp29_S2_5 | 512 | 5 | 0 | 0 | 0 (0) | 3 (5) |
| vdp29_S2_5 | 800 | 5 | 0 | 0 | 0 (0) | 3 (5) |
| vdp29_none | 512 | 5 | 0 | 0 | 0 (5) | 0 (5) |
| vdp29_none | 800 | 5 | 0 | 0 | 0 (5) | 0 (5) |
| vdp30_E1 | 512 | 5 | 0 | 0 | 0 (0) | 0 (0) |
| vdp30_E1 | 800 | 5 | 0 | 0 | 0 (0) | 0 (0) |
| vdp30_E3_10 | 512 | 5 | 0 | 0 | 0 (0) | 0 (5) |
| vdp30_E3_10 | 800 | 5 | 0 | 0 | 0 (0) | 0 (5) |
| vdp30_S2_5 | 512 | 5 | 0 | 0 | 0 (0) | 0 (5) |
| vdp30_S2_5 | 800 | 5 | 0 | 0 | 0 (0) | 0 (5) |
| vdp30_none | 512 | 5 | 0 | 0 | 0 (5) | 0 (5) |
| vdp30_none | 800 | 5 | 0 | 0 | 0 (5) | 0 (5) |
| vdp31_E1 | 512 | 5 | 0 | 0 | 0 (0) | 0 (0) |
| vdp31_E1 | 800 | 5 | 0 | 0 | 0 (0) | 0 (0) |
| vdp31_E3_10 | 512 | 5 | 5 | 0 | 0 (0) | 5 (5) |
| vdp31_E3_10 | 800 | 5 | 5 | 0 | 0 (0) | 5 (5) |
| vdp31_S2_5 | 512 | 5 | 1 | 0 | 0 (0) | 5 (5) |
| vdp31_S2_5 | 800 | 5 | 0 | 0 | 0 (0) | 5 (5) |
| vdp31_none | 512 | 5 | 0 | 0 | 0 (5) | 0 (5) |
| vdp31_none | 800 | 5 | 0 | 0 | 0 (5) | 0 (5) |

## Q3 ectopy dose-response

| outcome | windows | positive | OR per doubling of 1 + burden [95 % CI] | p | CHF OR (adjusted for burden) |
|---|---|---|---|---|---|
| **titration (P3)** | 396 | 199 | 2.43 [1.24, 4.78] | 0.0099 | 0.78 [0.34, 1.83], p = 0.57 |
| **K3 (P3)** | 396 | 0 | not estimable (fewer than 10 positives or negatives (or no burden variation)) | | |
| LLE alone | 396 | 126 | 1.33 [1.16, 1.52] | 3.1e-05 | 1.16 [0.56, 2.40], p = 0.7 |
| UPO alone | 396 | 8 | not estimable (fewer than 10 positives or negatives (or no burden variation)) | | |
| K1 | 396 | 3 | not estimable (fewer than 10 positives or negatives (or no burden variation)) | | |

| burden (beats) | windows | subjects | titration | LLE | UPO | K1 | K3 |
|---|---|---|---|---|---|---|---|
| 0 | 235 | 26 | 71 (30 %) | 55 (23 %) | 6 (3 %) | 2 (1 %) | 0 (0 %) |
| 1 | 42 | 19 | 25 (60 %) | 15 (36 %) | 1 (2 %) | 0 (0 %) | 0 (0 %) |
| 2-4 | 37 | 9 | 28 (76 %) | 10 (27 %) | 0 (0 %) | 0 (0 %) | 0 (0 %) |
| 5-15 | 47 | 10 | 46 (98 %) | 20 (43 %) | 0 (0 %) | 0 (0 %) | 0 (0 %) |
| >=16 | 35 | 7 | 29 (83 %) | 26 (74 %) | 1 (3 %) | 1 (3 %) | 0 (0 %) |

| paired comparison (3b) | subset | windows | raw positive | other arm positive | both | raw positives removed | mean subject-level change [95 % CI] |
|---|---|---|---|---|---|---|---|
| TIT_masked | all | 288 | 103 | 78 | 77 | 25.2 % | -14.0 % [-23.1, -6.7] |
| TIT_masked | n_masked_ge1 | 73 | 38 | 13 | 12 | 68.4 % | -33.2 % [-48.2, -19.4] |
| TIT_edited | all | 386 | 190 | 118 | 113 | 40.5 % | -19.7 % [-30.3, -9.8] |
| TIT_edited | n_masked_ge1 | 171 | 125 | 53 | 48 | 61.6 % | -33.0 % [-46.6, -19.9] |
| K1_edited | all | 386 | 3 | 2 | 2 | 33.3 % | -0.3 % [-0.8, 0.0] |
| K1_edited | n_masked_ge1 | 171 | 1 | 0 | 0 | 100.0 % | -0.3 % [-1.0, 0.0] |
| LLE_edited | all | 386 | 118 | 94 | 81 | 31.4 % | -6.1 % [-12.9, -0.5] |
| LLE_edited | n_masked_ge1 | 171 | 67 | 43 | 30 | 55.2 % | -11.8 % [-23.0, -1.8] |
| UPO_edited | all | 386 | 8 | 10 | 7 | 12.5 % | 1.3 % [-0.3, 3.8] |
| UPO_edited | n_masked_ge1 | 171 | 2 | 4 | 1 | 50.0 % | 1.5 % [-0.4, 4.2] |

Q3c: titration CHF OR unadjusted 3.32 [1.61, 6.84], p = 0.0012.

## Q4 robustness

| detector | windows | default + / 499 + | both | default only | 499 only | agreement | κ | McNemar p |
|---|---|---|---|---|---|---|---|---|
| K1 | 36 | 0 / 0 | 0 | 0 | 0 | 100.0 % | – | 1 |
| LLE | 36 | 7 / 7 | 7 | 0 | 0 | 100.0 % | 1.00 | 1 |
| UPO | 36 | 1 / 1 | 1 | 0 | 0 | 100.0 % | 1.00 | 1 |
| K3 | 36 | 0 / 0 | 0 | 0 | 0 | 100.0 % | – | 1 |

| detector | windows (all 3 lengths) | 256 | 512 | 1024 |
|---|---|---|---|---|
| K1 | 36 | 0.0 % | 0.0 % | 5.6 % |
| LLE | 36 | 16.7 % | 19.4 % | 38.9 % |
| UPO | 36 | 0.0 % | 2.8 % | 11.1 % |
| K3 | 36 | 0.0 % | 0.0 % | 2.8 % |

| detector | NSR night | NSR day | CHF night | CHF day |
|---|---|---|---|---|
| TIT | 29/49 | 51/167 | 30/44 | 89/136 |
| LLE | 8/49 | 43/167 | 19/44 | 56/136 |
| UPO | 3/49 | 3/167 | 1/44 | 1/136 |
| K1 | 1/49 | 1/167 | 1/44 | 0/136 |
| K3 | 0/49 | 0/167 | 0/44 | 0/136 |
| K4 | 0/49 | 0/167 | 0/44 | 0/136 |

Titration GEE night effect (adjusted for burden and group): OR 2.06 [1.15, 3.71], p = 0.016.

Q4d rates: {"K3": {"NSR": 0.0, "CHF": 0.0}, "K4": {"NSR": 0.0, "CHF": 0.0}, "K3RR": {"NSR": 0.0, "CHF": 0.011111111111111112}}; K3 vs K4 agreement {"both": 0, "a_only": 0, "b_only": 0, "neither": 396}; K3 vs K3RR {"both": 0, "a_only": 0, "b_only": 2, "neither": 394}.

## Q3 MIT-BIH (development, exploratory)

| outcome | windows | positive | OR per doubling of 1 + burden [95 % CI] | p | CHF OR (adjusted for burden) |
|---|---|---|---|---|---|
| **titration (P3)** | 89 | 37 | 1.23 [1.03, 1.46] | 0.022 | – |
| **K3 (P3)** | 89 | 0 | not estimable (fewer than 10 positives or negatives (or no burden variation)) | | |
| LLE alone | 89 | 31 | 1.40 [1.13, 1.72] | 0.0017 | – |
| UPO alone | 89 | 0 | not estimable (fewer than 10 positives or negatives (or no burden variation)) | | |
| K1 | 89 | 0 | not estimable (fewer than 10 positives or negatives (or no burden variation)) | | |

| burden (beats) | windows | subjects | titration | LLE | UPO | K1 | K3 |
|---|---|---|---|---|---|---|---|
| 0 | 17 | 11 | 1 (6 %) | 0 (0 %) | 0 (0 %) | 0 (0 %) | 0 (0 %) |
| 1 | 4 | 4 | 2 (50 %) | 1 (25 %) | 0 (0 %) | 0 (0 %) | 0 (0 %) |
| 2-4 | 6 | 5 | 2 (33 %) | 1 (17 %) | 0 (0 %) | 0 (0 %) | 0 (0 %) |
| 5-15 | 9 | 5 | 3 (33 %) | 3 (33 %) | 0 (0 %) | 0 (0 %) | 0 (0 %) |
| >=16 | 53 | 26 | 29 (55 %) | 26 (49 %) | 0 (0 %) | 0 (0 %) | 0 (0 %) |

| paired comparison (3b) | subset | windows | raw positive | other arm positive | both | raw positives removed | mean subject-level change [95 % CI] |
|---|---|---|---|---|---|---|---|
| TIT_masked | all | 23 | 5 | 1 | 1 | 80.0 % | -23.1 % [-46.2, -3.8] |
| TIT_masked | n_masked_ge1 | 6 | 4 | 0 | 0 | 100.0 % | -66.7 % [-100.0, -33.3] |
| TIT_edited | all | 56 | 19 | 3 | 2 | 89.5 % | -28.5 % [-44.6, -14.0] |
| TIT_edited | n_masked_ge1 | 39 | 18 | 2 | 1 | 94.4 % | -39.3 % [-58.7, -20.7] |
| K1_edited | all | 56 | 0 | 0 | 0 | – | 0.0 % [0.0, 0.0] |
| K1_edited | n_masked_ge1 | 39 | 0 | 0 | 0 | – | 0.0 % [0.0, 0.0] |
| LLE_edited | all | 56 | 12 | 3 | 2 | 83.3 % | -16.1 % [-29.0, -4.8] |
| LLE_edited | n_masked_ge1 | 39 | 12 | 3 | 2 | 83.3 % | -20.0 % [-36.0, -6.0] |
| UPO_edited | all | 56 | 0 | 1 | 0 | – | 3.2 % [0.0, 9.7] |
| UPO_edited | n_masked_ge1 | 39 | 0 | 1 | 0 | – | 4.0 % [0.0, 12.0] |
