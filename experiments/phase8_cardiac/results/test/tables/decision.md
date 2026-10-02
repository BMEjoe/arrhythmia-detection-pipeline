# Phase 8 primary decision

Passing candidates: ['c1_frozen512']
WINNER: **c1_frozen512**

# Phase 8 TEST results (decide.py)

## c1_frozen512 (window 512)

PASS: **True**; failed conditions: []

Pooled TEST chaotic at (vi): 66/990; at (iv): 10/330

| regime                         |   N |   k |   errors |   limit | within   |
|:-------------------------------|----:|----:|---------:|--------:|:---------|
| E1_bigeminy                    | 100 |   0 |        0 |       7 | True     |
| E2_trigeminy                   | 100 |   0 |        0 |       7 | True     |
| E3_couplets_10pct              | 100 |   1 |        0 |       7 | True     |
| E3_couplets_5pct               | 100 |   0 |        0 |       7 | True     |
| E4_runs                        | 100 |   0 |        0 |       7 | True     |
| E5_atrial_10pct                | 100 |   5 |        0 |       7 | True     |
| E5_atrial_5pct                 | 100 |   0 |        0 |       7 | True     |
| E6_ectopic10_trend             | 100 |   2 |        0 |       7 | True     |
| N1_linear_rr                   | 100 |   0 |        0 |       7 | True     |
| N2_power_law                   | 100 |   0 |        0 |       7 | True     |
| N3_linear_rr_trend             | 100 |   0 |        0 |       7 | True     |
| N4_linear_rr_step              | 100 |   0 |        0 |       7 | True     |
| N5_linear_rr_warped            | 100 |   0 |        0 |       7 | True     |
| N6_noisy_rsa                   | 100 |   0 |        0 |       7 | True     |
| S1_setar                       | 100 |   0 |        0 |       7 | True     |
| S2_ectopic_10pct               | 100 |   0 |        0 |       7 | True     |
| S2_ectopic_5pct                | 100 |   0 |        0 |       7 | True     |
| av_node:H=40.0                 | 100 |   0 |        0 |       7 | True     |
| av_node:H=45.0                 | 100 |   0 |        0 |       7 | True     |
| av_node:H=52.0                 | 100 |   0 |        0 |       7 | True     |
| av_node:H=55.0                 | 100 |   0 |        0 |       7 | True     |
| coupled_vdp:rho=10.0,omega=3.3 | 100 |   0 |        0 |       7 | True     |
| coupled_vdp:rho=2.0,omega=5.6  | 100 |   0 |        0 |       7 | True     |
| coupled_vdp:rho=4.0,omega=5.6  | 100 |   0 |        0 |       7 | True     |
| coupled_vdp:rho=5.45,omega=5.6 | 100 |   0 |        0 |       7 | True     |
| coupled_vdp:rho=6.0,omega=5.6  | 100 |   0 |        0 |       7 | True     |
| coupled_vdp:rho=9.6,omega=2.1  | 100 |   0 |        0 |       7 | True     |
| mackey_glass:tau=14.0          | 100 |   0 |        0 |       7 | True     |
| mackey_glass:tau=15.0          | 100 |   0 |        0 |       7 | True     |
| mackey_glass:tau=16.0          | 100 |   0 |        0 |       7 | True     |
| phase_reset:tau=0.38           | 100 |   0 |        0 |       7 | True     |
| phase_reset:tau=0.62           | 100 |   0 |        0 |       7 | True     |
| phase_reset:tau=0.66           | 100 |   0 |        0 |       7 | True     |
| phase_reset:tau=1.18           | 100 |   0 |        0 |       7 | True     |
| phase_reset:tau=1.74           | 100 |   0 |        0 |       7 | True     |

## c2_nlp_iaaft512 (window 512)

PASS: **False**; failed conditions: ['coupled_vdp:rho=2.0,omega=5.6', 'coupled_vdp:rho=5.45,omega=5.6']

Pooled TEST chaotic at (vi): 220/990; at (iv): 141/330

| regime                         |   N |   k |   errors |   limit | within   |
|:-------------------------------|----:|----:|---------:|--------:|:---------|
| E1_bigeminy                    | 100 |   0 |        0 |       7 | True     |
| E2_trigeminy                   | 100 |   0 |        0 |       7 | True     |
| E3_couplets_10pct              | 100 |   0 |        0 |       7 | True     |
| E3_couplets_5pct               | 100 |   2 |        0 |       7 | True     |
| E4_runs                        | 100 |   0 |        0 |       7 | True     |
| E5_atrial_10pct                | 100 |   0 |        0 |       7 | True     |
| E5_atrial_5pct                 | 100 |   0 |        0 |       7 | True     |
| E6_ectopic10_trend             | 100 |   0 |        0 |       7 | True     |
| N1_linear_rr                   | 100 |   0 |        0 |       7 | True     |
| N2_power_law                   | 100 |   0 |        0 |       7 | True     |
| N3_linear_rr_trend             | 100 |   0 |        0 |       7 | True     |
| N4_linear_rr_step              | 100 |   0 |        0 |       7 | True     |
| N5_linear_rr_warped            | 100 |   0 |        0 |       7 | True     |
| N6_noisy_rsa                   | 100 |   0 |        0 |       7 | True     |
| S1_setar                       | 100 |   0 |        0 |       7 | True     |
| S2_ectopic_10pct               | 100 |   0 |        0 |       7 | True     |
| S2_ectopic_5pct                | 100 |   0 |        0 |       7 | True     |
| av_node:H=40.0                 | 100 |   2 |        0 |       7 | True     |
| av_node:H=45.0                 | 100 |   0 |        0 |       7 | True     |
| av_node:H=52.0                 | 100 |   0 |        0 |       7 | True     |
| av_node:H=55.0                 | 100 |   0 |        0 |       7 | True     |
| coupled_vdp:rho=10.0,omega=3.3 | 100 |   0 |        0 |       7 | True     |
| coupled_vdp:rho=2.0,omega=5.6  | 100 |  24 |        0 |       7 | False    |
| coupled_vdp:rho=4.0,omega=5.6  | 100 |   0 |        0 |       7 | True     |
| coupled_vdp:rho=5.45,omega=5.6 | 100 |  19 |        0 |       7 | False    |
| coupled_vdp:rho=6.0,omega=5.6  | 100 |   0 |        0 |       7 | True     |
| coupled_vdp:rho=9.6,omega=2.1  | 100 |   0 |        0 |       7 | True     |
| mackey_glass:tau=14.0          | 100 |   0 |        0 |       7 | True     |
| mackey_glass:tau=15.0          | 100 |   0 |        0 |       7 | True     |
| mackey_glass:tau=16.0          | 100 |   0 |        0 |       7 | True     |
| phase_reset:tau=0.38           | 100 |   0 |        0 |       7 | True     |
| phase_reset:tau=0.62           | 100 |   0 |        0 |       7 | True     |
| phase_reset:tau=0.66           | 100 |   0 |        0 |       7 | True     |
| phase_reset:tau=1.18           | 100 |   0 |        0 |       7 | True     |
| phase_reset:tau=1.74           | 100 |   0 |        0 |       7 | True     |

## c3_nlp_ep256 (window 256)

PASS: **False**; failed conditions: ['coupled_vdp:rho=2.0,omega=5.6']

Pooled TEST chaotic at (vi): 148/990; at (iv): 73/330

| regime                         |   N |   k |   errors |   limit | within   |
|:-------------------------------|----:|----:|---------:|--------:|:---------|
| E1_bigeminy                    | 100 |   0 |        0 |       7 | True     |
| E2_trigeminy                   | 100 |   0 |        0 |       7 | True     |
| E3_couplets_10pct              | 100 |   0 |        0 |       7 | True     |
| E3_couplets_5pct               | 100 |   0 |        0 |       7 | True     |
| E4_runs                        | 100 |   0 |        0 |       7 | True     |
| E5_atrial_10pct                | 100 |   0 |        0 |       7 | True     |
| E5_atrial_5pct                 | 100 |   0 |        0 |       7 | True     |
| E6_ectopic10_trend             | 100 |   0 |        0 |       7 | True     |
| N1_linear_rr                   | 100 |   0 |        0 |       7 | True     |
| N2_power_law                   | 100 |   0 |        0 |       7 | True     |
| N3_linear_rr_trend             | 100 |   0 |        0 |       7 | True     |
| N4_linear_rr_step              | 100 |   0 |        0 |       7 | True     |
| N5_linear_rr_warped            | 100 |   0 |        0 |       7 | True     |
| N6_noisy_rsa                   | 100 |   0 |        0 |       7 | True     |
| S1_setar                       | 100 |   0 |        0 |       7 | True     |
| S2_ectopic_10pct               | 100 |   1 |        0 |       7 | True     |
| S2_ectopic_5pct                | 100 |   0 |        0 |       7 | True     |
| av_node:H=40.0                 | 100 |   0 |        0 |       7 | True     |
| av_node:H=45.0                 | 100 |   0 |        0 |       7 | True     |
| av_node:H=52.0                 | 100 |   0 |        0 |       7 | True     |
| av_node:H=55.0                 | 100 |   0 |        0 |       7 | True     |
| coupled_vdp:rho=10.0,omega=3.3 | 100 |   0 |        0 |       7 | True     |
| coupled_vdp:rho=2.0,omega=5.6  | 100 |  62 |        0 |       7 | False    |
| coupled_vdp:rho=4.0,omega=5.6  | 100 |   0 |        0 |       7 | True     |
| coupled_vdp:rho=5.45,omega=5.6 | 100 |   0 |        0 |       7 | True     |
| coupled_vdp:rho=6.0,omega=5.6  | 100 |   0 |        0 |       7 | True     |
| coupled_vdp:rho=9.6,omega=2.1  | 100 |   0 |        0 |       7 | True     |
| mackey_glass:tau=14.0          | 100 |   1 |        0 |       7 | True     |
| mackey_glass:tau=15.0          | 100 |   0 |        0 |       7 | True     |
| mackey_glass:tau=16.0          | 100 |   0 |        0 |       7 | True     |
| phase_reset:tau=0.38           | 100 |   0 |        0 |       7 | True     |
| phase_reset:tau=0.62           | 100 |   0 |        0 |       7 | True     |
| phase_reset:tau=0.66           | 100 |   0 |        0 |       7 | True     |
| phase_reset:tau=1.18           | 100 |   0 |        0 |       7 | True     |
| phase_reset:tau=1.74           | 100 |   0 |        0 |       7 | True     |

## c4_titration256 (window 256)

PASS: **False**; failed conditions: ['E1_bigeminy', 'E2_trigeminy', 'E3_couplets_10pct', 'E3_couplets_5pct', 'E4_runs', 'E5_atrial_10pct', 'E5_atrial_5pct', 'E6_ectopic10_trend', 'N4_linear_rr_step', 'N5_linear_rr_warped', 'S1_setar', 'S2_ectopic_10pct', 'S2_ectopic_5pct', 'av_node:H=40.0', 'av_node:H=45.0', 'av_node:H=52.0', 'coupled_vdp:rho=10.0,omega=3.3', 'coupled_vdp:rho=2.0,omega=5.6', 'coupled_vdp:rho=4.0,omega=5.6', 'coupled_vdp:rho=5.45,omega=5.6', 'coupled_vdp:rho=6.0,omega=5.6', 'coupled_vdp:rho=9.6,omega=2.1', 'mackey_glass:tau=14.0', 'mackey_glass:tau=15.0', 'mackey_glass:tau=16.0', 'phase_reset:tau=0.38', 'phase_reset:tau=0.62', 'phase_reset:tau=0.66', 'phase_reset:tau=1.18', 'phase_reset:tau=1.74']

Pooled TEST chaotic at (vi): 966/990; at (iv): 330/330

| regime                         |   N |   k |   errors |   limit | within   |
|:-------------------------------|----:|----:|---------:|--------:|:---------|
| E1_bigeminy                    | 100 |  85 |        0 |       7 | False    |
| E2_trigeminy                   | 100 | 100 |        0 |       7 | False    |
| E3_couplets_10pct              | 100 | 100 |        0 |       7 | False    |
| E3_couplets_5pct               | 100 | 100 |        0 |       7 | False    |
| E4_runs                        | 100 | 100 |        0 |       7 | False    |
| E5_atrial_10pct                | 100 |  58 |        0 |       7 | False    |
| E5_atrial_5pct                 | 100 |  89 |        0 |       7 | False    |
| E6_ectopic10_trend             | 100 | 100 |        0 |       7 | False    |
| N1_linear_rr                   | 100 |   3 |        0 |       7 | True     |
| N2_power_law                   | 100 |   2 |        0 |       7 | True     |
| N3_linear_rr_trend             | 100 |   7 |        0 |       7 | True     |
| N4_linear_rr_step              | 100 |  98 |        0 |       7 | False    |
| N5_linear_rr_warped            | 100 | 100 |        0 |       7 | False    |
| N6_noisy_rsa                   | 100 |   2 |        0 |       7 | True     |
| S1_setar                       | 100 |  43 |        0 |       7 | False    |
| S2_ectopic_10pct               | 100 | 100 |        0 |       7 | False    |
| S2_ectopic_5pct                | 100 | 100 |        0 |       7 | False    |
| av_node:H=40.0                 | 100 |  56 |        0 |       7 | False    |
| av_node:H=45.0                 | 100 |  85 |        0 |       7 | False    |
| av_node:H=52.0                 | 100 |  84 |        0 |       7 | False    |
| av_node:H=55.0                 | 100 |   0 |        0 |       7 | True     |
| coupled_vdp:rho=10.0,omega=3.3 | 100 | 100 |        0 |       7 | False    |
| coupled_vdp:rho=2.0,omega=5.6  | 100 |  72 |        0 |       7 | False    |
| coupled_vdp:rho=4.0,omega=5.6  | 100 | 100 |        0 |       7 | False    |
| coupled_vdp:rho=5.45,omega=5.6 | 100 | 100 |        0 |       7 | False    |
| coupled_vdp:rho=6.0,omega=5.6  | 100 | 100 |        0 |       7 | False    |
| coupled_vdp:rho=9.6,omega=2.1  | 100 | 100 |        0 |       7 | False    |
| mackey_glass:tau=14.0          | 100 |  60 |        0 |       7 | False    |
| mackey_glass:tau=15.0          | 100 | 100 |        0 |       7 | False    |
| mackey_glass:tau=16.0          | 100 | 100 |        0 |       7 | False    |
| phase_reset:tau=0.38           | 100 |  74 |        0 |       7 | False    |
| phase_reset:tau=0.62           | 100 |  56 |        0 |       7 | False    |
| phase_reset:tau=0.66           | 100 |  65 |        0 |       7 | False    |
| phase_reset:tau=1.18           | 100 | 100 |        0 |       7 | False    |
| phase_reset:tau=1.74           | 100 | 100 |        0 |       7 | False    |

## baseline_frozen256 (window 256)

PASS: **True**; failed conditions: []

Pooled TEST chaotic at (vi): 52/990; at (iv): 10/330

| regime                         |   N |   k |   errors |   limit | within   |
|:-------------------------------|----:|----:|---------:|--------:|:---------|
| E1_bigeminy                    | 100 |   0 |        0 |       7 | True     |
| E2_trigeminy                   | 100 |   0 |        0 |       7 | True     |
| E3_couplets_10pct              | 100 |   0 |        0 |       7 | True     |
| E3_couplets_5pct               | 100 |   0 |        0 |       7 | True     |
| E4_runs                        | 100 |   0 |        0 |       7 | True     |
| E5_atrial_10pct                | 100 |   1 |        0 |       7 | True     |
| E5_atrial_5pct                 | 100 |   0 |        0 |       7 | True     |
| E6_ectopic10_trend             | 100 |   3 |        0 |       7 | True     |
| N1_linear_rr                   | 100 |   0 |        0 |       7 | True     |
| N2_power_law                   | 100 |   0 |        0 |       7 | True     |
| N3_linear_rr_trend             | 100 |   0 |        0 |       7 | True     |
| N4_linear_rr_step              | 100 |   0 |        0 |       7 | True     |
| N5_linear_rr_warped            | 100 |   0 |        0 |       7 | True     |
| N6_noisy_rsa                   | 100 |   0 |        0 |       7 | True     |
| S1_setar                       | 100 |   0 |        0 |       7 | True     |
| S2_ectopic_10pct               | 100 |   2 |        0 |       7 | True     |
| S2_ectopic_5pct                | 100 |   1 |        0 |       7 | True     |
| av_node:H=40.0                 | 100 |   0 |        0 |       7 | True     |
| av_node:H=45.0                 | 100 |   0 |        0 |       7 | True     |
| av_node:H=52.0                 | 100 |   0 |        0 |       7 | True     |
| av_node:H=55.0                 | 100 |   0 |        0 |       7 | True     |
| coupled_vdp:rho=10.0,omega=3.3 | 100 |   0 |        0 |       7 | True     |
| coupled_vdp:rho=2.0,omega=5.6  | 100 |   0 |        0 |       7 | True     |
| coupled_vdp:rho=4.0,omega=5.6  | 100 |   0 |        0 |       7 | True     |
| coupled_vdp:rho=5.45,omega=5.6 | 100 |   0 |        0 |       7 | True     |
| coupled_vdp:rho=6.0,omega=5.6  | 100 |   0 |        0 |       7 | True     |
| coupled_vdp:rho=9.6,omega=2.1  | 100 |   0 |        0 |       7 | True     |
| mackey_glass:tau=14.0          | 100 |   0 |        0 |       7 | True     |
| mackey_glass:tau=15.0          | 100 |   0 |        0 |       7 | True     |
| mackey_glass:tau=16.0          | 100 |   0 |        0 |       7 | True     |
| phase_reset:tau=0.38           | 100 |   0 |        0 |       7 | True     |
| phase_reset:tau=0.62           | 100 |   0 |        0 |       7 | True     |
| phase_reset:tau=0.66           | 100 |   0 |        0 |       7 | True     |
| phase_reset:tau=1.18           | 100 |   0 |        0 |       7 | True     |
| phase_reset:tau=1.74           | 100 |   0 |        0 |       7 | True     |

