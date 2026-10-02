# Phase 8 close-out tables (from committed TEST results only)

Windows: {256: 20160, 512: 13440} method-window rows; errors: 0

## TEST CHAOTIC regimes by family and variant (chaos_primary)

|                       | c1_frozen512   | c2_nlp_iaaft512   | c3_nlp_ep256   | c4_titration256   | baseline_frozen256   |
|:----------------------|:---------------|:------------------|:---------------|:------------------|:---------------------|
| ('coupled_vdp', 'iv') | 4/120          | 39/120            | 38/120         | 120/120           | 6/120                |
| ('coupled_vdp', 'vi') | 29/360         | 109/360           | 29/360         | 354/360           | 35/360               |
| ('phase_reset', 'iv') | 6/210          | 102/210           | 35/210         | 210/210           | 4/210                |
| ('phase_reset', 'vi') | 37/630         | 111/630           | 119/630        | 612/630           | 17/630               |

## TEST CHAOTIC at (vi) by ectopy pattern

| variant   | c1_frozen512   | c2_nlp_iaaft512   | c3_nlp_ep256   | c4_titration256   | baseline_frozen256   |
|:----------|:---------------|:------------------|:---------------|:------------------|:---------------------|
| vi_E1     | 15/330         | 14/330            | 6/330          | 317/330           | 19/330               |
| vi_E3_10  | 26/330         | 118/330           | 65/330         | 327/330           | 17/330               |
| vi_S2_5   | 25/330         | 88/330            | 77/330         | 322/330           | 16/330               |

## TEST CHAOTIC by regime (iv + vi pooled, 120 windows each)

| regime                        | c1_frozen512   | c2_nlp_iaaft512   | c3_nlp_ep256   | c4_titration256   | baseline_frozen256   |
|:------------------------------|:---------------|:------------------|:---------------|:------------------|:---------------------|
| coupled_vdp:rho=2.0,omega=2.7 | 24/120         | 84/120            | 9/120          | 117/120           | 34/120               |
| coupled_vdp:rho=6.0,omega=3.3 | 6/120          | 37/120            | 4/120          | 120/120           | 2/120                |
| coupled_vdp:rho=6.0,omega=4.0 | 0/120          | 18/120            | 29/120         | 117/120           | 0/120                |
| coupled_vdp:rho=8.0,omega=3.3 | 3/120          | 9/120             | 25/120         | 120/120           | 5/120                |
| phase_reset:tau=0.58          | 1/120          | 47/120            | 64/120         | 120/120           | 4/120                |
| phase_reset:tau=0.6           | 37/120         | 90/120            | 90/120         | 119/120           | 14/120               |
| phase_reset:tau=1.14          | 0/120          | 2/120             | 0/120          | 120/120           | 0/120                |
| phase_reset:tau=1.16          | 5/120          | 11/120            | 0/120          | 110/120           | 3/120                |
| phase_reset:tau=1.2           | 0/120          | 33/120            | 0/120          | 114/120           | 0/120                |
| phase_reset:tau=1.58          | 0/120          | 0/120             | 0/120          | 120/120           | 0/120                |
| phase_reset:tau=1.6           | 0/120          | 30/120            | 0/120          | 119/120           | 0/120                |

## TEST CHAOTIC, secondary variants (chaos_second, 11 regimes x 10 seeds)

| variant    | c1_frozen512   | c2_nlp_iaaft512   | c3_nlp_ep256   | c4_titration256   | baseline_frozen256   |
|:-----------|:---------------|:------------------|:---------------|:------------------|:---------------------|
| i_clean    | 15/110         | 57/110            | 31/110         | 110/110           | 8/110                |
| ii_dyn_hi  | 28/110         | 55/110            | 27/110         | 90/110            | 6/110                |
| ii_dyn_lo  | 21/110         | 70/110            | 38/110         | 109/110           | 14/110               |
| iii_meas20 | 3/110          | 17/110            | 13/110         | 110/110           | 5/110                |
| iii_meas30 | 7/110          | 55/110            | 29/110         | 110/110           | 5/110                |
| v_E1       | 4/110          | 11/110            | 4/110          | 110/110           | 8/110                |
| v_E3_10    | 10/110         | 34/110            | 22/110         | 109/110           | 5/110                |
| v_S2_5     | 12/110         | 36/110            | 24/110         | 108/110           | 8/110                |

## DEV Mackey-Glass CHAOTIC at TEST seeds (dev_chaos)

| variant   | c1_frozen512   | c2_nlp_iaaft512   | c3_nlp_ep256   | c4_titration256   | baseline_frozen256   |
|:----------|:---------------|:------------------|:---------------|:------------------|:---------------------|
| iv_q      | 26/70          | 53/70             | 56/70          | 70/70             | 21/70                |
| vi_E1     | 0/70           | 0/70              | 0/70           | 64/70             | 0/70                 |
| vi_E3_10  | 18/70          | 29/70             | 29/70          | 70/70             | 16/70                |
| vi_S2_5   | 23/70          | 31/70             | 42/70          | 65/70             | 25/70                |

## NON-CHAOTIC regimes under ectopy + jitter (noncha_ect; secondary, not in the PASS rule)

|                              | c1_frozen512   | c2_nlp_iaaft512   | c3_nlp_ep256   | c4_titration256   | baseline_frozen256   |
|:-----------------------------|:---------------|:------------------|:---------------|:------------------|:---------------------|
| ('av_node', 'vi_E1')         | 0/40           | 0/40              | 0/40           | 17/40             | 0/40                 |
| ('av_node', 'vi_E3_10')      | 0/40           | 6/40              | 0/40           | 40/40             | 0/40                 |
| ('av_node', 'vi_S2_5')       | 0/40           | 3/40              | 0/40           | 35/40             | 0/40                 |
| ('coupled_vdp', 'vi_E1')     | 11/60          | 4/60              | 13/60          | 60/60             | 8/60                 |
| ('coupled_vdp', 'vi_E3_10')  | 7/60           | 18/60             | 10/60          | 60/60             | 5/60                 |
| ('coupled_vdp', 'vi_S2_5')   | 0/60           | 12/60             | 7/60           | 55/60             | 0/60                 |
| ('mackey_glass', 'vi_E1')    | 0/30           | 7/30              | 10/30          | 30/30             | 0/30                 |
| ('mackey_glass', 'vi_E3_10') | 0/30           | 30/30             | 28/30          | 30/30             | 0/30                 |
| ('mackey_glass', 'vi_S2_5')  | 0/30           | 22/30             | 26/30          | 29/30             | 0/30                 |
| ('phase_reset', 'vi_E1')     | 0/50           | 12/50             | 28/50          | 49/50             | 0/50                 |
| ('phase_reset', 'vi_E3_10')  | 3/50           | 27/50             | 22/50          | 50/50             | 4/50                 |
| ('phase_reset', 'vi_S2_5')   | 0/50           | 6/50              | 5/50           | 48/50             | 1/50                 |

## Phase 5 flows G2 / G3 (not out-of-sample)

| regime          | c1_frozen512   | c2_nlp_iaaft512   | c3_nlp_ep256   | c4_titration256   | baseline_frozen256   |
|:----------------|:---------------|:------------------|:---------------|:------------------|:---------------------|
| G2_rossler_flow | 0/100          | 49/100            | 100/100        | 100/100           | 0/100                |
| G3_mackey_glass | 0/100          | 59/100            | 100/100        | 100/100           | 0/100                |

## Pooled detections on PASS conditions

| group   | c1_frozen512   | c2_nlp_iaaft512   | c3_nlp_ep256   | c4_titration256   | baseline_frozen256   |
|:--------|:---------------|:------------------|:---------------|:------------------|:---------------------|
| noncha  | 0/1800         | 45/1800           | 63/1800        | 1452/1800         | 0/1800               |
| null    | 8/1700         | 2/1700            | 1/1700         | 1187/1700         | 7/1700               |

## Components of the frozen detector (C1 at 512, baseline at 256)

| grp              | method             |    N |   AND |   LLE |   UPO |
|:-----------------|:-------------------|-----:|------:|------:|------:|
| chaos_primary:iv | baseline_frozen256 |  330 |    10 |   143 |    19 |
| chaos_primary:iv | c1_frozen512       |  330 |    10 |   152 |    14 |
| chaos_primary:vi | baseline_frozen256 |  990 |    52 |   681 |    60 |
| chaos_primary:vi | c1_frozen512       |  990 |    66 |   740 |    74 |
| chaos_second     | baseline_frozen256 |  880 |    59 |   509 |    70 |
| chaos_second     | c1_frozen512       |  880 |   100 |   560 |   116 |
| dev_chaos        | baseline_frozen256 |  280 |    62 |   199 |    64 |
| dev_chaos        | c1_frozen512       |  280 |    67 |   219 |    67 |
| flows            | baseline_frozen256 |  200 |     0 |     0 |    38 |
| flows            | c1_frozen512       |  200 |     0 |     0 |    37 |
| noncha           | baseline_frozen256 | 1800 |     0 |    29 |   147 |
| noncha           | c1_frozen512       | 1800 |     0 |    23 |    94 |
| noncha_ect       | baseline_frozen256 |  540 |    18 |   240 |    33 |
| noncha_ect       | c1_frozen512       |  540 |    21 |   281 |    41 |
| null             | baseline_frozen256 | 1700 |     7 |   898 |    12 |
| null             | c1_frozen512       | 1700 |     8 |   929 |    11 |

## Median runtime per window (s, 4-worker load)

|   c1_frozen512 |   c2_nlp_iaaft512 |   c3_nlp_ep256 |   c4_titration256 |   baseline_frozen256 |
|---------------:|------------------:|---------------:|------------------:|---------------------:|
|           8.67 |              2.51 |           0.86 |              0.07 |                 3.71 |

