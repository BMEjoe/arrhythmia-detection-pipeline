# Phase 4 summary -- dev seeds (explore)

det = windows with >= 1 detected peak (method rule), k/N [Wilson 95 %]; levelC = detected and Level-C verified_unstable; AND / OR = combined with fp.lle_chaos_test on the same window; lle = lle_chaos_test alone; loc = median error of the nearest detected peak to the analytic on-attractor fixed point.

## baseline

| N | condition | det | levelC | AND | OR | lle | loc (n) | status | s/run |
|---:|---|---|---:|---:|---:|---:|---|---|---:|
| 256 | white_noise | 3/36 = 0.083 [0.029, 0.218] | 0 | 0 | 3 | 0 |  | ok: 32; embedding_not_saturated: 4 | 5.7 |
