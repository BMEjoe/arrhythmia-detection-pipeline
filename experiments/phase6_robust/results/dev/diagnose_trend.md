# Part C1: trend mechanism (development seeds 0-29, C2 UPO analysis)

Mode drift = |So-transform mode in the last third - mode in the first third| (s); trend change = |trend difference between the centres of the first and last thirds| (s).

| system | variant | Level-B windows | gate (UPO) windows | median peak height | median mode drift (s) | median trend change (s) | median drift / trend change |
|---|---|---|---|---|---|---|---|
| henon | clean | 30/30 | 30/30 | 0.830 | 0.0000 | 0.0430 | 0.00 |
| henon | trend | 8/30 | 8/30 | 0.181 | 0.0453 | 0.0430 | 1.11 |
| henon | trend_linear_detrended | 30/30 | 30/30 | 0.655 | 0.0017 | 0.0430 | 0.04 |
| logistic | clean | 30/30 | 30/30 | 0.585 | 0.0018 | 0.0378 | 0.05 |
| logistic | trend | 12/30 | 12/30 | 0.157 | 0.0449 | 0.0378 | 1.08 |
| logistic | trend_linear_detrended | 28/30 | 28/30 | 0.394 | 0.0055 | 0.0378 | 0.11 |
