# Step 0: analyze_segment path vs stored Phase 4 decisions

70 windows (14 Phase 4 conditions x test seeds 2000-2004, 256 samples). Differences: **0**.

| condition | UPO C2 detections (stored / analyze_segment) | lle_chaos_test detections (stored / analyze_segment) | identical windows |
|---|---|---|---|
| white_noise | 0 / 0 | 1 / 1 | 5/5 |
| ar1 | 0 / 0 | 0 / 0 | 5/5 |
| sinusoid | 0 / 0 | 0 / 0 | 5/5 |
| two_tone | 0 / 0 | 0 / 0 | 5/5 |
| logistic_p4 | 0 / 0 | 0 / 0 | 5/5 |
| logistic | 5 / 5 | 5 / 5 | 5/5 |
| henon | 5 / 5 | 5 / 5 | 5/5 |
| skewed_henon | 5 / 5 | 5 / 5 | 5/5 |
| logistic@30dB | 5 / 5 | 5 / 5 | 5/5 |
| logistic@20dB | 5 / 5 | 5 / 5 | 5/5 |
| logistic@10dB | 3 / 3 | 5 / 5 | 5/5 |
| henon@30dB | 5 / 5 | 5 / 5 | 5/5 |
| henon@20dB | 5 / 5 | 5 / 5 | 5/5 |
| henon@10dB | 3 / 3 | 5 / 5 | 5/5 |
