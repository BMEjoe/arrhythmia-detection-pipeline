# B3 waveform surrogate false positives (DEV, 2-min windows at 90 Hz)

The run was interrupted by a container restart after 8 of the planned 65 windows and was NOT resumed: every surrogate already exceeded the 7 % UNUSABLE limit on the strictly periodic null (B3 rule: > 7 % on ANY null => UNUSABLE), so the remaining windows could not change any verdict; the per-window cost (B5) is also recorded.

| null         | noise   |   windows | PPS h1   | PPS htau   | CS h1   | CS htau   | TS h1   | TS htau   |   median s/window |
|:-------------|:--------|----------:|:---------|:-----------|:--------|:----------|:--------|:----------|------------------:|
| periodic:0.8 | mix12   |         3 | 3/3      | 3/3        | 3/3     | 3/3       | 0/3     | 0/3       |               239 |
| periodic:0.8 | none    |         5 | 5/5      | 5/5        | 5/5     | 5/5       | 5/5     | 5/5       |               245 |
