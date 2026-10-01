# Phase 8 Part B: development diagnosis (DEV material only)

Material:
- the DEV model (Mackey–Glass maxima intervals; 7 CHAOTIC and 3 NON-CHAOTIC regimes);
- the Phase 5–6 null generators (`nulls.py`);
- Phase 5 positives P1, G2 and G3.

Seeds are development seeds (0–99). No TEST model series was generated for any
detector. Raw results are in `results/dev/*.jsonl` (`dev_runner.py`).

## D1. Frozen Phase 6 detector (`frozen256.jsonl`; 256 beats, seeds 0–4)

**On Mackey–Glass maxima intervals:**

| label | variant | AND | LLE alone | UPO alone |
|---|---|---|---|---|
| CHAOTIC | clean | 23 % | 66 % | 23 % |
| CHAOTIC | quantized + 30 dB | 20 % | 71 % | 20 % |
| CHAOTIC | (vi) S2 5 % | 37 % | 89 % | 37 % |
| CHAOTIC | (vi) bigeminy | 0 % | 37 % | 0 % |
| NON-CHAOTIC | clean | 0 % | 0 % | 33 % |

**Why it fails.**
- **The UPO component is the bottleneck.** The LLE test fires on most chaotic MG windows,
  but the instability-gated UPO test fires only at τ = 23–26 (AND 29–71 %). It passes
  in 0–11 % of windows at τ = 16.5–20.
- **It is also not specific on MG periodic orbits.** UPO alone fires on 33 % of clean
  non-chaotic MG windows (τ = 15). The AND stays at 0 only because the LLE test does not
  fire on them.
- **An embedding search over the LLE test would therefore not raise AND power.** The
  limiting component is UPO, and Phase 5 Section 6 already showed that UPO power falls at
  m = 3 and 4.

**Decision.** The frozen-detector candidate drops the embedding search and keeps only the
longer-window variant (512 beats).

**Nulls:** AND 0/5 on every condition except 1/5 on S2 10 %.

**Flows:** G2 and G3 remain 0/5.

## D2. Nonlinear prediction vs IAAFT surrogates (`nulls256`, `nulls512`, `mg256`, `mg512`)

**Statistic:** the robust normalized leave-one-out local-average prediction error (k = 5,
Theiler 5, horizon 1), embeddings m = 2–5, 39 IAAFT surrogates, z = (mean_sur − obs) /
SD_sur, combined as zmax over m.

**Sensitivity is high, including to flows.**
- zmax medians at 256 / 512 beats:

| series | zmax, 256 | zmax, 512 |
|---|---:|---:|
| P1 Hénon | 18 | 25 |
| G2 Rössler | 12 | 16 |
| G3 Mackey–Glass | 12 | 16 |
| MG chaos, clean / quantized | about 16 | |

- This is the first statistic in the project that sees the sampled flows of Phase 5.

**Ectopy dominates the signal.**
- Linear RR with ectopy reaches zmax 4–13 (95th percentiles):
  - S2 5 %: 8.7 (256) and 12.9 (512);
  - E3 5 %: 10.4 and 13.5.
- MG periodic orbits **with** ectopy reach 12–16, as high as chaos with ectopy (about 10).

**Why it fails without a threshold.**
- The compensatory pause after a premature beat (x_{i+1} ≈ 2 RR − x_i) is genuine nonlinear
  determinism, which the IAAFT null destroys.
- A test of "determinism beyond linear-Gaussian" therefore cannot separate chaos plus
  ectopy from noise plus ectopy. Specificity has to come from a threshold on the effect
  size calibrated on the ectopy nulls, at the cost of power on ectopic chaos.

**Bigeminy (E1)** removes the signal in both directions:
- null zmax ≤ 3.7;
- chaos with bigeminy about 2, i.e. 0 % power.

Half of all intervals are replaced, which destroys the underlying map.

## D3. Ectopy-preserving IAAFT (`methods.ep_iaaft`)

**Construction.** Intervals flagged by the pipeline's own robust local outlier rule are
replaced by interpolation, the IAAFT surrogate is drawn from the cleaned series, and the
flagged intervals are re-inserted at their positions.

**Effect.** It lowers the null tail where flagging is good (E3 and E4 at 512: 95th
percentiles from 12–13 down to about 4–9). It does not help S2 much, because only about
40 % of ectopy-affected intervals are flagged.

**Dev-calibrated power** (threshold above every null window):
- 256 beats, threshold 9.3:

| variant | MG chaos power |
|---|---|
| clean | 80 % |
| quantized | 80 % |
| (vi) S2 5 % | 63 % |
| (vi) E3 10 % | 43 % |
| (vi) bigeminy | 0 % |

- IAAFT version at 512 beats, threshold 14.2: 86 %, 86 %, 57 %, 77 % and 0 %.

## D4. Noise titration (Poon–Barahona, as described in PNAS 2001 and arXiv:nlin/0606032)

**Implementation.** Volterra–Wiener–Korenberg polynomial autoregression with κ ≤ 6 and
d ≤ 3, C(r) = log ε + r/N, and a nested F-test at 1 %. NL is found by bisection on added
white noise, and detection means NL > 0.

**Sensitivity:** MG chaos detected in 89–100 % of windows in every variant; P1, G2 and G3
in 100 %.

**False positives on the nulls, at 256 beats over 20 seeds:**

| null | rate |
|---|---|
| every ectopy null | 75–100 % |
| N4 step | 100 % |
| N5 static warp | 100 % |
| S1 SETAR | 45 % |
| N1, N3 | 10 % |
| MG periodic orbits | 87–100 % |

**Why it fails.** The VWK test detects any departure from a linear AR model: static
nonlinearity, steps, ectopy and limit cycles. NL > 0 is therefore not specific to chaos.
This agrees with the published criticism that titration measures nonlinearity rather than
chaos. **It will fail the PASS rule.** It is still carried as a candidate so that its
false-positive rate is measured on TEST seeds, as the task requires.

## D5. Consequences for the candidates (fixed in PREREGISTRATION.md)

1. **C1 `frozen512`:** the frozen detector at 512 beats. This tests whether longer windows
   relieve the UPO bottleneck; there are no new parameters.
2. **C2 `nlp_iaaft512`:** prediction zmax vs IAAFT at 512 beats, with a threshold
   calibrated on development nulls only.
3. **C3 `nlp_ep256`:** prediction zmax vs ectopy-preserving IAAFT at 256 beats, with a
   calibrated threshold.
4. **C4 `titration256`:** Poon–Barahona noise titration, NL > 0.

**Threshold calibration** (C2, C3). The threshold is the smallest z such that each of the
17 PASS null conditions has at most 2 % of development windows (100 seeds, all lengths
generated by `nulls.py`) at or above it, plus a safety margin of 0.5. MG non-chaotic
regimes at (iv) are checked too.

## D6. Calibrated thresholds (calibrate.py; results/dev/calibration.json)
- C2 nlp_iaaft512: Z = 16.92 (binding null E3 couplets 5 %, 100 dev seeds); dev MG power
  clean 74 %, (iv) 77 %, (vi) S2 46 %, E3 51 %, E1 0 %; DEV non-chaotic MG at (iv) max zmax 5.2.
- C3 nlp_ep256: Z = 9.31 (binding S2 5 %); dev MG power clean 80 %, (iv) 80 %, (vi) S2 63 %,
  E3 43 %, E1 0 %; DEV non-chaotic MG at (iv) max 4.1.
