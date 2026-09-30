# Phase 5 preregistration: frozen combined detector on RR-interval-like signals

Committed and pushed **before any test window is generated**. Nothing here,
nor the detector, the generators, `decide.py` or `analysis.py`, changes after
the test run starts. Any change goes in a dated amendment at the end of this
file, with the reason. `run_phase5.py --phase test` reuses the Phase 3 guard:
it refuses to start unless this file is committed, unmodified and pushed, and
lists the method below.

- method: `and_detector`

## 1. Detector (frozen; evaluated, not tuned)

- `config = fp.phase4_upo_config(replace(fp.CFG, lle_chaos_test=True))`.
- Every window goes through `fp.analyze_segment(rr, config)` on raw RR
  (`use_corrected_rr_for_dynamics = False`), the path the MIT-BIH experiment
  uses. The components are:
  - **LLE**: `out["lle_chaos_test"]["detected"]`. This is the Phase 3 winner:
    τ = 1, m = 2, Rosenstein slope of steps 1–5, 99 IAAFT surrogates,
    p ≤ 0.05;
  - **UPO**: `out["upo"]["instability_gate_detected"]`. This is the Phase 4
    winner C2: fixed m = 2, M = 15, 50 AAFT surrogates, median instability gate
    at δ = 0.2.
- **PRIMARY = LLE AND UPO.** Each component alone and LLE OR UPO are reported
  as secondary outcomes.
- **Step 0 check.** Before this preregistration, this path reproduced the
  stored Phase 4 C2 and `lle_chaos_test` decisions and p-values on 70/70
  windows (`results/equivalence_check.md`).
- **Errors.** If `analyze_segment` raises, the window counts as NOT detected
  by every rule and the error is counted per condition. The decisions the two
  components would have made (direct calls) are also recorded and reported as
  "AND if errors were not counted". Development seeds 0–9 gave 0 errors in 240
  windows.

## 2. Data (frozen at commit 5a7bbfa, before any detector run on Phase 5 signals)

All generator parameters are in `config.py`. The generators are in
`systems.py` and are described with sources in `README.md`. Each window is
256 RR intervals in seconds, drawn with:
- μ ~ U(0.65, 0.95) s;
- σ ~ U(0.035, 0.060) s;
- quantization to 1/360 s after all transformations, unless a condition is
  marked unquantized.

| id | condition | construction |
|---|---|---|
| N1 | `N1_linear_rr` | McSharry (2003) bimodal spectrum: LF 0.1 Hz and HF 0.25 Hz, SD 0.01 Hz, LF/HF ~ U(1.5, 2.0), random phases; a 256-beat segment of 1024 |
| N2 | `N2_power_law` | Gaussian 1/f^β, β ~ U(0.9, 1.1); a 256-beat segment of 4096 |
| N3 | `N3_linear_rr_trend` | N1 + centred linear trend, total change U(5, 10) % of μ, random sign |
| N4 | `N4_linear_rr_step` | N1 + centred abrupt step, U(5, 10) % of μ, random sign, at a beat in [0.3N, 0.7N] |
| N5 | `N5_linear_rr_warped` | N1 through z → exp(0.5 z), then rescaled |
| N6 | `N6_noisy_rsa` | respiratory sinusoid 0.2–0.3 Hz (50 % of variance) + white Gaussian noise |
| S1 | `S1_setar` | SETAR(2;1,1): slope 0.7 if x ≤ 0, −0.5 otherwise, N(0, 1) noise (stable fixed point) |
| S2 | `S2_ectopic_{2,5,10}pct` | N1 + 5 / 13 / 26 isolated premature beats (c ~ U(0.6, 0.8)) with full compensatory pause |
| S3 | `S3_linear_rr_unq`, `S3_linear_rr_trend_unq`, `S3_linear_rr_step_unq` | N1 / N3 / N4 unquantized (the same draws as the quantized twins) |
| P1, P2 | `P1_henon_rr`, `P2_logistic_rr` | Hénon x (1.4, 0.3) / logistic r = 4, rescaled |
| P3 | `P3_henon_rr_trend`, `P3_logistic_rr_trend` | P1 / P2 + the N3 trend distribution |
| P4 | `P4_{henon,logistic}_rr_{30,20}dB` | P1 / P2 + white measurement noise at 30 / 20 dB of σ², before quantization |
| G1 | `G1_lorenz_maxima` | successive maxima of Lorenz z (10, 28, 8/3), RK4 dt 0.005 |
| G2 | `G2_rossler_flow` | Rössler x (0.2, 0.2, 5.7), RK4 dt 0.01, sampled every 1.0 |
| G3 | `G3_mackey_glass` | Mackey–Glass τ = 17 (0.2, 0.1, n = 10), RK4 dt 0.1, sampled every 6.0 |

**Seeds.** Phase 2E scheme `SeedSequence([20260926, code, 256, seed])`, with a
base stream and a modifier stream (`config.py`).
- **Development:** 0–99 for N* and S*, 0–29 for P* and G*.
- **Test:** starts at 3000:
  - **N1–N6: 300 test seeds each (3000–3299);**
  - **every other condition: 200 (3000–3199).**

## 3. PRIMARY rule (`decide.py`, committed in 1affece before this file)

**PASS** iff the AND detector (m = 2, through `analyze_segment`) fires on at
most floor(0.07 N) windows on EACH of N1–N6 (quantized, 256 RR intervals).
With N = 300 the limit is **21/300** (7 %; equivalent to 7/100).
Otherwise **FAIL**, reported for every failing condition.

**Power** (reported, not pass/fail): AND detections on P1–P4 and G1–G3.

## 4. Secondary outcomes (report only)

- The LLE and UPO components alone, and OR, on every condition.
- S1, S2 (2, 5, 10 %) and S3. Their results decide how MIT-BIH detections are
  interpreted:
  - S2 bears on windows with ectopic beats;
  - S1 on noise-driven nonlinearity;
  - S3 on quantization.
- **m-sensitivity:** every condition at m = 3 and m = 4, with BOTH
  `lle_chaos_test_m` and `upo_fixed_dimension` set to m. Each component, AND
  and OR are reported.
  - These arms call `fp.lle_chaos_test(rr, cfg)` and
    `fp.run_upo_analysis(rr, cfg, rng=default_rng(cfg.random_seed))` directly.
    This is what `analyze_segment` does for these components.
  - The direct path matched `analyze_segment` at m = 2 on 48/48 development
    windows (`results/dev/m_path_check.json`).
- `analyze_segment` error counts, UPO status distribution and runtime.
- 95 % Wilson intervals throughout.
- **Step 3** (separate, report only; `phase4_m_check.py`): the Phase 4 C2 UPO
  configuration with `upo_fixed_dimension` = 3 and 4 on the Phase 4 test seeds
  2000–2149 at 256 samples. The conditions are white_noise, ar1, sinusoid,
  two_tone, henon, henon@30dB, henon@20dB and logistic@20dB, shown beside the
  stored m = 2 results.

## 5. Runtime budget (measured on development seeds, 4 workers on 4 cores)

- Development run `results/dev/runtime.jsonl`: 240 windows (24 conditions ×
  seeds 0–9), each with m = 2, 3 and 4. It took **496 s wall time, 2.07 s per
  window**.
- Mean CPU seconds per window: m = 2 path 2.56, m = 3 2.65, m = 4 2.96, total
  8.24.
- Test run: 6 × 300 + 18 × 200 = **5,400 windows × 2.07 s ≈ 11,200 s ≈
  3.1 h**.
- Step 3: 8 conditions × 150 seeds × 2 m ≈ 2,400 UPO runs, about 0.5 h.
- Total ≈ 3.6 h, under the 8-hour limit. The margin covers resuming after
  container restarts; the runner skips finished windows.
- Seeds were not reduced for any condition. N1–N6 were given more seeds
  (300) because they carry the PASS rule.

## 6. What was seen before this preregistration (development seeds only)

- Generators and their realism check were committed in 5a7bbfa before any
  detector run on Phase 5 signals.
- The detector was then run on development seeds 0–9 of every condition for
  the runtime budget (`results/dev/tables/runtime_summary.md`). Nothing was
  changed in response. The detector is frozen, and the generators were
  already committed. The results, at m = 2, as counts over 10:

  | condition | AND | LLE alone | UPO alone |
  |---|---|---|---|
  | N1–N6 | 0 each | 0–2 | 0–1 |
  | S2 ectopic 2 / 5 / 10 % | 0 / 0 / 2 | 10 / 10 / 10 | 0 / 0 / 2 |
  | P1, P2, P4 | 10 each | 10 each | 10 each |
  | P3 Hénon / logistic | 2 / 5 | 10 / 10 | 2 / 5 |
  | G1 | 10 | 10 | 10 |
  | G2 | 0 | 0 | 0 |
  | G3 | 0 | 0 | 2 |

## 7. Amendments

(none)
