# Phase 6 preregistration: detrending candidates (Part C) and the Phase 6 test run

This file is committed and pushed **before any test window (seed ≥ 4000) is
generated**. After that, nothing below changes, and neither do:
- the methods (`methods.py`);
- the generators (`systems.py`, `config.py`, `experiments/phase5_rr/`);
- the rule (`decide.py`) or the tables (`analysis.py`);
- `final_pipeline.py`.

Any change goes in a dated amendment at the end of this file. `run_phase6.py
--phase test` reuses the Phase 3 guard: it refuses to start unless this file is
committed, unmodified and pushed, and lists every method below.

## 1. Methods (`methods.py`; all through `fp.analyze_segment` on the RR series)

- method: `baseline_k`
  - BASELINE-K = `combined_chaos_config()` with
    `keep_upo_on_short_lle_embedding = True`.
  - Decision: `combined_chaos_detected` = `lle_chaos_test` AND the Phase 4 C2
    instability-gated UPO test, m = 2.
  - Reference only; not eligible to win.
- method: `d1_linear_g05`
  - D1 = BASELINE-K + `rr_detrend = "linear"` and
    `rr_detrend_min_trend_sd = 0.5`.
  - The least-squares line is subtracted, and the window mean added back, only
    when |total change of the line| ≥ 0.5 × SD of the residual.
- method: `d2_linear_g07`
  - D2 = D1 with the gate at 0.7.
- method: `d3_smoothprior300_g07`
  - D3 = BASELINE-K + `rr_detrend = "smoothness_priors"`,
    `rr_detrend_lambda = 300`, `rr_detrend_min_trend_sd = 0.7`.
  - This is the Tarvainen et al. (2002) detrend on the beat index. Its cut-off
    is about 0.013 cycles per beat.

**Where the detrend is applied.** It is applied to `rr_dynamics` at the start of
`analyze_segment`, before the production LLE, the UPO analysis and
`lle_chaos_test`. Every surrogate is therefore generated from the detrended
series.

**Errors.** A window where `analyze_segment` raises counts as not detected. A
UPO analysis ending in a failure status counts as UPO not detected
(`combined_chaos_detected`, fixed in 7b4f1c2).

**How the candidates were chosen.** The choice used development seeds only:
0–29 for P*/G1 and 0–49 for N1–N6/S2 (`results/dev/tables/tuning.md`).
- The ungated detrends were screened:
  - linear;
  - moving median with windows 31 / 61 / 101;
  - smoothness priors with λ 10 / 50 / 300.
- All of them restore trend power (P3 pooled: baseline 20/60; linear 58/60).
  All of them also cost power on untrended chaos. Worst was Lorenz maxima:
  29/30 → 21/30 with linear, and 4–11/30 with moving median or smoothness
  priors.
- The loss is caused by subtracting a random fitted trend from a stationary
  window. This smears the narrow So-transform peak (C1,
  `results/dev/diagnose_trend.md`).
- The **trend gate** removes the loss: on development windows the fitted-trend
  / residual-SD ratio is at most 0.7 on untrended chaos, and mostly between
  0.7 and 2.2 on P3.
- The gated variants were evaluated exactly from the stored decisions:

  | variant | pooled P3 | largest loss vs BASELINE-K on P1/P2/P4/G1 | development null and S2 AND |
  |---|---|---|---|
  | linear, gate 0.5 | 58/60 | 1/30 | ≤ 2/50 |
  | linear, gate 0.7 | 57/60 | 0 | ≤ 2/50 |
  | smoothness priors 300, gate 0.7 | 54/60 | 0 | ≤ 2/50 |
  | linear, gate 0.3 | 58/60 | 5/30 on G1 | – |
  | moving median 101, gate 0.7 | 50/60 | 0 | – |

  Linear gate 0.3 and moving median 101 gate 0.7 were rejected.
- There are three candidates, covering two gates for the linear detrend and
  one smoothness-priors detrend.

## 2. Data

- **Conditions.** All 37 conditions (`config.CONDITIONS`):
  - the 24 Phase 5 conditions (N1–N6, S1, S2 2/5/10 %, S3 × 3, P1, P2,
    P3 × 2, P4 × 4, G1–G3), with generators unchanged;
  - the 13 Part B conditions: E1 bigeminy, E2 trigeminy, E3 couplets 5/10 %,
    E4 runs, E5 atrial 5/10 %, E6 ectopic 10 % + trend, and the
    annotation-edited S2 5 %, S2 10 %, E1, E2 and E3 10 % (`README.md`).
- **Windows.** 256 RR intervals, with 1/360 s quantization (edited values not
  re-quantized).
- **Seeds.** Phase 2E scheme and Phase 5 paired design. **Test seeds**:
  - **N1–N6: 4000–4399 (400 each)**;
  - **every other condition: 4000–4299 (300 each)**;
  - 13,200 windows in all, each evaluated with all four methods.
- **Runner.** Methods with the same effective configuration on a window share
  one `analyze_segment` call (`methods.effective_key`):
  - a gated detrend below its gate is exactly BASELINE-K;
  - D1 and D2 are identical when both gates pass.

  This was checked against full runs on 39 development windows with 0
  differences (`results/dev/shortcut_check.json`). `computed_as` in each
  record shows the sharing.

## 3. Part A and Part B (report only, no rule)

- **Part A** was completed before this file (`results/partA/partA.md`):
  - the 3 Phase 5 windows that raised no longer error;
  - 0 of 1,800 N1–N6 decisions changed.
- **Part B**: AND, each component and OR, with 95 % Wilson intervals, for
  every ectopy condition, raw and annotation-edited, for all four methods.
  The primary report is BASELINE-K.

## 4. PRIMARY rule for Part C (`decide.py`, committed in 7026a12, unchanged since)

A candidate **PASSES** only if, on the test windows:
1. its AND detector fires on at most floor(0.07 N) windows on EACH of N1–N6
   (N = 400: limit **28**) AND on EACH S2 rate, 2, 5 and 10 % (N = 300:
   limit **21**); and
2. its AND detection on EACH of P1, P2, P4 Hénon 30 / 20 dB, P4 logistic
   30 / 20 dB and G1 is at most floor(0.03 N) windows below BASELINE-K on the
   same windows (N = 300: at most **9** below; the task's 6/200 scaled).

**WINNER**: among passing candidates, the highest pooled AND detection count
on P3 (Hénon + trend, logistic + trend; out of 600).

**Tie-breaks**, in order:
- fewer pooled AND detections on the Part B ectopy conditions (all 13 E*
  conditions, raw and edited; out of 3,900);
- the lower candidate index (D1 < D2 < D3).

If no candidate passes, there is no winner and nothing is adopted.

**Adoption (C5).** The winner becomes an opt-in option, which already exists as
`rr_detrend*`, default `"none"`. `combined_chaos_config` gains an optional
argument that selects it, and the default call stays unchanged. Tests follow the
ground rules.

## 5. Runtime budget (development, 4 workers on 4 cores)

- `results/dev/runtime.jsonl`: all 37 conditions × seeds 0–4 × 4 methods, with
  sharing. It took **201 s wall time for 185 windows**, 4.32 CPU-s per window
  on average.
- Trend conditions cost about 9.5 CPU-s per window (two extra analyses);
  others about 3.
- Test estimate: 36.4 CPU-s summed over the six N conditions × 400, plus
  123.4 CPU-s summed over the other 31 conditions × 300. That is
  ≈ 51,600 CPU-s ≈ 12,900 s wall ≈ **3.6 h**, under the 8 h limit.
- Seeds exceed the minimums (300 / 200). The margin allows resuming after
  container restarts; the runner skips finished windows.

## 6. What was seen before this preregistration (development seeds only)

- The C1 diagnosis (`results/dev/diagnose_trend.md`,
  `plots/C1_trend_mechanism.png`).
- The C2 tuning table (`results/dev/tables/tuning.md`).
- The shortcut check.
- The runtime run: every method on every condition at development seeds 0–4.
  Its detection counts were not used to choose anything; the candidates were
  fixed from the tuning table before it ran.
- **Discovered during development:** `combined_chaos_detected` raised on a UPO
  failure status (edited bigeminy windows are constant). This was fixed to
  count them as not detected (7b4f1c2).

## 7. Amendments

- **2026-09-30, before any test window was generated.** Arithmetic
  correction: 6 × 400 + 31 × 300 = **11,700** windows, not 13,200 as written
  in Section 2. The seeds, conditions, methods, rule and budget calculation
  are unchanged. The Section 5 estimate already used these counts.
