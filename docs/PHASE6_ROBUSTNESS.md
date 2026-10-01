# Phase 6: Robustness of the Combined Chaos Detector

Research-development record in the style of the Phase 5 report. Every statement
is tied to a results file, table, preregistration or test under
`experiments/phase6_robust/`. Paths are relative to that directory unless they
start with `docs/`, `tests/` or `final_pipeline.py`. Synthetic RR-like data
only. **No MIT-BIH record was loaded or analyzed.**

## 1. Provenance and design

| Item | Value |
|---|---|
| Branch | `claude/amazing-allen-p4ssx3`, restarted from `main` at `ea7e068` (Phases 2E–5) |
| Environment | Python 3.13.12, numpy 2.1.3, BLAS 1 thread, 4 workers; manifests in `results/{dev,test}/manifest.jsonl` |
| BASELINE-K | `combined_chaos_config()` + `keep_upo_on_short_lle_embedding = True`, through `fp.analyze_segment` on the RR series; decision `combined_chaos_detected` = `lle_chaos_test` AND Phase 4 C2 instability-gated UPO (m = 2) |
| Signals | Phase 5 generators unchanged: 256 RR intervals, 1/360 s quantization, Phase 2E seed scheme, paired design. The 13 Part B conditions (Section 3) are built on the N1 base |
| Development seeds | as Phase 5: 0–99 for N*, S*, E*; 0–29 for P*, G* |
| Test seeds | **N1–N6: 4000–4399 (400 each); all other 31 conditions: 4000–4299 (300 each)**. 11,700 windows, each evaluated with all four methods |
| Seed discipline | `PREREGISTRATION.md` committed and pushed in `3325e16` before any test window. A window-count amendment followed in `97a302b`, also before any test window. `decide.py` was committed in `7026a12`. `run_phase6.py --phase test` reuses the Phase 3 guard |

## 2. Part A: crash-fix check (`partA_crash_check.py`, `results/partA/partA.md`)

The check reruns BASELINE-K on Phase 5 test windows.

| Window | Phase 5 | BASELINE-K |
|---|---|---|
| P1 Hénon seed 3008 | raised `ValueError` (AND counted 0) | no error; LLE embedding too short flagged; LLE yes, UPO yes, **AND yes** |
| P1 Hénon seed 3154 | raised | no error; **AND yes** |
| P4 Hénon 30 dB seed 3014 | raised | no error; **AND yes** |

- **The flag affects only windows that previously raised.** On N1–N6 with
  Phase 5 test seeds 3000–3299, 1,800 windows, **0 of 1,800 windows changed
  any decision** (LLE, UPO or AND). The AND counts are identical: 1, 1, 0, 0,
  1, 0.
- In the Phase 6 test run the flag kept 7 more windows (6 P1, 1 P4 Hénon
  20 dB) that would have raised. **0 of 46,800 method evaluations raised.**

## 3. Part B: heavier ectopy (evaluation only; `README.md` for definitions and sources)

**Premature beats.** Every premature interval is c × the underlying sinus RR,
with c ~ U(0.6, 0.8) (as Phase 5 S2).

**Ventricular patterns** have a full compensatory pause: the sinus node is not
reset, and the interval count and total duration are preserved.
- **E1 bigeminy:** every other beat premature, 50 %.
- **E2 trigeminy:** every third beat premature, 33 %.
- **E3 couplets:** pairs of premature beats; 6 or 13 couplets per window
  (4.7 % / 10.2 % of beats).
- **E4 runs:** 1–3 runs per window, each 3–6 premature beats, with an absolute
  run cycle U(0.40, 0.55) s. Consecutive ventricular beats at more than
  100 bpm are ventricular tachycardia, per Al-Khatib et al. (2017).

**Other patterns:**
- **E5 atrial-type:** no compensatory pause. The sinus node is reset, and the
  post-ectopic interval is d·RR with d ~ U(1.0, 1.1). Rates 5 % and 10 %.
- **E6:** S2 10 % plus the N3 trend, reusing exactly the Phase 5 draws.

**Annotation-edited (NN) variants** of S2 5 %, S2 10 %, E1, E2 and E3 10 %:
- every interval touching an ectopic beat is replaced by linear interpolation
  (by beat index) between the neighbouring NN intervals;
- edges take the nearest NN value;
- editing is applied after quantization, and interpolated values are not
  re-quantized.

The fraction of genuine NN intervals is fixed by the pattern:

| edited condition | NN fraction |
|---|---|
| S2 5 % | 0.90 |
| S2 10 % | 0.80 |
| E3 10 % | 0.85 |
| E2 | 0.34 |
| E1 | 0.008 |

### 3.1 BASELINE-K on every ectopy condition (test seeds; detections / windows, 95 % Wilson interval)

| condition | AND | LLE alone | UPO alone | OR |
|---|---|---|---|---|
| *N1 reference (no ectopy)* | 0/400 (0–1 %) | 36/400 (7–12 %) | 0/400 (0–1 %) | 36/400 |
| S2 isolated 2 % | 1/300 (0–2 %) | **289/300 (94–98 %)** | 1/300 | 289/300 |
| S2 isolated 5 % | 5/300 (1–4 %) | **290/300 (94–98 %)** | 5/300 | 290/300 |
| S2 isolated 10 % | 7/300 (1–5 %) | **293/300 (95–99 %)** | 7/300 | 293/300 |
| E1 bigeminy | 0/300 (0–1 %) | 37/300 (9–17 %) | 0/300 | 37/300 |
| E2 trigeminy | 0/300 (0–1 %) | 103/300 (29–40 %) | 0/300 | 103/300 |
| E3 couplets 5 % | 1/300 (0–2 %) | **292/300 (95–99 %)** | 1/300 | 292/300 |
| E3 couplets 10 % | 3/300 (0–3 %) | **298/300 (98–100 %)** | 3/300 | 298/300 |
| E4 runs | 0/300 (0–1 %) | **275/300 (88–94 %)** | 0/300 | 275/300 |
| E5 atrial 5 % | 4/300 (1–3 %) | **269/300 (86–93 %)** | 4/300 | 269/300 |
| E5 atrial 10 % | **11/300 (2–6 %)** | **254/300 (80–88 %)** | 12/300 | 255/300 |
| E6 isolated 10 % + trend | **10/300 (2–6 %)** | **299/300 (98–100 %)** | 10/300 | 299/300 |

### 3.2 Raw vs annotation-edited (BASELINE-K; same seeds and ectopic draws)

| pattern | NN fraction | AND raw → edited | LLE alone raw → edited | UPO alone raw → edited |
|---|---|---|---|---|
| S2 isolated 5 % | 0.90 | 5 → **0** /300 | 290 → **23** /300 | 5 → 4 |
| S2 isolated 10 % | 0.80 | 7 → **2** /300 | 293 → **46** /300 | 7 → 4 |
| E3 couplets 10 % | 0.85 | 3 → **1** /300 | 298 → **27** /300 | 3 → **24** |
| E2 trigeminy | 0.34 | 0 → **13** /300 (4.3 %; 2.5–7.2 %) | 103 → **291** /300 | 0 → 13 |
| E1 bigeminy | 0.008 | 0 → 0 /300 | 37 → 0 /300 | 0 → 0 (analysis fails, see below) |

**Edited bigeminy is not analyzable.** Only the two edge intervals are NN, so
the edited window is a constant or a straight line. The UPO status is
`constant_data` in 144/300 windows and `tube_empty` or `too_few_in_tube` in
125/300. `lle_chaos_test` is undefined in every window. It counts as "not
detected", which is not evidence of anything.

The four methods gave identical ectopy AND counts, except D1 on E6 (9 vs 10).
The gated detrends rarely act on ectopy windows.

### 3.3 What Part B shows

1. **The LLE component alone cannot be interpreted on raw RR with ectopy.**
   - It fires in 80–100 % of windows with isolated beats, couplets, runs or
     atrial-type ectopy, even at 2 %.
   - The strictly regular patterns (bigeminy 12 %, trigeminy 34 %) trigger it
     less. There, the premature beats form a periodic, largely linear
     structure that the IAAFT surrogates reproduce.
2. **The AND detector stays at or below 3.7 % on every raw ectopy pattern**
   (worst: atrial 10 %, 11/300; isolated 10 % with trend, 10/300). It
   increases with the ectopic rate. The UPO gate supplies the specificity.
3. **Annotation editing works when most intervals are NN** (NN fraction
   ≥ 0.80: S2 5 %, S2 10 %, couplets 10 %):
   - the AND false-positive rate falls to 0–0.7 %;
   - LLE-alone firing falls from 97–99 % to 8–15 %, close to its 9 % rate on
     clean N1.

   Couplet editing raises UPO-alone firing (3 → 24 of 300), because
   interpolated stretches create locally smooth segments, but the AND stays at
   1/300.
4. **Editing fails for dense ectopy.**
   - Trigeminy (NN 0.34): editing *raises* AND false positives from 0 to
     13/300 (4.3 %). Two-thirds of the window becomes piecewise-linear
     interpolation, which both components mistake for structure.
   - Bigeminy (NN 0.008): the edited window is degenerate.

## 4. Part C: detrending for the drift power loss

### 4.1 C1: mechanism (development seeds 0–29; `diagnose_trend.py`, `results/dev/diagnose_trend.md`, `plots/C1_trend_mechanism.png`)

Paired windows (same map realization and RR scale) were run through the C2 UPO
analysis. The So-transform mode was located separately for transformed points
from the first, middle and last third of the window.

| system | variant | Level-B windows | UPO (gate) windows | median So-histogram peak height | median mode drift, first → last third | drift / trend change |
|---|---|---|---|---|---|---|
| Hénon | clean P1 | 30/30 | 30/30 | 0.830 | 0.000 s | 0.00 |
| Hénon | + trend (P3) | 8/30 | 8/30 | 0.181 | 0.045 s | **1.11** |
| Hénon | P3, linear detrend | 30/30 | 30/30 | 0.655 | 0.002 s | 0.04 |
| logistic | clean P2 | 30/30 | 30/30 | 0.585 | 0.002 s | 0.05 |
| logistic | + trend (P3) | 12/30 | 12/30 | 0.157 | 0.045 s | **1.08** |
| logistic | P3, linear detrend | 28/30 | 28/30 | 0.394 | 0.006 s | 0.11 |

**The mechanism is confirmed:**
- The trend moves the fixed point across the window. The So-transform mode
  drifts by about the trend's own change (ratio 1.1).
- The single sharp peak splits into one peak per third (see the plot), and the
  peak height falls to about a fifth.
- Level B then fails, and the UPO component misses. The LLE component is
  unaffected (Phase 5: 200/200).
- Removing the trend restores the peak.

### 4.2 C2: candidates (development seeds only; `results/dev/tables/tuning.md`)

**Implementation.** `PipelineConfig.rr_detrend` (default `"none"`) is applied
to `rr_dynamics` at the start of `analyze_segment`, so all surrogates come from
the detrended series; the window mean is added back. The options are:
- `"linear"`;
- `"moving_median"` (`rr_detrend_window`);
- `"smoothness_priors"` (Tarvainen et al. 2002, `rr_detrend_lambda`).

`rr_detrend_min_trend_sd` (default 0) applies the detrend only when the
least-squares trend's total change is at least that multiple of the residual
SD (`fp.linear_trend_ratio`).

**Screening** (P*/G1 seeds 0–29, N1–N6/S2 seeds 0–49):
- Every ungated detrend restores trend power: pooled P3 58/60 linear, 56/60
  smoothness priors λ = 300, 51/60 moving median 101, against 20/60 for
  BASELINE-K.
- Every ungated detrend also costs power on untrended chaos. G1 Lorenz maxima
  fell from 29 to 21/30 with linear and to 4–11/30 with the nonlinear filters.
- A detrend fits and subtracts a random trend from a stationary window, which
  smears the narrow So peak exactly as a real trend does.
- The fitted-trend / residual-SD ratio is at most about 0.7 on untrended chaos
  and mostly 0.7–2.2 on P3, so a gate on it removes the loss.
- Gated variants were evaluated exactly from the stored decisions (below the
  gate the window is analysed as BASELINE-K). Three candidates were chosen:

| id | method | pooled P3 (dev) | largest loss vs BASELINE-K (dev) |
|---|---|---|---|
| D1 | `d1_linear_g05`: linear, gate 0.5 | 58/60 | 1/30 |
| D2 | `d2_linear_g07`: linear, gate 0.7 | 57/60 | 0 |
| D3 | `d3_smoothprior300_g07`: smoothness priors λ = 300, gate 0.7 | 54/60 | 0 |

### 4.3 C3: preregistered rule (`PREREGISTRATION.md`, `decide.py`)

A candidate **PASSES** only if, on test windows, both hold:
1. Its AND fires on at most floor(0.07 N) windows on EACH of N1–N6 (≤ 28/400)
   and each S2 rate (≤ 21/300).
2. Its AND on EACH of P1, P2, P4 Hénon 30/20 dB, P4 logistic 30/20 dB and G1
   is at most floor(0.03 N) = 9/300 below BASELINE-K on the same windows.

**WINNER:** the highest pooled AND on P3.
**Tie-breaks:** fewer pooled AND on the 13 Part B conditions, then the
candidate index.

**Budget.** The development estimate was 3.6 h.
- Wall time was 3.5 h across three segments.
- A container restart at 5,990 windows was resumed with the same command.
- A JSON serialization crash at 10,972 windows was caused by the diagnostic
  `trend_ratio` being infinite on constant edited-bigeminy windows. It was
  fixed by recording null (one line), documented as amendment 2, and the run
  was resumed. No decision was affected.

### 4.4 C4: test results (`results/test/tables/test_primary.md`, `test_summary.md`)

**Specificity (AND detections):**

| method | N1 | N2 | N3 | N4 | N5 | N6 | S2 2 % | S2 5 % | S2 10 % |
|---|---|---|---|---|---|---|---|---|---|
| BASELINE-K | 0 | 0 | 1 | 1 | 2 | 0 | 1 | 5 | 7 |
| D1 | 0 | 0 | 0 | 0 | 2 | 0 | 1 | 5 | 7 |
| D2 | 0 | 0 | 0 | 0 | 2 | 0 | 1 | 5 | 7 |
| D3 | 0 | 0 | 0 | 0 | 2 | 0 | 1 | 5 | 7 |
| limit | 28/400 | 28 | 28 | 28 | 28 | 28 | 21/300 | 21 | 21 |

**Power on untrended positives (AND, out of 300; at most 9 below BASELINE-K):**

| method | P1 | P2 | P4 H 30 | P4 H 20 | P4 L 30 | P4 L 20 | G1 |
|---|---|---|---|---|---|---|---|
| BASELINE-K | 300 | 294 | 300 | 296 | 299 | 284 | 294 |
| D1 | 300 | 293 | 300 | 296 | 298 | 282 | **279 (FAIL, −15)** |
| D2 | 300 | 294 | 300 | 296 | 299 | 284 | 291 (−3) |
| D3 | 300 | 294 | 300 | 296 | 299 | 284 | 292 (−2) |

**Primary summary:**

| method | PASS | pooled P3 AND (Hénon + logistic, /600) | pooled Part B ectopy AND (/3,900) |
|---|---|---|---|
| BASELINE-K | (not eligible) | 195 (74 + 121) | 45 |
| D1 linear, gate 0.5 | **no** (G1) | 581 (300 + 281) | 44 |
| **D2 linear, gate 0.7** | **yes** | **571 (298 + 273)** | 45 |
| D3 smoothness priors 300, gate 0.7 | yes | 539 (295 + 244) | 45 |

**WINNER: D2** (`d2_linear_g07`). Pooled P3 AND power rises from 32.5 % to
95.2 %, with no false-positive increase on any null or ectopy condition.
- **Other conditions (D2 vs BASELINE-K), all out of 300:**
  - S1: AND 0 for both;
  - S3: 0/0/0 vs 0/0/2 (N1/N3/N4 unquantized);
  - G2 and G3: 0 for both.
- **Side effect.** Detrending makes the LLE component alone fire more often on
  the trend and step nulls (N3 23 → 40/400, N4 17 → 40/400). The AND is
  unaffected (0/400).
- **D1 failed** because its lower gate detrends more stationary Lorenz-maxima
  windows. This is the development finding, reproduced on the test seeds.

### 4.5 C5: adoption

- `final_pipeline.py` gains `PHASE6_DETREND = {"rr_detrend": "linear",
  "rr_detrend_min_trend_sd": 0.7}` and
  `combined_chaos_config(config=None, detrend=False)`.
  - With `detrend=True` it adds D2.
  - The default call is unchanged; this is tested.
  - The exact Phase 6 tested configuration is
    `combined_chaos_config(replace(CFG, keep_upo_on_short_lle_embedding=True), detrend=True)`.
- **Other `final_pipeline.py` changes in Phase 6**, all opt-in, with defaults
  unchanged:
  - `rr_detrend`, `rr_detrend_window`, `rr_detrend_lambda`,
    `rr_detrend_min_trend_sd`, `detrend_rr`, `linear_trend_ratio`;
  - **a fix to `combined_chaos_detected`.** A UPO analysis that ends in a
    failure status (for example `constant_data`) has no gate fields. It now
    counts as not detected instead of raising `ValueError`. It was found on
    edited-bigeminy development windows; Phase 5 had no such window.
- **Checks:**
  - `tests/test_phase6_detrend.py` (8 tests);
  - `tests/test_phase6_combined_detrend.py` (7 tests);
  - `tests/test_phase5_combined_chaos.py` (+1 test for the fix);
  - `adoption_check.py`: 111 stored test windows (seeds 4000–4002 × 37
    conditions) recomputed through the helper, **0 differences**;
  - ground rules (`experiments/phase3_lle/groundrule_check.sh`): pytest 400
    passed plus the pre-existing `test_hybrid_period1_stability...[1-prl_norm]`
    failure (401/401 with AVX-512 disabled); same-machine replicability 44/44;
    `run_phase2e --replicate` 13/44 with an identical file.

## 5. Recommended MIT-BIH analysis configuration

Based only on the evidence from Phases 3–6:

**Detector.** Use the AND decision:
`combined_chaos_detected(analyze_segment(rr, combined_chaos_config(replace(CFG, keep_upo_on_short_lle_embedding=True), detrend=True)))`.
- **Why AND and not a component:**
  - AND specificity held on every linear and nonlinear-stochastic null (Phases
    5–6: at most 2/400);
  - the LLE component alone fires on 80–100 % of ectopy windows (Part B, and
    Phase 5 S2);
  - the UPO component alone fires on 8 % of edited-couplet windows (Part B).
- **Why `keep_upo_on_short_lle_embedding = True`:** otherwise strongly chaotic
  windows are lost to an error (Part A; 7 more windows in the Phase 6 test
  run) without any change elsewhere.
- **Why `detrend=True` (D2):** it recovers power on drifting windows (pooled
  P3 from 32.5 % to 95.2 %). It costs at most 3/300 on untrended chaos and
  adds no false positives (Section 4.4).

**Input: annotation-edited NN intervals, not raw RR**, for windows whose
fraction of genuine NN intervals is **at least 0.80**. MIT-BIH has beat
annotations, so each ectopic interval and its following interval can be
replaced by linear interpolation as in Part B.
- On the tested patterns with NN fraction 0.80–0.90 (isolated 5 % and 10 %,
  couplets 10 %), editing lowered the AND false-positive rate from 1.7–2.3 %
  (raw) to 0–0.7 %.
- It returned the LLE component to near its clean-null rate (8–15 % vs 9 %),
  so the components remain interpretable.
- Raw RR gave the AND 1–3.7 % false positives purely from ectopy.
- **Windows with an NN fraction below 0.80 should be excluded** from chaos
  inference, or reported separately:
  - edited trigeminy (0.34) raised the AND to 4.3 %;
  - edited bigeminy (0.008) is degenerate.

  0.80 is the lowest NN fraction at which editing was shown to behave. Values
  between 0.34 and 0.80 were not tested.

**Detrending applies to the edited series.** The detrend acts on whatever
series is passed, so detrending an edited window is the same code path. On the
edited conditions the D2 AND counts were identical to BASELINE-K.

**Question being asked.** Editing changes the scientific question. The edited
series describes the sinus-rhythm dynamics between ectopic beats, not the
arrhythmic beat sequence itself. If the aim is to ask whether the *arrhythmic
sequence* is chaotic, raw RR is the object. In that case:
- AND detections must be compared against the ectopy-only false-positive rates
  in Section 3.1, up to 3.7 % (upper 95 % bound 6.4 %) for the patterns
  tested;
- they cannot be compared against zero.

## 6. Limitations

- **Power under editing was not measured.** No condition combined a chaotic
  signal with ectopy followed by editing. Editing replaces 10–20 % of intervals
  by interpolation and may reduce AND power on genuinely chaotic NN dynamics.
  This should be tested before relying on edited-NN negatives.
- **NN fractions between 0.34 and 0.80 were not tested.** The 0.80 threshold is
  the lowest tested value that behaved, not an estimated boundary.
- **The ectopy models are idealized:**
  - fixed coupling range;
  - exact full compensatory pauses for ventricular beats;
  - a simple reset model for atrial beats;
  - no interpolated PVCs, aberrancy, atrial fibrillation, or missed or extra
    beat detections.

  Real MIT-BIH annotations include non-ectopic abnormal labels (for example
  bundle-branch-block beats) that the NN edit would also remove.
- **The gate threshold is a single number,** chosen on development data from
  these generators. The ratio distribution of real RR windows is unknown.
  Windows with genuine slow non-linear drift (curved trends, gradual steps) are
  only partly handled by a linear detrend. Moving-median and smoothness-priors
  alternatives cost more power in this study.
- **Detrending side effect.** The LLE component alone fires more often on
  detrended trend and step nulls (up to 10 %), although the AND is unaffected.
- **Sampled flows** (Rössler, Mackey–Glass) remain undetected by every method,
  as in Phase 5.
- **Test-run interruptions.** A container restart and a serialization crash
  (amendment 2) interrupted the run. The 11,700 records are unique and complete,
  and each depends only on its seed. The shared-computation shortcut was
  verified to give identical decisions (39/39 development windows); the
  adopted helper reproduced 111/111 stored decisions.
- **Web sources.** Ectopy definitions were checked against web-search results,
  since publisher pages are blocked by the environment's proxy.
  physionet.org is still not reachable.
- **One machine.**

## 7. Reproducibility

```bash
python -m experiments.phase6_robust.partA_crash_check                       # Part A
python -m experiments.phase6_robust.diagnose_trend                          # C1 (dev)
python -m experiments.phase6_robust.run_phase6 --phase dev --tag tuning \
    --methods baseline_k,linear,movmed_31,movmed_61,movmed_101,sp_10,sp_50,sp_300 \
    --conditions P3_henon_rr_trend,P3_logistic_rr_trend,P1_henon_rr,P2_logistic_rr,P4_henon_rr_30dB,P4_henon_rr_20dB,P4_logistic_rr_30dB,P4_logistic_rr_20dB,G1_lorenz_maxima
python -m experiments.phase6_robust.run_phase6 --phase dev --tag tuning \
    --methods baseline_k,linear,movmed_31,movmed_61,movmed_101,sp_10,sp_50,sp_300 --seeds 0-49 \
    --conditions N1_linear_rr,N2_power_law,N3_linear_rr_trend,N4_linear_rr_step,N5_linear_rr_warped,N6_noisy_rsa,S2_ectopic_2pct,S2_ectopic_5pct,S2_ectopic_10pct
python -m experiments.phase6_robust.tune_summary                            # C2
python -m experiments.phase6_robust.shortcut_check
python -m experiments.phase6_robust.run_phase6 --phase test                 # guarded
python -m experiments.phase6_robust.analysis --phase test --file test
python -m experiments.phase6_robust.decide   --phase test --file test
python -m experiments.phase6_robust.adoption_check
python -m pytest tests
```

Commits:
- `7026a12` scaffold, Part B generators, C1, rule;
- `f574f98` Part A and the `rr_detrend` option;
- `7b4f1c2` tuning, gate option, helper fix, candidates;
- `3325e16` and `97a302b` preregistration;
- `37b1391` test results;
- `a981bfc` adoption.

## References

The ectopy definitions (LITFL, AMBOSS, ECGpedia and the 2017 AHA/ACC/HRS
guideline) and the smoothness-priors detrend (Tarvainen et al. 2002) are
referenced in `experiments/phase6_robust/README.md`. The HRV norms are in
`experiments/phase5_rr/README.md`.
