# Phase 5: RR-Interval Stress Test of the Frozen Combined Detector

Research-development record in the style of the Phase 3 and Phase 4 reports.
Every statement is tied to a results file, table, preregistration or test under
`experiments/phase5_rr/` (paths are relative to it unless they start with
`docs/`, `tests/` or `final_pipeline.py`). Synthetic RR-like data only.
**No MIT-BIH record was loaded or analyzed.**

This phase is an **evaluation, not a tuning round**. The detector's parameters
were fixed before any Phase 5 signal existed and were not changed in response to
any result.

## 1. Provenance and design

| Item | Value |
|---|---|
| Branch | `claude/amazing-allen-p4ssx3`. When Phase 5 started, PR #2 (Phase 4) was still open, so the branch is main (`858d616`) + Phase 4 (`c561e85`) + Phase 5 |
| Environment | Python 3.13.12, numpy 2.1.3 (`requirements.txt`), BLAS 1 thread, 4 workers; manifests in `results/{dev,test}/manifest.jsonl` |
| Detector | `fp.analyze_segment(rr, fp.phase4_upo_config(replace(fp.CFG, lle_chaos_test=True)))` on raw RR (`use_corrected_rr_for_dynamics = False`) |
| Components | LLE = `out["lle_chaos_test"]["detected"]` (Phase 3 winner); UPO = `out["upo"]["instability_gate_detected"]` (Phase 4 winner C2) |
| PRIMARY | LLE **AND** UPO. Each component alone and OR are secondary |
| Conditions | 24: nulls N1–N6 (PASS rule); secondary nulls S1–S3; positives P1–P4; generalization G1–G3 (Section 3) |
| Windows | 256 RR intervals (seconds), quantized to 1/360 s unless marked unquantized |
| Development seeds | Phase 2E seed numbers: 0–99 for N* and S*, 0–29 for P* and G* |
| Test seeds | **3000–3299 for N1–N6 (300 each); 3000–3199 for every other condition (200 each)**: 5,400 windows |
| Seed discipline | Generators frozen in `5a7bbfa` before any detector run on Phase 5 signals. `decide.py` in `1affece`. `PREREGISTRATION.md` committed and pushed in `03787bf` before any test window. `run_phase5.py --phase test` reuses the Phase 3 guard |
| m-sensitivity | Every test window was also run at m = 3 and m = 4, with `lle_chaos_test_m` and `upo_fixed_dimension` both set to m (report only) |

## 2. Step 0: the `analyze_segment` path reproduces Phase 4

**Setup.** Before anything else, 5 windows of each Phase 4 condition were run
through `analyze_segment` with the configuration above: 14 conditions × test
seeds 2000–2004 at 256 samples, 70 windows (`equivalence_check.py`,
`results/equivalence_check.md`).

**Result.** **0 differences** from the stored Phase 4 results:
- C2 UPO decisions;
- `lle_chaos_test` decisions and p-values (bitwise).

No window raised an error.

The m = 3 and m = 4 arms call `fp.lle_chaos_test` and `fp.run_upo_analysis`
directly, exactly as `analyze_segment` does. At m = 2 this direct path matched
`analyze_segment` on 48/48 development windows (`results/dev/m_path_check.json`).

## 3. Generators and realism check (Step 1)

All parameters are in `config.py`, the generators in `systems.py`. The sources
are in `README.md`, which also gives each construction in full.

**RR scale and pairing.** Every window draws a mean RR μ ~ U(0.65, 0.95) s and
an SD σ ~ U(0.035, 0.060) s. Conditions that share a base process and seed
share the realization, a paired design. For example, N1, N3, N4, S2 and S3
differ only in the modification or in quantization.

| id | condition | construction |
|---|---|---|
| N1 | linear_rr | McSharry et al. (2003) bimodal spectrum (LF 0.1 Hz / HF 0.25 Hz, SD 0.01 Hz), LF/HF ~ U(1.5, 2.0), random phases |
| N2 | power_law | Gaussian 1/f^β, β ~ U(0.9, 1.1) |
| N3 | linear_rr_trend | N1 + centred linear trend, total change U(5, 10) % of μ, random sign |
| N4 | linear_rr_step | N1 + centred abrupt step, U(5, 10) % of μ, at a beat in [0.3N, 0.7N] |
| N5 | linear_rr_warped | N1 → exp(0.5 z) → rescaled (skewness ≈ 1.5) |
| N6 | noisy_rsa | respiratory sinusoid 0.2–0.3 Hz (50 % of the variance) + white noise |
| S1 | setar | SETAR(2;1,1), slopes 0.7 / −0.5 (stable fixed point), noise-driven |
| S2 | ectopic 2 / 5 / 10 % | N1 + 5 / 13 / 26 premature beats (c ~ U(0.6, 0.8)) with full compensatory pause |
| S3 | *_unq | N1, N3, N4 unquantized, identical draws to the quantized twins |
| P1 / P2 | henon_rr / logistic_rr | Hénon x (1.4, 0.3) / logistic r = 4, rescaled |
| P3 | + trend | P1 / P2 + the N3 trend |
| P4 | 30 / 20 dB | P1 / P2 + white measurement noise, before quantization |
| G1 | lorenz_maxima | successive maxima of Lorenz z (10, 28, 8/3) |
| G2 | rossler_flow | Rössler x (0.2, 0.2, 5.7) sampled every 1.0 (about 6 samples per orbit) |
| G3 | mackey_glass | Mackey–Glass τ = 17 sampled every 6.0 |

**Realism** (`results/realism/realism_dev.md` and `plots/examples_*.png`;
development seeds; no detector run). The targets are:
- mean RR 0.6–1.0 s (task specification);
- LF/HF 1.5–2.0 and the Task Force (1996) 5-min norms;
- short-term SDNN 50 ± 16 ms and RMSSD 42 ± 15 ms (Nunan et al. 2010).

Measures use the Task Force bands on 4 Hz resampled RR with Welch spectra.

| group | mean RR in 0.6–1.0 s (every window) | SDNN, median (ms) | RMSSD, median (ms) | LF/HF, median | verdict |
|---|---|---|---|---|---|
| N1, N3–N5 | yes | 47.6–55.6 | 37–39 | 1.64–1.70 | all targets met |
| N2 | yes | 45.8 | 35.8 | 1.20 | met except the narrow Task Force LF/HF band (within Nunan) |
| N6 | yes | 47.0 | 60.2 | 0.11 | mean and SD only; HF-dominated by design |
| S1 | yes | 49.6 | 54.9 | 0.85 | met except the Task Force LF/HF band |
| S2 2 / 5 / 10 % | yes | 68 / 91 / 119 | 92 / 138 / 191 | 0.99 / 0.57 / 0.36 | mean only; ectopics inflate SDNN and RMSSD as in unedited recordings |
| P1–P4, G1 | yes | 45–51 | 62–74 | 0.24–0.58 | mean and SD only |
| G2 / G3 | yes | 48.0 / 48.4 | 50.5 / 39.6 | 0.03 / 0.34 | mean, SD and RMSSD |

**Reading the verdict:**
- The linear nulls are realistic in mean, variability and LF/HF balance.
- The chaotic positives are realistic in mean, SD and quantization only. A map
  or flow keeps its own spectrum when rescaled; Hénon and logistic are
  anticorrelated beat to beat.
- **physionet.org is not reachable** from this environment (proxy
  `CONNECT tunnel failed, response 403`, a policy denial), so the optional
  comparison with NSRDB summary statistics was not done.

## 4. Preregistered rule and budget (`PREREGISTRATION.md`)

**PASS** iff the AND detector (m = 2, through `analyze_segment`) fires on at
most floor(0.07 N) windows on EACH of N1–N6 (quantized, 256 intervals). With
N = 300 the limit is **21/300**.
- A window where `analyze_segment` raises counts as not detected.
- Power (AND on P1–P4, G1–G3) is reported, not pass/fail.

**Budget.** Development seeds 0–9 × 24 conditions took 496 s for 240 windows,
i.e. 2.07 s per window, each with m = 2, 3 and 4.
- Test: 5,400 windows ≈ 3.1 h, plus Step 3 ≈ 0.5 h.
- Actual test run: **11,421 s (3.2 h)** in one uninterrupted pass.
- No condition's seeds were reduced. N1–N6 were given 300 because they carry
  the rule.

**Seen before preregistration.** Detector results on development seeds 0–9 at
m = 2, over 10 windows (PREREGISTRATION.md Section 6):
- AND 0 on all of N1–N6;
- LLE 10/10 on every S2 rate;
- P3 AND 2 and 5;
- G2 and G3 AND 0.

Nothing was changed in response.

## 5. Test-seed results

### 5.1 Primary decision (m = 2, `analyze_segment`; `results/test/tables/test_primary.md`)

| condition | AND | limit | within | LLE alone | UPO alone | OR | errors |
|---|---|---|---|---|---|---|---|
| N1 linear_rr | **1/300** (0.1–1.9 %) | 21 | yes | 20/300 (6.7 %) | 4/300 | 23/300 | 0 |
| N2 power_law | **1/300** | 21 | yes | 8/300 (2.7 %) | 4/300 | 11/300 | 0 |
| N3 linear_rr_trend | **0/300** (0–1.3 %) | 21 | yes | 24/300 (8.0 %) | 2/300 | 26/300 | 0 |
| N4 linear_rr_step | **0/300** | 21 | yes | 14/300 (4.7 %) | 3/300 | 17/300 | 0 |
| N5 linear_rr_warped | **1/300** | 21 | yes | **53/300 (17.7 %)** | 6/300 | 58/300 | 0 |
| N6 noisy_rsa | **0/300** | 21 | yes | 17/300 (5.7 %) | 1/300 | 18/300 | 0 |

**Result: PASS.** The AND detector fired on 3 of 1,800 null windows in all.
The per-condition 95 % Wilson upper bound is 1.9 % (1/300) or 1.3 % (0/300).
The PASS rule is on the AND. Components alone and OR are secondary.

### 5.2 Every condition at m = 2 (detections / windows; 95 % Wilson intervals in `results/test/tables/test_summary.md`)

| condition | AND | LLE alone | UPO alone | OR | errors |
|---|---|---|---|---|---|
| S1 setar | 1/200 | **26/200 (13 %)** | 2/200 | 27/200 | 0 |
| S2 ectopic 2 % | 1/200 | **194/200 (97 %)** | 1/200 | 194/200 | 0 |
| S2 ectopic 5 % | 2/200 | **195/200 (98 %)** | 2/200 | 195/200 | 0 |
| S2 ectopic 10 % | **7/200 (3.5 %; 1.7–7.0 %)** | **197/200 (99 %)** | 7/200 | 197/200 | 0 |
| S3 N1 unquantized | 0/200 | 15/200 | 2/200 | 17/200 | 0 |
| S3 N3 unquantized | 0/200 | 25/200 | 2/200 | 27/200 | 0 |
| S3 N4 unquantized | 0/200 | 8/200 | 3/200 | 11/200 | 0 |
| P1 henon_rr | 198/200 | 198/200 | 198/200 | 198/200 | **2** |
| P2 logistic_rr | 194/200 | 200/200 | 194/200 | 200/200 | 0 |
| P3 henon_rr + trend | **41/200 (20.5 %)** | 200/200 | 41/200 | 200/200 | 0 |
| P3 logistic_rr + trend | **86/200 (43 %)** | 200/200 | 86/200 | 200/200 | 0 |
| P4 henon_rr 30 dB | 199/200 | 199/200 | 199/200 | 199/200 | **1** |
| P4 henon_rr 20 dB | 198/200 | 200/200 | 198/200 | 200/200 | 0 |
| P4 logistic_rr 30 dB | 200/200 | 200/200 | 200/200 | 200/200 | 0 |
| P4 logistic_rr 20 dB | 188/200 | 200/200 | 188/200 | 200/200 | 0 |
| G1 lorenz_maxima | 199/200 | 200/200 | 199/200 | 200/200 | 0 |
| G2 rossler_flow | **0/200** | 0/200 | 0/200 | 0/200 | 0 |
| G3 mackey_glass | **0/200** | 0/200 | 62/200 | 62/200 | 0 |

### 5.3 Why components fire or miss (descriptive; `results/test/tables/component_breakdown.md`)

- **UPO misses under a trend** (P3) are Level-B failures, not gate failures:
  - Hénon + trend: 44/200 windows have a Level-B peak and 41 pass the gate;
  - logistic + trend: 86 and 86.

  A trend of 5–10 % of μ is 40–90 ms across the window, comparable to σ. It
  moves the fixed point across the window and smears the So-transform
  histogram peak.
- **The LLE test cannot detect the sampled flows** G2 and G3. Its p-value is
  1.00 in every window: the divergence slope over steps 1–5 of the smooth
  signal is below every IAAFT surrogate's, and the test is one-sided.
  - On G2, UPO finds a Level-B peak in 200/200 windows, but none passes the
    instability gate.
  - On G3, UPO finds one in 106 windows and passes the gate in 62.
- **Ectopic windows.** The LLE p-value is 0.01 (the minimum with 99
  surrogates) in 97–99 % of S2 windows at every rate. UPO fires in 1, 2 and 7
  of 200 windows, rising with the rate.
- **Quantization** (S3 against its quantized twin, the same 200 seeds): LLE
  fires 14 vs 15 (N1), 18 vs 25 (N3) and 8 vs 8 (N4) times; UPO 3 vs 2, 1 vs 2
  and 2 vs 3. There is no systematic effect.
- **Errors.** 3 of 5,400 windows raised `ValueError: Insufficient embedded
  points after reconstruction` in `analyze_segment`: P1 seeds 3008 and 3154,
  and P4 Hénon 30 dB seed 3014.
  - The production LLE embedding (TDMI τ × Cao m; median τ 10–11, m 10 on
    Hénon) left fewer than 50 points, and `keep_upo_on_short_lle_embedding` is
    False in the tested configuration.
  - Both components would have detected all three windows. The preregistered
    AND counts them as not detected.
  - No null window raised.

## 6. m-sensitivity of both components (report only)

Counts at m = 2 / 3 / 4, with both `lle_chaos_test_m` and
`upo_fixed_dimension` set to m (`results/test/tables/test_summary.md`).

| condition | N | AND | LLE alone | UPO alone |
|---|---|---|---|---|
| N1 linear_rr | 300 | 1 / 0 / 0 | 20 / 16 / 6 | 4 / 3 / 0 |
| N2 power_law | 300 | 1 / 0 / 0 | 8 / 10 / 11 | 4 / 2 / 2 |
| N3 trend | 300 | 0 / 1 / 0 | 24 / 40 / 24 | 2 / 6 / 0 |
| N4 step | 300 | 0 / 0 / 0 | 14 / 6 / 3 | 3 / 7 / 0 |
| N5 warped | 300 | 1 / 1 / 0 | 53 / 13 / 2 | 6 / 13 / 1 |
| N6 noisy_rsa | 300 | 0 / 0 / 0 | 17 / 7 / 7 | 1 / 0 / 0 |
| S1 setar | 200 | 1 / 0 / 0 | 26 / 23 / 38 | 2 / 0 / 2 |
| S2 ectopic 2 % | 200 | 1 / 3 / 0 | 194 / 158 / 119 | 1 / 3 / 0 |
| S2 ectopic 5 % | 200 | 2 / 2 / 2 | 195 / 172 / 157 | 2 / 2 / 2 |
| S2 ectopic 10 % | 200 | 7 / 9 / 12 | 197 / 194 / 197 | 7 / 9 / 12 |
| S3 N1 / N3 / N4 unq | 200 | 0 / 0 / 0; 0 / 2 / 0; 0 / 0 / 0 | 15 / 13 / 4; 25 / 32 / 10; 8 / 7 / 2 | 2 / 4 / 0; 2 / 3 / 0; 3 / 4 / 0 |
| P1 henon_rr | 200 | 198 / 197 / 187 | 198 / 197 / 197 | 198 / 200 / 190 |
| P2 logistic_rr | 200 | 194 / 194 / 185 | 200 / 200 / 200 | 194 / 194 / 185 |
| P3 henon + trend | 200 | 41 / 6 / 1 | 200 / 199 / 200 | 41 / 6 / 1 |
| P3 logistic + trend | 200 | 86 / 42 / 23 | 200 / 200 / 200 | 86 / 42 / 23 |
| P4 henon 30 / 20 dB | 200 | 199 / 198 / 180; 198 / 195 / 169 | 199–200 at every m | 199 / 198 / 180; 198 / 195 / 169 |
| P4 logistic 30 / 20 dB | 200 | 200 / 200 / 190; 188 / 194 / 159 | 200 at every m | 200 / 200 / 190; 188 / 194 / 159 |
| G1 lorenz_maxima | 200 | 199 / 147 / 101 | 200 / 200 / 200 | 199 / 147 / 101 |
| G2 rossler_flow | 200 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 |
| G3 mackey_glass | 200 | 0 / 0 / 0 | 0 / 0 / 0 | 62 / 45 / 184 |

**What changes with m:**
- **AND specificity does not depend on m.** It stays at most 1/300 on every
  N1–N6 condition at m = 3 and 4.
- **UPO power falls at m = 4:**
  - P4 logistic 20 dB drops from 188 to 159;
  - Lorenz maxima drops from 199 through 147 to 101;
  - trend power collapses (41 → 6 → 1).
- **The LLE component's false-positive rates move with m in both directions:**
  - N5 warped falls 53 → 13 → 2;
  - N1 falls 20 → 16 → 6;
  - N3 peaks at 40 at m = 3;
  - S1 rises to 38 at m = 4.
- **The flows are not rescued by m.** The LLE test stays at 0/200 on Rössler
  and Mackey–Glass at every m. UPO on Mackey–Glass rises to 184/200 at m = 4,
  but the AND stays at 0 because the LLE component never fires. The fixed
  m = 2 does not generalize to these sampled flows, and neither does m = 3 or 4
  for the LLE test.
- **Ectopic false positives grow with m.** AND on 10 % ectopics rises
  7 → 9 → 12.

## 7. Phase 4 m-check (Step 3; report only)

`phase4_m_check.py`, `results/phase4_m_check.md`. The Phase 4 C2 UPO
detector (`fp.phase4_upo_config()`: M = 15, 50 AAFT surrogates, median gate at
δ = 0.2) was rerun with only `upo_fixed_dimension` changed. It used the Phase 4
test seeds 2000–2149 at 256 samples. The m = 2 column is the stored Phase 4
result.

| condition | m = 2 (stored) | m = 3 | m = 4 |
|---|---|---|---|
| white_noise | 0/150 | 0/150 | 0/150 |
| ar1 | 0/150 | 0/150 | 0/150 |
| sinusoid | 0/150 | 0/150 | 0/150 |
| two_tone | 1/150 (0.1–3.7 %) | **7/150 (2.3–9.3 %)** | 0/150 |
| henon | 150/150 | 150/150 | 150/150 |
| henon@30dB | 150/150 | 150/150 | 148/150 |
| henon@20dB | 143/150 | 149/150 | 140/150 |
| logistic@20dB | 142/150 | 148/150 | **125/150** |

**Findings:**
- **m = 3.** Noisy power is slightly higher: Hénon 20 dB goes from 143 to 149
  and logistic 20 dB from 142 to 148. Two-tone false positives rise to 7/150
  (4.7 %). That is within the Phase 4 7 % limit (10/150), but seven times the
  m = 2 count.
- **m = 4.** Specificity is perfect on these controls, but noisy logistic power
  falls to 125/150.
- **m = 2 is not dominated at 256 samples.** Phase 4 fixed m = 2 on
  development data without testing m = 3 and 4 on test seeds, and this check
  fills that gap. m = 3 trades a small power gain for a real two-tone
  specificity cost. Statuses: 25 (m = 3) and 20 (m = 4) of 1,200 windows were
  `no_peaks`; the rest were `ok`.

## 8. What S1–S3 imply for interpreting MIT-BIH detections

1. **Do not interpret `lle_chaos_test` alone on windows with ectopic beats.**
   - With isolated premature beats at only 2 % of beats (5 in 256), the LLE
     test fires in 97 % of otherwise linear-Gaussian windows. At 5 % and 10 %
     it fires in 98–99 %.
   - An LLE-only "chaos" detection in an abnormal MIT-BIH window is therefore
     expected from the ectopics themselves and says nothing about
     deterministic chaos.
   - It also cannot separate abnormal from normal windows by dynamics: the
     arrhythmia's beat-timing irregularity alone triggers it.
2. **The AND detector is largely protected from isolated ectopics, but not
   completely.**
   - It fires in 0.5 %, 1 % and **3.5 % (95 % CI 1.7–7.0 %)** of windows at
     2 %, 5 % and 10 % ectopics, against 0.3 % on the clean N1 null.
   - It grows with the ectopic rate, and with m (12/200 at m = 4 at 10 %).
   - The 10 % rate is the pipeline's abnormal-window threshold
     (`CLASSIFIER_ABNORMAL_FRACTION`). Abnormal MIT-BIH windows may contain far
     more ectopics (bigeminy is 50 %) or runs of them, which were **not**
     tested here. The false-positive rate there is unknown and may be higher.
   - An AND detection rate in abnormal MIT-BIH windows should be compared
     against a few percent expected from isolated ectopics, not against 0.
3. **Non-chaotic nonlinearity triggers the LLE test** (S1 SETAR: 13 %; N5
   static warp: 17.7 %). The IAAFT null is linear-Gaussian, so rejecting it
   means "not linear-Gaussian", not "chaotic". The UPO gate supplies the
   specificity to chaos; the AND stayed at 1/200 on SETAR.
4. **Quantization at 1/360 s is not a concern** for either component (S3).
5. **Nonstationarity costs power, not specificity.** Trends and steps did not
   raise the AND false-positive rate (0/300), but a trend comparable to the RR
   SD cuts UPO power on chaotic signals to 20–43 %. A negative AND on a
   drifting MIT-BIH window is weak evidence against chaos.
6. **Continuous-time chaos is not detected.** A negative result means "no
   evidence of low-dimensional map-like chaos", not "no chaos". The tested
   sampled flows are Rössler at about 6 samples per orbit and Mackey–Glass at
   step 6; both were missed at every m.
7. **Windows can be lost to errors.** `analyze_segment` raised on 3 of 5,400
   windows, all strongly chaotic, when the production LLE embedding was too
   short. On MIT-BIH such windows would be missing, not negative.
   `keep_upo_on_short_lle_embedding` (Phase 3 A4) prevents the loss, but it
   was not part of the tested configuration.

## 9. Adoption (PASS branch of Step 5)

Two opt-in helpers were added to `final_pipeline.py`. The defaults are
unchanged: `CFG` is not modified and no existing function's behaviour changes.
- **`combined_chaos_config(config=None)`** returns
  `phase4_upo_config(replace(config or CFG, lle_chaos_test=True))`. This is the
  tested configuration. It is to be used with `analyze_segment` on raw RR.
- **`combined_chaos_detected(segment_result)`** returns the AND decision:
  `lle_chaos_test["detected"] and upo["instability_gate_detected"]`. It raises
  `ValueError` if either component is missing.

Checks:

- `tests/test_phase5_combined_chaos.py` (7 tests):
  - the helper returns exactly the tested configuration;
  - `CFG` defaults are unchanged;
  - missing components raise;
  - the stored decisions of 5 test windows are reproduced (N1, S2 10 %, P1,
    P3 Hénon, G3).
- `adoption_check.py`: 120 stored test windows (seeds 3000–3004 × 24
  conditions) recomputed through
  `combined_chaos_detected(analyze_segment(rr, combined_chaos_config()))`, with
  **0 differences** in LLE, UPO and AND (`results/test/adoption_check.json`).
- Ground rules (`experiments/phase3_lle/groundrule_check.sh`):
  - pytest 384 passed plus the pre-existing
    `test_hybrid_period1_stability...[1-prl_norm]` failure (385/385 with numpy
    AVX-512 dispatch disabled);
  - same-machine replicability 44/44 bitwise identical;
  - `run_phase2e --replicate` 13/44, with `replicability.json` identical to the
    pre-change file.

## 10. Limitations

- **Synthetic RR is not RR.** The nulls match published norms for mean,
  variability and LF/HF. They are still idealized:
  - Gaussian;
  - stationary within the window except where a trend or step is added;
  - with no respiration–heart-rate coupling, no beat-detection jitter and no
    missed or extra beats.

  The positives are chaotic maps and flows rescaled to RR statistics, not
  physiological dynamics.
- **The NSRDB comparison was not possible**: physionet.org is blocked here.
  The literature values were checked against web-search results, since the
  full texts were also blocked.
- **Ectopic model.** Only isolated premature beats with a full compensatory
  pause, at most 10 % of beats and at least 3 beats apart, were tested. Not
  tested:
  - couplets, runs, bigeminy or trigeminy;
  - interpolated beats or non-compensatory (atrial) pauses;
  - atrial fibrillation;
  - missed or extra beat detections.

  These are exactly what abnormal MIT-BIH windows contain.
- **Warping result.** The LLE test's 17.7 % false positives on N5 cannot be
  attributed to the static nonlinearity alone: a warped *unquantized* condition
  was not in the design, and the IAAFT's imperfect spectrum match under strong
  skew is a plausible cause. The AND held (1/300).
- **The PASS margin is large, but N is finite.** The worst per-condition
  upper bound is 1.9 %. The rule used 300 windows per null. This establishes
  specificity for these generators only.
- **Generalization parameters were choices.** G2's sampling step (1.0) and
  G3's step (6.0) were predeclared but not derived from a physiological
  argument. Coarser or finer sampling may change the flow results.
- **Detector frozen by design.** Known weaknesses (trend sensitivity of Level B;
  the LLE test's one-sided failure on smooth flows; LLE sensitivity to
  ectopics) were documented, not fixed.
- **Machine dependence.** Everything was run on one machine, with the ground
  rules as in Phases 3 and 4.
- **Branch state.** Phase 5 was developed on the branch backing the open
  Phase 4 PR #2, so that PR now also contains the Phase 5 commits.

## 11. Reproducibility

```bash
python -m experiments.phase5_rr.equivalence_check                      # step 0
python -m experiments.phase5_rr.realism                                # step 1 (dev)
python -m experiments.phase5_rr.run_phase5 --phase dev --seeds 0-9 --tag runtime
python -m experiments.phase5_rr.m_path_check
python -m experiments.phase5_rr.run_phase5 --phase test                # guarded
python -m experiments.phase5_rr.analysis --phase test --file test
python -m experiments.phase5_rr.decide   --phase test --file test
python -m experiments.phase5_rr.report_extras
python -m experiments.phase5_rr.phase4_m_check                         # step 3
python -m pytest tests
```

Commits:
- `cdd6a57` step 0;
- `5a7bbfa` generators;
- `1affece` runner and rule;
- `03787bf` preregistration;
- `579c30d` test results;
- `9c26ac2` step 3;
- the helpers and this report in the final commits.
