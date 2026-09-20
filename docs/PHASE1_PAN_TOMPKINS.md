# Phase 1: Pan-Tompkins Detector Correction and Validation

Research-development record. Factual and technical; every statement below is
tied to code or tests present in commit `72eea3f`.

## 1. Baseline

| Item | Value |
|---|---|
| Baseline repository commit | `9383a31` (Add current arrhythmia detection pipeline) |
| Phase 1 commit | `72eea3f` (Correct Pan-Tompkins implementation and add validation tests) |
| Files changed in Phase 1 | `final_pipeline.py`, `tests/synth_ecg.py`, `tests/test_pt_detector.py`, `tests/test_pt_filters.py` |

Phase 1 focused specifically on the Pan-Tompkins (PT) QRS detector and its
validation (`_pan_tompkins_filter_chain_200hz`,
`_pan_tompkins_processing_delay_200hz`, `_RRTracker`, `_pan_tompkins_decide`,
`detect_r_peaks`, and the related `pt_*` fields of `PipelineConfig`).

The following were intentionally **not** modified in Phase 1: downstream RR
extraction, RR correction, TDMI, Cao embedding dimension, Takens embedding,
Rosenstein LLE, UPO detection, surrogate testing, TDA, classifier logic,
labels, and experimental design.

## 2. Pan-Tompkins corrections

Reference for "baseline" below is `9383a31`.

1. **Low-pass filter recurrence / gain normalization.** The baseline divided
   the entire recurrence by 36, which is not the intended filter (its DC gain
   is not 1; `test_lp_baseline_bug_is_real`). The production filter is now the
   FIR-equivalent 11-tap triangle `conv(ones(6), ones(6)) / 36`, with unit DC
   gain.
2. **High-pass filter.** Implemented from the intended normalized transfer
   function `H(z) = z^-16 - (1/32)(1 - z^-32)/(1 - z^-1)`, realized as a 32-tap
   FIR (`-1/32` on every tap, `+1` added at tap 16). The literal reading of the
   printed equation places a pole at `z = -1` that does not cancel the
   `(1 - z^-32)` zeros; this is checked in
   `test_hp_printed_eq6_read_literally_has_pole_at_minus_one`.
3. **Derivative filter.** Baseline coefficients were incorrect
   (`test_derivative_baseline_coefficients_were_wrong`). Now
   `d[n] = (1/8)(x[n] + 2x[n-1] - 2x[n-3] - x[n-4])`, i.e. the paper's
   `(1/8)(-z^-2 - 2z^-1 + 2z + z^2)` with a two-sample causal delay. The `1/T`
   scale factor is omitted (T = 1 sample); thresholds are data-adaptive.
4. **Startup / edge handling.** The record is extended on the left by 128
   samples of edge replication of `x[0]`, filtered from zero state, and the
   padded outputs are discarded. No negative-index reads occur and no ECG
   sample is altered. Resampling of non-200 Hz input uses `padtype="line"` to
   avoid a step transient from DC offset.
5. **Full 200 ms candidate refractory spacing.** Candidate peaks of the
   integrated signal are now found with `distance = refractory` (baseline used
   `refractory // 2`). Candidates inside the refractory window of the last
   accepted beat are rejected.
6. **Functional searchback.** Rejected sub-threshold candidates are retained
   until the next accepted QRS; once `RR_MISSED` elapses, the maximal reserved
   peak above `0.5 * THRESHOLD_I1` and outside the refractory period is
   accepted. Repeats while further missed intervals remain. A tail check runs
   at the end of the record.
7. **RR_AVERAGE1 / RR_AVERAGE2 tracking** (`_RRTracker`). AVERAGE1 is the mean
   of the last 8 RRs regardless of value; AVERAGE2 is the mean of the last 8
   RRs inside the limits. A single previous RR is not used.
8. **RR_LOW / RR_HIGH / RR_MISSED.** Limits are `0.92`, `1.16`, and `1.66`
   times RR_AVERAGE2 (`pt_rr_low_factor`, `pt_rr_high_factor`,
   `pt_searchback_seconds`). The rhythm is treated as irregular only when 8 RRs
   exist and not all lie within [LOW, HIGH]; in that case only THRESHOLD_I1 is
   halved.
9. **Searchback SPKI update.** Searchback-accepted peaks update SPKI with
   weight 0.25 (`SPKI = 0.25*h + 0.75*SPKI`); regular detections use 0.125.
10. **T-wave slope window.** For candidates within 360 ms of the previous QRS,
    the maximum absolute derivative is taken over the integration window
    ending at the candidate (not a single sample) and compared with 0.5 times
    that of the previous QRS. A premature true QRS inside 360 ms is not
    rejected on this basis if its slope is sufficient.
11. **Integration-window configuration.** Window length is derived from
    `pt_integration_seconds` (0.150 s, 30 samples at 200 Hz) and passed
    consistently through the filter chain, decision logic, and delay estimate.
12. **Nominal processing-delay handling.** The delay is computed as the
    centroid of the squared-derivative impulse response plus `(N-1)/2` for the
    integrator, rounded (37.486 -> 37 samples at N = 30). See Section 4.

## 3. Validation

The Phase 1 test suite contained **73 tests, all passing** (re-confirmed on the
committed tree with `python -m pytest tests -q`: `73 passed`). Coverage:

- **Filter-response tests** (`tests/test_pt_filters.py`): low-pass unit gain and
  triangle impulse response, DC gain and cutoff; high-pass DC rejection and
  low-frequency rejection; combined band-pass response; derivative impulse,
  ramp, and sinusoid response; delay measurement.
- **FIR / transfer-function equivalence tests**: FIR vs. recursion for the
  low-pass (impulse) and high-pass (impulse and random input).
- **Chain/startup tests**: no negative-index wraparound, no startup transient
  for constant input, equivalence to unpadded steady-state output beyond the
  transient, input preserved.
- **Synthetic ECG detector tests** (`tests/synth_ecg.py`,
  `tests/test_pt_detector.py`): clean beats, baseline wander/noise/polarity,
  rate ranges, irregular RR, premature beats, pauses, broad QRS, noisy beats,
  short records.
- **Sampling-rate tests**: clean detection across multiple sampling rates.
- **DC-offset robustness tests**: multiple offsets, and offset invariance of
  the result.
- **Searchback tests**: weak-beat recovery, quarter weighting with trace,
  no invented beats for a genuinely missing beat, persistence of rejected
  candidates until the next QRS.
- **RR-tracking tests**: AVERAGE1 vs. AVERAGE2 divergence, limit factors,
  8-beat requirement for irregularity, scaling with rhythm, first-threshold-only
  halving.
- **T-wave discrimination tests**: tall T waves not detected as QRS, slope
  window coverage, slope rule, premature QRS within 360 ms retained.
- **Refractory / noise tests**: refractory spacing of detections and full
  refractory spacing of candidates.

**MIT-BIH validation was NOT performed in Phase 1.** All detector evidence is
from synthetic signals.

## 4. Methodological qualifications

- The implementation is **Pan-Tompkins-style**. It is not claimed to be a
  literal line-by-line reproduction of the 1985 algorithm.
- The production implementation uses the **intended FIR-equivalent filter
  forms**, not a literal transcription of the printed recurrences (see the
  high-pass discussion above).
- The **~37-sample processing delay is a nominal compensation**, not an exact
  universal delay. The integrator acts on the squared signal, so the integrated
  peak position depends on QRS morphology (about 30-45 samples on the synthetic
  beats examined per the code comment). Very broad QRS complexes may fall
  outside the refinement window.
- The **±80 ms raw-ECG fiducial refinement** (largest deviation from the local
  median in the original ECG) is an implementation choice. The reported
  fiducial is this refined location and must not be described as the intrinsic
  Pan-Tompkins timing output. Refined detections closer than the refractory
  period are merged by keeping the stronger excursion.
- **Startup edge replication** (128 samples of `x[0]`) is an implementation
  choice, not a claim about the original paper.
- **Final decision channel.** The detector uses the **integrated channel** for
  the final accept/reject decision. Filtered-channel SPKF/NPKF quantities are
  maintained but are not used as a mandatory second gating channel, and the
  filtered-channel threshold halving in irregular rhythm (paper Eq. 23) is not
  implemented. The source describes both integrated and band-pass-filtered
  signals in its detection logic, so this is a **documented methodological
  discrepancy requiring separate evaluation**, not exact reproduction.
- Phase 1 has **not** established full Pan-Tompkins fidelity on MIT-BIH.

## 5. Scope boundary

Phase 1 deliberately did **not** attempt to resolve:

- UPO methodology
- LLE methodology
- TDMI / Cao methodological questions
- RR correction methodology
- UPO coverage interpretation
- TDA
- classifier design
- MIT-BIH experimental validation

## 6. Reproducibility

- Test files are in `tests/` (`test_pt_filters.py`, `test_pt_detector.py`,
  `synth_ecg.py`).
- The Phase 1 tests passed (73 of 73).
- Code and tests are committed in `72eea3f`.
- Run with: `python -m pytest tests -q`.
- MIT-BIH validation remains a future validation step.

## 7. Phase 2: UPO reconstruction and validation

Phase 2 will independently address:

- local-linear Jacobian correction
- fixed-point transformation verification
- period-p cyclic reduction
- candidate residual calculation
- candidate stability evaluation at the candidate itself
- candidate extraction / peak-picking validation
- UPO coverage validation
- end-to-end validation on logistic, Hénon, periodic, noisy-chaotic, and
  nonchaotic/noise systems
