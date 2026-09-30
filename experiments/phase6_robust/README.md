# Phase 6: robustness of the combined chaos detector

Evaluation and preregistered selection. There are three parts:
- crash fix (Part A);
- heavier ectopy with annotation-edited NN variants (Part B);
- detrending candidates for the drift power loss (Part C).

No MIT-BIH record is loaded in this phase.

**BASELINE-K** = `combined_chaos_config()` with
`keep_upo_on_short_lle_embedding = True` (`methods.py`), run through
`fp.analyze_segment` on the RR series. The decision is `combined_chaos_detected`:
`lle_chaos_test` AND the C2 instability-gated UPO test, m = 2.

**Signals.** The Phase 5 generators are reused unchanged
(`experiments/phase5_rr/`: 256 RR intervals, 1/360 s quantization, Phase 2E seed
scheme, paired base and modifier streams). The Part B conditions are added in
`systems.py`, on the **N1 linear_rr base**: the same realization and RR scale
as N1 for a given seed.
- Development seeds: 0–99 for null-type conditions (N*, S*, E*), 0–29 for P*
  and G*.
- **Test seeds start at 4000.**

## Part B: ectopy patterns

Every premature (ectopic) interval is c × the underlying sinus RR, with
c ~ U(0.6, 0.8). This coupling range is the same as Phase 5 S2.

**Ventricular-type patterns (E1–E4, E6)** assume a **full compensatory pause**.
The ectopic beat does not reset the sinus node, so the next conducted sinus beat
keeps its scheduled time:
- a group of k premature beats replaces k scheduled sinus beats;
- the interval after the group is the remainder up to that scheduled sinus beat;
- the interval count and total duration are preserved (checked).

| id | condition | definition used | source |
|---|---|---|---|
| E1 | `E1_bigeminy` | every other beat premature: premature at interval s, s+2, … with s ∈ {1, 2} (50 % of beats) | bigeminy = each premature complex alternates with a normal complex, i.e. every second beat premature [LITFL PVC; AMBOSS PVC] |
| E2 | `E2_trigeminy` | every third beat premature: premature at s, s+3, … with s ∈ {1, 2, 3} (33 %) | trigeminy = each premature complex follows two normal complexes [LITFL PVC; AMBOSS PVC] |
| E3 | `E3_couplets_{5,10}pct` | pairs of consecutive premature beats; round(rate·256/2) = 6 or 13 couplets (4.7 % / 10.2 % of beats); both premature intervals c·RR with independent c | couplet = two consecutive premature ventricular complexes [LITFL PVC; AMBOSS PVC] |
| E4 | `E4_runs` | U{1, 2, 3} runs per window, each of U{3, …, 6} consecutive premature beats; first coupling c·RR, then an absolute run cycle U(0.40, 0.55) s (109–150 bpm) | ≥ 3 consecutive ventricular complexes at > 100 bpm (cycle < 600 ms) is ventricular tachycardia; non-sustained if it terminates spontaneously [Al-Khatib 2017] |
| E5 | `E5_atrial_{5,10}pct` | isolated premature beats **without** a full compensatory pause: premature interval c·RR_i, then the sinus node is reset and the post-ectopic interval is d·RR_{i+1}, d ~ U(1.0, 1.1). Later beats shift earlier (a non-compensatory pause). round(rate·256) = 13 / 26 ectopics | premature atrial complexes usually enter and reset the sinus node, giving a less-than-fully-compensatory pause [ECGpedia; LITFL PAC] |
| E6 | `E6_ectopic10_trend` | N3's trend draw, then S2 10 %'s isolated-ectopic draw (the Phase 5 modifier streams) | as S2 and N3 |

- **Placement.** Groups are placed one at a time at uniformly drawn starts,
  within intervals 1 … 254, at least 3 normal intervals apart. Bigeminy and
  trigeminy fill the whole window.
- **Compensatory pause.** The next sinus beat after a ventricular premature
  beat arrives about two normal cycles after the preceding sinus beat
  [LITFL PVC; AMBOSS PVC].

### Annotation-edited (NN) variants

These are `S2_ectopic_{5,10}pct_edited`, `E1_bigeminy_edited`,
`E2_trigeminy_edited` and `E3_couplets_10pct_edited`. They mimic analysing
annotation-edited normal-to-normal intervals, the standard HRV practice
[Task Force 1996].
- Every interval that begins or ends at an ectopic beat (premature, couplet or
  run, and compensatory) is removed.
- It is replaced by linear interpolation (by beat index) between the nearest NN
  intervals before and after its block. A block at a window edge takes the
  nearest NN value.
- Editing is applied to the quantized raw series of the same seed. Interpolated
  values are not re-quantized.

**Bigeminy consequence.** In bigeminy only the two edge intervals are NN, so the
edited window is by construction a near-straight line between them. This is
reported, not patched. The NN fractions (development seed 0) are:
- S2 5 % / 10 %: 90 % / 80 %;
- E1: 0.8 %;
- E2: 34 %;
- E3 10 %: 85 %.

## Part C: detrending

`fp.detrend_rr` is selected by `PipelineConfig.rr_detrend`, default `"none"`.
It is applied inside `analyze_segment` to `rr_dynamics`, before the LLE, the
UPO analysis and `lle_chaos_test`, so all surrogates come from the detrended
series. The window mean is added back.
- **Linear:** subtract the least-squares line.
- **Moving median:** subtract a centred running median of `rr_detrend_window`
  beats, with reflected edges.
- **Smoothness priors** [Tarvainen 2002]:
  z_stat = (I − (I + λ² D₂ᵀD₂)⁻¹) z, over beat index, with λ =
  `rr_detrend_lambda`. Tarvainen et al. use λ = 500 on 4 Hz resampled series
  (cut-off about 0.04 Hz); here the series is beat-indexed, so λ is tuned on
  development seeds.

`diagnose_trend.py` checks the mechanism (C1) and `results/dev/tables/tuning.md`
records the tuning (C2).

## References

- [LITFL PVC] Premature Ventricular Complex (PVC), LITFL ECG Library.
  https://litfl.com/premature-ventricular-complex-pvc-ecg-library/
- [AMBOSS PVC] Premature ventricular contractions, AMBOSS.
  https://www.amboss.com/us/knowledge/premature-ventricular-contractions
- [Al-Khatib 2017] Al-Khatib SM, et al. 2017 AHA/ACC/HRS Guideline for
  Management of Patients With Ventricular Arrhythmias and the Prevention of
  Sudden Cardiac Death. *Circulation* 138:e272–e391, 2018 (published online
  2017). https://www.ahajournals.org/doi/10.1161/CIR.0000000000000549
- [ECGpedia] Ectopic Complexes, ECGpedia.
  https://en.ecgpedia.org/wiki/Ectopic_Complexes
- [LITFL PAC] Premature Atrial Complex (PAC), LITFL ECG Library.
  https://litfl.com/premature-atrial-complex-pac/
- [Task Force 1996] Heart rate variability: standards of measurement,
  physiological interpretation and clinical use. *Circulation*
  93(5):1043–1065, 1996.
- [Tarvainen 2002] Tarvainen MP, Ranta-aho PO, Karjalainen PA. An advanced
  detrending method with application to HRV analysis. *IEEE Trans Biomed Eng*
  49(2):172–175, 2002. https://doi.org/10.1109/10.979357

The definitions were checked against web-search results. The full texts of
these pages were not fetched: the environment's proxy blocks most publisher
sites.
