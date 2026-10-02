# Phase 9 development diagnosis (DEV material only)

All material in this file uses DEVELOPMENT families and DEVELOPMENT seeds (9000–9499):
- the 17 Phase 5–6 nulls;
- Mackey–Glass (Phase 8 DEV family);
- the Phase 5 positives;
- the DEV-only circle-map controls (`devmaps.py`).

No TEST family (phase_reset, coupled_vdp, av_node, KTz) and no TEST seed (≥ 9500) was used.
Prior knowledge from Phase 8 TEST is disclosed in HANDOFF.md.

## 1. What was removed before the diagnosis

- **Waveform arm.**
  - PPS, CS and TS reject a strictly periodic ECG 5/5, so all three are UNUSABLE under B3.
  - PPS and CS also fail their own Rössler verification (10/10 rejections of periodic + white
    noise).
  - Cost is about 200 s per window (B5).
- **T-wave and QT features.**
  - The NeuroKit2 DWT delineator fails QTDB verification.
  - Fixed-window RTp / Ta fail at nstdb 12 dB (T apex SD 80 ms).
- Consequence: the candidates work on the **beat (RR) level**, with beat annotations.

## 2. Statistics bank (`dev_diag.py`, `diag_summary.py`, `results/dev/diag_summary.md`)

**What it computes.** Each 512-interval window is compared with 19 ectopy-preserving IAAFT
surrogates (Phase 8 `ep_iaaft`) on:
- NLP error, h = 1–5, m = 2–5;
- Sugihara–May ρ(h);
- SDLE (4 shells), FSLE (14 levels), ε-entropy;
- PE / CECP, RQA DET;
- NLP after GHKSS.

**Coverage.** The run was stopped after 379 of the 1,820 planned windows to free the CPUs for
calibration. The groups were interleaved (random task order), so every group is represented,
with 14–65 windows each. The pattern below was already unambiguous at that point.

**Findings** (group medians; full tables in `results/dev/diag_summary.md`):

1. **Every determinism statistic separates deterministic from linear-stochastic, but not chaotic
   from non-chaotic.** Median z against EP-IAAFT, chaotic MG / non-chaotic MG / circle-map QP /
   locked:

   | statistic | chaotic MG | non-chaotic MG | QP | locked |
   |---|---|---|---|---|
   | NLP (m 3, h 1) | −8.6 | −10.0 | −11.2 | −14.4 |
   | ε-entropy (m 2, 0.5 SD) | −5.6 | −7.1 | −9.3 | −14.2 |
   | PE complexity C | +14.0 | +2.9 | +6.1 | +4.2 |
   | DET | +7.6 | +4.6 | +8.8 | +7.7 |

   This is the Phase 8 failure mode, quasi-periodic coupled vdP, reproduced on development data
   (the circle map).
2. **Ectopy in non-chaotic deterministic rhythms looks like chaos even with the RR-rule mask.**
   - The Phase 8 C3 statistic (NLP vs EP-IAAFT with the pipeline's RR outlier mask) has a
     median zmax of about 20 for non-chaotic MG + E3_10, QP + E3_10 and locked + E3_10.
   - Error growth on the RR-rule-masked series is 0.66 for non-chaotic MG + E3_10 (smoke test),
     so the outlier rule misses ectopy-related intervals.
3. **The forecast-error growth separates chaos from periodic / QP when ectopy is absent or
   annotated.** With no ectopy, E(h 5) − E(h 2) at m = 3 has these medians:

   | chaotic MG | non-chaotic MG | QP | locked | nulls |
   |---|---|---|---|---|
   | 0.37–0.40 | 0.00 | 0.006 | 0.00 | ≈ 0 or negative |

   The nulls fail the determinism gate.
4. **FSLE at δ = 0.32 SD (level 10)** also separates (medians; nulls fail the determinism gate):

   | chaotic MG | Phase 5 positives | non-chaotic MG | QP | locked | nulls |
   |---|---|---|---|---|---|
   | 0.39–0.79 | 0.30 | 0.11–0.19 | 0.15 | 0.18 | 0.33–0.38 |

5. **Measures that did not help:**
   - SDLE at 512 beats: chaotic large-shell λ 0.016–0.079 vs QP −0.005 and locked 0.024;
     overlapping, few pairs.
   - GHKSS: |z| is lower after noise reduction, as METHODS.md predicted for short series.
   - RQA DET and PE/CECP: they respond to any determinism (finding 1).
6. **The null size of the EP-IAAFT NLP test is liberal.** Nulls have median z ≈ −2 to −3 at
   h = 1, so thresholds must be calibrated (Phase 8 recipe), not taken as z = 1.645.

**Annotation-masked smoke test** (`candidates9.py`; DEV seed 9400; nstdb mix 12, detected
beats, jitter):

| window | zmax | G |
|---|---|---|
| MG τ = 23 + E3_10 | 12.9 | 0.72 |
| MG τ = 15 + E3_10 (non-chaotic) | 9.7 | 0.012 |
| QP + S2_5 | 11.3 | 0.003 |
| N1 null | 1.7 | −0.54 |
| MG τ = 23 + E1 (bigeminy) | not analysable: every vector touches an ectopic interval | — |

## 3. Candidates (≤ 6; `candidates9.py`)

| # | name | idea | null control |
|---|---|---|---|
| 1 | `k1_frozen512` | frozen Phase 6 detector (Phase 8 C1), baseline | LLE IAAFT + UPO surrogates |
| 2 | `k2_mnlp_ann` | annotation-masked NLP determinism only (zmax ≥ Z_DET). The sensitive reference, expected to fail on QP / ectopic non-chaotic | masked IAAFT |
| 3 | `k3_mnlp_growth_ann` | k2 AND forecast-error growth G ≥ G_MIN (chaos signature) | masked IAAFT + growth gate calibrated on DEV non-chaotic |
| 4 | `k4_mnlp_growth_rr` | as k3 with the RR-rule mask (no annotations needed) | as k3 |
| 5 | `k5_mnlp_fsle_ann` | k2 AND FSLE(0.32 SD) ≥ F_MIN (Part D measure as the chaos gate) | as k3 |

**Thresholds.** Calibrated on DEV only (`calibrate9.py`, rules in its docstring):
- Z: the smallest value with ≤ 2 % of each DEV null condition at or above it, + 0.5.
- Gates: the smallest value with ≤ 2 % of each DEV null and DEV non-chaotic condition
  passing both gates, + 0.05.

**Expected limitations, stated before TEST:**
- **Bigeminy (E1).** Nearly every interval is ectopy-related, so the annotation-masked
  candidates cannot analyse such windows. They count as not detected, which caps power on the
  E1 third of the primary pool.
- **KTz (TEST morphology family).** Its RR is the pacing input, so no RR candidate can detect
  its chaos. All KTz CHAOTIC windows are expected to be missed by every candidate.
- **Annotations.** The `ann` candidates assume beat labels of annotation quality, as in the
  PhysioNet databases that Part F would use. On the synthetic data they are the generator's
  beat types matched to the detected beats; unmatched detections are labelled 'X' and masked.
