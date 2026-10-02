# Phase 8: Chaos Detection on Cardiac-Model RR Intervals (close-out)

Research-development record in the style of the Phase 5–7 reports. Every
statement is tied to a file committed under `experiments/phase8_cardiac/`.
Paths are relative to that directory unless they start with `docs/`,
`experiments/phase9_waveform/` or `final_pipeline.py`.

**Status: Phase 8 was STOPPED by the user after its preregistered TEST
decision.** This close-out was written at the start of Phase 9, from results
already committed in Phase 8. **No Phase 8 run was repeated or extended for
it.** The only new code is a read-only tabulation of the committed TEST records,
`experiments/phase9_waveform/phase8_closeout/tables.py`; its output is in
`experiments/phase9_waveform/phase8_closeout/phase8_closeout_tables.md`.

**Not done, because the phase was stopped:**
- **Part D adoption.** The winner `c1_frozen512` was **not** added to
  `final_pipeline.py` as an opt-in option. No adoption tests were written.
  `final_pipeline.py` is unchanged; its SHA-256 prefix is `27bf93808ab45487`,
  as frozen in the preregistration.
- **Part E (real data).**
  - The separate Part E preregistration was **not** written.
  - nsr2db and chf2db were **not** downloaded or analysed.
  - The planned exploratory run on the Phase 7 MIT-BIH windows was **not** done.
- **Part F, the planned report steps.** No further secondary analysis or plot
  planned for the Phase 8 report was produced. This document replaces that
  report with a summary of what exists.

## 1. Provenance and design

| Item | Value |
|---|---|
| Branch | `claude/amazing-cray-9jbx6p`, from `main` at `1b4b7a0` (Phases 2E–6), after Phase 7 (head `4283743`) |
| Environment | Python 3.13.14, numpy 2.1.3, scipy 1.18.1, scikit-learn 1.9.1, wfdb 4.3.1, numba 0.68.0 (model integration only). BLAS 1 thread, 4 workers |
| Literature access (checked 2026-10-01) | **Reachable:** arXiv, PMC / Europe PMC, PLOS supplements, physionet.org. **Blocked:** APS, AIP, Elsevier, physiology.org, MDPI, Springer. Open-access supplements are in `sources/` |
| Goal | A detector for cardiac-relevant chaos (flows, delay systems, cardiac models) with a false-positive rate near zero on realistic non-chaotic RR, including ectopy |
| Discipline | Model-level split fixed in `HANDOFF.md` before any detector run on any Phase 8 model. Every candidate designed and calibrated on DEVELOPMENT material only. `PREREGISTRATION.md` committed and pushed alone in `991ede3`, before any TEST seed. `run_test.py` reuses the Phase 3 guard. **No amendment** |
| TEST run | 13,440 windows (6,720 at 256 beats, 6,720 at 512), each evaluated by every method of its length. **0 analysis errors.** Resumed once after a container restart at 12,542 windows (`7a05f91`). Each record depends only on its key, and the files hold no duplicate task IDs |

## 2. Model library (`MODELS.md`, `models/`)

Each model was implemented from an accessible published source and verified by
reproducing a published result before use. Unverifiable models were dropped.

### 2.1 Retained models

| family | source | beat events | verification (published → this implementation) |
|---|---|---|---|
| **Phase-resetting map**, hiPSC-CM PRC, aggregate A (Glass–Shrier lineage) | Diagne, Bury, …, Glass et al., *PLoS Comput Biol* 2026, PMC12900431: main-text Eq. 1; S1 Text Eqs. 1–3 (piecewise PRC); S1 Table (A = 60, φr = 0.71, B = 50, S = 0.80) | the model's own action potentials (native beat-to-beat) | τ = 0.66: stable period 3 → 3 distinct phases, LE −0.853. τ = 0.38, 0.53, 1.74: periodic → 5, 2, 8 phases, LE < 0. τ = 1.2: chaotic → 999/1,000 distinct phases, **LE +0.135 ± 0.002**. **Verified** |
| **AV-nodal conduction** (Sun, Amellal, Glass, Billette 1995) | restated in Zhao & Schaeffer, arXiv:math/0609106, Eq. 51 (primary *J Theor Biol* paper blocked). Amin 33 ms, τrec 70 ms, τfat 30,000 ms, γ 0.3 ms | His activations; interval H + A (native, seconds) | border-collision period doubling at H_bif = 56.9078 ms → A* = 130 crossing at **H = 56.90784 ms**, period 2 below and period 1 above. **Verified** |
| **Mackey–Glass** delayed feedback | Glass & Mackey, *Scholarpedia* 5(3):6908, Eq. 1; β = 0.2, γ = 0.1, n = 10 (Farmer 1982 set, as Phase 5 G3) | parabolically interpolated maxima of x(t), rescaled to mean 0.8 s, SD 0.05 s | LE at τ = 50: published 5.46e-3 (Galerkin), 5.76e-3, 5.83e-3 (arXiv:1810.01016 Table 1) → **5.65–5.90e-3**. Period-doubling route: period 1 at τ ≤ 13, period 2 at τ = 14, chaos from τ ≈ 17. **Verified** |
| **Three coupled modified van der Pol** oscillators (SA, AV, HP), delayed unidirectional coupling, forcing F_SA = ρ sin ωt | da Silva Lima, Savi, Bessa, *Sci Rep* 14:23446, 2024, PMC11458860: Eqs. 1–5, Fig. 1b parameter table, β_t = 0.1048 | upward midpoint crossings of u_HP (ventricular activations), real time = model time × β_t | normal rhythm ≈ 90 bpm → **89.4 bpm**. Pathological "high-frequency, highly dispersed" at four published (ρ, ω) → 137–256 bpm, RR CV 0.34–0.44. **Verified** |

**Rescaling** (`MODELS.md` Section 6):
- phase-resetting map: each window to mean 0.8 s, keeping its coefficient of
  variation;
- Mackey–Glass: mean 0.8 s, SD 0.05 s;
- AV node and coupled vdP: native seconds.

### 2.2 Dropped models

| model | reason |
|---|---|
| Seidel–Herzel baroreflex | **Not verifiable.** Primary paper blocked; restatement arXiv:q-bio/0603016. Six readings of ambiguous details gave a regular T of 0.93–0.98 s instead of the published 0.76 s, and none reproduced the published oscillations at θvNa = 3 s. Code kept in `models/seidel_herzel.py` for the record |
| DeBoer beat-to-beat model, incl. chaotic respiration | no accessible parameter values |
| AV node during atrial fibrillation | no accessible equations |
| Gois–Savi, restated in PMC9938421 | its Table 1 (even after fixing three evident typos) gives a stable equilibrium, not its own Table 3 chaos |
| Modulated parasystole (Diagne et al. Eq. 2) | the ectopic PRC's functional form is not given in the accessible sources |

## 3. Model-level split (`HANDOFF.md`, fixed 2026-10-01 before any detector run on any Phase 8 model)

| family | chaotic regimes? | split |
|---|---|---|
| Mackey–Glass | yes (τ ≥ 17) | **DEVELOPMENT** |
| phase-resetting map | yes | **TEST** |
| coupled modified vdP | yes (forced) | **TEST** |
| AV node | **no** (LE ≤ 0 at every valid H) | **TEST**, non-chaotic regimes only |

- Development positives could also use the Phase 5 generators P1, P2, G1–G3.
  Phase 5 G2 / G3 results are therefore **not out-of-sample**.
- Before the preregistration, only the following were computed for TEST models:
  - ground-truth LEs;
  - verification;
  - realism statistics and plots, with dedicated realism seeds 900000–900009.

## 4. Ground-truth Lyapunov exponents (`ground_truth.py`, `results/ground_truth/`)

**Scan.** 426 tasks: every grid point × 3 independent initial conditions or
seeds. Methods:
- phase-resetting map: mean log|f′(φ)| along the orbit (1-D map, analytic
  derivative), 200,000 iterations;
- AV node: QR (Benettin) on the exact 2-D Jacobian, 100,000 beats;
- Mackey–Glass: RK4 with the tangent delay equation, dt = 0.05, 2 × 10⁵ time
  units;
- coupled vdP: tangent integration of the delay system with periodic
  renormalization (Benettin), dt = 0.001, 8,000 model units.

Standard errors come from batch means.

**Labelling rule** (predeclared in `ground_truth.py`, applied mechanically):

| label | condition |
|---|---|
| **CHAOTIC** | λ − 2.576 SE > 0 for the combined estimate AND for every initial condition, AND converged: full-run and half-run estimates within 25 % |
| **NON-CHAOTIC** | λ + 2.576 SE ≤ tol. For maps tol = 0. For flows and delay systems tol = 2.576 SE, plus λ ≤ 0.1 × the family's smallest CHAOTIC λ |
| otherwise | AMBIGUOUS, discarded |

**36 regimes frozen** in `results/ground_truth/regimes.json` (`9a8035d`).

| family | split | CHAOTIC regimes (λ) | NON-CHAOTIC regimes (λ) |
|---|---|---|---|
| Mackey–Glass (per model time unit) | DEV | τ = 16.5 (0.0018), 17 (0.0052), 18 (0.0073), 20 (0.0076), 23 (0.0095), 26 (0.0095), 30 (0.0070) | τ = 14, 15, 16 (λ between −8.5e-6 and −2.3e-6) |
| phase-resetting map (per stimulus) | TEST | τ = 0.58 (0.081), 0.60 (0.158), 1.14 (0.035), 1.16 (0.185), 1.20 (0.138), 1.58 (0.081), 1.60 (0.158) | τ = 0.38 (−0.045), 0.62 (−0.725), 0.66 (−0.853), 1.18 (−0.209), 1.74 (−0.192) |
| coupled vdP (per model time unit) | TEST | (ρ, ω) = (2, 2.7) (0.085), (6, 3.3) (0.079), (6, 4.0) (0.101), (8, 3.3) (0.075) | (5.45, 5.6) (−0.0002), (9.6, 2.1) (−0.081), (2, 5.6) (−0.0001), (4, 5.6) (+0.0001 ± 0.001), (10, 3.3) (−0.017), (6, 5.6) (−0.049) |
| AV node (per beat) | TEST | none | H = 40, 45, 52, 55 ms (λ between −0.0067 and −0.0053) |

Totals: 18 CHAOTIC regimes (7 DEV, 11 TEST) and 18 NON-CHAOTIC (3 DEV, 15 TEST).

**Variants** (`series.py`):

| code | content |
|---|---|
| (i) | clean |
| (ii) | dynamical noise, two levels |
| (iii) | measurement noise, 30 / 20 dB |
| (iv) | 30 dB + quantized to 1/360 s |
| (v) | (iv) + Phase 6 ectopy: S2 5 %, E1 bigeminy, E3 couplets 10 % |
| (vi) | (v) + the Phase 7 measured V-beat R-peak jitter |

The realism table is `results/realism/realism.md`; plots are
`plots/A_examples_*.png`.

**Realism caveats:**
- the AV node runs at about 0.17–0.19 s per beat (a native His rhythm, not a
  sinus rate);
- the chaotic coupled-vdP regimes are tachycardic (mean 0.30–0.56 s) with very
  large dispersion.

## 5. Development diagnosis (`DIAGNOSIS.md`; DEV material only, seeds 0–999)

| item | finding |
|---|---|
| **D1. Frozen Phase 6 detector** (AND of the LLE IAAFT test and the instability-gated UPO test) | The UPO component is the bottleneck. LLE alone fires on 66–89 % of chaotic MG windows, but UPO passes only at τ = 23–26. AND power: 23 % clean, 0 % with bigeminy. UPO alone fires on 33 % of clean non-chaotic MG windows; the AND stays specific only because LLE does not. Flows G2 / G3: 0/5 |
| **D2. Nonlinear prediction vs IAAFT** (robust local-average prediction error, m = 2–5, 39 surrogates, zmax) | First statistic in the project that sees sampled flows (G2 / G3 zmax about 12–16). But ectopy dominates: linear RR with ectopy reaches zmax 9–13, and MG periodic orbits with ectopy 12–16. The compensatory pause is genuine nonlinear determinism that the IAAFT null destroys. Bigeminy removes the signal |
| **D3. Ectopy-preserving IAAFT** (flagged intervals re-inserted) | Lowers the null tail for E3 / E4, but not S2: only about 40 % of ectopy-affected intervals are flagged |
| **D4. Poon–Barahona noise titration** (VWK, κ ≤ 6, d ≤ 3) | Detects MG chaos in 89–100 % of windows, but also 75–100 % of every ectopy null, N4, N5, and 87–100 % of MG periodic orbits. It measures nonlinearity, not chaos. Carried as a candidate only so that its false-positive rate is measured on TEST |
| **D6. Calibrated thresholds** (smallest z with ≤ 2 % of dev windows on each of the 17 PASS nulls, + 0.5) | Z_C2 = 16.92 (binding: E3 couplets 5 %); Z_C3 = 9.31 (binding: S2 5 %) |

## 6. Preregistered rule (`PREREGISTRATION.md`, `991ede3`)

**Candidates:**

| candidate | definition |
|---|---|
| **C1 `c1_frozen512`** | the frozen Phase 6 recommended detector on 512-beat windows, with no new parameters |
| **C2 `c2_nlp_iaaft512`** | nonlinear prediction zmax vs IAAFT, 512 beats, Z ≥ 16.92 |
| **C3 `c3_nlp_ep256`** | nonlinear prediction zmax vs ectopy-preserving IAAFT, 256 beats, Z ≥ 9.31 |
| **C4 `c4_titration256`** | noise limit NL > 0, 256 beats |
| baseline `baseline_frozen256` | the frozen detector at 256 beats; reported, never eligible |

**TEST seeds** 5000–5099, never used before.

**PASS** iff detections ≤ floor(0.07 × 100) = **7/100** on EACH of:
- the 17 Phase 5–6 nulls: N1–N6, S1, S2 5 / 10 %, E1, E2, E3 5 / 10 %, E4,
  E5 5 / 10 %, E6;
- EACH of the 18 NON-CHAOTIC regimes at (iv), DEV and TEST,

at the candidate's window length. An analysis error counts as not detected.

**WINNER:** the passing candidate with the highest pooled detection on the TEST
CHAOTIC regimes at (vi), out of 990 windows. Tie-breaks:
1. pooled detection at (iv), out of 330;
2. candidate index.

## 7. TEST decision (`results/test/tables/decision.md`, `816c0db`; `decide.py` unchanged)

### 7.1 Primary outcome

| candidate | window | PASS | conditions failed | TEST chaotic at (vi) | TEST chaotic at (iv) |
|---|---:|---|---:|---:|---:|
| **C1 `c1_frozen512`** | 512 | **yes** | 0 | **66/990 (6.7 %)** | **10/330 (3.0 %)** |
| C2 `c2_nlp_iaaft512` | 512 | no | 2 | 220/990 | 141/330 |
| C3 `c3_nlp_ep256` | 256 | no | 1 | 148/990 | 73/330 |
| C4 `c4_titration256` | 256 | no | 30 | 966/990 | 330/330 |
| baseline `baseline_frozen256` | 256 | yes (not eligible) | 0 | 52/990 (5.3 %) | 10/330 (3.0 %) |

**WINNER: `c1_frozen512`, the only passing candidate.** Its power is low:
- 66/990 TEST chaotic windows detected at the most realistic variant (vi);
- 10/330 at (iv).

**Why the others fail:**
- **C2 and C3** fail on non-chaotic coupled-vdP regimes at ω = 5.6, the
  quasi-periodic-looking forced rhythms:
  - C2: (ρ, ω) = (2, 5.6) 24/100 and (5.45, 5.6) 19/100;
  - C3: (2, 5.6) 62/100.

  Both pass every Phase 5–6 null.
- **C4 (noise titration) failed specificity on every ectopy pattern:**
  - E1 85, E2 100, E3 5 % 100, E3 10 % 100, E4 100, E5 5 % 89, E5 10 % 58,
    E6 100, S2 5 % 100, S2 10 % 100 per 100 windows;
  - also N4, N5, S1;
  - and 17 of the 18 non-chaotic regimes.

  This is the TEST confirmation of DIAGNOSIS D4.

### 7.2 Full specificity table (detections per 100 TEST windows; ✗ = above the limit of 7)

| condition | C1 | C2 | C3 | C4 | baseline |
|:--|--:|--:|--:|--:|--:|
| E1_bigeminy | 0 | 0 | 0 | 85 ✗ | 0 |
| E2_trigeminy | 0 | 0 | 0 | 100 ✗ | 0 |
| E3_couplets_10pct | 1 | 0 | 0 | 100 ✗ | 0 |
| E3_couplets_5pct | 0 | 2 | 0 | 100 ✗ | 0 |
| E4_runs | 0 | 0 | 0 | 100 ✗ | 0 |
| E5_atrial_10pct | 5 | 0 | 0 | 58 ✗ | 1 |
| E5_atrial_5pct | 0 | 0 | 0 | 89 ✗ | 0 |
| E6_ectopic10_trend | 2 | 0 | 0 | 100 ✗ | 3 |
| N1_linear_rr | 0 | 0 | 0 | 3 | 0 |
| N2_power_law | 0 | 0 | 0 | 2 | 0 |
| N3_linear_rr_trend | 0 | 0 | 0 | 7 | 0 |
| N4_linear_rr_step | 0 | 0 | 0 | 98 ✗ | 0 |
| N5_linear_rr_warped | 0 | 0 | 0 | 100 ✗ | 0 |
| N6_noisy_rsa | 0 | 0 | 0 | 2 | 0 |
| S1_setar | 0 | 0 | 0 | 43 ✗ | 0 |
| S2_ectopic_10pct | 0 | 0 | 1 | 100 ✗ | 2 |
| S2_ectopic_5pct | 0 | 0 | 0 | 100 ✗ | 1 |
| av_node H = 40 | 0 | 2 | 0 | 56 ✗ | 0 |
| av_node H = 45 | 0 | 0 | 0 | 85 ✗ | 0 |
| av_node H = 52 | 0 | 0 | 0 | 84 ✗ | 0 |
| av_node H = 55 | 0 | 0 | 0 | 0 | 0 |
| coupled_vdp (10, 3.3) | 0 | 0 | 0 | 100 ✗ | 0 |
| coupled_vdp (2, 5.6) | 0 | **24 ✗** | **62 ✗** | 72 ✗ | 0 |
| coupled_vdp (4, 5.6) | 0 | 0 | 0 | 100 ✗ | 0 |
| coupled_vdp (5.45, 5.6) | 0 | **19 ✗** | 0 | 100 ✗ | 0 |
| coupled_vdp (6, 5.6) | 0 | 0 | 0 | 100 ✗ | 0 |
| coupled_vdp (9.6, 2.1) | 0 | 0 | 0 | 100 ✗ | 0 |
| mackey_glass τ = 14 (DEV) | 0 | 0 | 1 | 60 ✗ | 0 |
| mackey_glass τ = 15 (DEV) | 0 | 0 | 0 | 100 ✗ | 0 |
| mackey_glass τ = 16 (DEV) | 0 | 0 | 0 | 100 ✗ | 0 |
| phase_reset τ = 0.38 | 0 | 0 | 0 | 74 ✗ | 0 |
| phase_reset τ = 0.62 | 0 | 0 | 0 | 56 ✗ | 0 |
| phase_reset τ = 0.66 | 0 | 0 | 0 | 65 ✗ | 0 |
| phase_reset τ = 1.18 | 0 | 0 | 0 | 100 ✗ | 0 |
| phase_reset τ = 1.74 | 0 | 0 | 0 | 100 ✗ | 0 |
| **pooled, 17 nulls (1,700)** | 8 | 2 | 1 | 1,187 | 7 |
| **pooled, 18 non-chaotic regimes (1,800)** | 0 | 45 | 63 | 1,452 | 0 |

### 7.3 Secondary outcomes (descriptive; tabulated from the committed TEST records)

**TEST CHAOTIC by family:**

| family, variant | C1 | C2 | C3 | C4 | baseline |
|---|---|---|---|---|---|
| coupled vdP, (iv) | 4/120 | 39/120 | 38/120 | 120/120 | 6/120 |
| coupled vdP, (vi) | 29/360 | 109/360 | 29/360 | 354/360 | 35/360 |
| phase-resetting, (iv) | 6/210 | 102/210 | 35/210 | 210/210 | 4/210 |
| phase-resetting, (vi) | 37/630 | 111/630 | 119/630 | 612/630 | 17/630 |

**C1's detections are concentrated in two regimes:**
- phase-resetting τ = 0.60 (37/120);
- coupled vdP (2, 2.7) (24/120).

**Five of the 11 TEST chaotic regimes were never detected by C1:**
phase-resetting τ = 1.14, 1.20, 1.58, 1.60 and vdP (6, 4.0) (0/120 each).

**TEST CHAOTIC at (vi) by ectopy pattern:**

| pattern | C1 | baseline |
|---|---|---|
| S2 5 % | 25/330 | 16/330 |
| E1 bigeminy | 15/330 | 19/330 |
| E3 couplets 10 % | 26/330 | 17/330 |

**Secondary variants** (C1, 110 windows each):

| variant | C1 |
|---|---|
| (i) clean | 15 |
| (ii) dynamical noise, low / high | 21 / 28 |
| (iii) measurement noise, 30 / 20 dB | 7 / 3 |
| (v) S2 / E1 / E3 | 12 / 4 / 10 |

Measurement noise at 20 dB costs most of the remaining power.

**Components of the frozen detector:**

| component | C1, TEST chaotic (vi) | C1, Phase 5–6 nulls | role |
|---|---|---|---|
| LLE | 740/990 | 929/1,700 | fires on most windows, chaotic or not |
| UPO gate | 74/990 | 11/1,700 | supplies the specificity, and limits power, as in Phases 6–7 |

**DEV Mackey–Glass chaos at TEST seeds** (C1):

| variant | C1 |
|---|---|
| (iv) | 26/70 |
| (vi) S2 | 23/70 |
| (vi) E3 | 18/70 |
| (vi) E1 | 0/70 |

**Phase 5 flows** (100 each):

| method | G2 Rössler | G3 Mackey–Glass |
|---|---|---|
| C1 | 0 | 0 |
| baseline | 0 | 0 |
| C2 | 49 | 59 |
| C3 | 100 | 100 |
| C4 | 100 | 100 |

These are not out-of-sample; they entered the development diagnosis.

**Non-chaotic regimes with ectopy and jitter** (secondary, outside the PASS
rule; 10 seeds per regime and pattern):
- C1 fired in 21/540 windows overall;
- including **11/60 in non-chaotic coupled-vdP regimes with bigeminy** and
  7/60 with couplets;
- the baseline gives 18/540, with 8/60 in coupled-vdP bigeminy.

**C1's specificity under ectopy is therefore weaker than the PASS table
suggests.** The PASS rule tested non-chaotic regimes only at (iv).

**Median runtime** per window under 4-worker load:

| method | s |
|---|---:|
| C1 | 8.7 |
| C2 | 2.5 |
| C3 | 0.9 |
| C4 | 0.07 |
| baseline | 3.7 |

## 8. Interpretation

- **The trade-off was not broken.** From RR intervals alone, every candidate
  that was sensitive to the TEST chaotic regimes failed specificity:
  - C2 and C3 on forced quasi-periodic coupled-vdP rhythms;
  - C4 on everything, including every ectopy pattern.

  The only specific candidate, the frozen detector, detected 6.7 % of TEST
  chaotic windows at (vi).
- **Increasing the window from 256 to 512 beats** raised C1's (vi) power only
  from 52 to 66 of 990. The UPO component remains the bottleneck (DIAGNOSIS
  D1).
- **The PASS rule did not cover non-chaotic regimes with ectopy.** There C1
  reached 11/60 on one family. A future rule should include them. Phase 9
  does: its rule covers ectopy and nstdb-noise versions of every null.

## 9. Limitations

- **The phase was stopped:**
  - no adoption;
  - no real-data check;
  - no planned report analyses.

  The winner's behaviour on real RR (nsr2db / chf2db) is unknown.
- **Four TEST families were evaluated in aggregate.** The phase-resetting map,
  coupled vdP and AV node have now been seen at TEST seeds 5000–5099. Any later
  phase that reuses them must use new seeds and disclose this.
- **Model realism is limited:**
  - the AV node is not a sinus rhythm;
  - the chaotic coupled-vdP regimes are tachycardic and highly dispersed;
  - Mackey–Glass has no physiological beat definition (maxima intervals, with
    predeclared rescaling).
- **Restated sources.** The AV node uses a restatement (the primary paper is
  blocked), and the coupled-vdP parameter table was read from a figure image.
- **C1 specificity under ectopy** was measured only as a secondary outcome, on
  10 seeds per regime and pattern (Section 7.3).
- **One machine, one run,** resumed once after a container restart.

## 10. Reproducibility

```bash
python -m experiments.phase8_cardiac.ground_truth --scan && python -m experiments.phase8_cardiac.ground_truth --select
python -m experiments.phase8_cardiac.realism
python -m experiments.phase8_cardiac.dev_runner ...           # development (see DIAGNOSIS.md)
python -m experiments.phase8_cardiac.calibrate
python -m experiments.phase8_cardiac.run_test --workers 4     # guarded by PREREGISTRATION.md
python -m experiments.phase8_cardiac.decide
python -m experiments.phase9_waveform.phase8_closeout.tables  # this close-out's secondary tables (read-only)
```

Commits:
- `6afdae8` scaffold and source survey;
- `bba105a` models and verification;
- `6118115` `MODELS.md`;
- `9a8035d` ground truth;
- `72f709c`–`360ea37` Part B development, diagnosis and calibration;
- `991ede3` preregistration;
- `be4d6cc`–`4f7d27f` TEST run;
- `816c0db` decision.
