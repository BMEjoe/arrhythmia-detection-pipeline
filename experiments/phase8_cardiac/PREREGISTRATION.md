# Phase 8 preregistration: candidate chaos detectors for cardiac-relevant chaos

This file is committed and pushed **before any TEST seed is generated and before any
detector is run on any TEST-model series**.

**What exists at that point:**
- for the TEST models (phase-resetting map, coupled vdP, AV node): the implementations,
  their verification against published results, the ground-truth Lyapunov exponents from
  the model equations, and the realism statistics and plots with realism seeds
  900000–900009. No detector output on any of them.
- all candidate design and calibration used DEVELOPMENT material only: the DEV model
  Mackey–Glass, the Phase 5–6 null generators, and the Phase 5 positives, with development
  seeds 0–999 (`DIAGNOSIS.md`, `results/dev/`).

**What is frozen after this file.** Nothing below changes, and neither do:
- `models/`, `series.py`, `nulls.py`, `methods.py`;
- `run_test.py`, `decide.py`;
- `results/ground_truth/regimes.json`;
- `final_pipeline.py`.

Any unavoidable change goes in a dated amendment at the end, saying whether any TEST
output had been seen. `run_test.py` reuses the Phase 3 guard: it refuses to start unless
this file is committed, unmodified, pushed, and lists every method below.


**Frozen code** (git `360ea37`; SHA-256 prefixes):

| file | SHA-256 prefix |
|---|---|
| `methods.py` | 5e56f35bea69c36b |
| `series.py` | fcd0ee3d3bbe7719 |
| `nulls.py` | 01cd54fc4d45aa1c |
| `run_test.py` | 0778a84020d29ed9 |
| `decide.py` | 56b4c1939d21b898 |
| `calibrate.py` | 28980ce47075745b |
| `models/av_node.py` | caf89ecdad6bf116 |
| `models/coupled_vdp.py` | 35754c2595f68258 |
| `models/mackey_glass.py` | d2f92763c950800a |
| `models/phase_reset.py` | 2f038fdd51438438 |
| `results/ground_truth/regimes.json` | 3226b4359a7cee3d |
| `final_pipeline.py` | 27bf93808ab45487 |

## 1. Model library and split (`MODELS.md`, `HANDOFF.md`)

| family | source | split | CHAOTIC regimes | NON-CHAOTIC regimes |
|---|---|---|---|---|
| Mackey–Glass maxima intervals | Glass & Mackey, Scholarpedia Eq. 1 | DEV | τ = 16.5, 17, 18, 20, 23, 26, 30 | τ = 14, 15, 16 |
| phase-resetting map, hiPSC-CM PRC (aggregate A) | Diagne et al., PLoS Comput Biol 2026, Eq. 1 | TEST | τ = 0.58, 0.60, 1.14, 1.16, 1.20, 1.58, 1.60 | τ = 0.38, 0.62, 0.66, 1.18, 1.74 |
| coupled modified vdP (SA, AV, HP), delayed coupling | da Silva Lima et al., Sci Rep 2024, Eqs. 1–5 | TEST | (ρ, ω) = (2, 2.7), (6, 3.3), (6, 4.0), (8, 3.3) | (5.45, 5.6), (9.6, 2.1), (2, 5.6), (4, 5.6), (10, 3.3), (6, 5.6) |
| AV node | Sun et al. via Zhao & Schaeffer, Eq. 51 | TEST | none | H = 40, 45, 52, 55 ms |

Labels come from the ground-truth LE rule in `ground_truth.py`, frozen in
`results/ground_truth/regimes.json`.

**Variants** (`series.py`):
- (i) clean;
- (ii) dynamical noise, low / high;
- (iii) measurement noise, 30 / 20 dB;
- (iv) 30 dB + quantized to 1/360 s;
- (v) (iv) + ectopy: S2 5 %, E1 bigeminy, E3 couplets 10 %;
- (vi) (v) + Phase 7 measured V-beat timing jitter.

## 2. Candidates (`methods.py`; every candidate has a null test)

- method: `c1_frozen512`
  - **C1** = the frozen Phase 6 recommended detector on 512-beat windows:
    `combined_chaos_detected(analyze_segment(rr, combined_chaos_config(replace(CFG, keep_upo_on_short_lle_embedding=True), detrend=True)))`.
  - Null control: the LLE IAAFT test (99 surrogates) AND the UPO surrogate test plus
    instability gate.
  - No embedding search. DIAGNOSIS D1 shows the UPO component, not the LLE embedding, limits
    power.
- method: `c2_nlp_iaaft512`
  - **C2** = nonlinear-prediction determinism test on 512-beat windows.
    - Statistic: robust normalized leave-one-out local-average prediction error (k = 5
      nearest neighbours, Theiler 5, horizon 1, lag-1 embedding).
    - For m = 2, 3, 4, 5: z_m = (mean − observed) / SD over 39 IAAFT surrogates
      (`fp.iaaft_surrogate`, data-seeded RNG).
    - Detected iff zmax = max_m z_m ≥ **Z_C2 = 16.92**.
- method: `c3_nlp_ep256`
  - **C3** = the same statistic on 256-beat windows against **ectopy-preserving IAAFT**
    surrogates (`methods.ep_iaaft`).
    - Intervals flagged by `fp.correct_rr_intervals` (CFG) are re-inserted at their
      positions in IAAFT surrogates of the cleaned series.
    - Detected iff zmax ≥ **Z_C3 = 9.31**.
- method: `c4_titration256`
  - **C4** = Poon–Barahona noise titration on 256-beat windows, as described in PNAS 2001
    and arXiv:nlin/0606032 (`methods.vwk_nonlinearity`, `noise_titration`).
    - VWK polynomial autoregression with κ = 1–6 and d = 1–3; C(r) = log ε(r) + r/N.
    - Nonlinear iff C_nl < C_lin and the nested F-test gives p < 0.01.
    - Detected iff the noise limit NL > 0.
- method: `baseline_frozen256` (reference, never eligible)
  - The frozen detector at 256 beats (Phases 6–7) on the same TEST signals at 256 beats.

**Threshold calibration** (`calibrate.py`, `results/dev/calibration.json`; development data
only). Z is the smallest value such that each of the 17 PASS null conditions has at most
2 % of its development windows at or above Z, plus 0.5:
- C2: 512 beats, seeds 0–99;
- C3: 256 beats, seeds 0–99.

The DEV non-chaotic Mackey–Glass regimes at (iv) were checked to lie below Z.

## 3. TEST seeds and windows (`run_test.py`)

TEST seeds are 5000–5099, never used before. Each window is generated once per length,
and every method of that length evaluates it. Lengths:
- 256: C3, C4, baseline;
- 512: C1, C2.

| group | content | seeds | windows per length |
|---|---|---|---:|
| null | the 17 PASS nulls: N1–N6, S1, S2 5 % and 10 %, E1, E2, E3 5 % and 10 %, E4, E5 5 % and 10 %, E6 (all quantized, `nulls.py` at the length) | 5000–5099 | 1,700 |
| noncha | every NON-CHAOTIC regime (15 TEST + 3 DEV) at (iv) | 5000–5099 | 1,800 |
| chaos_primary | 11 TEST CHAOTIC regimes × (iv, vi_S2_5, vi_E1, vi_E3_10) | 5000–5029 | 1,320 |
| chaos_second | TEST CHAOTIC × (i, ii lo, ii hi, iii 30, iii 20, v_S2_5, v_E1, v_E3_10) | 5000–5009 | 880 |
| dev_chaos | 7 DEV CHAOTIC × (iv, vi ×3) | 5000–5009 | 280 |
| noncha_ect | NON-CHAOTIC × (vi ×3) | 5000–5009 | 540 |
| flows | Phase 5 G2 Rössler, G3 Mackey–Glass | 5000–5099 | 200 |

## 4. Primary rule (`decide.py`)

**PASS** iff the detections are at most floor(0.07 N) = **7/100** on EACH of the 17 null
conditions and EACH of the 18 NON-CHAOTIC regimes at (iv), at the candidate's window
length. An analysis error counts as not detected.

**WINNER:** among passing candidates, the highest pooled detection count on the TEST
CHAOTIC regimes at (vi), i.e. vi_S2_5 + vi_E1 + vi_E3_10, 990 windows.

**Tie-breaks:**
1. higher pooled detection at (iv), 330 windows;
2. candidate index.

**If no candidate passes, that is the result.**

## 5. Secondary outcomes (report only)

- Detection by model, regime, window length and variant.
- The baseline frozen detector on the same signals at 256 beats.
- Phase 5 G2 / G3 detection.
  - Note: G2 and G3 entered the development diagnosis, so they are not out-of-sample.
- Non-chaotic regimes under ectopy.
- The DEV Mackey–Glass chaos at test seeds.
- Component decisions (LLE, UPO) for C1 and the baseline.

## 6. Runtime and adoption

**Runtime.** Measured on development windows under 4-worker load:
- median s per window: frozen 512 = 9.6, C2 512 = 3.9, C3 256 = 1.2, C4 256 = 0.1, frozen 256 (baseline) = 4.1, plus series generation (coupled vdP about 1-3 s)
- 13,440 windows in all. The estimate is about 6720 x 14.5 s (512) + 6720 x 6 s (256) = 1.4e5 CPU-s, i.e. about 9.5 h on 4 workers, within the ~10 h budget.

**Adoption.** If there is a winner, it is added to `final_pipeline.py` as an opt-in option
(default off) with tests. `groundrule_check.sh` must pass.

**Part E** (only if there is a winner) has its own preregistration, committed before any
nsr2db / chf2db download.

## Amendments

(none)
