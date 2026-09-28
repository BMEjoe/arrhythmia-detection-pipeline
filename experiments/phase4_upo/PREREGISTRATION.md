# Phase 4 preregistration: UPO detector specificity and noise robustness

Committed and pushed **before any test-seed run**. After that the candidates,
parameters and decision rule are frozen; any change goes in a dated amendment at
the end of this file with the reason. `run_phase4.py --phase test` reuses the
Phase 3 guard and refuses to start unless this file is committed, unmodified and
pushed, and lists every method with a `method:` line.

This design follows `docs/PHASE4_UPO_PROPOSAL.md` Section 4, with the user's
changes taking precedence: at most 3 candidates plus the baseline; at least one
with the P1-a instability gate; robust aggregates for any hybrid-stability
gate; a runtime budget first; the primary rule below; the combined detector
with `lle_chaos_test` as a secondary outcome.

## 1. Data

- **Systems** (Phase 2E generators and seed scheme; `systems.py`):
  - `white_noise`, `ar1` (φ = 0.8): nulls;
  - `sinusoid` (period 7.3), `two_tone`, `logistic_p4` (r = 3.5): periodic /
    quasi-periodic controls. `two_tone` is new:
    x_k = sin(2π k / 7.3 + φ1) + 0.6 sin(2π k φ_g / 7.3 + φ2), with φ_g the
    golden ratio and phases ~ U(0, 2π). It uses system code 9 under the Phase 2E
    SeedSequence scheme;
  - `logistic` (r = 4), `henon`, `skewed_henon` (PRE Eq. 31): chaotic;
  - `logistic` and `henon` plus white Gaussian observation noise at 30, 20 and
    10 dB.

  That is 14 conditions.
- **Window lengths:** 256 (primary) and 512.
- **Test seeds: 2000–2149 (150 per condition and length)**, i.e. 4,200 windows,
  with `SeedSequence([20260926, system_code, window_length, seed])`. They start
  at 2000, not 1000, so that no window was seen in any earlier phase (Phase 3
  used 1000–1099).
- **Development seeds** (Phase 2E seeds: nulls 0–99 at 256 and 0–49 at 512,
  all other conditions 0–29) were used for all exploration
  (`results/dev/explore.jsonl`, `results/dev/dev512.jsonl`).

## 2. Methods (all parameters fixed)

Every method is a decision rule applied to one saved detector RUN
(`detector.run`, which calls the pipeline's own functions). Common settings:
- period 1 only (production `so_periods = (1,)`);
- τ = 1 (`one_sample` map);
- So transform with `prl_norm`, R = 100, κ = 3, 30 histogram bins, 10 % tube;
- **50 surrogates** (production `so_surrogate_count`, AAFT, J < 0.05 per peak);
- Level-C verification with production gates;
- `random_seed` = CFG default.

A window is a **detection** (for candidates, a "verified detection") when at
least one period-1 peak is Level B (surrogate-significant) **and** passes the
method's gate.

- method: `baseline`
  - The production detector, unchanged: `fp.run_upo_analysis` (Cao at lag 1,
    M = 7 Jacobian neighbours). Detection ⇔ ≥ 1 Level-B peak.
  - Reported, **not eligible** to win.
- method: `c1_prod_extgate`
  - C1: the same production run, plus the **P1-a gate on the candidate-centred
    extension monodromy**: leading multiplier modulus ≥ 1.2 (δ = 0.2).
- method: `c2_m2M15_mediangate`
  - C2: fixed lag-1 embedding dimension m = 2 (P3-b), M = 15 Jacobian
    neighbours (P3-a, larger neighbourhood);
  - **P1-a gate on the hybrid PRL/PRE source-stability estimate with the
    element-wise MEDIAN** of the peak's member Jacobians: leading modulus ≥ 1.2.
- method: `c3_m2M7_trimgate`
  - C3: fixed m = 2, production M = 7;
  - **P1-a gate on the hybrid estimate with the 10 % TRIMMED MEAN** (per
    element) of the member Jacobians: leading modulus ≥ 1.2.

Code: `methods.py` and `detector.py` at the commit that adds this file. The
runner records their SHA-256 in `results/test/manifest.jsonl`. Neither file is
edited between that commit and the end of the test evaluation. No candidate
uses the arithmetic-mean hybrid estimate. That estimate is recorded only as a
diagnostic (`src_mean`).

## 3. PRIMARY decision rule (fixed in advance)

All on **test seeds at 256 samples**, N = 150 windows per condition.

1. **PASS.** A candidate passes only if its detection count on EACH of
   `white_noise`, `ar1`, `sinusoid` and `two_tone` is **at most 10/150**
   (= floor(0.07 × 150); 7 %).
2. **WINNER.** Among passing candidates, the winner has the highest pooled
   verified-detection count on `henon@30dB` plus `henon@20dB` (k out of 300).
3. **Tie-breaks,** in order:
   - (a) the lower median localization error on clean `henon`: per window with
     a detection, |nearest detected peak location − 0.63135 (analytic
     on-attractor fixed point)|, then the median over those windows; a
     candidate with no clean-Hénon detection ranks last;
   - (b) the lower pooled false-positive count over the four PASS systems;
   - (c) the lower candidate index (C1 < C2 < C3).
4. If no candidate passes, there is no winner, and this is reported plainly.
5. **Adoption.** The winner, if any, is added to `final_pipeline.py` as an opt-in
   `PipelineConfig` option, default off, which reproduces the method exactly
   (unit tests). `experiments/phase3_lle/groundrule_check.sh` must pass.

The rule is implemented in `decide.py`, committed before this file (dry run on
development seeds only: `results/dev/tables/explore_primary.md`).

## 4. Secondary outcomes (every method, including the baseline)

- Detection counts at 256 and 512 for every condition, with Wilson intervals:
  logistic r = 3.5 rejection, 30 / 20 / 10 dB detection, skewed Hénon, and
  the 512 results.
- Level-C verified rate: detected and `verified_unstable`.
- Localization error (median; nearest detected peak to the analytic fixed
  point) for logistic (0.75), Hénon and skewed Hénon (0.93019).
- Logistic r = 3.5 status distribution.
- **Combined detector** on the same windows, reported both ways:
  `fp.lle_chaos_test` detected AND the method's detection, and
  `lle_chaos_test` detected OR the method's detection. `lle_chaos_test` alone is
  also reported.
- Runtime per run.

The proposal's "bin-mean vs KDE mode" localization comparison is **not**
included: the saved runs record the pipeline's bin-mean locations only. This is
stated as a deviation. The `cao_e1_undefined` diagnostic for logistic r = 3.5
is likewise not recomputed. The baseline's status (`embedding_not_saturated`)
is reported instead, and Phase 3 A5 already established that all logistic
r = 3.5 windows of that status are duplicate-vector cases.

## 5. Development evidence (DEVELOPMENT seeds only)

256 samples, 50 surrogates (`results/dev/tables/explore_gates_256.txt`,
`explore_summary.md`, `explore_primary.md`); detections per condition
(nulls /100, others /30):

| method | WN | AR(1) | sinusoid | two-tone | Hénon | Hénon 30 | Hénon 20 | clean Hénon loc. error |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| baseline | 7 | 11 | 21 | 15 | 29 | 11 | 8 | 0.068 |
| C1 | 0 | 2 | 0 | 0 | 29 | 11 | 7 | 0.068 |
| C2 | 0 | 1 | 0 | 0 | 30 | 30 | 30 | 0.009 |
| C3 | 0 | 3 | 0 | 0 | 30 | 30 | 27 | 0.010 |

512 samples (`results/dev/tables/dev512_summary.md`; nulls /50):

| method | WN | AR(1) | sinusoid | two-tone | Hénon 30 | Hénon 20 |
|---|---:|---:|---:|---:|---:|---:|
| baseline | 5 | 2 | 19 | 21 | 18 | 8 |
| C1 | 1 | 0 | 0 | 0 | 18 | 7 |
| C2 | 0 | 0 | 0 | 0 | 30 | 30 |
| C3 | 2 | 2 | 0 | **10** | 30 | 28 |

Why these candidates:
- The exploration covered 5 detector runs × 13 gates (extension modulus, median-
  and trimmed-mean hybrid, each at δ ∈ {0.1, 0.2, 0.3, 0.5}; Level-C verified).
- The P1-a gate at δ = 0.2 removed every sinusoid and two-tone detection of the
  production run.
- Fixing m = 2 restored noisy-Hénon power and cut the localization error from
  0.068 to 0.010.
- C1 isolates the gate alone. C2 and C3 differ in Jacobian neighbourhood
  (15 vs 7) and in the robust aggregate (median vs trimmed mean).
- Fixed m = 3 raised two-tone detections (22/30 unfiltered). Cao m with M = 15
  did not change noisy-Hénon power.
- Known from development before the test run: C3 flags two-tone in 10/30
  windows at 512 samples. This affects only secondary outcomes.

## 6. Runtime budget (measured on development runs, 4 workers on 4 cores)

Mean seconds per window, under the same 4-worker load:

| run | 256 | 512 |
|---|---:|---:|
| production (`mcao_M7_S50`: baseline, C1) | 4.95 | 8.67 |
| fixed m 2, M 15 (`m2_M15_S50`: C2) | 2.95 | 5.10 |
| fixed m 2, M 7 (`m2_M7_S50`: C3) | 2.94 | 5.08 |
| `lle_chaos_test` (combined detector) | 1.16 | 2.78 |
| **total per window** | **12.0** | **21.6** |

- One test seed covers 14 conditions × (12.0 + 21.6) s ≈ 470 CPU-seconds.
- 150 seeds ≈ 70,560 CPU-seconds ≈ 17,640 s wall time on 4 workers
  ≈ **4.9 hours**, about 60 % headroom under the 8-hour limit.
- 200 seeds (≈ 6.5 h) would leave too little headroom for overruns.
- More surrogates would raise the cost roughly in proportion. 50 is the
  production count (nominal size of J(W) < 0.05: 3/51 per window,
  `docs/PHASE2E_SYNTHETIC_VALIDATION.md` 2.7), so the baseline stays the
  unchanged production detector.

## 7. Amendments

(none)
