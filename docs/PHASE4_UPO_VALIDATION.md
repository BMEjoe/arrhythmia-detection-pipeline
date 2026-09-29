# Phase 4: UPO Detector Specificity and Noise Robustness

Research-development record in the style of the Phase 2E and Phase 3 reports.
Every statement is tied to a results file, table, preregistration or test under
`experiments/phase4_upo/` (paths are relative to it unless they start with
`docs/`, `tests/` or `final_pipeline.py`). Synthetic data only; no MIT-BIH or
real ECG data. The design is `docs/PHASE4_UPO_PROPOSAL.md` Section 4, with the
user's changes taking precedence (see `HANDOFF.md`).

## 1. Provenance and design

| Item | Value |
|---|---|
| Branch | `claude/amazing-allen-p4ssx3`, restarted from `main` at `858d616` |
| Environment | Python 3.13.12, numpy 2.1.3 (`requirements.txt`), BLAS 1 thread; manifests `results/{dev,test}/manifest.jsonl` |
| Systems (14 conditions) | white noise, AR(1) φ = 0.8 (nulls); sinusoid (period 7.3), **two-tone** quasi-periodic (new: sin(2πk/7.3 + φ1) + 0.6 sin(2πkφ_g/7.3 + φ2), golden-ratio frequency ratio, system code 9), logistic r = 3.5 (periodic controls); logistic r = 4, Hénon, skewed Hénon (chaotic); logistic and Hénon at 30 / 20 / 10 dB (`systems.py`) |
| Windows | 256 (primary) and 512 |
| Development seeds | Phase 2E seeds: nulls 0–99 (256) / 0–49 (512), other conditions 0–29 |
| Test seeds | **2000–2149** (150 per condition and length; 4,200 windows). They start at 2000 so that no window was seen in an earlier phase (Phase 3 used 1000–1099) |
| Seed discipline | `PREREGISTRATION.md` committed and pushed in `a6fe566` before any test window. `decide.py` was committed earlier, in `4ea3b98`. `run_phase4.py --phase test` reuses the Phase 3 guard |
| Detector settings | Period 1 only (production `so_periods = (1,)`); 50 AAFT surrogates (production); Level-C verification with production gates |

A **detection** is a window with at least one period-1 peak that is Level B
(surrogate J < 0.05) and passes the method's gate. For the baseline there is no
gate. Rates are detections / all windows, with 95 % Wilson intervals.

Each window's detector runs are saved once (`results/test/test.jsonl`, with
per-peak features). Every method is a fixed rule on those saved runs
(`methods.py`, `detector.py`). `fp.lle_chaos_test` (the Phase 3 winner) is
computed on the same window for the combined detector.

## 2. Development work (development seeds only)

- **Exploration** (`results/dev/explore.jsonl`,
  `results/dev/tables/explore_gates_256.txt`): 560 windows at 256 samples ×
  5 detector runs:
  - production (Cao m, M = 7);
  - fixed m = 2 with M = 7 and M = 15;
  - fixed m = 3;
  - Cao m with M = 15.

  Each run was evaluated under 13 gates: extension-monodromy, median-hybrid and
  trimmed-mean-hybrid leading modulus ≥ 1 + δ for δ ∈ {0.1, 0.2, 0.3, 0.5},
  and Level-C verified. The findings:
  - on the production run, the P1-a gate at δ = 0.2 removes every sinusoid
    and two-tone detection (21/30 → 0/30 and 15/30 → 0/30);
  - fixing m = 2 restores noisy-Hénon power (30 + 20 dB: 19/60 → 57–60/60) and
    cuts the clean-Hénon localization error from 0.068 to 0.010;
  - fixed m = 3 raises two-tone detections (22/30 before gating);
  - M = 15 at Cao m did not raise noisy-Hénon power (18/60).
- **512 samples** (`results/dev/dev512.jsonl`,
  `results/dev/tables/dev512_summary.md`): the three shortlisted runs. Known
  before the test run: C3 flags two-tone in 10/30 windows at 512.
- **Hybrid stability aggregates.** Following the user's rule (Phase 3 A2),
  every hybrid-stability gate uses the element-wise median or the 10 % trimmed
  mean of the member Jacobians, never the arithmetic mean.

## 3. Candidates and preregistered rule (`PREREGISTRATION.md`)

| id | method | run | gate (P1-a, δ = 0.2) |
|---|---|---|---|
| — | `baseline` (not eligible) | production: Cao m at lag 1, M = 7 | none (Level B) |
| C1 | `c1_prod_extgate` | production | candidate-centred extension monodromy leading modulus ≥ 1.2 |
| C2 | `c2_m2M15_mediangate` | fixed m = 2, M = 15 | hybrid source stability, **element-wise median**, ≥ 1.2 |
| C3 | `c3_m2M7_trimgate` | fixed m = 2, M = 7 | hybrid source stability, **10 % trimmed mean**, ≥ 1.2 |

**Primary rule** (test seeds, 256 samples, N = 150):
- **PASS** only if detections ≤ **10/150** (7 %) on EACH of white noise,
  AR(1), sinusoid and two-tone.
- **WINNER** = most pooled verified detections on Hénon at 30 dB + 20 dB.
- **Tie-breaks:** the median clean-Hénon localization error, then fewer pooled
  false positives, then candidate index.

**Budget** (Section 6 of the preregistration), measured on development runs
under 4-worker load:
- all runs plus `lle_chaos_test` take 12.0 s per 256-sample window and 21.6 s
  per 512-sample window;
- 150 seeds × 14 conditions ≈ 70,600 CPU-seconds ≈ 4.9 hours on 4 workers
  (limit about 8 hours);
- 50 surrogates is the production count.

## 4. Test-seed results

Sources: `results/test/test.jsonl` (4,200 unique windows; 0 duplicates or torn
records), `results/test/tables/test_summary.{md,csv}` and
`results/test/tables/test_primary.md` (`decide.py`).

### 4.1 Primary decision (256 samples)

| method | white noise | AR(1) | sinusoid | two-tone | PASS | Hénon 30 dB | Hénon 20 dB | pooled (of 300) | clean Hénon loc. error |
|---|---|---|---|---|---|---|---|---|---|
| baseline | 12/150 [0.046, 0.135] | 7/150 [0.023, 0.093] | **89/150** [0.513, 0.669] | **67/150** [0.369, 0.527] | no (not eligible) | 70/150 | 40/150 | 110 | 0.0679 |
| C1 | 0/150 [0, 0.025] | 3/150 [0.007, 0.057] | 0/150 | 0/150 | yes | 65/150 | 34/150 | 99 | 0.0679 |
| **C2** | 0/150 [0, 0.025] | 0/150 [0, 0.025] | 0/150 | 1/150 [0.001, 0.037] | yes | **150/150** [0.975, 1] | **143/150** [0.907, 0.977] | **293** | 0.0102 |
| C3 | 4/150 [0.010, 0.067] | 3/150 [0.007, 0.057] | 0/150 | 1/150 | yes | 149/150 | 140/150 | 289 | 0.0100 |

**All three candidates pass. Winner: C2** (`c2_m2M15_mediangate`) with 293/300,
against C3 289/300 (a margin of 4 windows, 1.3 percentage points) and C1 99/300.
No tie-break was needed. Against the baseline:
- noisy-Hénon detection rises from 110/300 (0.37) to 293/300 (0.98);
- the baseline's periodic-control false positives (sinusoid 89/150, two-tone
  67/150) fall to 0/150 and 1/150;
- the clean-Hénon localization error falls from 0.068 to 0.010.

### 4.2 Detections per condition (secondary outcomes), out of 150

| method | N | WN | AR1 | sinus. | 2-tone | log r3.5 | log | Hén | skew Hén | log 30 | log 20 | log 10 | Hén 30 | Hén 20 | Hén 10 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| baseline | 256 | 12 | 7 | 89 | 67 | 0 | 150 | 140 | 133 | 128 | 130 | 106 | 70 | 40 | 9 |
| C1 | 256 | 0 | 3 | 0 | 0 | 0 | 150 | 140 | 132 | 59 | 41 | 20 | 65 | 34 | 8 |
| C2 | 256 | 0 | 0 | 0 | 1 | 0 | 150 | 150 | 150 | 150 | 142 | 48 | 150 | 143 | 74 |
| C3 | 256 | 4 | 3 | 0 | 1 | 0 | 150 | 150 | 150 | 150 | 139 | 53 | 149 | 140 | 56 |
| baseline | 512 | 8 | 12 | 99 | 103 | 0 | 150 | 149 | 149 | 139 | 131 | 141 | 75 | 45 | 13 |
| C1 | 512 | 0 | 1 | 0 | 2 | 0 | 150 | 149 | 148 | 97 | 77 | 46 | 73 | 40 | 10 |
| C2 | 512 | 0 | 1 | 0 | 2 | 0 | 150 | 150 | 150 | 150 | 150 | 89 | 150 | 148 | 83 |
| C3 | 512 | 1 | 1 | **4** | **29** | 0 | 150 | 150 | 150 | 150 | 149 | 92 | 150 | 144 | 98 |

- **C1 (gate only)** fixes specificity but loses noisy-logistic power (256,
  20 dB: 130 → 41). It does not help noisy Hénon.
- **C3** is specific at 256. At 512 it flags two-tone 29/150 and sinusoid
  4/150, confirming the development warning; C2 flags 2/150 and 0/150.
- **At 10 dB** every method is weak on Hénon (C2 74/150 and 83/150). Logistic
  at 10 dB is detected best by the baseline (106 and 141). The fixed-m methods
  detect 48–92.
- **Logistic r = 3.5:** 0 detections for every method. Its status is
  `embedding_not_saturated` in 150/150 windows (production run, Cao) and
  `no_valid_transforms` in 150/150 windows (fixed m = 2). With 4 distinct
  values, no local Jacobian can be estimated.

### 4.3 Level-C verified detections (detected and `verified_unstable`), out of 150

| method | N | logistic | Hénon | skewed Hénon | log 30 dB | Hén 30 dB | Hén 20 dB | nulls / periodic |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| baseline, C1 | 256 | 1 | 1 | 3 | 0 | 0 | 1 | 0 |
| C2 | 256 | 149 | 30 | 141 | 79 | 16 | 2 | 0 |
| C3 | 256 | 148 | 31 | 140 | 56 | 18 | 2 | 0 |
| baseline, C1 | 512 | 0 | 0 | 0 | 2 | 1 | 0 | 0 |
| C2 | 512 | 150 | 78 | 150 | 57 | 11 | 0 | 0 |

At production dimension almost nothing passes Level C (Problem 2 of the
proposal). At m = 2, logistic and skewed Hénon verify in 94–100 % of windows,
and Hénon in 20–52 %. Under noise verification stays rare. No null or periodic
window verifies under any method.

### 4.4 Localization error (median, nearest detected peak to the analytic fixed point)

| method | N | logistic (0.75) | Hénon (0.6314) | skewed Hénon (0.9302) | Hénon 30 dB | Hénon 20 dB |
|---|---:|---:|---:|---:|---:|---:|
| baseline | 256 | 0.0077 | 0.0679 | 0.0157 | 0.0592 | 0.1056 |
| C2 | 256 | 0.0009 | 0.0102 | 0.0117 | 0.0268 | 0.0218 |
| C3 | 256 | 0.0008 | 0.0100 | 0.0114 | 0.0273 | 0.0229 |
| baseline | 512 | 0.0071 | 0.0716 | 0.0125 | 0.0780 | 0.1915 |
| C2 | 512 | 0.0008 | 0.0102 | 0.0120 | 0.0279 | 0.0203 |

C1 matches the baseline because it uses the same run. The baseline's Hénon
error (about 0.07) is the bin-edge quantization documented in Phase 2E 3.3.

### 4.5 Combined binary detector with `lle_chaos_test` (same windows)

`lle_chaos_test` alone (Phase 3 winner):
- at 256: 11/150 white noise and 8/150 AR(1); 0 on sinusoid, two-tone and
  logistic r = 3.5; 150/150 on every chaotic condition including 10 dB;
- at 512: 4/150 and 4/150.

| method | rule | N | WN | AR1 | sinus. | 2-tone | log 10 dB | Hén 30 dB | Hén 20 dB | Hén 10 dB |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| baseline | AND | 256 | **0** | **0** | **0** | **0** | 106 | 70 | 40 | 9 |
| baseline | OR | 256 | 23 | 15 | 89 | 67 | 150 | 150 | 150 | 150 |
| C1 | AND | 256 | 0 | 0 | 0 | 0 | 20 | 65 | 34 | 8 |
| C1 | OR | 256 | 11 | 11 | 0 | 0 | 150 | 150 | 150 | 150 |
| **C2** | **AND** | 256 | **0** | **0** | **0** | **0** | 48 | 150 | 143 | 74 |
| C2 | OR | 256 | 11 | 8 | 0 | 1 | 150 | 150 | 150 | 150 |
| C3 | AND | 256 | 0 | 0 | 0 | 0 | 53 | 149 | 140 | 56 |
| C3 | OR | 256 | 15 | 11 | 0 | 1 | 150 | 150 | 150 | 150 |
| C2 | AND | 512 | 0 | 0 | 0 | 0 | 89 | 150 | 148 | 83 |
| C2 | OR | 512 | 4 | 5 | 0 | 2 | 150 | 150 | 150 | 150 |

- **AND:** 0 false positives on all four controls for every method at both
  lengths (also 0 on logistic r = 3.5). Chaotic detection equals the UPO
  method's own, because `lle_chaos_test` detects every chaotic window. C2 AND
  keeps 150/150 and 143/150 at 30 / 20 dB.
- **OR:** detects every chaotic window, but its null false positives are
  essentially those of `lle_chaos_test` (C2: 11/150 and 8/150 at 256, i.e.
  7.3 % and 5.3 %). For the baseline OR it adds the periodic false positives.

### 4.6 Runtime (median s per run, clean Hénon windows, 4 workers)

| run | 256 | 512 |
|---|---:|---:|
| production (baseline, C1) | 2.2 | 5.0 |
| m 2, M 15 (C2) | 1.8 | 4.1 |
| m 2, M 7 (C3) | 1.7 | 4.0 |

The final segment of the test run (4,096 windows after two resumes; Section 7)
took 11,773 s wall time.

## 5. Sensitivity of the Phase 3 winner to m (report only)

Source: `m_sensitivity.py`, `results/m_sensitivity.jsonl`,
`results/tables/m_sensitivity.md`. `fp.lle_chaos_test` was run with
`lle_chaos_test_m` = 2, 3, 4 on the **Phase 3 test seeds 1000–1099 at 256**.
The m = 2 arm reproduces the stored Phase 3 C1 results in 700/700 windows.

| condition | m = 2 | m = 3 | m = 4 | median LLE (bias) m = 2 / 3 / 4 |
|---|---|---|---|---|
| white noise | 4/100 [0.016, 0.098] | 7/100 [0.034, 0.137] | 9/100 [0.048, 0.162] | 0.122 / 0.182 / 0.195 |
| AR(1) | 3/100 | 1/100 | 1/100 | 0.204 / 0.242 / 0.242 |
| sinusoid | 0/100 | 0/100 | 0/100 | 0.019 / 0.038 / 0.041 |
| logistic | 100/100 | 100/100 | 100/100 | 0.692 (−0.001) / 0.684 (−0.009) / 0.659 (−0.034) |
| Hénon | 100/100 | 100/100 | 100/100 | 0.415 (−0.005) / 0.415 (−0.005) / 0.412 (−0.008) |
| logistic 20 dB | 100/100 | 100/100 | 100/100 | 0.481 (−0.213) / 0.491 (−0.202) / 0.462 (−0.231) |
| Hénon 20 dB | 100/100 | 100/100 | 100/100 | 0.352 (−0.068) / 0.346 (−0.073) / 0.320 (−0.099) |

- **Detection does not depend on the fixed m = 2.** Chaotic power is 100/100
  and sinusoid 0/100 at every m.
- **The false-positive rate shifts with m:** white noise 4 → 7 → 9 per 100 and
  AR(1) 3 → 1 → 1. All intervals contain 0.05, but the white-noise trend is
  upward.
- Clean bias grows at m = 4 (logistic −0.034).
- Nothing was re-selected. m = 2 remains the preregistered Phase 3 choice.

## 6. Adoption of the winner

C2 is in `final_pipeline.py` as opt-in options, **all default off**:
- `upo_fixed_dimension` (None = Cao, production);
- `upo_instability_gate`, `upo_instability_gate_delta = 0.2`,
  `upo_instability_gate_aggregate = "median"` (or `"trim10"`);
- new functions `robust_source_stability` and `phase4_upo_config()`, which
  returns m = 2, `so_jacobian_neighbors = 15`, significance on, median gate at
  δ = 0.2.

With the gate on, `run_upo_analysis` adds `instability_gated_uop_candidates` and
`instability_gate_detected`. Levels A, B, C and every existing field are
unchanged. Checks:
- `tests/test_phase4_upo_gate.py` (10 tests): defaults off and no new fields
  by default; the preregistered values; robustness of the median / trimmed-mean
  aggregate to one ill-conditioned member (the arithmetic mean is dominated);
  bitwise agreement with `experiments/phase4_upo` C2 on five signal types;
- recomputing 140 stored test-seed C2 results (5 per condition and length)
  through `fp.run_upo_analysis(x, fp.phase4_upo_config())`: 0 mismatches in
  decision and peak locations;
- ground rules (`experiments/phase3_lle/groundrule_check.sh`):
  - pytest 377 passed plus the pre-existing stability failure (378/378 with
    numpy AVX-512 dispatch disabled);
  - same-machine replicability 44/44 bitwise identical;
  - `run_phase2e --replicate` unchanged (13/44 on this machine, identical file).

## 7. Limitations

- **Synthetic maps are not RR intervals.** The fixed m = 2 suits the one- and
  two-dimensional maps tested here and was chosen on development data from such
  maps. RR-interval dynamics may need a different dimension. The fixed-m choice
  replaces a data-driven estimate (Cao) whose overestimation was the problem,
  but it is not itself validated for physiological data.
- **Primary rule decided on power, with a small margin.** C2 beat C3 by 4
  windows (293 vs 289 of 300). Both use m = 2; they differ in Jacobian
  neighbourhood (15 vs 7) and aggregate (median vs trimmed mean). The data do
  not separate those two factors.
- **Specificity at 512 was not part of the primary rule.** C3 fails the 7 %
  threshold on two-tone at 512 (29/150). C2 stays at 2/150. A candidate chosen
  at 256 is not guaranteed to hold at other lengths.
- **10 dB noise.** C2 detects Hénon at 10 dB in only 74/150 (256) and 83/150
  (512) windows, and noisy logistic at 10 dB less often than the baseline
  (48 vs 106 at 256). The fixed low dimension trades power at high noise on the
  one-dimensional logistic map for specificity and noisy-Hénon power.
- **Level C remains rare under noise** (C2: 16/150 at Hénon 30 dB, 2/150 at
  20 dB). "Verified detection" in the primary rule means Level B plus the
  instability gate, not Level-C verification.
- **Stable periodic orbits with exact duplicates** (logistic r = 3.5) are
  reported as `no_valid_transforms` under fixed m = 2, not as "no UPO". RR data
  quantized at 1/360 s may behave similarly (Phase 2E 3.7).
- **Surrogate null.** Level B still tests against AAFT (linear-Gaussian)
  surrogates. The gate removes neutral (periodic / quasi-periodic) peaks by
  stability, not by a periodic-orbit null. Proposal P1-b / P1-c were not tested.
- **Deviations from the proposal**, stated in the preregistration:
  - no bin-mean vs KDE-mode localization comparison;
  - `cao_e1_undefined` not recomputed (Phase 3 A5 covers logistic r = 3.5);
  - test seeds start at 2000.
- **Combined detector.** AND gave 0 false positives on all controls. With
  150 windows per control, the 95 % upper bound per control is still 2.5 %
  (Wilson).
- **Run interruptions.** Container restarts killed the test run twice (after 63
  and 104 windows). It was resumed with the same preregistered configuration;
  finished windows are skipped. The 4,200 results are unique, with no torn
  records. Each window's result depends only on its data seed, not on run
  order.
- **Machine dependence.** One machine. The robust aggregates were chosen partly
  because the arithmetic-mean estimate is machine-dependent (Phase 2E 5.3).
  Other ill-conditioned steps (local least squares) can still differ at the ulp
  level across machines.

## 8. Reproducibility

```bash
python -m experiments.phase4_upo.run_phase4 --phase dev --runs cao:7:50,2:7:50,2:15:50,3:7:50,cao:15:50 --tag explore --windows 256
python -m experiments.phase4_upo.explore_summary
python -m experiments.phase4_upo.run_phase4 --phase dev --runs cao:7:50,2:15:50,2:7:50 --tag dev512 --windows 512
python -m experiments.phase4_upo.run_phase4 --phase test --methods all        # guarded
python -m experiments.phase4_upo.analysis --phase test --file test
python -m experiments.phase4_upo.decide   --phase test --file test
python -m experiments.phase4_upo.m_sensitivity && python -m experiments.phase4_upo.m_sensitivity --summary
python -m pytest tests
```

Commits:
- `1695cfb` scaffold;
- `55778ac` m-sensitivity;
- `7b8f6e1` exploration;
- `4ea3b98` candidates and `decide.py`;
- `a6fe566` preregistration;
- `c897263` test results;
- `02a58ed` adoption.
