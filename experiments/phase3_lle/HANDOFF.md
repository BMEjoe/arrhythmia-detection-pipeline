# Phase 3 handoff (living document; update after every step)

Read this first if you are resuming. Branch: `claude/amazing-allen-p4ssx3`.
Never merge into main, never open a pull request. Commit and push after every step.

## Ground rules
- `final_pipeline.py`: only new opt-in `PipelineConfig` options whose defaults
  reproduce current behaviour exactly. After any change: full `pytest tests`
  must pass (except the pre-existing stability failure documented under A2), and
  `python -m experiments.phase2e.run_phase2e --replicate` must stay identical.
  If either fails, revert and report it.
- Seed split for Part B: DEVELOPMENT = Phase 2E core seeds; TEST = new seeds
  (not looked at until B4). PREREGISTRATION.md must be committed and pushed
  before any test-seed run.

## Environment
- Python 3.13 venv at `/root/venv313` (created with
  `uv venv -p /usr/bin/python3.13 /root/venv313`, then
  `uv pip install -p /root/venv313/bin/python -r requirements.txt`).
- Run everything with `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1`.

## Status
| Step | State | Commit |
|---|---|---|
| A1 Tables B/C aggregation fix | done | 5e22d04 |
| A2 requirements.txt, pytest on 3.13 | done | 88f344e |
| A3 posthoc lle.json | done | b415607 |
| A4 opt-in keep UPO on short LLE embedding | done | 4aac419 |
| A5 cao_e1_undefined diagnostic | done | 7f53a25 |
| B1 baseline (dev seeds) | done | fd3c751 |
| B2 diagnosis (dev seeds) | done | 655233a |
| B3 candidates + PREREGISTRATION.md | done (preregistration pushed before any test run) | 25ec207 |
| B4 test-seed evaluation | running (started after 25ec207; decide.py at 7c4a496) | |
| B5 adopt winner (opt-in) | pending | |
| B6 docs/PHASE3_LLE_VALIDATION.md | draft sections 1-4 committed; 5+ after B4 | |
| C docs/PHASE4_UPO_PROPOSAL.md | done | a14da95 |

## Notes
- A1: rerunning Phase 2E analysis under a newer library stack also rewrites
  plots, Table E alignment and 1-ulp float digits in two CSVs; those were
  restored (not content changes). Only Tables B/C are committed.
- A2: pytest under Python 3.13: 334 passed, 1 failed (stability test [1-prl_norm]).
  Cause: one ill-conditioned member Jacobian (norm 4,078) dominates the
  arithmetic mean in source_period1_stability. The value depends on the CPU
  floating-point path: with
  `NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL AVX512_SPR"`
  the suite passes 335/335. Default env: 334/1 (pre-existing, not caused by us).
- Replicability: `run_phase2e --replicate` gives 13/44 in this container even
  with the unmodified pipeline (different machine; 1e-15 differences amplified by
  ill-conditioned quantities). Therefore the ground-rule check for
  final_pipeline.py changes is done as:
  1. before the change: `python -m experiments.phase3_lle.replicate_check --save /tmp/claude-0/rep_before.jsonl`
     (regenerate it from the pre-change commit if /tmp was lost);
  2. after: `--compare` must give 44/44, AND `run_phase2e --replicate` must still
     give 13/44 with a replicability.json identical to the pre-change one
     (saved at /tmp/claude-0/replicability_env_before.json); then
     `git checkout experiments/phase2e/results/replicability.json` (keep the
     original-machine file committed).
  3. pytest: default env must give exactly the same 1 failure; with AVX-512
     disabled it must give 335/335 (plus any new tests).
- A3: lle.json regenerated via posthoc_diagnostics.lle_diagnostics() only (198 s).
  Production LLE path is 3.8-6.4x (logistic) / 2.3-3.2x (Henon) too low;
  tau=1 + true m gives 0.69 / 0.40. White noise at tau=1, m=2 gives 0.42-0.58,
  so a surrogate/null test is needed, not only bias correction (input to B2/B3).
- A4: `PipelineConfig.keep_upo_on_short_lle_embedding=False`. Tests in
  tests/test_lle_short_embedding.py (3). Ground rules checked with
  experiments/phase3_lle/groundrule_check.sh:
  default pytest 337 passed + 1 pre-existing failure; AVX-512 off 338/338;
  same-env replicability 44/44; run_phase2e --replicate 13/44, file identical.
- A5: `cao_method(..., return_diagnostics=False)` and
  `PipelineConfig.upo_report_cao_diagnostics=False`. Tests in
  tests/test_cao_e1_diagnostics.py (7). Ground rules: pytest 344 + 1 pre-existing
  failure; AVX-512 off 345/345; same-env 44/44; --replicate 13/44 identical.

## Part B notes
- Runner: `python -m experiments.phase3_lle.run_phase3 --phase dev --methods baseline --workers 4`
  (resumes; log in /tmp/claude-0/b1.log). Test phase is refused until
  PREREGISTRATION.md is committed + pushed and lists the method.
- B2 (results/dev/diagnosis.jsonl, results/dev/tables/B2_diagnosis.md,
  plots/B2_divergence_curves_N*.png). Effect sizes at 256 (median bias):
  * delay: V0 (TDMI tau*, m*) -> V1 (tau 1, m*): logistic -0.571 -> -0.400,
    Henon -0.301 -> -0.141.
  * dimension: V1 -> V2 (tau 1, Cao lag-1 m ~4): -0.022 / -0.046; V3 (true m): ~0.
    Under noise Cao lag-1 m rises to 7-9 and V2 bias returns (logistic 20 dB -0.44).
  * fit region: production curves (V0) saturate by k95 = 2 in 30/30 chaotic
    windows, i.e. the 0..5 fit spans the plateau; at tau 1 / true m saturation is
    k95 = 9 (logistic), 12 (Henon), so 0..5 is inside the linear region.
  * noise floor: at tau 1 the k=0->1 jump is 2.0 (WN), 1.7 (AR1), 3.9 (logistic
    20 dB) vs 0.69 clean; it inflates the 0..5 slope (WN 0.42, logistic 20 dB 0.90).
    Fitting from k=1 removes it (WN 0.09, AR1 0.16, logistic 20 dB 0.31,
    Henon 20 dB 0.32; clean unchanged 0.696 / 0.414).
- B3 tuning: experiments/phase3_lle/tune.py (dev, 256, 39 surrogates), writes
  results/dev/tuning.jsonl; `--summary` prints rates per grid combination.
  Vectorized divergence_curve in estimators.py reproduces fp.rosenstein_lle's
  curve to 7e-16.
- B3 tuning result (results/dev/tuning.jsonl; `python -m experiments.phase3_lle.tune --summary`):
  tau 1 + fit from k=1 gives dev power 1.00 at clean/30/20/10 dB for most LLE
  combos; sinusoid and logistic_p4 never flagged; IAAFT keeps AR(1) closer to
  nominal than AAFT; 0-1 test has power <= 0.07. Chosen candidates (estimators.py,
  99 IAAFT surrogates, window-specific RNG): c1_rosenstein_m2_iaaft,
  c2_kantz_m3_sat_iaaft, c3_eps_m2_iaaft, c4_zero_one_iaaft.
  Dev run: `run_phase3 --phase dev --methods c1_...,c2_...,c3_...,c4_... --workers 3`
  (log /tmp/claude-0/b3dev.log). NEXT: analysis --phase dev, write
  PREREGISTRATION.md, commit + push it, only then `--phase test --methods all`.
- B1 (results/dev/baseline.jsonl, results/dev/tables/B1_baseline.md): existing
  Rosenstein + 199 AAFT, p <= 0.05. Dev 256: FPR WN 7/100, AR1 4/100, sinusoid
  0/30, logistic_p4 0/30 (LLE undefined); power logistic 14/30, Henon 12/30,
  20 dB 19/30 and 14/30; bias -0.571 / -0.301. 512: power 30/30 and 27/30.
  Median runtime 14 s (256) / 32 s (512) under 4-way contention.
- B3 done: PREREGISTRATION.md committed and pushed. From here on do NOT edit
  estimators.py or the rule; changes only as dated amendments in the file.
  NEXT (B4): `python -m experiments.phase3_lle.run_phase3 --phase test --methods all --workers 4`
  (resumable; log /tmp/claude-0/b4.log), then `analysis --phase test`.
