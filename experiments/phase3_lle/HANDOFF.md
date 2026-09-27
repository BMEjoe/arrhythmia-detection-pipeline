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
| A2 requirements.txt, pytest on 3.13 | done | (this commit) |
| A3 posthoc lle.json | pending | |
| A4 opt-in keep UPO on short LLE embedding | pending | |
| A5 cao_e1_undefined diagnostic | pending | |
| B1 baseline (dev seeds) | pending | |
| B2 diagnosis (dev seeds) | pending | |
| B3 candidates + PREREGISTRATION.md | pending | |
| B4 test-seed evaluation | pending | |
| B5 adopt winner (opt-in) | pending | |
| B6 docs/PHASE3_LLE_VALIDATION.md | pending | |
| C docs/PHASE4_UPO_PROPOSAL.md | pending | |

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
