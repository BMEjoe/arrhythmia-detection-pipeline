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
| A1 Tables B/C aggregation fix | done | (this commit) |
| A2 requirements.txt, pytest on 3.13 | pending | |
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
- Pytest under Python 3.13 (numpy 2.1.3, scipy 1.18.1): 334 passed, 1 failed
  (`test_hybrid_period1_stability_matches_analytic_henon_multipliers[1-prl_norm]`),
  same as under 3.11. To be diagnosed in A2.
