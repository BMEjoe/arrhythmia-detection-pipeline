# Phase 6 handoff (living document; update after every step)

Read this first when resuming. Branch `claude/amazing-allen-p4ssx3`, restarted
from main at `ea7e068` (PR #2 merged; main has Phases 2E-5). Never merge into
main, never open a pull request. Commit and push after every step.

## Task (user, Phase 6) in brief
BASELINE-K = combined_chaos_config() + keep_upo_on_short_lle_embedding=True,
analyze_segment on raw RR. Phase 5 generators, seed scheme, 256 intervals,
1/360 s quantization; dev seeds as Phase 5; TEST seeds from 4000. No MIT-BIH.
- Part A (report): BASELINE-K on the 3 Phase 5 crashed windows (P1 3008, 3154;
  P4 henon 30 dB 3014) + N1-N6 seeds 3000-3299; AND counts must equal Phase 5.
- Part B (report, no rule): E1 bigeminy, E2 trigeminy, E3 couplets 5/10 %,
  E4 runs (1-3 runs of 3-6), E5 atrial (no compensatory pause) 5/10 %,
  E6 S2-10 % + N3 trend; annotation-edited variants of S2 5/10 %, E1, E2,
  E3 10 % (linear interpolation of ectopic + compensatory intervals). AND,
  components, OR with Wilson CIs. Purpose: raw RR vs edited NN for MIT-BIH.
- Part C: C1 dev diagnosis of trend mechanism; C2 <= 3 detrending candidates
  as opt-in PipelineConfig option (default "none"), tuned on dev only;
  C3 PREREGISTRATION.md pushed before any test seed (>= 300 seeds N1-N6,
  >= 200 others, <= ~8 h on 4 workers); rule: PASS iff AND <= floor(0.07N) on
  each N1-N6 and each S2 rate AND power on P1, P2, P4 (4), G1 no more than
  6/200 below BASELINE-K; WINNER max pooled P3 AND; tie fewer pooled Part B
  ectopy AND, then index. C4 test; C5 adopt winner opt-in + extend
  combined_chaos_config with an optional argument (default unchanged) + tests.
- Part D: docs/PHASE6_ROBUSTNESS.md incl. "Recommended MIT-BIH analysis
  configuration" (raw RR vs edited NN; detrending), limitations.

## Ground rules (as Phases 3-5)
- final_pipeline.py: opt-in options only, defaults unchanged; after any change
  `bash experiments/phase3_lle/groundrule_check.sh` (default pytest: all pass
  except the pre-existing stability failure; AVX-512 off: all pass; same-env
  44/44 vs /tmp/claude-0/rep_before.jsonl; --replicate 13/44 identical to
  /tmp/claude-0/replicability_env_before.json).
- /root/venv313, OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1. Long runs die when
  the container idles between turns: babysit in-turn; runners are resumable.
- physionet.org not reachable (proxy 403).

## Status
| Step | State | Commit |
|---|---|---|
| Part A crash-fix check | done: 3 crashed windows no longer error (all AND-detected); N1-N6 0/1800 decision differences | (this commit) |
| Part B generators (systems.py, config.py, README.md) | written, smoke-tested | 7026a12 |
| final_pipeline.py rr_detrend option (none/linear/moving_median/smoothness_priors) + tests/test_phase6_detrend.py | done; groundrule_check passes (391+1 pre-existing; 392/392 AVX-512 off; 44/44; 13/44 identical) | (this commit) |
| C1 diagnosis (diagnose_trend.py, dev 0-29) | done: So-mode drifts by ~1.1x the trend change; peak height 0.83 -> 0.18 (Henon); linear detrend restores Level B 30/30 | 7026a12 |
| runner / analysis / decide / methods | written | 7026a12 |

## Notes
- NEXT: C2 dev tuning: run_phase6 --phase dev --tag tuning with TUNING methods
  on P3, P1, P2, P4, G1 (dev 0-29) and N1-N6, S2 (dev 0-49); choose <= 3
  candidates -> methods.PREREGISTERED; Part B dev baseline run; runtime budget;
  PREREGISTRATION.md.
