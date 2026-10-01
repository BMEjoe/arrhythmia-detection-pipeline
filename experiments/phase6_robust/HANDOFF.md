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
| Part A crash-fix check | done: 3 crashed windows no longer error (all AND-detected); N1-N6 0/1800 decision differences | f574f98 |
| Part B generators (systems.py, config.py, README.md) | written, smoke-tested | 7026a12 |
| final_pipeline.py rr_detrend option (none/linear/moving_median/smoothness_priors) + tests/test_phase6_detrend.py | done; groundrule_check passes | f574f98 |
| C2 dev tuning (8 configs; results/dev/tables/tuning.md); gate option rr_detrend_min_trend_sd; combined_chaos_detected fix (UPO failure status -> not detected, was ValueError; found on edited bigeminy = constant_data) | done; groundrule_check passes (393+1; 394/394; 44/44; 13/44) | 7b4f1c2 |
| Candidates D1 linear gate 0.5, D2 linear gate 0.7, D3 smoothness priors lambda 300 gate 0.7 (methods.PREREGISTERED); runner exact-sharing shortcut verified 39/39 dev windows (results/dev/shortcut_check.json) | done | 7b4f1c2 |
| Dev runtime (seeds 0-4, 37 conditions, 4 methods): 201 s / 185 windows | done | 7b4f1c2 |
| C3 PREREGISTRATION.md (+ amendment: window count 11,700) | done, pushed before any test window | 3325e16, 97a302b |
| C4 test run + decide | done: 11,700 windows (resumed after a container restart at 5,990 and after a JSON inf crash at 10,972, amendment 2); WINNER d2_linear_g07 (P3 571/600; D1 fails G1 279 vs 294; D3 passes 539) | 37b1391 |
| C5 adopt D2: combined_chaos_config(config=None, detrend=False), PHASE6_DETREND; tests/test_phase6_combined_detrend.py; adoption_check 111/111; groundrule_check passes (400+1; 401/401; 44/44; 13/44) | done | (this commit) |
| C1 diagnosis (diagnose_trend.py, dev 0-29) | done: So-mode drifts by ~1.1x the trend change; peak height 0.83 -> 0.18 (Henon); linear detrend restores Level B 30/30 | 7026a12 |
| runner / analysis / decide / methods | written | 7026a12 |

## Notes
- PREREGISTRATION.md pushed: do NOT edit methods.py, systems.py, config.py,
  decide.py, analysis.py, run_phase6.py or final_pipeline.py until the test
  evaluation is done (amendments only).
- NEXT: test run `PYTHONDONTWRITEBYTECODE=1 nohup /root/venv313/bin/python -m experiments.phase6_robust.run_phase6 --phase test --workers 4 > /tmp/claude-0/p6_test.log 2>&1 &`
  (resumable; ~3.6 h; 13,200 windows), then analysis/decide --phase test --file test.
- Test run started after 97a302b. Container restart at 5,990/11,700 windows
  (all lines intact); resumed with the same command (log /tmp/claude-0/p6_test2.log).
- NEXT: C5 adopt D2: combined_chaos_config(config=None, detrend=False); detrend=True
  -> rr_detrend='linear', rr_detrend_min_trend_sd=0.7; tested config =
  combined_chaos_config(replace(CFG, keep_upo_on_short_lle_embedding=True), detrend=True);
  tests + adoption check + groundrule_check.sh; then docs/PHASE6_ROBUSTNESS.md.
| D docs/PHASE6_ROBUSTNESS.md | done | (this commit) |
- ALL PHASE 6 STEPS DONE. Nothing left to resume.
