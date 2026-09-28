# Phase 4 handoff (living document; update after every step)

Read this first when resuming. Branch `claude/amazing-allen-p4ssx3` (restarted
from main at 858d616 after PR #1 was merged). Never merge into main, never open
a pull request. Commit and push after every step.

## Ground rules (identical to Phase 3)
- `final_pipeline.py`: only new opt-in `PipelineConfig` options, defaults
  unchanged. After any change run `experiments/phase3_lle/groundrule_check.sh`:
  default pytest = all pass except the pre-existing stability failure;
  AVX-512 disabled = all pass; same-machine replicability 44/44
  (baseline `/tmp/claude-0/rep_before.jsonl`; regenerate from commit 7f53a25
  with `python -m experiments.phase3_lle.replicate_check --save` if lost);
  `run_phase2e --replicate` 13/44 with a file identical to
  `/tmp/claude-0/replicability_env_before.json`.
- Seeds: DEV = Phase 2E seeds; TEST = 2000+ (config.TEST_SEED_START), count
  fixed by the preregistration budget; runner refuses `--phase test` until
  PREREGISTRATION.md is committed, unmodified, pushed and lists the methods.

## Environment
- `/root/venv313` (Python 3.13, `requirements.txt`); run with
  `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1`.

## User's changes to the proposal (take precedence)
1. structure as phase3_lle, runner guard reused;
2. runtime budget measured on dev first; surrogates + test seeds (>= 50 per
   condition) so the test run takes about 8 h or less; calculation in the
   preregistration;
3. <= 3 candidates + unchanged baseline; >= 1 with the P1-a instability gate;
   hybrid stability only with a robust aggregate (median / trimmed mean);
4. PRIMARY: test seeds, 256: PASS iff detection count on EACH of white_noise,
   ar1, sinusoid, two_tone <= 7 % (stated as a count); WINNER = highest pooled
   verified detections on henon@30dB + henon@20dB; tie-break median
   localization error on clean Henon, then fewer pooled FPs, then index;
5. secondary: everything in the proposal + combined detector with
   lle_chaos_test (AND and OR);
6. m-sensitivity of lle_chaos_test (m = 3, 4) on Phase 3 test seeds at 256
   (report only);
7. adopt winner opt-in + tests; docs/PHASE4_UPO_VALIDATION.md.

## Status
| Step | State | Commit |
|---|---|---|
| Scaffold (config, systems incl. two_tone, detector, runner, methods) | done | 1695cfb |
| Dev exploration at 256 (5 run configs) | done: results/dev/explore.jsonl; summary via explore_summary | 7b8f6e1 |
| Dev 512 for shortlisted runs (cao:7, 2:15, 2:7; 50 surr) | done: results/dev/dev512.jsonl | (this commit) |
| m-sensitivity check | done (results/tables/m_sensitivity.md) | (this commit) |
| Candidates + budget + PREREGISTRATION.md | done: pushed before any test run; 150 test seeds (2000-2149), 50 surrogates | (this commit) |
| Test run + decide | pending | |
| Adopt winner (opt-in) | pending | |
| docs/PHASE4_UPO_VALIDATION.md | pending | |

## Notes
- Timing (single process, 50 surrogates, production run): 3-5 s per 256
  window, 7-11 s per 512 window.
- Phase 2E dev records (core/noise, 256): extension-modulus gate at delta 0.2
  removes all sinusoid Level-B windows (21/30 -> 0/30), WN 7/100 -> 0,
  AR1 11/100 -> 2, keeps clean chaos; Henon 30 dB stays 11/30 (power problem).
- detector.run with a fixed m equal to Cao's m reproduces run_upo_analysis
  bitwise (checked on 3 windows).
- m-sensitivity (Phase 3 test seeds, 256): m=2 reproduces stored C1 700/700.
  Detections/100 m=2/3/4: WN 4/7/9, AR1 3/1/1, sinusoid 0/0/0, logistic,
  henon, logistic@20, henon@20 all 100/100/100. Clean bias logistic
  -0.001/-0.009/-0.034, henon -0.005/-0.005/-0.008.
- Dev exploration (256) shortlist: C1 production + extension-modulus gate >= 1.2
  (WN0 AR2 sin0 2tone0; Henon30+20 18/60); C2 fixed m=2, M=15 + median hybrid
  >= 1.2 (WN0 AR1 sin0 2t0; 60/60; loc 0.009); C3 fixed m=2, M=7 + trimmed-mean
  hybrid >= 1.2 (WN0 AR3 sin0 2t0; 57/60). Baseline: WN7 AR11 sin21 2t15; 19/60.
- PREREGISTRATION.md pushed. Do NOT edit methods.py, detector.py, decide.py or
  the rule; changes only as dated amendments. NEXT: `python -m
  experiments.phase4_upo.run_phase4 --phase test --methods all --workers 4`
  (resumable; log /tmp/claude-0/p4_test.log; expected ~4.9 h), then
  `analysis --phase test --file test` and `decide --phase test --file test`.
