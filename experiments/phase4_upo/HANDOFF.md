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
| Scaffold (config, systems incl. two_tone, detector, runner, methods) | done | (this commit) |
| Dev exploration at 256 (5 run configs) | running: results/dev/explore.jsonl, log /tmp/claude-0/p4_explore.log | |
| m-sensitivity check | pending | |
| Candidates + budget + PREREGISTRATION.md | pending | |
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
