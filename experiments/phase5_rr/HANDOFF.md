# Phase 5 handoff (living document; update after every step)

Read this first when resuming. Branch `claude/amazing-allen-p4ssx3`. NOTE: the
task said main contains Phase 4, but at the start of Phase 5 PR #2 (Phase 4)
was still OPEN and main was at 858d616 (Phase 3). This branch = main + Phase 4
(c561e85) + Phase 5, so Phase 5 commits also appear in the open PR #2.
Never merge into main, never open a pull request. Commit and push after every step.

## Task (user, Phase 5) in brief
Stress test of the FROZEN combined detector on RR-like signals; evaluation,
not tuning. Detector: `fp.analyze_segment(rr, fp.phase4_upo_config(replace(CFG,
lle_chaos_test=True)))`; LLE = out["lle_chaos_test"]["detected"], UPO =
out["upo"]["instability_gate_detected"]; PRIMARY = AND; secondary each alone, OR.
No MIT-BIH data. Steps: 0 equivalence check; 1 generators + realism (commit
before prereg); 2 preregistration (>= 100 test seeds from 3000, <= ~8 h on 4
workers; cut secondary/generalization seeds first, never N1-N6); PASS iff AND
fires <= 7/100 (scaled) on EACH of N1-N6 quantized at 256; 3 Phase 4 m-check
(C2 with upo_fixed_dimension 3, 4 on Phase 4 test seeds 2000-2149, 8
conditions); 4 test run; 5 FAIL -> dev-only diagnosis + docs/PHASE6_PROPOSAL.md
(no fixes), PASS -> combined_chaos_config() / combined_chaos_detected() opt-in
helpers + tests; 6 docs/PHASE5_RR_STRESS_TEST.md.

## Ground rules (as Phases 3 and 4)
- `final_pipeline.py`: only opt-in options, defaults unchanged; then
  `experiments/phase3_lle/groundrule_check.sh` (pytest default = all pass
  except the pre-existing stability failure; AVX-512 off all pass; same-env
  44/44 vs /tmp/claude-0/rep_before.jsonl; --replicate 13/44 identical to
  /tmp/claude-0/replicability_env_before.json).
- Environment: /root/venv313, `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1`.
- Long runs die when the container idles between turns: babysit in-turn;
  runners are resumable (skip finished task ids).
- physionet.org: NOT reachable (proxy 403 CONNECT, policy denial, checked
  2026-09-29). pypi.org reachable; wfdb not installed.

## Status
| Step | State | Commit |
|---|---|---|
| 0 equivalence check (70 windows) | done: 0 differences (decisions and p identical) | (this commit) |
