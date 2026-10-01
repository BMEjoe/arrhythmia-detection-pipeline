# Phase 8 handoff (living document; update after every step)

Read this first when resuming. Branch `claude/amazing-cray-9jbx6p` (Phases 2E-7 on it,
Phase 7 head 4283743). Never merge into main, never open a pull request. Commit and
push after every step.

## Task (user, Phase 8) in brief
Build/evaluate a chaos detector for cardiac-relevant chaos (flows, delay systems,
cardiac models) with near-zero FP on realistic non-chaotic RR incl. ectopy.
- Part A: cardiac model library from accessible published sources (MODELS.md: source,
  equation numbers, parameter table), each verified by reproducing a published result;
  >= 4 families, >= 2 with native beat-to-beat intervals. CHAOTIC / NON-CHAOTIC regimes
  by ground-truth largest Lyapunov exponent from the equations (Benettin / delay-valid,
  convergence + uncertainty); ambiguous regimes discarded. RR series 256/512/1024;
  variants (i) clean (ii) dynamical noise x2 levels (iii) measurement noise 30/20 dB
  (iv) quantized 1/360 s (v) (iv)+Phase 6 ectopy (S2 5 %, E1, E3 10 %) (vi) (v)+Phase 7
  V jitter. Realism table + example plots. MODEL-LEVEL SPLIT fixed in this file before
  any detector tuning: >= 2 chaotic families entirely TEST.
- Part B: <= 4 candidate detectors with surrogate/null control, development models +
  Phase 5-6 nulls only; diagnose successes/failures.
- Part C: PREREGISTRATION.md before any TEST seed: PASS iff <= floor(0.07 N) on each of
  N1-N6, S1, S2 5/10 %, E1-E6, and every NON-CHAOTIC regime (test+dev), quantized, at the
  candidate's window length; WINNER max pooled detection on TEST chaotic regimes at (vi);
  tie: TEST (iv), then index. Secondary: by model/regime/length/variant; frozen Phase 6
  detector baseline; Phase 5 G2/G3. Budget <= ~10 h on 4 workers.
- Part D: run, report all, adopt winner opt-in (default off) + tests.
- Part E (only if winner): separate prereg BEFORE downloading nsr2db / chf2db; R1 rate
  per group (subject-cluster bootstrap) vs synthetic FP; R2 group difference; then
  EXPLORATORY Phase 7 MIT-BIH windows.
- Part F: docs/PHASE8_CARDIAC_CHAOS.md.
- MIT-BIH never used for tuning/selection.

## Ground rules (as Phases 3-7)
- final_pipeline.py: opt-in options only, defaults unchanged; then
  `bash experiments/phase3_lle/groundrule_check.sh` (needs /tmp/claude-0/rep_before.jsonl
  and /tmp/claude-0/replicability_env_before.json created BEFORE the change, see Phase 3
  HANDOFF).
- Env: /root/venv313 (Python 3.13.14, numpy 2.1.3, scipy 1.18.1, sklearn 1.9.1, wfdb 4.3.1)
  + numba 0.68.0 (Phase 8 model integration only). OMP_NUM_THREADS=1.
- Literature access (checked 2026-10-01): arXiv, PMC/Europe PMC (full-text XML with MathML),
  PLOS supplements, physionet.org reachable; APS, AIP, Elsevier/ScienceDirect, physiology.org,
  MDPI, Springer (bot challenge) BLOCKED. Downloaded sources in sources/ (open-access only).

## Source survey (Part A)
| family | accessible source | status |
|---|---|---|
| Seidel-Herzel baroreflex | Dudkowska & Makowiec, arXiv q-bio/0603016 (restates Physica D 115:145, 1998, which is blocked) | to implement |
| DeBoer beat-to-beat (+ chaotic respiration) | AJP 253:H680 blocked; only simplified restatement without parameters (arXiv physics/0503053) | DROP (no accessible parameters) |
| phase-resetting / circle map (Glass, Shrier) | Diagne, Bury, ..., Glass et al. PLoS Comput Biol 2026, PMC12900431, Eq. 1 + S1 Text/S1 Table (PRC fit); parasystole Eq. 2 / S2 Text Eq. 4,7 | to implement |
| AV node (Sun, Amellal, Glass, Billette 1995) | Zhao & Schaeffer arXiv math/0609106 Eq. 51 (restatement with parameters) | to implement |
| AV node during AF | no accessible equations | DROP |
| modified van der Pol, coupled SA-AV-HP with delays | single oscillator arXiv 2202.00058 (Eq. 2; 2-D autonomous, cannot be chaotic); three-coupled version restated in PMC9938421 (Eq. 1, Table 1, LE-confirmed chaos windows) | to implement / verify (restatement has typos) |
| Mackey-Glass | Phase 5 G3 | to verify against an open source |

## MODEL-LEVEL SPLIT (fixed 2026-10-01, before ANY detector run on any Phase 8 model)
| family | module | chaotic regimes? | split |
|---|---|---|---|
| Mackey-Glass delay feedback (Glass & Mackey Scholarpedia Eq. 1) | models/mackey_glass.py | yes (tau >= 17) | **DEV** |
| phase-resetting map, hiPSC-CM PRC (Diagne et al. 2026 Eq. 1, S1 Text/Table) | models/phase_reset.py | yes (e.g. tau = 1.2) | **TEST** |
| three coupled modified vdP with delays (da Silva Lima, Savi, Bessa 2024 Eqs. 1-5, Fig 1b) | models/coupled_vdp.py | yes (forced) | **TEST** |
| AV node (Sun et al. 1995 via Zhao & Schaeffer Eq. 51) | models/av_node.py | NO (LE <= 0 everywhere valid) | **TEST** (non-chaotic regimes only) |
Development positives may also use the Phase 5 synthetic chaotic generators (P1, P2, G1-G3);
Phase 5 G2/G3 results are then NOT out-of-sample (noted for the report).
TEST models: only ground-truth LE (model equations), verification and realism/plots with
dedicated non-test seeds are computed before the preregistration; NO detector is run on
any TEST-model series before PREREGISTRATION.md is pushed.

Dropped: Seidel-Herzel (could not reproduce the restating source's Figs 4-5 with 6 readings
of ambiguous details; MODELS.md), DeBoer (no accessible parameters), Gois-Savi restatement
PMC9938421 (restated Table 1 gives a stable equilibrium, not its own Table 3), AV node in AF
(no equations), modulated parasystole (PRC form with m, p, r, N not given).

## Status
| Step | State | Commit |
|---|---|---|
| source survey | done | 6afdae8 |
| models implemented + verified: phase_reset (tau 0.66 period-3, 1.2 chaotic LE +0.135), av_node (H_bif 56.9078 exact; LE ~ -0.006 everywhere), mackey_glass (LE tau 50 = 0.0057-0.0059 vs published 0.0055-0.0058), coupled_vdp (normal 89.4 bpm vs ~90; pathological 137-256 bpm dispersed); Seidel-Herzel DROPPED | done | (this commit) |
