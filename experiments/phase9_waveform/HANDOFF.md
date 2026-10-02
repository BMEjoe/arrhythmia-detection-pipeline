# Phase 9 handoff (living document; update after every step)

Read this first when resuming. Branch `claude/bold-dirac-fc69kl`. Never merge into main,
never open a pull request. Commit and push after every step.

**Branch base.** The task said "start from main, which contains Phases 2E-7 and a
partially completed Phase 8". origin/main actually stops at Phase 6 (`1b4b7a0`); Phases
7-8 were on `claude/amazing-cray-9jbx6p` (head `816c0db`), a descendant of main. This
branch was fast-forwarded to `816c0db` (no merge into main, main untouched).

## Task (user, Phase 9) in brief
ECG waveform analysis with appropriate surrogates + noise-robust chaos measures.
- Part 0: docs/PHASE8_CARDIAC_CHAOS.md from committed Phase 8 results only (no runs).
- Part A: ECGSYN (McSharry 2003) driven by EXTERNAL RR (verify shape + beat times);
  RR from Phase 5-6 nulls (N1-N6, S1, S2 5/10, E1, E2, E3 5/10, E4, E5 5/10, E6) with
  predeclared ectopic morphology, and every Phase 8 regime (keep labels and split).
  A2 APD restitution (morphology) family with ground-truth LE for the exact RR input;
  >= 1 morphology family in TEST. A3 360 Hz + MIT-BIH resolution, nstdb noise (bw, ma,
  em) at predeclared SNRs, Phase 7 V jitter. A4 windows (waveform 2/5 min; beats
  256/512/1024). Realism table + plots. SPLIT recorded here BEFORE any measure is run on
  any Phase 9 series.
- Part B: waveform embedding (AMI delay, predeclared dimension), PPS (Small 2001, incl.
  noise-radius rule), cycle shuffle (Theiler), twin surrogates (if verifiable). B3 FP on
  waveform nulls BEFORE use; FP > 7 % on any null => UNUSABLE. B4 NLP + Part D measures.
  B5 compute check (drop waveform arm if over budget).
- Part C: published delineator (verified), per-beat RR, QT/RT, QRSd, QRS amp, T amp;
  delineation accuracy; multivariate embedding; multivariate IAAFT (Schreiber-Schmitz).
- Part D: SDLE, FSLE, (eps,tau)-entropy, PE + complexity-entropy plane, RQA DET, LP noise
  reduction (optional preprocessing); verify each; detection rule each; apply to RR (full
  Phase 8 library), multivariate features, waveform.
- Part E: <= 6 candidates on DEV only; DIAGNOSIS.md; PREREGISTRATION.md pushed BEFORE any
  TEST seed; PASS iff <= floor(0.07N) on EACH null (Phase 5-6 nulls + ectopy at input
  level, their nstdb versions, every NON-CHAOTIC regime of every model and morphology
  family, dev and test); WINNER max pooled detection on TEST CHAOTIC at most realistic
  variant; tie fewer pooled null detections, then index. Baseline c1_frozen512 on same RR
  windows. Budget <= ~16 h on 4 workers (reduce secondary seeds first, never nulls).
  Adopt winner opt-in (default off) with tests.
- Part F (only if a winner): separate prereg BEFORE downloading nsrdb + chfdb; W1 rate per
  group (subject-cluster bootstrap) vs synthetic FP; W2 group difference; then EXPLORATORY
  Phase 7 MIT-BIH.
- Part G: docs/PHASE9_WAVEFORM_NOISE_ROBUST.md.
- Seeds: DEV from 9000, TEST from 9500. Never reuse earlier seeds. Never use TEST families
  or TEST seeds for development. MIT-BIH only EXPLORATORY.

## Ground rules (as Phases 3-8)
- final_pipeline.py: opt-in options only, defaults unchanged; then
  `bash experiments/phase3_lle/groundrule_check.sh` (needs /tmp/claude-0/rep_before.jsonl
  and /tmp/claude-0/replicability_env_before.json created BEFORE the change:
  `python -m experiments.phase3_lle.replicate_check --save /tmp/claude-0/rep_before.jsonl`;
  `python -m experiments.phase2e.run_phase2e --replicate --workers 4` then copy
  experiments/phase2e/results/replicability.json to /tmp/claude-0/replicability_env_before.json
  and `git checkout experiments/phase2e/results/replicability.json`).
- New code in experiments/phase9_waveform/.

## Environment (rebuilt 2026-10-02 in a fresh container)
- `uv venv -p /usr/bin/python3.13 /root/venv313 && uv pip install -p /root/venv313/bin/python -r requirements.txt wfdb==4.3.1 numba==0.68.0`
- Python 3.13.14, numpy 2.1.3, scipy 1.18.1, scikit-learn 1.9.1, pandas 3.0.6,
  matplotlib 3.11.2, wfdb 4.3.1, numba 0.68.0. OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1.
- pytest baseline (2026-10-02): default env **455 passed, 1 failed**
  (pre-existing `test_upo_stability.py::...[1-prl_norm]`); AVX-512 disabled **456 passed**.
  Same pattern as recorded in Phases 3-6 (Phase 7 added the extra tests).
- Machine: 4 CPUs, 15 GB RAM.

## Status
| Step | State | Commit |
|---|---|---|
| env rebuild + pytest baseline | done (455+1 / 456) | (this commit) |
| Part 0 docs/PHASE8_CARDIAC_CHAOS.md (+ read-only tables phase8_closeout/) | done | (this commit) |

## Notes
