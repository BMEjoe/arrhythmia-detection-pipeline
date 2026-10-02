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

## PHASE 9 SPLIT (recorded 2026-10-02, BEFORE any detector or chaos measure was run on any Phase 9 series)

Only generator verification (ECGSYN, KTz, LR1), ground-truth Lyapunov exponents (KTz),
nstdb download/checksums and the Part A realism table/plots (RR/HRV statistics, amplitudes,
beat-detection counts; no chaos measure, no detector) had been run on Phase 9 material at
this point.

| family | source | labels | split |
|---|---|---|---|
| Phase 5-6 nulls (17: N1-N6, S1, S2 5/10, E1, E2, E3 5/10, E4, E5 5/10, E6) | Phase 5/6 generators via families.null_window | NULL | DEV seeds 9000-9499; TEST seeds >= 9500 |
| Mackey-Glass maxima intervals (Phase 8) | Phase 8 | 7 CHAOTIC, 3 NON-CHAOTIC (Phase 8 regimes.json) | **DEV** |
| Phase 5 positives P1 Henon RR, P2 logistic RR, G1 Lorenz maxima, G2 Rossler, G3 Mackey-Glass | Phase 5 generators | chaotic (development positives only; not in any TEST pool) | **DEV** |
| phase-resetting map (Phase 8) | Phase 8 | 7 CHAOTIC, 5 NON-CHAOTIC | **TEST** |
| coupled modified vdP (Phase 8) | Phase 8 | 4 CHAOTIC, 6 NON-CHAOTIC | **TEST** |
| AV node (Phase 8) | Phase 8 | 4 NON-CHAOTIC | **TEST** |
| KTz paced-cell map -> APD -> QT / T wave (A2 morphology family) | morph/ktz.py, ground_truth_ktz.py | regimes = (P_nom, input) pairs: 11 CHAOTIC (92,none) (92,S2_5) (92,E3_10) (100,none) (146,none) (146,S2_5) (146,E1) (230,none) (232,none) (260,none) (278,none); 31 NON-CHAOTIC (92,E1) (100,E1) (100,E3_10) (230,S2_5) (230,E1) (230,E3_10) (232,E1) (232,E3_10) (260,E1) (278,E1) (278,E3_10) and P_nom = 120, 150, 200, 250, 300 with every input; 6 AMBIGUOUS discarded | **TEST** (only verified morphology family; **development had no morphology family**) |

Morphology families: the modified Luo-Rudy EAD family (Tran et al. 2009) was DROPPED (not
verifiable; METHODS.md 3.2), so only one morphology family exists and it is entirely TEST.

Disclosure (for the report): the Phase 8 TEST families (phase_reset, coupled_vdp, av_node)
were evaluated in aggregate in Phase 8 with TEST seeds 5000-5099. Known before Phase 9 was
designed: the frozen Phase 6 detector detected 52-66/990 TEST chaotic windows at (vi) and was
specific on every null and non-chaotic regime at (iv) but fired 11/60 on non-chaotic coupled vdP
with bigeminy; nonlinear prediction vs (ectopy-preserving) IAAFT was sensitive but failed on the
forced quasi-periodic coupled-vdP regimes (rho, omega) = (2, 5.6) and (5.45, 5.6); noise titration
failed on every ectopy pattern and 17/18 non-chaotic regimes; the UPO component limits the frozen
detector's power.  Phase 9 uses new seeds (>= 9500) for every TEST family and never uses a TEST
family or a TEST seed for development.

Seeds: DEV windows 9000-9499; TEST windows >= 9500; KTz ground truth 990000-990002; realism
980000-980004.  nstdb noise: DEV windows use nstdb signal 0, TEST windows signal 1.

Windows (A4, predeclared): beat-level windows of 256, 512 and 1024 beats; waveform windows of
2 and 5 minutes cut from the start of a beat-level record (512 or 1024 beats).  Any waveform
downsampling rate is predeclared in Part B and checked on a verification subset.

Variants (library.py): clean (true beats, no noise); ectopy (models: S2_5 / E1 / E3_10; KTz:
part of the regime; nulls: built in); nstdb noise <spec><snr> with spec in mix, bw, ma, em and
SNR 24 / 12 / 6 dB (nst definitions); beats true or detected (Pan-Tompkins); jitter (Phase 7 V
offsets at detected V beats).  Most realistic variant = ectopy + nstdb mix 12 dB + detected beats
+ jitter.

## Status
| Step | State | Commit |
|---|---|---|
| env rebuild + pytest baseline | done (455+1 / 456) | (this commit) |
| Part 0 docs/PHASE8_CARDIAC_CHAOS.md (+ read-only tables phase8_closeout/) | done | cdaf022 |
| A1 ECGSYN external-RR + verification vs ecgsyn.c / paper | done | c634d05 |
| A1 ectopic morphology (morphology.py), RR sources with beat types (families.py; nulls bit-identical to Phase 8), window builder (library.py), nstdb + nst SNR (noise.py) | done | (this commit) |
| A2 KTz verified (LE 0.00344 at P=92 vs published ~0.0035); KTz ground truth for exact inputs (ktz_labels.json); LR1-EAD DROPPED (not verifiable) | done | (this commit) |
| A realism table + plots (realism.py, results/realism, plots/A_*) | done | (this commit) |
| PHASE 9 SPLIT recorded (above) before any measure/detector on Phase 9 series | done | (this commit) |
| D1-D6 measures implemented + verified (measures.py, verify_measures.py, METHODS 5) | done | (see git log) |
| B1-B2 waveform embedding + PPS / CS / TS (surrogates_wave.py); TS verified (sync test); PPS/CS verification run verify_surrogates.py -> results/verification/surrogates.json | TS done; PPS/CS run in progress | |
| B3 waveform FP (fp_waveform.py; summary results/dev/fp_waveform_summary.md): periodic ECG rejected 5/5 (PPS, CS, TS) -> all UNUSABLE; B5 ~200 s/window -> **waveform arm DROPPED** | done | |
| C1 NK2 DWT delineator vs QTDB (verify_delineation.py): FAILED -> dropped; QT, QRSd dropped | done | |
| C1' fixed-window RTp / Ta (beatfeat.py, verify_beatfeat.py): FAIL at mix 12 (T apex SD 80 ms) -> not usable at the primary variant; KTz chaos (T wave only) invisible to verified features | done | |
| C2 multivariate IAAFT verified (verify_multivariate.py) | done | |
| DEV-only circle-map controls (devmaps.py; QP + locked, NON-CHAOTIC) added before any measure ran | done | |
| E dev diagnosis bank (dev_diag.py; stopped at 379/1,820, pattern clear) + DIAGNOSIS.md | done | |
| E candidates9.py (k1-k5) + DEV calibration (calibrate9.py, 6,450 windows) -> thresholds frozen | done | |
| E PREREGISTRATION.md final, pushed BEFORE any TEST seed | done | |
| E TEST run (12,290 windows, 8.3 h) + decide9: PASS k3, k4, k5; WINNER k3_mnlp_growth_ann (90/1,320; 0/10,100); k1 fails (non-chaotic vdP + ectopy), k2 fails | done | 801256a |
| E adoption: final_pipeline.masked_growth_chaos_test (standalone opt-in; analyze_segment unchanged), tests/test_phase9_masked_growth.py; groundrule_check passed | done | 8939ab7 |
| F PREREGISTRATION_F.md pushed before download (89d4ba9; 1 amendment: adopted form); nsrdb 1/180, chfdb 0/150; W1 no, W2 no | done | |
| F EXPLORATORY MIT-BIH: 1/190 | done | |
| G docs/PHASE9_WAVEFORM_NOISE_ROBUST.md | done | |

**ALL PHASE 9 STEPS DONE.**

## Notes
- Background jobs (restart if the container restarted; both resumable / rerunnable):
  `OMP_NUM_THREADS=1 nohup /root/venv313/bin/python -m experiments.phase9_waveform.dev_diag --workers 3 &`
  `OMP_NUM_THREADS=1 nohup /root/venv313/bin/python -m experiments.phase9_waveform.verify_surrogates --workers 4 &`
- Never `pkill -f <pattern>` with the plain pattern (kills the calling shell); use `pkill -f "phase9_waveform[.]dev_diag"`.
- KTz windows: RR is the pacing input (constant or the ectopy pattern); its chaos is in the T wave only.
