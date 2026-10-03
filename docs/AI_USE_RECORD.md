# Record of work done by Claude Code (from the git history)

Purpose: a factual record, taken from the repository's git history, of what was implemented, run and written by
Claude Code (Anthropic's coding agent), to support the author's disclosure of AI use. It does not describe the
author's role; the author will do that.

**How this was compiled (2026-10-03).** From the full history of this repository (207 commits before the final
refinement), with `git log` author fields, commit messages, `Co-Authored-By` trailers and the files each commit
added or changed. Claude Code committed under the identity `Claude <noreply@anthropic.com>`; those commits are
listed below by phase with their hash ranges (inclusive, in history order). Each phase was carried out in a
Claude Code session on a `claude/...` branch and later merged into `main` by the repository owner (merge commits
`858d616`, `ea7e068`, `1b4b7a0` on GitHub; later phases reached `main` without merge commits). Every phase
report states the instructions it followed in its handoff file (`experiments/*/HANDOFF.md`).

**Commits not authored by Claude Code.** 11 commits carry the repository owner's identities (BMEjoe / Brooks):
`77c2996`, `4c71f61`, `338770d`, `9383a31`, `72eea3f`, `4d9f2d6`, `dd43fd0`, `88835d5` and the three merge commits.
One of them, `88835d5` ("Complete Phase 2E analysis, diagnostics, and validation outputs", 2026-09-27, 79 files:
`experiments/phase2e/` code, results, tables and plots), carries the trailer
`Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. The others carry no Claude trailer. This record does not
attribute the content of those commits.

## Summary by phase

| phase | Claude Code commits | range | dates |
|---|---:|---|---|
| Phase 2E report and Phase 3 (incl. Parts A and C) | 18 | `689c0ab` .. `42a062f` | 2026-09-27 to 09-28 |
| Phase 4 | 21 | `1695cfb` .. `c561e85` | 2026-09-28 to 09-29 |
| Phase 5 | 13 | `cdd6a57` .. `f910d4d` | 2026-09-29 to 09-30 |
| Phase 6 | 17 | `7026a12` .. `ae2cfb6` | 2026-09-30 |
| Phase 7 | 22 | `8f23a40` .. `4283743` | 2026-10-01 |
| Phase 8 | 43 | `6afdae8` .. `816c0db` | 2026-10-01 to 10-02 |
| Phase 9 (incl. Phase 8 close-out) | 28 | `cdaf022` .. `71a4f7c` | 2026-10-02 |
| Phase 10 | 34 | `0a98b75` .. `1147f76` | 2026-10-03 |
| Final refinement (this phase) | all commits after `1147f76` on `claude/serene-sagan-2qgi8e` | `4655009` .. tag `paper-v1` | 2026-10-03 |

Many commits in Phases 4-10 are checkpoints of running jobs ("test-run checkpoint", "spike-in progress").

## Phase 2E report and Phase 3 (`689c0ab` .. `42a062f`)

- **Written:** `docs/PHASE2E_SYNTHETIC_VALIDATION.md` (report on the Phase 2E outputs committed in `88835d5`),
  `docs/PHASE3_LLE_VALIDATION.md`, `docs/PHASE4_UPO_PROPOSAL.md` (proposal, no code), `experiments/phase3_lle/`
  README, HANDOFF and `PREREGISTRATION.md`.
- **Implemented:** Phase 3 modules `estimators.py` (baseline LLE test and candidates C1-C4: Rosenstein slope at
  τ = 1, Kantz-style divergence, ε-excluded variant, 0-1 test; IAAFT surrogates), `diagnose.py`, `tune.py`,
  `run_phase3.py` (with the preregistration guard reused by later runners), `analysis.py`, `decide.py`,
  `replicate_check.py`, `groundrule_check.sh`; `requirements.txt`; tests `test_lle_chaos_test.py`,
  `test_lle_short_embedding.py`, `test_cao_e1_diagnostics.py`, `test_phase3_estimators.py`.
- **Changed existing code:** `final_pipeline.py`: opt-in `keep_upo_on_short_lle_embedding` (A4), opt-in Cao E1
  diagnostics (A5), opt-in `lle_chaos_test` = the preregistered winner C1 (B5); `experiments/phase2e/analysis.py`:
  fix of the Tables B/C pooling of surrogate counts (A1); regenerated `posthoc` LLE file (A3).
- **Run:** diagnosis and tuning on development seeds; the preregistered test-seed evaluation (12,000 windows);
  the decision rule; the diagnosis of the pre-existing AVX-512-dependent test failure (A2).

## Phase 4 (`1695cfb` .. `c561e85`)

- **Implemented:** `experiments/phase4_upo/` (`systems.py`, `detector.py`, `methods.py`, `run_phase4.py`,
  `analysis.py`, `decide.py`, `explore_summary.py`, `m_sensitivity.py`, `config.py`); test `test_phase4_upo_gate.py`;
  opt-in instability-gated UPO options in `final_pipeline.py` (winner C2).
- **Run:** development exploration (256 and 512), the preregistered test run (4,200 windows; two resumes after
  container restarts), decision, m-sensitivity.
- **Written:** `PREREGISTRATION.md`, HANDOFF, `docs/PHASE4_UPO_VALIDATION.md`.

## Phase 5 (`cdd6a57` .. `f910d4d`)

- **Implemented:** `experiments/phase5_rr/` (RR-interval generators `systems.py`, `realism.py`,
  `equivalence_check.py`, `run_phase5.py`, `analysis.py`, `decide.py`, `detector.py`, `m_path_check.py`,
  `phase4_m_check.py`, `adoption_check.py`, `report_extras.py`, `config.py`); opt-in `combined_chaos_config()` /
  `combined_chaos_detected()` in `final_pipeline.py`; test `test_phase5_combined_chaos.py`.
- **Run:** the Phase 4 equivalence check, realism checks, development runtime run, the preregistered test run
  (5,400 windows), decision, m-sensitivity.
- **Written:** README (generator sources), `PREREGISTRATION.md`, HANDOFF, `docs/PHASE5_RR_STRESS_TEST.md`.

## Phase 6 (`7026a12` .. `ae2cfb6`)

- **Implemented:** `experiments/phase6_robust/` (ectopy generators, `methods.py`, `detector.py`,
  `diagnose_trend.py`, `partA_crash_check.py`, `shortcut_check.py`, `tune_summary.py`, `run_phase6.py`,
  `analysis.py`, `decide.py`, `adoption_check.py`); in `final_pipeline.py` the opt-in `rr_detrend` option, the
  crash fix for short LLE embeddings, the fix of `combined_chaos_detected` on UPO failure status and the adopted
  D2 detrend; tests `test_phase6_detrend.py`, `test_phase6_combined_detrend.py`.
- **Run:** Part A crash check, trend diagnosis and tuning (development), the preregistered test run (11,700
  windows; one serialization fix recorded as an amendment), decision.
- **Written:** README, `PREREGISTRATION.md` (and its two amendments), HANDOFF, `docs/PHASE6_ROBUSTNESS.md`.

## Phase 7 (`8f23a40` .. `4283743`)

- **Implemented:** `experiments/phase7_mitbih/` (`data.py`, `qc.py`, `qc_tables.py`, `detector.py`,
  `run_phase7.py`, `analysis.py`, `exploratory.py`, `audit.py`, `spike_in.py`); tests `test_phase7_analysis.py`,
  `test_phase7_detector.py`, `test_phase7_mitbih_helpers.py`.
- **Run:** MIT-BIH download with checksum verification (`DATA_MANIFEST.json`), R-peak quality control, the
  preregistered detector run (891 windows), preregistered and exploratory analyses, code audit, spike-in controls.
- **Written:** `QC.md`, `PREREGISTRATION.md`, HANDOFF, `docs/PHASE7_MITBIH_RESULTS.md`.

## Phase 8 (`6afdae8` .. `816c0db`)

- **Implemented:** `experiments/phase8_cardiac/` cardiac models (`models/phase_reset.py`, `av_node.py`,
  `mackey_glass.py`, `coupled_vdp.py`, `seidel_herzel.py` (dropped)), `ground_truth.py` (Lyapunov-exponent scans),
  `series.py`, `nulls.py`, `realism.py`, `methods.py` (candidates incl. an approximate noise titration),
  `dev_runner.py`, `calibrate.py`, `run_test.py`, `decide.py`.
- **Run:** ground-truth scans, development diagnosis and calibration, the preregistered TEST run (13,440 windows
  over two lengths), decision.
- **Written:** `MODELS.md`, `DIAGNOSIS.md`, `PREREGISTRATION.md`, HANDOFF. (The Phase 8 report was written in
  Phase 9.)

## Phase 9 (`cdaf022` .. `71a4f7c`)

- **Written first:** `docs/PHASE8_CARDIAC_CHAOS.md` (Phase 8 close-out from committed results) with
  `phase8_closeout/tables.py`.
- **Implemented:** `experiments/phase9_waveform/`: ECGSYN driven by external beat times (`ecgsyn.py`, with a build
  of the authors' C code in `tools/ecgsyn_c/` for verification), morphology families (`morph/ktz.py`, `lr1.py`
  (dropped), scans), `morphology.py`, `noise.py`, `library.py`, `families.py`, `regimes9.py`, `devmaps.py`,
  waveform surrogates (`surrogates_wave.py`), multivariate IAAFT (`multivariate.py`), beat features
  (`beatfeat.py`), six noise-robust measures (`measures.py`), candidates k1-k5 (`candidates9.py`),
  `calibrate9.py`, `run_test9.py`, `decide9.py`, `partf.py`, verification scripts (`verify_*.py`), diagnosis and
  realism scripts; the opt-in `final_pipeline.masked_growth_chaos_test` (K3) with `tests/test_phase9_masked_growth.py`.
- **Run:** verifications against published behaviour, development diagnosis and calibration, the preregistered
  TEST run (12,290 windows), decision, Part F on nsrdb / chfdb (with its own preregistration), exploratory
  MIT-BIH run, ground-rule check.
- **Written:** `METHODS.md`, `DIAGNOSIS.md`, `PREREGISTRATION.md`, `PREREGISTRATION_F.md`, HANDOFF,
  `docs/PHASE9_WAVEFORM_NOISE_ROBUST.md`.

## Phase 10 (`0a98b75` .. `1147f76`)

- **Implemented:** `experiments/phase10_final/`: noise titration from the accessible sources (`titration10.py`)
  and its verification (`verify_titration.py`), data layer with checksums (`data10.py`), `harness10.py`,
  `run10.py` (guarded), `overlap10.py` (subject-overlap checks and the amended rule), `analysis10.py`,
  `tables10.py`, `figures10.py`.
- **Run:** titration verification, development runs (partly; spike-in and robustness runs stopped at session
  end), confirmatory download of nsr2db / chf2db after the preregistration, overlap checks, the confirmatory run
  (988 windows, 9,536 segments, 8,200 synthetic windows, 1,328 spike-in tasks, 331 robustness windows), analyses,
  sensitivity analysis.
- **Written:** `METHODS.md`, `PREREGISTRATION.md` (draft, final, amendments 1-2), HANDOFF, report drafts,
  `docs/PHASE10_FINAL_RESULTS.md`.

## Final refinement (this phase; commits from `4655009`, tagged `paper-v1`)

- **Implemented:** `experiments/final_refinement/verify_all.py` (independent recomputation of 339 reported
  numbers), `errata_phase9_decision.py`, `prereg_integrity.py`, `make_figures.py`, `download_data.py`;
  `reproduce.sh`; `requirements-lock.txt`; `CITATION.cff`.
- **Run:** pytest baseline and ground-rule check, verification, regeneration of every phase's tables from stored
  results, figure generation.
- **Written:** `VERIFICATION.md`, `ERRATA.md` (and the marked corrections in the Phase 9 and 10 reports),
  `PREREG_INTEGRITY.md`, `KNOWN_ISSUES.md`, `FIGURE_NOTES.md`, the tables in `experiments/final_refinement/tables/`,
  `README.md`, and in `docs/`: `RESULTS_INDEX.md`, `CONCLUSIONS_CHECKLIST.md`, `SOURCES_USED.md`, this file,
  `FINAL_SUMMARY.md`. No figure caption and no paper text was written.

To reproduce this record: `git log --author=Claude --format='%h %ad %s' --date=short` and
`git log --format='%h %an %s%n%b' --author=BMEjoe --author=Brooks`.
