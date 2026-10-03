# arrhythmia-detection-pipeline

Code, preregistrations, results and reports of a ten-phase study that asks whether deterministic chaos can be
detected in human heart-rate (RR-interval) data, and how ectopic beats, noise and detector design produce false
detections. The work runs from a reconstructed R-peak and unstable-periodic-orbit pipeline (Phases 1-2E) through
validated detectors on synthetic and cardiac-model data (Phases 3-9) to a final confirmatory analysis on two
PhysioNet databases (Phase 10). Every confirmatory analysis was preregistered and committed before its data were
analysed (`experiments/final_refinement/PREREG_INTEGRITY.md`).

Release `paper-v1` adds independent verification, errata, publication figures and tables, and this
reproducibility package (`docs/FINAL_SUMMARY.md`).

## Repository structure

| path | contents |
|---|---|
| `final_pipeline.py` | the pipeline: Pan-Tompkins R peaks, delay embedding, Rosenstein LLE, So et al. UPO detection; opt-in detectors adopted by later phases (`lle_chaos_test`, gated UPO, `combined_chaos_*`, `rr_detrend`, `masked_growth_chaos_test`) |
| `tests/` | pytest suite (463 tests) |
| `experiments/phase*/` | one directory per phase: code, `PREREGISTRATION*.md`, `HANDOFF.md` (living notes), sources/methods files, `results/` (raw JSONL and tables) |
| `experiments/final_refinement/` | verification (`verify_all.py`, `VERIFICATION.md`), `ERRATA.md`, `KNOWN_ISSUES.md`, `PREREG_INTEGRITY.md`, figures (`make_figures.py`, `figures/`, `tables/`, `FIGURE_NOTES.md`), `download_data.py` |
| `docs/` | phase reports and the records for writing: `RESULTS_INDEX.md`, `CONCLUSIONS_CHECKLIST.md`, `SOURCES_USED.md`, `AI_USE_RECORD.md`, `FINAL_SUMMARY.md` |
| `reproduce.sh` | single entry point (below) |
| `requirements-lock.txt` | pinned environment for all phases |

## Phases: reports and preregistrations

| phase | question | report | preregistration |
|---|---|---|---|
| 1 | Pan-Tompkins R-peak detector | `docs/PHASE1_PAN_TOMPKINS.md` | – |
| 2E | frozen UPO pipeline on synthetic systems | `docs/PHASE2E_SYNTHETIC_VALIDATION.md` | – |
| 3 | LLE as a size-controlled chaos test | `docs/PHASE3_LLE_VALIDATION.md` | `experiments/phase3_lle/PREREGISTRATION.md` |
| 4 | UPO detector specificity (proposal: `docs/PHASE4_UPO_PROPOSAL.md`) | `docs/PHASE4_UPO_VALIDATION.md` | `experiments/phase4_upo/PREREGISTRATION.md` |
| 5 | combined detector on RR-like series | `docs/PHASE5_RR_STRESS_TEST.md` | `experiments/phase5_rr/PREREGISTRATION.md` |
| 6 | ectopy and drift robustness | `docs/PHASE6_ROBUSTNESS.md` | `experiments/phase6_robust/PREREGISTRATION.md` |
| 7 | MIT-BIH Arrhythmia Database | `docs/PHASE7_MITBIH_RESULTS.md` | `experiments/phase7_mitbih/PREREGISTRATION.md` |
| 8 | cardiac-model RR intervals | `docs/PHASE8_CARDIAC_CHAOS.md` | `experiments/phase8_cardiac/PREREGISTRATION.md` |
| 9 | ECG waveform and noise-robust measures; nsrdb vs chfdb | `docs/PHASE9_WAVEFORM_NOISE_ROBUST.md` | `experiments/phase9_waveform/PREREGISTRATION.md`, `PREREGISTRATION_F.md` |
| 10 | detection limits, faithful noise titration, ectopy dose-response (nsr2db, chf2db) | `docs/PHASE10_FINAL_RESULTS.md` | `experiments/phase10_final/PREREGISTRATION.md` |

## Environment

Python 3.13 (runs used 3.13.12-3.13.15; the final refinement 3.13.14), single-threaded BLAS.

```bash
uv venv -p python3.13 .venv
uv pip install -p .venv/bin/python --no-deps -r requirements-lock.txt   # --no-deps: see KNOWN_ISSUES.md 8
export PY=.venv/bin/python OMP_NUM_THREADS=1
```

## How to reproduce

```bash
bash reproduce.sh fast            # verification, every phase's tables and decisions, figures (~15 min, no data)
bash reproduce.sh tests           # pytest (~5 min; one known platform-dependent failure, KNOWN_ISSUES.md 1)
bash reproduce.sh data            # download the 7 PhysioNet databases, SHA-256 verified
bash reproduce.sh slow list       # full reruns per phase: commands and expected runtimes (~55-65 h in total)
bash reproduce.sh slow 10         # e.g. rerun Phase 10 from scratch in a separate clone
```

`fast` regenerates everything from the stored raw results and checks the outputs against the committed files:
all are byte-identical except Phase 2E's tables and plots, which differ only in formatting
(KNOWN_ISSUES.md 7). `slow` reruns a phase's runners in a separate clone; full reruns on another machine
reproduce the counts approximately, not bit for bit (KNOWN_ISSUES.md 2). Run every command from the repository
root.

## Verification, errata and known issues

- `experiments/final_refinement/VERIFICATION.md`: independent recomputation of the key reported numbers from the
  raw result files.
- `experiments/final_refinement/ERRATA.md`: every discrepancy found, with evidence; corrected passages in the
  reports are marked `[Erratum En]`. No erratum changes a result, decision or conclusion; one
  conclusion-wording issue (E14) was identified and left for the author's decision.
- `experiments/final_refinement/KNOWN_ISSUES.md`: the AVX-512-dependent unit test, cross-machine replicability,
  the frozen Phase 9 decision script, partial development runs, GEE working correlations, the neurokit2/pandas pin.

## Data and licences

No PhysioNet file is in this repository; `reproduce.sh data` (`experiments/final_refinement/download_data.py`)
downloads them and checks every file against PhysioNet's `SHA256SUMS.txt` and the pinned manifests in
`experiments/final_refinement/data_manifests/`. All seven databases are distributed by PhysioNet under the
**Open Data Commons Attribution License v1.0** (as stated on each database page, checked 2026-10-03), which
requires attribution:

| database (version 1.0.0) | DOI | used in | publication PhysioNet asks users to cite |
|---|---|---|---|
| MIT-BIH Arrhythmia Database | 10.13026/C2F305 | Phases 7, 9, 10 | Moody GB, Mark RG. IEEE Eng Med Biol 20(3):45-50 (2001) |
| MIT-BIH Normal Sinus Rhythm Database | 10.13026/C2NK5R | Phases 9, 10 | – |
| BIDMC Congestive Heart Failure Database | 10.13026/C29G60 | Phases 9, 10 | Baim DS et al. J Am Coll Cardiol 7(3):661-670 (1986) |
| Normal Sinus Rhythm RR Interval Database | 10.13026/C2S881 | Phase 10 | – |
| Congestive Heart Failure RR Interval Database | 10.13026/C2F598 | Phase 10 | – |
| MIT-BIH Noise Stress Test Database | 10.13026/C2HS3T | Phase 9 | Moody GB, Muldrow WE, Mark RG. Computers in Cardiology 11:381-384 (1984) |
| QT Database | 10.13026/C24K53 | Phase 9 | Laguna P et al. Computers in Cardiology 24:673-676 (1997) |

Every database also requires the standard PhysioNet citation: Goldberger AL et al. PhysioBank, PhysioToolkit,
and PhysioNet: components of a new research resource for complex physiologic signals. Circulation
101(23):e215-e220 (2000). One open-access article (Wu et al. 2009, PLoS ONE, CC-BY) and five PLOS supplementary files used to
implement methods are stored under `experiments/*/sources/`. Other methods and papers used are listed in
`docs/SOURCES_USED.md`.

## Citing and archiving

Citation metadata: `CITATION.cff` (author Brooks Bezanson). No software licence file has been chosen yet; add one
before archiving if the code is to be reusable by others.

To archive a release on Zenodo (done by the author):
1. Log in to zenodo.org with GitHub and enable the repository under *Account → GitHub*.
2. Create a GitHub release from tag `paper-v1` (merge the branch into `main` first if the release should be
   cut from `main`). Zenodo archives the release and mints a DOI.
3. Add the DOI to `CITATION.cff` (`doi:` field) and to this README in a follow-up commit.
