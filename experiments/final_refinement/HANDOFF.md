# Final refinement handoff (living document; update after every step)

Read this first when resuming. Branch `claude/serene-sagan-2qgi8e`. Never merge into main,
never open a pull request. Commit and push after every step.

**Branch base.** `origin/main` at `1147f76` (Phases 2E–10). Confirmed on 2026-10-03 that main
contains `docs/PHASE10_FINAL_RESULTS.md`. The session branch started at the same commit.
The clone was shallow; it was unshallowed (`git fetch --unshallow origin`) to read the full
history (207 commits) for the preregistration and AI-use records.

## Task (user, FINAL REFINEMENT PHASE) in brief
NO NEW SCIENCE: no new question, method, data download for new analysis, parameter change or
rerun with different choices. Allowed: independent verification, error correction (via the
error policy), figures, tables, documentation, packaging. Anything that needs new analysis goes
to "Open items" in `docs/FINAL_SUMMARY.md`.

- **Error policy.** A verification error is never fixed silently. Document it in
  `experiments/final_refinement/ERRATA.md` (what, evidence, affected numbers and reports, does a
  conclusion change). Fix only if the fix restores the preregistered analysis exactly; rerun only
  the affected computation deterministically; update every report/table that used it. If a fix
  would change a conclusion: STOP and tell the user before proceeding.
- **Part 1.** `verify_all.py` recomputes key numbers from raw result files with independent code;
  output `VERIFICATION.md` (reported, recomputed, source, MATCH/MISMATCH) + cross-report checks.
- **Part 2.** pytest + `groundrule_check.sh`; `KNOWN_ISSUES.md` (AVX-512 test); preregistration
  integrity vs git history with prereg commit vs first confirmatory-result commit; no restricted raw
  data committed.
- **Part 3.** `make_figures.py`: F1–F6 (+ supplementary), vector PDF + 300-dpi PNG from result
  files only; tables CSV + MD in `tables/`; `FIGURE_NOTES.md` (what is plotted, source, n,
  statistics, files read). NO captions, NO paper text.
- **Part 4.** `reproduce.sh` (fast/slow), PhysioNet download scripts with checksums, pinned
  environment, `CITATION.cff` (Brooks Bezanson), Zenodo notes, README rewrite.
- **Part 5.** `docs/RESULTS_INDEX.md`, `docs/CONCLUSIONS_CHECKLIST.md`, `docs/AI_USE_RECORD.md`
  (Claude Code's work only; not the author's role), `docs/SOURCES_USED.md` (only sources actually
  used, "to be read and verified by the author before citation").
- **Finish.** `docs/FINAL_SUMMARY.md`; tag the final commit `paper-v1` on this branch.

## Environment (rebuilt 2026-10-03 in a fresh container)
- Same recipe as Phase 10 (`experiments/phase10_final/requirements-phase10.txt`):
  `uv venv -p /usr/bin/python3.13 /root/venv313`; `uv pip install -r requirements.txt wfdb==4.3.1
  numba==0.68.0`; `neurokit2==0.2.13 ordpy==1.2.3`; re-pin `pandas==3.0.6`; statsmodels==0.15.0
  installed after the pytest baseline.
- Python 3.13.14. `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1`. 4 CPUs, 15 GB RAM.
- pytest baseline (2026-10-03, default env): **462 passed, 1 failed** (pre-existing
  `test_upo_stability.py::test_hybrid_period1_stability_matches_analytic_henon_multipliers[1-prl_norm]`).
  Identical to Phase 10.

## Status
| Step | State | Commit |
|---|---|---|
| main contains Phase 10 report; branch = main `1147f76`; env rebuilt; pytest default 462 + 1 | done | 4655009 |
| Part 1 verify_all.py + VERIFICATION.md: 339 checks, 327 MATCH, 12 MISMATCH | done | 6def4d3 |
| ERRATA.md E1-E12 (none changes a conclusion; later E13 stale dev analysis regenerated, E14 Conclusion 8 NOT corrected: author decision); Phase 9 decision re-tabulated (errata_phase9_decision.py); Phase 10 report text corrected, marked [Erratum En] | done | 3581197 |
| Part 2: pytest default 462+1, AVX-512 off 463; groundrule replicate_check 44/44; prereg_integrity (9/9 intact); KNOWN_ISSUES.md; no raw PhysioNet file in history | done | 408fd9c |
| analysis10 --phase conf rerun regenerates analysis.json byte-identical (317 s) | done | (recorded) |
| Part 3 make_figures.py: F1-F6, S1-S8 (PDF + 300-dpi PNG), tables T1-T9 + S tables, FIGURE_NOTES.md; byte-stable on rerun | done | a7351e3 |
| groundrule_check.sh complete: 462+1 / 463 / replicate 44/44 / replicability.json identical | done | (recorded) |
| E13: Phase 10 dev analysis.json/tables.md regenerated (were stale vs committed dev robust.jsonl) | done | a4c56bc |
| Part 4: reproduce.sh (fast verified: 781 s, all byte-identical except 17 known Phase 2E formatting files), download_data.py + data_manifests (licences ODC-By 1.0 checked on PhysioNet), requirements-lock.txt (--no-deps), CITATION.cff (validated), README | done | 0b3c338 + 3015a21 |
| E14 (Conclusion 8 wording vs MG maxima-interval detections, check X-21): documented, NOT corrected | done | 3015a21 |
| Part 5: docs/RESULTS_INDEX.md, CONCLUSIONS_CHECKLIST.md, AI_USE_RECORD.md, SOURCES_USED.md | done | 3015a21 |
| Provenance: pinned mitdb manifest = Phase 7 DATA_MANIFEST hashes (144/144) | done | (recorded) |
| FINISH: docs/FINAL_SUMMARY.md (340 checks: 327 MATCH, 13 MISMATCH; E1-E13 corrected, E14 for the author; figures, tables, reproduction, open items) | done | 851ba36 |
| Tag `paper-v1`: created locally on the final commit; `git push origin paper-v1` refused by the session's git proxy (HTTP 403; branch pushes allowed, tag pushes not). The author pushes or creates the tag | blocked | - |

## Resume
- Environment: see above (statsmodels 0.15.0 installed after the pytest baseline, as in Phase 10).
- Final refinement complete. Waiting on the author: the tag `paper-v1` on the last commit of this branch (tag push
  refused from this session), E14 (Conclusion 8 wording), licence choice, open items in docs/FINAL_SUMMARY.md 5.
  Not merged into main; no pull request.
- `bash reproduce.sh fast` rewrites Phase 2E tables/plots with formatting-only differences; restore them with
  `git checkout -- experiments/phase2e` before committing.
- Re-run everything from stored results: `python -m experiments.final_refinement.verify_all`,
  `python -m experiments.final_refinement.errata_phase9_decision`, `python -m experiments.final_refinement.prereg_integrity`,
  `python -m experiments.final_refinement.make_figures` (all < 1 min, run from the repository root).
- groundrule_check.sh needs before-files: `replicate_check --save <file>` and a copy of
  `experiments/phase2e/results/replicability.json` after `run_phase2e --replicate` (then `git checkout` that file);
  pass them as REP_BEFORE / REPJSON_BEFORE.
