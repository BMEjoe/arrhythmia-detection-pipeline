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
| main contains Phase 10 report; branch = main `1147f76`; env rebuilt; pytest default 462 + 1 | done | (this commit) |

## Resume
- Environment: see above. Then `bash experiments/final_refinement/reproduce.sh fast` (once written).
