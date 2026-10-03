# Phase 10 handoff (living document; update after every step)

Read this first when resuming. Branch `claude/blissful-brahmagupta-th4tdh`. Never merge into
main, never open a pull request. Commit and push after every step.

**Branch base.** origin/main at `71a4f7c` (Phases 2E-9). Confirmed on 2026-10-03 that main
contains `docs/PHASE9_WAVEFORM_NOISE_ROBUST.md` and `final_pipeline.masked_growth_chaos_test`
(line 2309). The session branch started at the same commit.

## Task (user, Phase 10, FINAL QUESTION PHASE) in brief
Answer only Questions 1-4; new ideas go to "Future work".
- **Detectors (frozen, no retuning).**
  - K1 = Phase 9 baseline `k1_frozen512` = Phase 6 combined detector
    (`keep_upo_on_short_lle_embedding=True`, detrend D2), 512 intervals
    (`experiments.phase7_mitbih.detector.evaluate`). Components LLE alone, UPO alone secondary.
  - K3 = `final_pipeline.masked_growth_chaos_test` with frozen `PipelineConfig.masked_growth_*`.
- **Data split.** DEVELOPMENT = MIT-BIH, nsrdb, chfdb (already analysed). CONFIRMATORY = nsr2db +
  chf2db: NOT downloaded until PREREGISTRATION.md is committed + pushed. After download: subject
  overlap check vs nsrdb/chfdb (docs + record metadata), exclude overlaps; check annotation
  sampling rate; if coarser than 1/360 s, quantize every Q1 spike-in / synthetic check at it.
- **Q1** detection limits: additive spike-in c(n) (Henon, logistic over a lambda grid computed from
  the equations; Phase 8 coupled-vdP + phase-reset chaotic regimes) at chaos fractions f; Phase 7
  replacement design secondary; p(f, lambda, family) with 95 % CIs, predeclared smooth curve;
  exclusion bound pi_upper(f, lambda) from real detection rate (subject-cluster bootstrap) and FP.
- **Q2** faithful Poon-Barahona noise titration (VWK, information criterion, noise increments,
  noise limit) from published sources; VERIFY on published examples (stop titration claims if it
  fails); study data/preprocessing of key positive studies; titration by group on confirmatory
  data; again after K3-style masking and Phase 7 editing; synthetic FP checks (Phase 5-6 nulls +
  ectopy; Phase 8-9 non-chaotic coupled vdP with ectopy).
- **Q3** ectopy dose-response (GEE or mixed logistic, predeclared): 3a burden effect for titration,
  LLE, UPO, K1, K3; 3b masking/editing removes positives (paired, cluster bootstrap); 3c group
  difference in titration adjusted for burden.
- **Q4** robustness (secondary): 499 surrogates; 256/512/1024 intervals; night vs day; K3 with
  annotation labels vs K4-style RR-rule masking.
- **Step P** PREREGISTRATION.md (grids, formulas, models, window sampling, exclusions, tests,
  outputs; primary P1 exclusion curves K3 + K1; P2 titration rate by group raw + change after
  masking; P3 burden effect for titration and K3; runtime <= ~20 h on 4 workers; reduce spike-in
  grid before real windows). Runner guard. Later changes only as dated amendments.
- **Report** docs/PHASE10_FINAL_RESULTS.md incl. "Conclusions supported by Phases 2E-10" and
  "Future work".

## Ground rules (as Phases 3-9)
- final_pipeline.py: opt-in options only, defaults unchanged; then
  `bash experiments/phase3_lle/groundrule_check.sh` (needs the before-files, see Phase 9 HANDOFF).
- New code in experiments/phase10_final/. Methods from primary published sources only; sources
  in METHODS.md.
- Seeds: Phase 10 synthetic seeds start at 10000 (never reuse earlier seeds).

## Environment (rebuilt 2026-10-03 in a fresh container)
- Recipe in `requirements-phase10.txt` (= Phase 9 environment).
- Python 3.13.14, numpy 2.1.3, scipy 1.18.1, scikit-learn 1.9.1, pandas 3.0.6,
  matplotlib 3.11.2, wfdb 4.3.1, numba 0.68.0. OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1.
- pytest baseline (2026-10-03): default env **462 passed, 1 failed** (pre-existing
  `test_upo_stability.py::test_hybrid_period1_stability_matches_analytic_henon_multipliers[1-prl_norm]`);
  AVX-512 disabled **463 passed**. Identical to Phase 9.
- Machine: 4 CPUs, 15 GB RAM.
- Literature access (2026-10-03): PMC / NCBI eutils, PLOS, nature.com landing pages, physionet.org
  reachable; AIP (Chaos) 403.

## Status
| Step | State | Commit |
|---|---|---|
| prerequisites on main confirmed; env rebuilt; pytest baseline 462+1 / 463 | done | (this commit) |

## Notes
- Never `pkill -f <pattern>` with the plain pattern (kills the calling shell).
