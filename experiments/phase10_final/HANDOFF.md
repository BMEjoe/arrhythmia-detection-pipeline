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
| prerequisites on main confirmed; env rebuilt; pytest baseline 462+1 / 463 | done | 5c0c5b1 |
| Q2a titration10.py from accessible sources (BP96/PM97 paywalled); readings tried; criteria predeclared | done | e5131dc |
| Q2b verify_titration.py: **VERIFIED** (V1-V6 pass; METHODS 2.6) | done | (this commit) |
| Q2c sources: Wu09 + PM97 used nsrdb/chfdb (dev); 12-min (~800 beat) segments; NN without interpolation | done | 7c988b7 |
| dev data (mitdb, nsrdb, chfdb; SHA OK) + harness10 + run10 (guarded) + analysis10 + overlap10 + figures10 | done | afd8296 |
| DEV run real (485), synth (410), seg (33 subjects) done; spike 96/528 and robust 3/132 done, PAUSED (resume after conf) | partial | a74aa99, 05e86f4 |
| PREREGISTRATION.md pushed BEFORE download | done | 8fc7bd9 |
| conf download (guarded): nsr2db 54, chf2db 29, SHA OK; 128 Hz; no header start times | done | |
| overlap10: original rule flagged 13 chance pairs -> Amendment 1 (time-consistent rule, validated) -> no exclusion | done | 4795cf6 |
| CONF run complete (real 988, seg, synth, spike, robust) | done | f979faa |
| conf analysis + sensitivity (flagged excluded) + tables + figures; Amendment 2 (GEE independence fallback) | done | 965a586 |
| docs/PHASE10_FINAL_RESULTS.md assembled (sections 1-14) | done | (this commit) |
| DEV spike (96/528) + robust (36/132): stopped at session end (resumable); partial, reported as such in report 9.1 | partial | |

## Resume
- Phase 10 is COMPLETE except the optional dev spike/robust completion. To finish it: rerun the dev command below, then `analysis10 --phase dev`, `tables10 --phase dev`, `figures10 --phase dev`, and update report 9.1.
- CONF run is resumable (same command; finished ids skipped):
  `cd /home/user/arrhythmia-detection-pipeline && nohup env OMP_NUM_THREADS=1 PYTHONPATH=. /root/venv313/bin/python -m experiments.phase10_final.run10 --phase conf --part all --workers 4 >> experiments/phase10_final/results/conf/run_conf.log 2>&1 &`
  (conf data: `run10 --phase conf --download` (guarded); exclusions.json already written by `overlap10 --amended`).
- After conf: `analysis10 --phase conf`, `analysis10 --phase conf --sensitivity-flagged`, `figures10 --phase conf`;
  then resume the DEV run (spike, robust) for the Q4 / Q1 development comparison.
- Never `pkill -f` a pattern that also appears later in the same shell command (it kills that shell).
- Dev run is resumable: re-run the same command; finished task ids are skipped:
  `cd /home/user/arrhythmia-detection-pipeline && nohup env OMP_NUM_THREADS=1 PYTHONPATH=. /root/venv313/bin/python -m experiments.phase10_final.run10 --phase dev --part real,synth,seg,spike,robust --workers 4 >> experiments/phase10_final/results/dev/run_dev.log 2>&1 &`
- Dev data: `python -m experiments.phase10_final.data10 --download dev` (gitignored dev_data/).
- Then `python -m experiments.phase10_final.analysis10 --phase dev` and `figures10 --phase dev`.
- Dev-run note: the dev seg/synth parts were started with code before two small additions (seg `t_start`,
  synth `TIT_raw_q128`); rerun synth for dev if the q128 arm is needed.

## Notes
- Never `pkill -f <pattern>` with the plain pattern (kills the calling shell).
