# Phase 2E handoff — state as of 2026-09-26 ~21:25 UTC

Read this file first if you are resuming Phase 2E in a new session.

## Rules that still apply
- The task is Phase 2E synthetic validation of the FROZEN So et al. UPO pipeline.
- No MIT-BIH or real ECG data.
- Do not modify `final_pipeline.py` unless a genuine bug (category F) is demonstrated.
- Do not change any predeclared grid in `config.py`.
- Do not commit and do not push. The user has not authorized either.
- Baseline: HEAD `dd43fd0`, clean tree, `pytest -q` → **335 passed in 48.07 s**. Phase 1
  definitions (PT chain, RR extraction/correction, TDMI, Cao, Takens, Rosenstein, AAFT) were
  verified byte-identical to `4d9f2d6`.

## How to resume
```bash
git clone <repo> arrhythmia-detection-pipeline && cd arrhythmia-detection-pipeline
git checkout dd43fd0                       # or main, if it still points there
unzip /path/phase2e_backup.zip             # restores experiments/ (code + results)
python -m experiments.phase2e.run_phase2e --experiments all --workers 2   # resumes; skips finished tasks
python -m experiments.phase2e.run_phase2e --replicate --workers 2         # Step 25
python -m experiments.phase2e.analysis                                     # tables, audits, plots
```
Then write `docs/PHASE2E_SYNTHETIC_VALIDATION.md` (Step 31 report) and print the executive
summary. The user asked for this deliverable; do not commit it.

## Progress when this file was written
- 1,988 of 3,643 tasks finished. No harness errors and no `analyze_segment` exceptions.
- Complete: `henon_loc`, `henon_loc_oracleJ`, the 512-sample core windows, and the
  100-surrogate arm (mostly).
- Remaining: 256- and 128-sample core windows, noise, tau_step, surrogate-20, algo_seed,
  oracle_dim, fp_sens, stability, flag_independence.
- `analysis.py` has been tested on the localization data only. The other tables may need
  small fixes once all the data exist.
- Report still to be written: `docs/PHASE2E_SYNTHETIC_VALIDATION.md`. The README links to it.

## Preliminary findings (verify against final data before reporting)
1. **Hénon ~0.01 localization discrepancy = finite resolution, not a bug.**
   - The bin-mean peak location error is 0.001–0.04, always < 0.27 bin widths, and does not
     change monotonically with bin count.
   - It is identical with the exact analytic Jacobian substituted (`henon_loc_oracleJ`).
   - A histogram-free KDE mode of the same tube points is within ~1e-4 (standard Hénon)
     and ~8e-4 (skewed Hénon) of the analytic fixed point.
   - Cause: the pipeline's location estimator (mean of tube points in the peak cell)
     under histogram quantization.
2. **Off-attractor fixed points.** The standard Hénon off-attractor fixed point (−1.131) was
   detected in 0 of 3 realizations at the production point. The skewed-Hénon one (−1.903) was
   detected in 3 of 3.
3. **Period-2 Level-A peak count depends on resolution.** Minimal-period-2 peaks grow from
   ~4 at 15 bins to ~56–78 at 90 bins. The true 2-cycle remains the nearest peak (error
   0.0005–0.03). Many period-2 source peaks are resolution artifacts.
4. **Cao at lag 1 overestimates dimension.** It chose d = 4–6 for Hénon (true 2) and 4 for
   logistic (true 1) at 256–2048 samples. Preserve this; do not force the dimension (the
   `oracle_dim` arm is diagnostic only).
5. **Status semantics.** `embedding_not_saturated` is also returned when Cao E1 is undefined
   because of exact duplicate vectors:
   - logistic r = 3.5 (4-cycle): 30 of 30 at 512
   - RR-quantized (1/360 s) sinusoid in `fp_sens`
   The label conflates "undefined" with "not saturated". This matters for MIT-BIH, where RR
   values are quantized at 1/360 s. Report it as a limitation, not a bug.
   Some white-noise and AR(1) windows also return `embedding_not_saturated` (about 14 so far).
6. **Rosenstein LLE on the sinusoid (period 7.3, exactly 73-periodic).** Nearest-neighbour
   distances are about 1e-14, so the fit is at the floating-point floor (slope ≈ 0.008,
   R² ≈ 0.03). One 128-sample run gave +0.24. Check this in Table H. White noise also gives
   positive LLE (0.01–0.08), so positive LLE is not evidence of chaos.
7. **Significance math.** With n = 50 surrogates, the J < 0.05 rule has nominal size 3/51 =
   0.0588 under exchangeability. W_i also uses a surrogate mean that includes surrogate i,
   a small anti-conservative asymmetry that is source-faithful. Compare these with the
   empirical null false-positive rate.
8. **Shared random streams.** Production seeds the So R draws and the surrogates from
   `config.random_seed`, so every window shares them. The `algo_seed` arm tests this.
