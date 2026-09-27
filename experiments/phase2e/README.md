# Phase 2E — synthetic end-to-end validation of the frozen So et al. UPO pipeline

Synthetic validation only. No MIT-BIH or real ECG data is used, and nothing in
`final_pipeline.py` is modified. The findings are in
[`docs/PHASE2E_SYNTHETIC_VALIDATION.md`](../../docs/PHASE2E_SYNTHETIC_VALIDATION.md).

## Reproduce

From the repository root:

```bash
python -m experiments.phase2e.run_phase2e --experiments all --workers 2   # ~3643 tasks
python -m experiments.phase2e.run_phase2e --replicate --workers 2         # Step 25 rerun check
python -m experiments.phase2e.analysis                                     # tables, audits, plots
```

The runner can be resumed. Tasks already present in `results/<experiment>.jsonl`
are skipped unless `--no-resume` is given. BLAS is pinned to one thread
(`OMP_NUM_THREADS=1`), so reruns are bitwise reproducible.
`results/manifest.json` records the git HEAD, the working-tree status, the SHA-256
of `final_pipeline.py` and of every Phase 2E source file, and the library versions.

## Files

| File | Content |
|---|---|
| `config.py` | Every predeclared parameter: systems, seeds, window lengths, SNR levels, surrogate counts, localization grid |
| `systems.py` | Generators and analytic references (Hénon fixed points and 2-cycle, skewed-Hénon PRE Table I values, Jacobians) |
| `metrics.py` | Result-row contract, built from `analyze_segment` / `run_upo_analysis` / `upo_summary_features` |
| `run_phase2e.py` | Experiment runner (multiprocessing, checkpointed JSONL) |
| `analysis.py` | Tables A–I, audits, plots |
| `results/*.jsonl` | Raw output: one line per task (task spec, full row, per-peak detail, surrogate W, histograms) |
| `results/*.csv` | Flat copy of every row. Failed windows are kept with their status and NaN fields |
| `results/tables/` | Tables A–I and supplementary tables (`.md` + `.csv`) |
| `results/audits.json` | Significance recomputation, status consistency, stability, flag independence, replicability |
| `plots/` | 300-dpi PNG figures |

## Experimental arms

| Arm | What varies | Design |
|---|---|---|
| `core` | 7 systems × {128, 256, 512} | one_sample, 50 surrogates. Nulls: 50/100/50 seeds; deterministic systems: 30 seeds; constant: 3 |
| `noise` | logistic, Hénon × SNR {30, 20, 10, 5, 0} dB | 256 samples, seeds 0–29, **paired** with the clean core windows |
| `tau_step` | 7 systems | 256 samples, **paired** with the core seeds (PROJECT EXTENSION map) |
| `surrogate_sens` | 20 and 100 surrogates | 256 samples, **paired** with the core seeds (50 = core) |
| `algo_seed` | `random_seed` = CFG + 1000 + seed | null systems, 256 samples, 50 seeds |
| `oracle_dim` | fixed d = 1 (logistic), d = 2 (Hénon) | **Diagnostic only**, never substituted for production results |
| `fp_sens` | nextafter; two RR-quantization routes | detection without significance |
| `henon_loc` | bins × tube/slab percentile × range policy × period | standard and skewed Hénon, N = 1024, d = 2, M = 2, 3 realizations |
| `henon_loc_oracleJ` | analytic Jacobian substituted | **Diagnostic only**: separates Jacobian-estimation error from transform/histogram error |
| `stability` | source stability disabled at runtime | checks that the stability diagnostic never changes Level A |
| `flag_independence` | `compute_surrogates` False/True | checks that UPO significance does not depend on the LLE surrogate flag |

### Systems (observable passed to the pipeline)

| System | Definition |
|---|---|
| constant | `c ~ U(0.5, 1.5)` |
| white_noise | `N(0, 1)` iid |
| ar1 | `x_n = 0.8 x_{n-1} + e_n`, `e_n ~ N(0,1)`, stationary start, 1000-sample transient |
| sinusoid | `sin(2πn/7.3 + φ0)`, `φ0 ~ U(0, 2π)`. The sequence repeats exactly every 73 samples |
| logistic_p4 | logistic `r = 3.5` (stable 4-cycle), `x0 ~ U(0.05, 0.95)`, 1000-iteration transient |
| logistic | logistic `r = 4`, `x0 ~ U(0.05, 0.95)`, 1000-iteration transient |
| henon | `x' = 1 − 1.4x² + y`, `y' = 0.3x`, `(x0, y0) ~ U(−0.1, 0.1)²`, 1000-iteration transient; observable `x` |
| skewed_henon | PRE Eq. 31 (`b1 = 0.965`, `b2 = 0.25`), localization study only |
| noisy_* | clean window + white Gaussian noise, `σ² = var(clean)/10^(SNR/10)` (observational noise) |

Seeds: the data for `(system, window_length, seed)` come from
`SeedSequence([20260926, system_code, window_length, seed])`. Each window is its
own trajectory. Noise uses a separate child stream per `(seed, SNR)`. Every row
records `seed` and the derived `data_seed`.

**Documented dependence.** The production pipeline seeds the So random matrices
and the surrogate phases from `config.random_seed`, so those random streams are
the same for every window. The `algo_seed` arm checks that the null
false-positive rate does not depend on that fixed choice.
