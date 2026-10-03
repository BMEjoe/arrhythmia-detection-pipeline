# Known issues

## 1. One platform-dependent unit test fails on AVX-512 machines (pre-existing; not changed)

`tests/test_upo_stability.py::test_hybrid_period1_stability_matches_analytic_henon_multipliers[1-prl_norm]`

- **Symptom.** Default environment: 462 passed, 1 failed. With numpy's AVX-512 dispatch disabled
  (`NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL AVX512_SPR"`):
  463 passed. Identical in Phases 3–10 and in this final refinement (2026-10-03, Python 3.13.14, numpy 2.1.3,
  scipy 1.18.1, Intel Xeon with AVX-512).
- **Cause** (diagnosed in Phase 3, A2; `docs/PHASE2E_SYNTHETIC_VALIDATION.md` 5.3). For an off-attractor Hénon
  peak, `source_period1_stability` averages about 108 member Jacobians with an arithmetic mean. With AVX-512
  arithmetic one member is an ill-conditioned M = 2 neighbour fit (Frobenius norm 4,078; next largest 22.9),
  which moves the leading modulus from about 2.0 to 17.24 and fails the test's tolerance. Without AVX-512 the
  ulp-level arithmetic differs, the member set has 106 points, the largest norm is 93.4, and the test passes.
- **Scope.** The failing quantity is the arithmetic-mean ("prl_norm") stability estimate of the original UPO
  pipeline. The adopted UPO gate (Phase 4 C2) uses the element-wise median, which is robust to the outlier
  (2.011 in the same case). No preregistered result uses the arithmetic-mean estimate.
- **Handling.** Neither the test nor the pipeline was changed. `groundrule_check.sh` runs the suite both ways.

## 2. Bitwise replicability across machines

Reruns on one machine are bitwise identical (`replicate_check`: 44/44). Against the stored Phase 2E records,
made on a different machine, 13/44 tasks are bitwise identical; the others differ at ~1e-15, amplified only in
ill-conditioned quantities (Phase 2E report 5.3). Surrogate random streams of the LLE and K3 tests are seeded by
a CRC of the input bytes, so a last-bit change in the RR series changes the surrogates and can move p-values and
z statistics (Phase 7 report 12.1: decisions unchanged in 10/10 audited windows, p-values moved by up to 0.5).
A full rerun on another machine should therefore reproduce the counts approximately, not bit for bit. The
figures, tables and verification in this repository are computed from the stored results and are exact.

## 3. Frozen Phase 9 decision script groups KTz regimes by name (erratum E1)

`experiments/phase9_waveform/decide9.py` (frozen by SHA-256 in the Phase 9 preregistration) tabulates 81
conditions instead of the preregistered 101. Rerunning it reproduces the 81-condition table. The corrected
tabulation is `python -m experiments.final_refinement.errata_phase9_decision`. Only k2's failed-condition count
differs (16 instead of 12); no decision changes (`ERRATA.md`).

## 4. Partial development runs (Phase 10)

The Phase 10 development spike-in (96 of 528 tasks) and robustness (36 of 132 tasks) runs were stopped at the
end of Phase 10 and are reported as partial (`docs/PHASE10_FINAL_RESULTS.md` 9.1). They are resumable but were
not completed here, because completing them would be new analysis. The development `analysis.json` was written
before the Amendment 2 fallback code existed and carries no `working_correlation` field; every development GEE
converged with the exchangeable working correlation (verified in `VERIFICATION.md`, 10-D-2).

## 5. Exchangeable GEE for the titration outcome

The preregistered exchangeable GEE diverges for the confirmatory titration P3 model (and the segment model);
the independence working correlation was used (Amendment 2). An independent implementation also diverges
(`VERIFICATION.md`, 10-P3-A2). Models whose exchangeable fit converged kept it, so the adjusted and unadjusted
titration ORs in PHASE10 Section 7 use different working correlations (`ERRATA.md` O3).

## 6. Data are not in the repository

No raw PhysioNet file is committed (checked over the full git history: no `.dat`, `.hea`, `.atr`, `.ecg`, `.q1c`
file was ever added). Results, manifests and checksums are committed; the data are downloaded with
`experiments/final_refinement/download_data.py` (SHA-256 verified). Committed derived files from MIT-BIH are
detected R-peak sample indices (`phase7_mitbih/results/qc/detected_peaks.npz`), window boundaries
(`results/qc/*.csv`) and detector timing offsets (`results/spike_in/v_offsets_samples.npy`).
