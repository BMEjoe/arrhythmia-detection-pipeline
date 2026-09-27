# Phase 4 proposal: UPO detector issues (diagnosis and candidate fixes)

Proposal only. No UPO code was changed. Every statement is tied to a
Phase 2E results file, table or post-hoc diagnostic
(`experiments/phase2e/results/...`), or to a Phase 3 diagnostic
(`experiments/phase3_lle/results/dev/...`). Per-peak numbers not already in a
table were computed from the committed `core.jsonl` / `noise.jsonl` detail records.
The proposed validation (Section 4) uses the same development/test seed split as
Phase 3.

Terminology as in `docs/PHASE2E_SYNTHETIC_VALIDATION.md`:
- **Level A** = So et al. source peaks;
- **Level B** = Level A with surrogate J < 0.05;
- **Level C** = PROJECT verification gates (residual ≤ 0.05, R² ≥ 0.9,
  support ≥ 2) plus candidate-centred monodromy.

---

## Problem 1: Level B flags the sinusoid in 63–70 % of windows (256 / 512)

### Evidence

- Sinusoid (period 7.3 samples, a linear rotation with a neutral fixed point at
  the centre): at least one Level-B peak in 21/30 (256) and 19/30 (512) windows,
  against 7/100 and 5/50 for white noise (Table A, Table B).
- Level-B peaks are **not at the fixed point**: 0/42 (256) and 0/38 (512) lie
  within 0.05 of the centre. The median |location| is 0.201 (IQR 0.200–0.202)
  on a unit-amplitude signal (`table_peak_locations.md`; `core.jsonl` peak detail).
- The peaks are **neutral, not unstable**. The Level-C candidate-centred
  monodromy gives a leading multiplier modulus of **1.000 (IQR 1.000–1.000)**
  for all Level-B sinusoid peaks at 256 and 512. For Level-B chaotic peaks it is
  3.33 (IQR 3.23–3.55, Hénon 256), 3.67 (Hénon 512), and 9.1 / 15.1 (logistic
  256 / 512) (`core.jsonl`, `extension_multiplier_moduli`). The hybrid source
  estimate gives 0.989 (IQR 0.977–1.028) at 256 for the sinusoid.
- Level C gates 0, 1 and 3 sinusoid peaks, and verified_unstable is 0 at every
  length (Table I). The failing gate is the residual (Table I, `C failed gates`).
- The So test uses AAFT ("Gaussian-scaled phase shuffle") surrogates at the
  observed lag and dimension (`final_pipeline.py::assess_so_significance`). Its
  null hypothesis is "linearly correlated Gaussian noise, monotonically
  rescaled".

### Likely cause

1. **Wrong null for this question.** A strictly periodic signal is not
   linear-stochastic, so rejecting the AAFT null is correct behaviour for the
   test. But it says nothing about unstable periodic orbits. Level B
   answers "is there more density near some transformed point than in a linear
   Gaussian surrogate", not "is there an unstable orbit".
2. **Degenerate geometry.** The sinusoid embeds (Cao m = 6 in 18/30 windows at
   256) as a closed one-dimensional curve in 6-D, so local Jacobian fits are
   rank-deficient. The transformed points then accumulate at reproducible
   off-centre locations (|loc| ≈ 0.20). Their existence is well established;
   their stability is neutral.
3. **Stability is not part of Level B.** The pipeline already computes a
   stability diagnostic per peak, but by design it never removes a peak
   (`detect_so_fixed_points` docstring).

### Candidate fixes

- **P1-a Instability gate on Level B (recommended first).** Require a Level-B
  peak to have a leading multiplier modulus ≥ 1 + δ, with a predeclared δ
  (for example 0.2), from the candidate-centred monodromy (Level C machinery)
  or a robust (median) source estimate. In the evidence above, sinusoid peaks sit
  at 1.000 and chaotic peaks at ≥ 3.2. Risk: the hybrid source estimate is
  outlier-sensitive (the arithmetic mean over member Jacobians; see Phase 2E
  report 5.3). An element-wise median or trimmed mean should be evaluated
  alongside it.
- **P1-b Periodic-orbit-preserving null.** Add a second surrogate test against a
  "periodic orbit plus noise" null: pseudo-periodic surrogates (Small, Yu &
  Harrison 2001) or cycle-shuffled surrogates. A peak is then reported only if
  it is significant against both the linear-stochastic and the periodic null.
- **P1-c Periodicity pre-screen with its own status.** Label windows with a
  dominant spectral line (for example spectral flatness below a threshold, or a
  recurrence-based return-time test) as `periodic`, as a new status. The
  existing statuses stay unchanged, following the Phase 3 A5 pattern of adding
  fields rather than altering statuses.

---

## Problem 2: Cao overestimates the dimension, so every chaotic peak fails Level C

### Evidence

- Cao at lag 1 chooses d = 4 for Hénon (true 2) in 25/30 (256) and 26/30 (512)
  windows, and d = 4 for logistic (true 1) in 30/30 (Table A; Phase 2E report 2.4).
- Level C gates 0 chaotic peaks at production dimension. The residual gate
  fails for every one (Table I).
- `posthoc/levelc.json` (3 seeds, N = 512, peak nearest the analytic fixed point):
  - Hénon at Cao m = 4: residual 0.115–0.689, affine-fit condition number 1,957
    (seed 0), candidate off the local 2-plane by 0.022.
  - Hénon at forced m = 2: residual 0.027–0.050, condition number 54, off-plane
    ≈ 1e-17. All gates pass in 3/3.
  - Logistic at m = 4: residual 6.9–11.9. At m = 1: 0.063–0.073, still above
    0.05. At m = 3: 0.81–1.51.
- With the oracle dimension, `verified_unstable` peaks appear at logistic x = 0
  (30 windows) and x ≈ 0.75 (1), and at Hénon x ≈ 0.63 (7)
  (`oracle_dim.jsonl`; `table_oracle_dimension.md`).
- Phase 3 B2 (`results/dev/tables/B2_diagnosis.md`): Cao at lag 1 gives m = 4
  for both maps, rising to 7–9 under 30–10 dB noise.

### Likely cause

For a 1-D or 2-D attractor embedded in 4-D, the attractor occupies a thin
subset. The Level-C local affine model is fitted with 30 neighbours on that
subset (`verify_neighbors`), which is ill-conditioned in the transverse
directions. The candidate is a diagonal point (c, …, c) that lies **off** the
local manifold. The model therefore extrapolates off-manifold, and the
one-step residual at the candidate is large even when the candidate location is
accurate. For Hénon the candidate error is 0.007–0.010 at m = 2, and at m = 4 it
is 0.013–0.016 in 2 of 3 seeds (−0.071 in the third).

### Candidate fixes

- **P2-a Tangent-space (rank-regularized) Level-C fit.** Fit the local affine
  model in the principal subspace of the neighbourhood, with the rank chosen by
  a singular-value threshold or from a local PCA dimension estimate. Evaluate the
  residual at the candidate projected onto that subspace. This directly targets
  the 1,957 vs 54 condition-number gap.
- **P2-b Dimension choice that does not overestimate.** Keep Cao's E1 but choose
  the smallest m whose E1 exceeds a fixed fraction (for example 0.9) of its
  plateau, instead of requiring two consecutive small changes. Alternatively use
  false-nearest-neighbours with a noise-aware threshold. Level C would then be
  run at that m. This must be validated for noise, where Cao inflates m further.
- **P2-c Dimension-sweep consensus for Level C only.** Run Level C at m = 1..m_Cao
  and accept a peak as verified if it passes at any m ≤ m_Cao, with location
  agreement across dimensions. The multiple testing is controlled by
  predeclaring the sweep. This leaves source detection (Levels A/B) unchanged.

---

## Problem 3: Hénon Level-B detection collapses at 30 dB

### Evidence

- Level-B detection for Hénon falls from 29/30 (clean) to 11/30 (30 dB) and
  8/30 (20 dB). Logistic stays at 29/30 and 27/30 (Table D).
- Hénon at 256 samples, from `core.jsonl` / `noise.jsonl`:

  | | clean | 30 dB | 20 dB |
  |---|---:|---:|---:|
  | median Cao m | 4 | 7 | 8 |
  | median r_J | 6.07 | 2.10 | 1.53 |
  | median error of the nearest Level-A peak to 0.631 | 0.068 | 0.154 | 0.222 |

- Median source coverage rises from 0.107 to 0.215 (30 dB) and 0.655 (20 dB).
  The transformed-point density spreads (Table D).
- Phase 3 B2 (D4): the mean log nearest-neighbour distance at τ = 1, m = 2 rises
  from −4.19 (clean) to −3.68 (30 dB) for Hénon. Neighbour distances are
  therefore set by the noise at 30 dB. For logistic (1-D, clean −7.09) the jump
  is larger, but its fixed point is strongly repelling (multiplier −2) and
  isolated, so its histogram peak survives.

### Likely cause

Two mechanisms act together. Noise raises the Cao dimension (4 → 7), and in
higher dimension the So transform has more poorly-determined directions. Noise
also enters the local Jacobian estimates directly: they use
`so_jacobian_neighbors = 7` nearest spatial neighbours whose separations are now
at the noise scale. The transformed points therefore scatter and the peak's
excess over the surrogate mean shrinks (r_J 6.07 → 2.10). Hénon's saddle
(multipliers −1.92 and 0.16) has a weaker stable direction contribution than
logistic's 1-D repeller, which may explain why it degrades first. This last
point is an **untested hypothesis**.

### Candidate fixes

- **P3-a Noise-scale Jacobian neighbourhoods.** Estimate the Jacobians from
  neighbours beyond an estimated noise scale, or from more neighbours with ridge
  regularization. The noise scale can be taken from the nearest-neighbour
  distance distribution, as in the Phase 3 `eps` rule. This is an opt-in
  replacement for the fixed M = 7.
- **P3-b Decouple the UPO dimension from noisy Cao.** Use the P2-b dimension rule,
  or a fixed low dimension for RR-like data, validated against the oracle
  dimension on noisy data. The existing `oracle_dim` arm was run only on clean
  data.
- **P3-c Nonlinear noise reduction before UPO detection.** Apply local
  projective noise reduction (Grassberger et al. 1993; Kantz & Schreiber) to
  the embedded data, with parameters fixed in advance. It must be applied
  identically to the surrogates, or the surrogate test becomes invalid.

---

## 4. Proposed validation design (same development/test split as Phase 3)

- **Seeds.** DEVELOPMENT = Phase 2E core/noise seeds (0–99 nulls at 256,
  0–49 at 512, 0–29 otherwise), used for diagnosis and tuning. TEST = 1000–1099
  per condition, same SeedSequence scheme, run only after a pushed
  `PREREGISTRATION.md` (reuse `experiments/phase3_lle/run_phase3.py`'s guard).
- **Systems.**
  - nulls: white noise, AR(1);
  - periodic controls: sinusoid, logistic r = 3.5, plus a quasi-periodic
    two-tone signal as a new harder control;
  - chaotic: logistic, Hénon, skewed Hénon (PRE Table I references);
  - noisy logistic / Hénon at 30, 20, 10 dB.

  Windows 256 (primary) and 512.
- **Candidates.** At most 4, each a combination of the fixes above with every
  parameter fixed, plus the unchanged production detector as baseline.
- **Primary rule (fixed in advance).** On test seeds at 256, a candidate PASSES
  if the Level-B (or gated Level-B) rate on white noise, AR(1) **and the
  sinusoid** each has a 95 % Wilson interval containing or below α. Among
  passing candidates, the WINNER has the highest rate of verified detections
  (Level B plus the candidate's gate) on Hénon at 30 dB. Tie-break: lowest
  median localization error to the analytic fixed point on clean Hénon.
- **Secondary outcomes.**
  - Level-C verified rate at production dimension;
  - localization error (bin-mean vs KDE mode);
  - logistic r = 3.5 status distribution, with `cao_e1_undefined`;
  - 20 / 10 dB detection and 512 results;
  - runtime.
- **Ground rules.** Every change goes into `final_pipeline.py` as an opt-in
  `PipelineConfig` option, default off. The full pytest suite and the
  same-machine replicability check
  (`experiments/phase3_lle/groundrule_check.sh`) must pass unchanged.
- **Budget note.** Surrogate UPO significance costs about 20–40 s per 256-sample
  window at 50 surrogates (`stability` / `flag_independence` runtimes in
  `audits.json`). A full test grid of about 1,500 windows per candidate is
  therefore several CPU-hours per candidate and needs to be estimated before the
  surrogate count is preregistered.
