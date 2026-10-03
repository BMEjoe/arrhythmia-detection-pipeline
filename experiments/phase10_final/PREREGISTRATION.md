# Phase 10 preregistration: detection limits, faithful noise titration, ectopy dose-response, robustness

This file is committed and pushed **before nsr2db or chf2db is downloaded**. No file of either
database (other than the PhysioNet documentation pages, read to plan the data handling) has been seen in
this project. The runner (`run10.py --phase conf`, including `--download`) refuses to start unless this
file is committed, unmodified and pushed (Phase 3 guard). Any later change is a dated amendment at the end
of this file stating whether any confirmatory output had been seen.

All code that produces the confirmatory results is frozen with this file: `data10.py`, `harness10.py`,
`titration10.py`, `run10.py`, `overlap10.py`, `analysis10.py` (commit of this file). Development (MIT-BIH,
nsrdb, chfdb) results are in `results/dev/` and are reported as DEVELOPMENT only.

Methods listed for the runner guard:
- method: `K1`
- method: `K3`
- method: `K4`
- method: `K3RR`
- method: `TIT`

## 0a. Development results seen before this file (DEVELOPMENT data only; reported as such)
- Titration verification: VERIFIED (METHODS.md 2.6).
- nsrdb + chfdb, 12 windows of 512 per subject (396 windows): TIT raw positive 199, LLE 126, UPO 8, K1 3,
  K3 0, K4 0, K3RR 1; TIT GEE OR per doubling of 1 + burden 2.43 (1.24-4.78); LLE 1.33 (1.16-1.52); masking
  removed 68 % of TIT positives among windows with masked intervals (analysable in both arms).
- 12-min segments (Wu-style), 33 subjects: TIT raw DR NSR 0.48, CHF 0.81; masking changed CHF by -0.58 and
  NSR by -0.02; Wu-preprocessed DR NSR 0.46, CHF 0.36; NSR night 0.70 vs day 0.41.
- Synthetic (5 development seeds per condition): TIT positive on N5 (static warp), couplets, runs, forced
  vdP (5.45, 5.6) without ectopy; 0 on the linear nulls and SETAR.
- Spike-ins on the first development subjects (reduced f grid): K3 detected Henon / logistic / vdP at
  f >= 0.5-0.9 in some base windows; phase-reset never; K1 only at f = 0.9 for some maps.
None of these changed a detector, the titration implementation, or any parameter.

## 0. What the analyses can and cannot show

- An **additive chaotic component** (Question 1) is a model of *hidden chaos* added to real human RR
  variability. It is **not a model of a chaotic heart**: the real variability stays, and the chaotic part
  is an independent low-dimensional signal. The exclusion bounds therefore refer only to this model.
- The titration indicator detects nonlinear predictability; a positive noise limit is "chaos" only under
  the source's own assumptions (autonomous, stationary, no ectopy artefacts). Phases 5-9 found it positive
  on ectopy, static nonlinearity and non-chaotic deterministic rhythms.
- Real-data detections by K1 / K3 are not proof of chaos (Phase 9 Part F limitations).

## 1. Data

- CONFIRMATORY: PhysioNet nsr2db 1.0.0 (Normal Sinus Rhythm RR Interval Database, 54 records) and chf2db
  1.0.0 (Congestive Heart Failure RR Interval Database, 29 records). Files: RECORDS, `*.hea`, `*.ecg`
  (beat annotations), SHA-256 checked against PhysioNet's SHA256SUMS.txt (`data10.download`). One record
  per subject (documentation). Group NSR = nsr2db, CHF = chf2db.
- Documentation (read before this file): nsr2db = 30 men (28.5-76 y) and 24 women (58-73 y), contributed by
  Washington University (P. Stein) and Columbia-Presbyterian (R. Goldsmith); chf2db = CHF NYHA I-III, 34-79 y,
  contributed by Columbia-Presbyterian; both "digitized at 128 samples per second, beat annotations obtained
  by automated analysis with manual review and correction". The development nsrdb (18 subjects, 20-50 y) and
  chfdb (15 subjects, NYHA 3-4) come from Beth Israel Hospital, Boston.
- **Checks before any detector run** (`overlap10.py`, output `results/conf/data_checks.json`):
  1. **Resolution:** the annotation sampling rate of every record. If it is coarser than 1/360 s (expected
     1/128 s), every Question 1 spike-in is quantized at that record's resolution (the code always quantizes
     at the record's own annotation rate).
  2. **Metadata:** header fields (fs, length, start time, comments) tabulated next to nsrdb / chfdb.
  3. **Subject overlap** with nsrdb / chfdb: documentation (above) + metadata (age / sex where available) +
     RR-sequence matching (5 probes of 500 intervals per confirmatory record matched against every
     development record at every lag; a probe matches iff the median absolute interval difference at the
     best-correlation lag is < 2/128 s; a pair overlaps iff >= 2 probes match). Every confirmatory subject
     with an overlapping pair, or identified as the same person from documentation / metadata, is
     **excluded** from every confirmatory analysis (`results/conf/exclusions.json`).

## 2. Beats, labels, windows, segments

- **Beats** are annotations with a wfdb beat symbol (N L R B A a J S V r F e j n E / f Q ?); others ignored.
- **Label** 'N' iff symbol in {N, L, R, B, e, j}; any other beat symbol is non-normal; an interval < 0.25 s
  or > 2.5 s also marks its ending beat non-normal (Phase 9 Part F rule).
- **Ectopy-related interval (mask)**: interval k (beat k -> k + 1) iff beat k or beat k + 1 is non-normal
  (K3 rule = Phase 7 editing rule).
- **Ectopy burden** of a window: number of its 513 beats whose annotation symbol is not in
  {N, L, R, B, e, j} (symbol-based).
- **Raw RR**: all intervals between consecutive beats, ectopy included.
- **Windows (Questions 1, 3, 4):** per subject 12 windows of 512 intervals. Window j = 1..12 nominally starts
  at the first beat at or after T0 + (j - 0.5)/12 (T_end - T0) (T0, T_end = first and last beat). If its
  512 intervals contain one outside [0.25, 3.0] s (the pipeline's hard RR limits: data gaps, spurious
  annotations), the start moves to the beat after the last such interval, repeatedly; the window must end
  before the next window's nominal start, else it is missing (`harness10.window_starts`). Windows of 256 and
  1024 intervals (Q4b) use the same rule with their own length.
- **Clock time:** window start clock time = header start time + elapsed time. Night = clock hour in
  [00:00, 05:00). If a record has no start time, its windows are labelled by elapsed time with an assumed
  start of 09:54 (the median start time of the 33 development long-term records, nsrdb + chfdb); such
  results are labelled an approximation and also reported excluding those records.
- **12-minute segments (Question 2):** consecutive 12-min segments from the first beat (Wu et al. 2009);
  intervals outside [0.25, 3.0] s are removed in every arm; a segment is analysable iff the remaining
  intervals cover >= 90 % of its 12 min and >= 400 remain.
- **Edited series** (Phase 6 `edit_nn`): masked intervals replaced by linear interpolation between NN
  intervals; eligible iff NN fraction >= 0.80 (Phases 6-7), else "not eligible".

## 3. Detectors and titration (frozen)

- **K1** = Phase 9 `k1_frozen512` = `experiments.phase7_mitbih.detector.evaluate(rr)` (Phase 6 combined AND
  detector, `keep_upo_on_short_lle_embedding = True`, D2 detrend; 99 LLE + 50 UPO surrogates) on raw RR.
  Components LLE alone (`lle_detected`) and UPO alone (`upo_detected`) are secondary. K1 edited = K1 on the
  edited window.
- **K3** = `final_pipeline.masked_growth_chaos_test(rr, labels, CFG)` (frozen `masked_growth_*`; not
  analysable counts as not detected).
- **K4** = Phase 9 candidate `k4_mnlp_growth_rr` (RR-rule mask, frozen Z = 12.1603, G = 0.85); **K3RR** = K3
  with labels from the RR rule (the beat ending an RR-rule-flagged interval marked non-normal). Q4d only.
- **TIT** = verified noise titration (`titration10.py`, METHODS.md 2; kappa <= 6, d <= 3, linear memory <= 83,
  variance-ratio F-test at 1 %, 10 noise realisations, bisection). Positive iff NL > 0. Arms: raw; masked
  (regression rows whose span contains a masked interval are not used; >= 300 rows else not analysable);
  edited (eligible windows only); Wu (12-min segments only: masked intervals eliminated without
  interpolation, NN intervals concatenated, >= 400 NN intervals).
- RNG streams: K1, K3 as frozen (seeded by a CRC of the data); titration seeded by (20261010, CRC of data).

## 4. Question 1: detection limits

**Base windows:** one per subject: window j = 1 + (s mod 12) where s is the subject's 0-based index in the
database's RECORDS file (the next available j cyclically if missing).

**Chaotic component c(n)** (standardised to mean 0, SD 1 over the window), 16 parameter sets:

| family | parameters | Lyapunov exponent (computed from the equations / Phase 8 ground truth) |
|---|---|---|
| Henon x (b = 0.3) | a = 1.08, 1.14, 1.22, 1.40 | 0.136, 0.244, 0.303, 0.419 per iteration (beat) |
| logistic | r = 3.58, 3.65, 3.88, 4.00 | 0.105, 0.255, 0.464, 0.693 per iteration (beat) |
| coupled vdP (Phase 8 chaotic regimes) | (rho, omega) = (2, 2.7), (6, 3.3), (6, 4.0), (8, 3.3) | 0.085, 0.079, 0.101, 0.075 per model time unit |
| phase-resetting map (Phase 8) | tau = 1.14, 0.58, 1.20, 1.16 | 0.035, 0.081, 0.138, 0.185 per stimulus |

Henon / logistic exponents: mean log |Jacobian growth| over 4e5 iterations after 2,000 transient; each
parameter is robustly chaotic (exponent > 0 at +/- 0.001 and +/- 0.002). Initial conditions and the
Phase 8 model RNG are keyed by (database, subject index, window, family, parameter); 2,000 transient
iterations for the maps; Phase 8 models as in Phase 8 (`series.model_rr`).

**Additive design (PRIMARY):** x'_k = x_k + s c_k for every interval k of the raw window, with
s^2 = f / (1 - f) * var_NN, var_NN = variance of the window's non-masked intervals, so that
f = var(chaos part) / (var(chaos part) + var(real NN part)); f in {0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9}. Beat
times t0 + cumsum(x') are rounded to the record's annotation resolution (1/fs). Labels unchanged (K3 masks
the real ectopy positions; K1 sees the real ectopy, as each detector handles it on real data).
**Replacement design (secondary, Phase 7):** the window replaced by c rescaled to its NN mean and SD (SD
reduced if an interval would fall below 0.25 s), the real ectopic beats inserted at their real positions
with the Phase 6 Part B model (Phase 7 `spike_in.insert_ectopy`), quantized to 1/fs.

**Detection probability** p(f, lambda, family) per detector (K3, K1; LLE alone and UPO alone secondary):
empirical rate with Wilson 95 % CI at each design f, and a smooth curve per (family, parameter):
**Firth-penalised logistic regression logit p = b0 + b1 f** on the additive design (finite under
separation), with pointwise 95 % bands and a one-sided 95 % lower band from 2,000 subject-cluster bootstrap
resamples, evaluated on f = 0.05, 0.06, ..., 0.90.

**Real-data detection rate bound:** for each detector, on all confirmatory windows (both groups pooled),
U = max(subject-cluster bootstrap 95th percentile of the mean subject-level detection rate (10,000
resamples), one-sided 95 % Clopper-Pearson upper bound of detected / all windows; not analysable = not
detected, as in the spike-ins). The CP term keeps the
bound informative when there are no detections (Phase 7 lesson: a zero-width bootstrap interval is not
precision).

**Exclusion bound (PRIMARY formula; model-free):** at each design f,

    pi_upper(f, lambda) = min(1, max(0, U - L) / p_emp(f, lambda)),   L = 0,

with p_emp = detected / base windows at that f (pi_upper = 1 when p_emp = 0): the smallest prevalence of
windows containing chaos at (f, lambda) that would be expected to produce more detections than the 95 %
bound allows (expected rate >= pi p). L = 0 because a lower bound on the synthetic false-positive rate is
not a valid lower bound on the real false-positive rate (real data need not contain the ectopy-heavy
synthetic conditions). The empirical p is used, not the fitted curve, because the Firth penalty adds
pseudo-detections: in development, families with zero detections at every f had fitted p of 0.04-0.07,
which would manufacture exclusions from data containing no detection.
**Secondary:** (i) p replaced by its Wilson 95 % lower bound; (ii) the Firth curve p_hat on the fine grid
(and its one-sided 95 % lower band), descriptive only; (iii) L = one-sided 95 % Clopper-Pearson lower bound
of the detector's pooled Phase 9 TEST false-positive rate over the 101 PASS conditions (K3 0/10,100 -> L = 0;
K1 97/10,100); (iv) per group (NSR, CHF).

**P1 (PRIMARY OUTCOME 1):** for K3 and for K1, the exclusion curves pi_upper(f) per (family, parameter) at
the design f values, and for each the smallest design f with pi_upper < 0.05 and with pi_upper < 0.20 at
that f and at every larger design f (or "no such f up to f = 0.9"). Per family, also the largest of these
over its parameters ("for every lambda tested").

## 5. Question 2: faithful replication of noise titration

**Verification (done before this file; METHODS.md 2.5-2.6):** VERIFIED. If the implementation is changed,
titration claims are void.

**P2 (PRIMARY OUTCOME 2), 12-min segments (Wu et al. 2009's segmentation, which they state follows Poon &
Merrill 1997):**
- per subject, the detection rate DR = fraction of analysable segments with TIT raw positive; per group
  the mean of the subject DRs with a subject-cluster bootstrap 95 % CI (10,000);
- group difference CHF - NSR: bootstrap 95 % CI (subjects resampled within groups) and a two-sided
  permutation test of group labels over subjects (10,000); "group difference" iff the CI excludes 0;
- **change after masking**: per subject DR(masked) - DR(raw) on the segments analysable in the masked arm;
  per group (and pooled) the mean change with subject-cluster bootstrap 95 % CI; "masking changes the
  positive rate" iff the CI excludes 0. Also the fraction of raw-positive segments that are negative after
  masking.
- The masked arm needs ~84 consecutive clean intervals per regression row (linear memory up to 83), so it is
  not analysable in ectopy-dense segments; the primary change uses analysable segments only, and the
  same change with non-analysable masked segments counted as negative is reported as secondary, with the
  analysable fractions by group.
- Prediction from Phases 5-9: masking lowers the CHF raw DR (ectopy-driven positives) more than the NSR DR.
- **Secondary:** edited arm (paired change, eligible segments); Wu-preprocessed arm (DR by group: the
  study-faithful replication, "noise-limit positives by group"); mean NL among positive segments (raw, Wu);
  night / day DR; analysable / eligible fractions by group; GEE of segment positivity on burden and group
  (as Section 6).
- **Q2e synthetic checks** (secondary; seeds 10100-10199, 100 windows per condition, 512 and 800 intervals):
  TIT raw / masked / edited (true beat types as labels) on the 17 Phase 5-6 PASS nulls (incl. their ectopy
  patterns; `phase9_waveform.families.null_window`) and the 6 non-chaotic coupled-vdP regimes with ectopy
  none / S2_5 / E1 / E3_10 (`families.model_window`, quantized 1/360 s); raw also after re-quantization to
  1/128 s. Reported as false-positive counts per condition.

## 6. Question 3: ectopy dose-response (confirmatory windows of 512 intervals)

- Outcomes per window: TIT raw positive, LLE alone, UPO alone, K1, K3 (not analysable = negative).
- **Model (predeclared): GEE**, binomial family, logit link, exchangeable working correlation within
  subject, robust (sandwich) standard errors (statsmodels 0.15.0):
  logit P(positive) = b0 + b1 log2(1 + burden) + b2 CHF.
- **P3 (PRIMARY OUTCOME 3):** b1 (odds ratio per doubling of 1 + burden, 95 % Wald CI, two-sided p) for TIT
  and for K3. "Positives track ectopy" iff the CI of b1 excludes 0 with b1 > 0. If a method has < 10
  positive (or < 10 negative) windows the model is not fitted ("not estimable"): the burden table
  (0, 1, 2-4, 5-15, >= 16 beats; counts, Wilson CIs) and the rate difference burden >= 1 vs 0 (subject-
  cluster bootstrap CI) are reported instead, and "positives track ectopy" cannot be claimed.
- Prediction from Phases 5-9: TIT and LLE-alone positives increase with burden; K3 positives do not.
  Reported either way. The prediction for K3 is CONTRADICTED iff K3's b1 is estimable with CI excluding 0
  and b1 > 0; if K3 has < 10 positives it is reported as "not contradicted; K3 positives too rare for a
  dose-response" with the burden table (Wilson bounds per burden bin).
- **3a secondary:** the same model for LLE alone, UPO alone, K1; unadjusted model (burden only); K3 restricted
  to analysable windows.
- **3b paired** (windows usable in both arms; all windows and windows with >= 1 masked interval): TIT raw vs
  masked, TIT raw vs edited, K1 / LLE / UPO raw vs edited: counts positive in each arm and in both, fraction
  of raw positives removed, mean subject-level paired difference with subject-cluster bootstrap 95 % CI.
- **3c:** b2 for TIT (group effect adjusted for burden) vs the unadjusted model logit P = b0 + b2 CHF.

## 7. Question 4: robustness (secondary; supports no new confirmatory claim)

Confirmatory and development data.
- **4a surrogates:** windows j in {1, 4, 7, 10}: K1 with 499 LLE and 499 UPO surrogates, K3 with 499 surrogates
  (`PipelineConfig` fields `lle_chaos_test_surrogates`, `so_surrogate_count`, `masked_growth_surrogates`;
  no change to final_pipeline.py); agreement with the default counts (2 x 2 table, agreement, Cohen's kappa,
  exact McNemar).
- **4b window length:** the same windows at 256 and 1024 intervals (own clean-start rule) vs 512: rates of K1,
  LLE, UPO, K3 by group.
- **4c time of day:** night vs day window counts and rates by detector and group (512); GEE of TIT on burden,
  group and night; 12-min segment DR night vs day.
- **4d beat labels:** K3 (annotation labels) vs K4 and K3RR (RR rule, no annotations): rates by group,
  agreement, burden table.

## 8. Window sampling and budget

- Confirmatory tasks (83 subjects before exclusions): real 996 windows; 12-min segments ~120 per subject;
  spike 83 x 16 tasks (7 f + replacement each); robust 332 windows; synthetic 41 conditions x 2 lengths x 100.
- **Runtime budget** (development timings under 4-worker load; projected wall time on 4 workers):

  | part | dev mean per task | confirmatory tasks | projected wall |
  |---|---|---|---|
  | real (512 windows, all detectors and arms) | 15.4 s | 996 | 1.1 h |
  | seg (12-min segments, 4 titration arms) | 345 s per subject | 83 | 2.0 h |
  | synth (Q2e) | 4.8 s | 8,200 | 2.7 h |
  | spike (Q1; 8 series per task) | ~45 s for 5 series -> ~78 s | 1,328 | 7.2 h |
  | robust (Q4a/4b) | ~100 s | 332 | 2.3 h |
  | **total** | | | **~15.3 h** |
- If the projected total exceeds 20 h, the spike-in grid is thinned first (drop f = 0.05 and 0.7, then the
  second and third parameter of each family), never the real-data windows; any thinning is recorded as an
  amendment.
- Order of the run: real, seg, synth, spike, robust (`run10.py --phase conf --part all`), after
  `overlap10.py`.

## 9. Outputs

`results/conf/{data_checks.json, exclusions.json, real.jsonl, seg.jsonl, synth.jsonl, spike.jsonl,
robust.jsonl}`; `results/conf/analysis/analysis.json` (`analysis10.py --phase conf`); tables and figures
for `docs/PHASE10_FINAL_RESULTS.md`.

## Amendments

### Amendment 1 (2026-10-03, after the download and the `overlap10.py` checks, BEFORE any detector, titration, HRV or other analysis of confirmatory data)

**Seen at this point:** the SHA-256 checks; annotation sampling rates; header fields; the RR-matching
statistics of `overlap10.py` (`results/conf/data_checks.json`) and the diagnostic below. No confirmatory
output of any detector or titration existed.

1. **Overlap rule.** The preregistered RR-matching rule flagged 13 pairs (10 chf2db subjects: chf203, 205,
   206, 209, 210, 214, 215, 216, 224, 228). All are chance matches, not duplicates:
   - three chf2db records "match" two different chfdb subjects each (impossible for one person);
   - in every flagged pair the matched positions of different probes are mutually inconsistent in time
     (e.g. chf205's probes taken at hours 14 and 18 both matched the same ~3-min stretch of chfdb chf09);
   - the matching probes have very low variability (SD 0.006-0.05 s), so a median difference of 2 samples at
     1/128 s is reachable by chance among ~1e5 lags; the rule did not account for this;
   - header metadata contradicts identity (e.g. chf205 is 39 y, M; it "matched" chf09, 63 y, F, and chf13,
     61 y, M); documentation: different institutions.

   **Amended rule** (`overlap10.consistent_overlap`, `--amended`): a pair overlaps iff >= 3 of the 5 probes match
   (same median-difference criterion) AND share one time offset within 120 s. Validated on development data
   before being applied: re-annotated duplicates of 4 development records (shifted 37 min, re-quantized at
   1/128 s with +/- 1 sample jitter) were all detected (5/5 consistent probes), and no distinct pair was flagged
   (<= 1 probe). Applied to all 2,739 confirmatory x development pairs: the largest number of time-consistent
   matching probes in any pair is 1; **no subject overlaps; no exclusion** (`results/conf/exclusions.json`,
   `data_checks_amended.json`).
   **Added secondary sensitivity analysis:** P1-P3 are repeated with the 10 originally flagged chf2db subjects
   excluded (`analysis10.py --phase conf --sensitivity-flagged`).
2. **Start times.** No nsr2db / chf2db header has a start time. As preregistered, every confirmatory window
   and segment is placed in clock time with the assumed 09:54 start; ALL confirmatory night / day results are
   therefore approximations (the preregistered "excluding those records" comparison is empty).
3. **Resolution.** Every confirmatory annotation file is at 128 Hz (1/128 s, coarser than 1/360 s); the Question
   1 spike-ins are quantized at 1/128 s as preregistered (automatic in the code).
4. Headers carry age, sex (chf2db: 21 unknown) and NYHA class (chf2db: I-III), reported in the data table.

### Amendment 2 (2026-10-03, after the confirmatory real-window detector run; Question 3 model fitting)

**Seen at this point:** the confirmatory 512-interval window results (`results/conf/real.jsonl`) and the
first run of `analysis10.py` on them, in which the preregistered exchangeable GEE for the titration outcome
(P3) returned non-finite estimates: statsmodels' iteration overflowed and the dependence parameter became
NaN (also when started from the independence estimates). The same model converged for LLE alone, UPO alone
and K1 (dependence parameters 0.085, -0.025, 0.002). While diagnosing, the independence-working-correlation
fit for titration was seen (burden coefficient 0.836, robust SE 0.166).

**Change (technical fallback, applied uniformly to every GEE in `analysis10.py`):** if the exchangeable fit
gives a non-finite coefficient or standard error, the same model is refitted with the independence working
correlation (still a GEE with subject clusters and robust sandwich standard errors, which remain valid
under any within-subject correlation; it is less efficient). Each result records the working correlation
used. No other analysis, threshold or decision rule changed.

