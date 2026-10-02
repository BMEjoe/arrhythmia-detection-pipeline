# Phase 9 Part F preregistration: the winner on real long-term RR (nsrdb vs chfdb)

This file is committed and pushed **before nsrdb or chfdb is downloaded**. No data from
either database has been seen in this project.

**Detector.** The Phase 9 winner `k3_mnlp_growth_ann`, adopted as
`final_pipeline.masked_growth_chaos_test`:
- turned on with `PipelineConfig.masked_growth_test`, default off;
- uses the preregistered thresholds Z = 9.4521 and G = 0.25, 39 surrogates;
- runs on 512-interval windows.

The detector code, the thresholds and this analysis are frozen once this file is pushed.

**What a detection means.** The window's beat sequence, with its ectopy-related intervals
masked, is:
1. more predictable than masked IAAFT surrogates; and
2. has a forecast error that keeps growing with horizon.

On synthetic data this rule was specific: 0 of 10,100 TEST null and non-chaotic windows
detected. Several properties of real data were NOT in the synthetic nulls, so a detection on
real data is NOT a proof of chaos:
- nonstationarity over minutes;
- sleep/wake transitions;
- non-Gaussian nonlinear stochastic regulation;
- atrial fibrillation;
- annotation errors.

## 1. Data

- PhysioNet nsrdb 1.0.0: MIT-BIH Normal Sinus Rhythm Database, 18 subjects.
- PhysioNet chfdb 1.0.0: BIDMC Congestive Heart Failure Database, 15 subjects.
- Each database is downloaded with `wfdb` into a gitignored directory. SHA-256 sums are
  checked against PhysioNet's SHA256SUMS.txt.
- Beat annotations used:
  - nsrdb: `atr`;
  - chfdb: `ecg` (the database's beat annotation file).
- If an annotation extension differs from the one stated here, the one PhysioNet documents for
  that database is used, and an amendment records it before any detector run.
- The unit is the **subject**; one record per subject in both databases.

## 2. Beats, labels and windows (predeclared)

**Beats** are the annotations whose symbol is a beat symbol: wfdb beat codes N L R B A a J S
V r F e j n E / f Q ?. Non-beat annotations (rhythm, noise, comments) are ignored.

**Label:**
- 'N' iff the symbol is N, L, R, B, e or j (normal / bundle-branch / escape beats of
  supraventricular origin with a normal conduction sequence);
- every other beat symbol is non-normal.
- An interval shorter than 0.25 s or longer than 2.5 s marks its ending beat non-normal (a
  missed or spurious annotation).

**Windows:** per subject, up to 10 windows of 512 consecutive intervals. Window j
(j = 1 … 10) starts at the first beat at or after j hours from the start of the record. A
window that does not fit inside the record, or that overlaps the previous window, is skipped.

**Decision:** `masked_growth_chaos_test(rr, labels)["detected"]`. Not analysable counts as
not detected. The fraction of analysable windows is reported.

## 3. Statistics (subject-level)

- **W1 (primary).** Per group (NSR, CHF):
  - the subject-level detection rate r̄ = mean over subjects of each subject's fraction of
    detected windows (subjects with 0 windows are excluded);
  - 95 % CI by subject-cluster bootstrap (10,000 resamples of subjects, seed 20261007,
    percentile interval).
  - Compared with the synthetic false-positive rate: the winner's pooled TEST rate on all
    PASS conditions was 0/10,100, with one-sided 95 % Clopper–Pearson upper bound 0.0003.
  - **"Exceeds synthetic FP"** iff the lower 95 % bound > 0.0003.
- **W2.** Difference r̄_CHF − r̄_NSR:
  - 95 % CI by bootstrap, resampling subjects within each group (10,000, seed 20261008);
  - two-sided permutation test of group labels over subjects (10,000 permutations, seed
    20261009).
  - **"Group difference supported"** iff the CI excludes 0.
- **Secondary:**
  - window-level counts;
  - analysable fraction and masked-interval fraction per group;
  - z and G distributions;
  - per-subject table.

## 4. Then (EXPLORATORY, decided now but not confirmatory)

Phase 7 MIT-BIH (48 records, subjects as in Phase 7):
- the same detector on 512-interval windows of the annotation beat series
  (`experiments.phase7_mitbih.data`, annotation-time source, label 'N' as Phase 7's
  effective label);
- non-overlapping windows from the record start;
- detection rate by subject, for records with and without abnormal rhythm (Phase 7 window
  classes).

MIT-BIH influenced nothing in Phase 9 and is reported as exploratory only.

## Amendments

**2026-10-02 (before any Part F detector run; data downloaded, only beat counts and window
positions of one record per database looked at).** The adopted form changed during the
ground-rule check. `analyze_segment` must keep a label-free signature
(`tests/test_upo_feature_contract.py`), so the winner is adopted as the standalone opt-in
function `final_pipeline.masked_growth_chaos_test(rr, labels, config)`. Its parameters live
in the `PipelineConfig.masked_growth_*` fields, and it is never called by default. There is
no `masked_growth_test` flag. The test itself is unchanged: equivalence with `candidates9.py`
is tested. Part F calls this function with `CFG`.
