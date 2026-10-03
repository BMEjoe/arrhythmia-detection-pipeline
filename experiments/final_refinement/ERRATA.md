# Errata found in the final refinement

Every MISMATCH in `VERIFICATION.md` (produced by `verify_all.py`) is listed here with its evidence. The
error policy of this phase applies: nothing is fixed silently; a fix is made only if it restores the
preregistered analysis (or the analysis output) exactly; only the affected computation is rerun,
deterministically; every report and table that used the wrong number is updated; and a fix that would change a
conclusion is not made without the author's decision.

**Outcome: 14 findings. E1-E13 are corrected; none of them changes any result, conclusion, decision, winner or
preregistered outcome. E14 (the wording of Phase 10 Conclusion 8) is NOT corrected: correcting it would change a
conclusion, so it is left for the author's decision.**
- E1 is a tabulation deviation from a preregistered decision rule (Phase 9). Corrected by a deterministic
  re-tabulation; the decision is unchanged.
- E2–E12 are report-text errors in `docs/PHASE10_FINAL_RESULTS.md`: a miscount, a rounding slip, and summary
  sentences that misstate the tables beneath them. The tables, analysis outputs and raw results are correct.
  The text was corrected in place; each corrected passage carries the marker `[Erratum En]`.
- E13 is a stale derived file: the development `analysis.json` / `tables.md` of Phase 10 predated the last
  development robustness records. Regenerated deterministically; no reported number used it.

| id | where | was | is | evidence (VERIFICATION.md) | conclusion changed? |
|---|---|---|---|---|---|
| E1 | Phase 9 TEST decision (`decide9.py` output; `docs/PHASE9_WAVEFORM_NOISE_ROBUST.md` 6) | k2 fails **12** conditions (81 conditions tabulated) | k2 fails **16** of the **101** preregistered conditions | 9-1.k2_mnlp_ann, 9-2 | no |
| E2 | PHASE10 3 (Windows row) | "6 subjects with 10–11" | 7 subjects (six with 11 windows, one with 10) | 10-1 | no |
| E3 | PHASE10 8.5 (chf218 window 12) | z_max 16.2 | 16.1 (16.147) | 10-K3t-chf218.12 | no |
| E4 | PHASE10 12 #13 | "Masking or editing ectopy removed 49–66 % of segment positives" | masking removed 49–66 %, editing 44–55 % | T-1 | no |
| E5 | PHASE10 12 #5 | LLE alone "fired in 94–99 %" of windows with 2–10 % isolated ectopy | 96–99 % (94 % is a Wilson lower bound, not a rate) | T-2 | no |
| E6 | PHASE10 1 (Q1) | "K3 never saw weak maps (λ ≤ 0.26 per beat) … no bound … at any f" | true for logistic λ ≤ 0.26 and Hénon λ = 0.14; Hénon a = 1.14 (λ 0.244) has a π < 0.20 bound from f = 0.7 | T-3, 10-Q1-K3.henon.a1.14 | no |
| E7 | PHASE10 1 (Q1) and 5 (Findings) | "K1 reached only π < 0.20 from f = 0.5–0.7, and only for the strong maps" | K1 also reached π < 0.05 at f = 0.9 (Hénon a = 1.40, logistic r = 3.88, 4.00) and π < 0.20 for Hénon a = 1.14 (f = 0.7) and phase-reset τ = 1.14 (f = 0.9) | T-4, 10-Q1-K1.* | no |
| E8 | PHASE10 1 (Q4) and 8.1 | "499 surrogates changed ≤ 3 decisions per detector" | up to 6 window decisions (LLE alone: 3 lost, 3 gained); K1 1, UPO alone 2, K3 0 | T-5, 10-Q4a-* | no |
| E9 | PHASE10 5 (Findings, first bullet) | stronger maps (λ ≥ 0.30) excluded at π < 0.05 for "f ≥ 0.5" | f ≥ 0.5, except logistic r = 4.00 (f ≥ 0.7) | T-6, 10-Q1-K3.logistic.r4.00 | no |
| E10 | PHASE10 7 (3b text) | "removes most titration and LLE positives" | most titration positives (57.5 % masked, 66.0 % edited), about half of LLE positives (47.1 %) | T-7, 10-P3b-* | no |
| E11 | PHASE10 8.2 | "Three subjects recur from the 512-interval detections" | five of the six subjects recur (nsr029, nsr049, chf211, chf216, chf218) | T-8 | no |
| E12 | PHASE10 6 | development nsrdb "lost only 2–6 % of positives" | 4–6 % (3.5–6.1 % over the masked, edited and Wu arms) | T-9 | no |
| E13 | `phase10_final/results/dev/analysis/{analysis.json, tables.md}` | Q4 development robustness computed on 3 windows | 36 windows (the committed `dev/robust.jsonl`) | reproduce check (rerun of `analysis10 --phase dev`) | no |
| E14 (**not corrected**) | PHASE10 12 #8 (Conclusion 8) | "None of the project's detectors sees continuous-flow chaos sampled as RR intervals at realistic lengths" | holds for flows sampled at fixed time steps; Mackey-Glass chaos as maxima intervals was detected (K1 26/70, K3 81/140) | X-21 | **would change Conclusion 8**: author's decision |

## E1 — Phase 9 TEST decision tabulated over 81 instead of the 101 preregistered conditions

**What was wrong.** The Phase 9 preregistration (`experiments/phase9_waveform/PREREGISTRATION.md` Sections 3–4)
defines PASS as ≤ 7/100 in EACH of 101 conditions: 17 + 17 nulls, 18 Phase 8 non-chaotic regimes at the input
level, and 18 Phase 8 + 31 KTz non-chaotic regimes at the realistic level. A KTz regime is a (pacing period P,
input pattern) pair (`results/ground_truth/ktz_labels.json`: e.g. P = 92 with input none, S2_5, E1). The TEST
windows name a KTz regime by P only (`ktz:P=120`), with the input pattern in the separate `ectopy` field. The
frozen `decide9.py` groups by `group + "|" + name`, so it pooled the 31 KTz non-chaotic regimes into 11 groups of
100–400 windows with limits floor(0.07 N) = 7–28, and tabulated 81 conditions.

**Evidence.** `verify_all.py` 9-1 / 9-2; `errata_phase9_decision.py` → `errata/phase9_decision_101.md`.
Under the preregistered 101 conditions:

| candidate | PASS (101) | failed (101) | PASS (decide9.py, 81) | failed (81) | primary chaotic | PASS-condition detections |
|---|---|---|---|---|---|---|
| k1_frozen512 | no | 4 | no | 4 | 74/1,320 | 97/10,100 |
| k2_mnlp_ann | no | **16** | no | 12 | 304/1,320 | 607/10,100 |
| k3_mnlp_growth_ann | **yes** | 0 | yes | 0 | **90/1,320** | **0/10,100** |
| k4_mnlp_growth_rr | yes | 0 | yes | 0 | 2/1,320 | 1/10,100 |
| k5_mnlp_fsle_ann | yes | 0 | yes | 0 | 0/1,320 | 3/10,100 |

The four extra k2 failures are KTz non-chaotic regimes without ectopy: P = 120 (8/100), 200 (12/100),
250 (10/100), 300 (9/100). Pooled over their P group they were within the pooled limit.

**Affected.** `docs/PHASE9_WAVEFORM_NOISE_ROBUST.md` Section 6 (k2 row: "12"); `results/test/tables/decision9.md`
and `by_condition9.csv` (81-condition tabulation). No other report quotes the k2 failure count.

**Conclusion changed?** No. Every candidate's PASS/FAIL is the same, the winner is still k3_mnlp_growth_ann (the
only quantity used later), the pooled counts are identical, and k3's specificity (0/10,100) is unaffected.

**Fix (restores the preregistered rule exactly).** `errata_phase9_decision.py` re-tabulates the decision from
the stored per-candidate detection flags (the same inputs as `decide9.py`) with the preregistered condition key
(KTz: group|name|input). `decide9.py` is not modified: its SHA-256 is frozen in the Phase 9 preregistration.
The Phase 9 report row was corrected (marked `[Erratum E1]`), and an erratum note was appended to
`decision9.md`; `by_condition9.csv` is kept as the frozen script produced it, and the corrected per-condition
table is `errata/by_condition9_101.csv`.

## E2 — number of subjects with fewer than 12 windows

**Was:** "988 windows of 512 intervals (49 + 27 subjects with all 12, 6 subjects with 10–11 …)". **Is:** 7
subjects with 10–11 windows (chf207, chf212, nsr016, nsr036, nsr041, nsr045 with 11; nsr015 with 10):
49 × 12 + 27 × 12 + 6 × 11 + 10 = 988. Source: `results/conf/real.jsonl`. Affects PHASE10 Section 3 only.
No analysis used the count. No conclusion changes.

## E3 — z_max of chf218 window 12

**Was:** 16.2. **Is:** 16.1 (stored 16.1466; `results/conf/real.jsonl`, `K3.zmax`). A rounding slip in the
descriptive table of PHASE10 Section 8.5. No conclusion changes.

## E4 — range of positives removed (Conclusion 13)

**Was:** "Masking or editing ectopy removed 49–66 % of segment positives (Phase 10, P2)." **Is:** masking removed
49 % (NSR) and 66 % (CHF), editing 44 % and 55 % (PHASE10 Section 6 table: 49.2, 66.5, 44.5, 55.2). The sentence
quoted the masked range for both arms. Conclusion 13 (titration on heart-rate data mostly detects ectopy) is
unchanged.

## E5 — LLE-alone firing rate on isolated ectopy (Conclusion 5)

**Was:** "It fired in 94–99 % of linear-Gaussian windows with 2–10 % isolated ectopic beats (Phases 5–6)."
**Is:** 96–99 % (Phase 5 S2: 194, 195, 197 / 200; Phase 6 S2: 289, 290, 293 / 300). 94 % is the lowest Phase 6
Wilson lower bound. Conclusion 5 is unchanged.

## E6 — "weak maps" in the Q1 summary

**Was (Section 1):** "K3 never saw weak maps (λ ≤ 0.26 per beat) or any phase-resetting regime, so no bound is
possible for them at any f." **Is:** no bound exists for logistic r = 3.58, 3.65 (λ ≤ 0.26), Hénon a = 1.08
(λ 0.14) and every phase-reset regime; but Hénon a = 1.14 (λ 0.244 ≤ 0.26) was detected in 8–9/83 windows at
f ≥ 0.7 and has a π < 0.20 bound from f = 0.7 (Section 5.1 table, which is correct; Section 5 "Findings" also
states it correctly). No conclusion changes.

## E7 — K1 bounds in the Q1 summary

**Was (Section 1):** "K1 reached only π < 0.20 from f = 0.5–0.7, and only for the strong maps." **(Section 5
Findings):** "π < 0.20 only from f = 0.5–0.7 for the strong maps, and no bound at all for coupled vdP."
**Is:** K1 reached π < 0.20 from f = 0.5–0.7 for Hénon a ≥ 1.14 and logistic r ≥ 3.88, and from f = 0.9 for
phase-reset τ = 1.14; it reached π < 0.05 only at f = 0.9 (Hénon a = 1.40, logistic r = 3.88 and 4.00); no bound
for coupled vdP (Section 5.2 table, which is correct). "K1 gives weaker bounds than K3" is unchanged.

## E8 — decisions changed by 499 surrogates

**Was (Section 1):** "499 surrogates changed ≤ 3 decisions per detector in 331 windows." **(Section 8.1):**
"Surrogate Monte Carlo error changes at most 3 decisions per detector." **Is:** the Section 8.1 table (correct)
shows LLE alone 3 default-only and 3 499-only windows, i.e. 6 changed window decisions (net 0); K1 1, UPO alone 2,
K3 0. Agreement 98.2–100 %, κ 0.91–1.00 are unchanged. No conclusion changes.

## E9 — exclusion fraction for the stronger maps

**Was (Section 5 Findings):** "For the stronger Hénon and logistic maps (λ ≥ 0.30 per beat) and one coupled-vdP
regime, chaos accounting for ≥ 50 % of beat-to-beat variance (f ≥ 0.5) can be present in at most 5 % of real
windows." **Is:** f ≥ 0.5 for Hénon a = 1.22, 1.40, logistic r = 3.88 and vdP (8, 3.3), but f ≥ 0.7 for logistic
r = 4.00 (Section 5.1 table, correct; Section 1 gives the correct range 0.5–0.7). No conclusion changes.

## E10 — "most … LLE positives" removed (Q3b)

**Was (Section 7, 3b):** "Masking or editing removes most titration and LLE positives in windows that contain
ectopy." **Is:** titration 57.5 % (masked) and 66.0 % (edited); LLE alone 47.1 % (edited), i.e. about half (3b
table, correct). No conclusion changes.

## E11 — subjects recurring among the 1,024-interval K3 detections

**Was (Section 8.2):** "Three subjects recur from the 512-interval detections (nsr029, nsr049, chf216; chf211 and
chf218 also reappear)." **Is:** the nine 1,024-interval detections come from nsr025, nsr029, nsr049, nsr051,
chf202, chf211, chf216 (two windows) and chf218; five of the six subjects with a 512-interval detection recur
(nsr029, nsr049, chf211, chf216, chf218); chf211 window 7 is detected at 256, 512 and 1,024 intervals. Source:
`results/conf/robust.jsonl`, `real.jsonl`. No conclusion changes.

## E12 — development nsrdb positives removed

**Was (Section 6):** "In development, nsrdb (20–50 y, almost no ectopy) lost only 2–6 % of positives." **Is:**
3.5–6.1 % of raw positives (Wu 3.5 %, edited 3.6 %, masked 6.1 %; `results/dev/analysis/analysis.json`), i.e.
4–6 %. No conclusion changes.

## E13 — stale development analysis output (Phase 10)

**What was wrong.** The committed `experiments/phase10_final/results/dev/analysis/analysis.json` (and the
`tables.md` generated from it) was written when the development robustness run had 3 windows. The development
run continued in the background and the committed `results/dev/robust.jsonl` holds 36 windows (commits `a3defa0`,
`1147f76`). Rerunning `analysis10 --phase dev` on the committed raw files changes only the Q4 development block
(surrogate-count agreement and window-length rates, now on 36 windows) and adds the `working_correlation` field
written by the current code (every development GEE converged with the exchangeable working correlation, as
before; their estimates are unchanged). The development spike-in, real-window and segment results are unchanged.

**Affected.** Only these two development files. `docs/PHASE10_FINAL_RESULTS.md` 9.1 reports the development
robustness run only as "36 of 132 tasks", which is correct; no development Q4 number appears in any report.

**Fix.** `analysis10 --phase dev` and `tables10 --phase dev` were rerun (deterministic; same frozen code) and the
regenerated files committed. **Conclusion changed?** No.

## E14 — Conclusion 8 is broader than its evidence (NOT corrected; author's decision required)

**What the report says (PHASE10 Section 12, Conclusion 8).** "None of the project's detectors sees
continuous-flow chaos sampled as RR intervals at realistic lengths. Rössler and Mackey–Glass were 0 % for K1 at
every m (Phases 5, 8). K3 needs strong low-dimensional predictability that decays within 5 beats (Phase 9)."

**Evidence.** The cited results are flows sampled at fixed time steps (Phase 5 G2 Rössler, G3 Mackey-Glass
τ = 17 sampled every 6.0; the same flows in Phase 8): K1-type AND 0/200 at m = 2, 3, 4 (checks 5-6, 8-5). But
Phases 8 and 9 also report Mackey-Glass chaos represented as intervals between successive maxima (the Phase 8
model's predeclared beat definition; development regimes τ = 16.5-30 evaluated at TEST seeds), and there the
detectors did detect it (check X-21, recomputed from the raw files):
- Phase 8, K1 (C1, 512 intervals): 26/70 at the input variant, 23/70 with isolated ectopy, 18/70 with couplets,
  0/70 with bigeminy (`docs/PHASE8_CARDIAC_CHAOS.md` 7.3);
- Phase 9, realistic level: K3 81/140 and K1 42/140 with isolated ectopy or couplets, 0/70 with bigeminy
  (`docs/PHASE9_WAVEFORM_NOISE_ROBUST.md` 6); by regime, K3 20/20 at τ = 18 and τ = 20, 0-1/20 at τ = 16.5-17.

**What is affected.** Only the wording of Conclusion 8 (and its mention in `docs/CONCLUSIONS_CHECKLIST.md`). No
number, table or other conclusion depends on it. These Mackey-Glass results were secondary (development regimes,
not out-of-sample for Phase 8, see the Phase 8 report), which may be why the conclusion did not cite them.

**Why it is not corrected.** A correction (for example restricting the statement to flows sampled at fixed time
steps and to weak, near-threshold Mackey-Glass chaos) changes what the conclusion claims. Under the error policy
of this phase such a fix is not made without the author's decision. The report text is unchanged.

## Observations that are not errors

- **O1. Monte Carlo term of the real-rate bound U (P1).** U = max(bootstrap 95th percentile, Clopper–Pearson).
  For K3 the bootstrap percentile takes one of two values, 0.0141 or 0.0151, depending on the seed (8 detections
  in 6 subjects). The seeded analysis drew 0.0141, so the CP term (0.0146) decided. With 0.0151, every f_min of
  every K3 and K1 row is unchanged (check 10-Q1-robust, 64/64). The analysis is seeded and deterministic; nothing
  was changed.
- **O2. Amendment 2 corroborated.** An independent implementation of the exchangeable GEE also fails to converge
  for the titration P3 model (check 10-P3-A2), as Amendment 2 states for statsmodels.
- **O3. Working correlations differ between titration models.** Under Amendment 2's rule, models whose
  exchangeable fit converged kept it. In PHASE10 Section 7, the adjusted titration ORs (2.31; CHF 0.77) use the
  independence working correlation, while the unadjusted burden OR (3.85) and the unadjusted CHF OR (2.60), the
  night model (Section 8.3) and the development models use the exchangeable one. This follows the amended
  rule, but the adjusted and unadjusted burden ORs are therefore not a like-for-like comparison. (While
  debugging the verifier, an independence fit of the unadjusted burden model was computed by mistake and gave
  2.21; it is not a reported result. A like-for-like comparison is listed as an open item.)
- **O4. V4 margin.** The titration V4 counts (78/78 chaotic, 18/18 periodic) use the predeclared ±0.02 margin on
  the Lyapunov exponent (METHODS.md 2.5). Without the margin the scan gives 80/80 and 21/21.
- **O5. Phase 9 Part F and Phase 10 use the stored detection flags.** For K3 and K4 the flags were re-derived from
  the stored z_max and G with the frozen thresholds; they agree in every window (checks 9-0, 10-0).
