# Phase 9 methods: sources, implementations and verification

Every method is implemented from its primary published source where that source is
accessible. Accessibility was checked on 2026-10-02:
- **Reachable:** arXiv; NCBI E-utilities full text for PMC articles (the PMC web pages are
  behind a reCAPTCHA; Europe PMC full-text returns HTTP 500 for several NIH manuscripts);
  PLOS; Frontiers; eScholarship metadata; physionet.org (WFDB, data, ECGSYN paper HTML);
  raw.githubusercontent.com; PyPI.
- **Blocked:** APS (link.aps.org), IEEE Xplore, AHA journals, Elsevier/ScienceDirect, Royal
  Society.

A method is used only after a published behaviour was reproduced (Section "Verification" of
each entry). Unverifiable methods are dropped and listed in Section 9.

## 1. Waveform generator: ECGSYN (A1) — VERIFIED

**Source (primary, open).** McSharry PE, Clifford GD, Tarassenko L, Smith LA, "A dynamical
model for generating synthetic electrocardiogram signals", *IEEE TBME* 50(3):289–294 (2003).
- The authors' HTML version is on PhysioNet (physionet.org/content/ecgsyn/1.0.0/paper/):
  Eqs. (1)–(4) and Table I.
- Reference implementation: the authors' `ecgsyn.c` (PhysioNet ECGSYN 1.0.0, GPL).

**Implementation.** `ecgsyn.py`.
- The equations are Eq. (1)–(4).
- The limit-cycle phase is driven by an **external beat-time series**: θ advances by 2π
  between consecutive input beats. On the limit cycle θ̇ = ω exactly, so this is Eq. (1) with
  ω(t) = 2π/T(t) piecewise constant per beat. The R event falls exactly on each input beat.
- z is integrated with RK4 at 720 Hz and every 2nd sample is kept (as `ecgsyn.c`).
- Heart-rate adjustment of the b_i and θ_i: as `ecgsyn.c`.
- Baseline z0 = 0.005 sin(2π 0.25 t + φ): the `ecgsyn.c` value, equal to the paper's
  0.15 mV after scaling. φ ~ U(0, 2π) per window.
- **Fixed mV map** from the 60-bpm reference beat to [−0.4, 1.2] mV. `ecgsyn.c` instead
  rescales each record by its own min/max.
- **Per-beat extrema parameters** (the only extension): beat k owns the half-cycle before
  and after its R event.

**Verification** (`verify_ecgsyn.py`, `results/verification/ecgsyn.json`,
`plots/A1_ecgsyn_verification.png`):

| check | published / reference | this implementation |
|---|---|---|
| vs `ecgsyn.c`, constant 60 / 75 bpm (C compiled with stand-in FFT and RNG for the two non-redistributable Numerical Recipes files; `tools/ecgsyn_c/`) | — | max \|Δ\| 1.0e-5 / 5.3e-5 mV, r = 0.99999999994 / 0.999999997 |
| vs `ecgsyn.c`, 60 ± 5 and 90 ± 3 bpm (beat times reconstructed from the C RR process) | — | r = 0.99998 / 0.999999; max \|Δ\| 0.027 / 0.004 mV (the C code changes ω mid-beat on its sample grid) |
| Table I extrema times, 60 bpm | P −0.2, Q −0.05, R 0, S 0.05, T 0.3 s | P −0.168, Q −0.046, R 0.000, S 0.045, T 0.247 s |
| Fig. 9: RS amplitude vs RR | strong positive correlation | r = 0.9997 |
| Fig. 10: QT vs RR | linear | R² = 0.977 |
| beat times reproduce the input RR | Fig. 7 | R-peak offset −0.25 ms (constant); RR error max 1.4 ms (< 1 sample) for a Phase 8 chaotic RR window with couplets |

The Table I times are approximate. The published angles (−π/3, −π/12, 0, π/12, π/2) at
1 s per beat give −0.167, −0.042, 0, 0.042 and 0.25 s, which is what the implementation
produces.

**Verified.**

## 2. Ectopic morphology (A1) — predeclared

Implemented in `morphology.py`; constants fixed before any detector or measure was run on a
Phase 9 series.

**Ventricular ectopic beats** ('V': Phase 5 S2; Phase 6 E1–E4 and E6; Phase 8 ectopy
variants).
- Basis: LITFL "Premature Ventricular Complex": "Broad QRS complex (≥ 120 ms) with abnormal
  morphology", "Discordant ST segment and T wave changes", "Usually followed by a full
  compensatory pause".
- No P wave.
- QRS widths and Q/S angles ×1.5. The model QRS (Q − 2σ to S + 2σ) is ~131 ms normal and
  ~197 ms for the PVC.
- QRS amplitudes ×0.68, giving R ≈ 1.53× normal.
- T inverted (×−0.50, width ×1.5), giving |T| ≈ 1.00× normal. The amplitude constants were
  calibrated numerically at 75 bpm.

**Atrial ectopic beats** ('A', E5).
- Basis: LITFL "Premature Atrial Complex": abnormal P, "an inverted P wave" from low atrial
  foci, normal QRS.
- Implementation: P inverted, QRST normal.

## 3. Repolarization (APD) morphology families (A2)

### 3.1 KTz paced-cell map — VERIFIED, TEST family

**Source (primary, open).** Gall WM et al., arXiv:2202.12406.
- Single-cell logistic KTz map, Eqs. (1)–(3); Jacobian Eq. (7).
- Parameters (Sect. 2–3.1): K = 0.6, Ie = 0, δ = λ = 0.001, T = 0.154, xr = −0.48.
- Stimulus 0.1 for 10 time steps every P steps.
- APD between zero crossings of x.

**Implementation.** `morph/ktz.py`: QR (Eckmann–Ruelle) Lyapunov exponent at every step.

**Verification** (published Fig. 3b): largest LE λ_L ≈ 0.0035 per time step at P = 92, the
maximum; chaotic insets at P = 92 and 232.

| P | this implementation |
|---|---|
| 92 | λ = 0.00344 (0.00346 ± 0.00005 over 3 initial conditions); chaotic |
| 232 | λ = +0.0020; chaotic |
| 250, 300 | period-1, λ < 0 |

Full scan of P = 80–300: `results/verification/ktz_scan.jsonl`. **Verified.**

**ECG mapping** (`morphology.py`, predeclared):
- Time step = 0.8 s / P_nom (paced at 0.8 s).
- T event moved by APD_k − APD_ref; T amplitude × APD_k / APD_ref.
- Basis: T-wave alternans is the ECG manifestation of APD alternans (Qu et al., *Phys Rep*
  2014, PMC4175480).

**Ground truth.** `ground_truth_ktz.py`, `results/ground_truth/ktz_labels.json`:
- The LE of the map is computed for the exact input: constant pacing, or pacing with Phase 8
  ectopy S2_5, E1 or E3_10, which alters the pacing intervals.
- 3 realizations, 8,000 beats each; the Phase 8 map rule (Section 3.2 of HANDOFF).
- Constant pacing: 7 CHAOTIC (P = 92, 100, 146, 230, 232, 260, 278) and 5 NON-CHAOTIC
  (P = 120, 150, 200, 250, 300).
- Ectopy in the input destroys most of the narrow chaotic windows. CHAOTIC remain only at
  (92, S2_5), (92, E3_10), (146, S2_5) and (146, E1). 8 pairs are AMBIGUOUS and discarded.

### 3.2 Modified Luo–Rudy I with EAD chaos — DROPPED (not verifiable)

**Source.** Tran DX et al., *PRL* 102:258103 (2009), NIH manuscript PMC2726623.
- LR1 with E_si = 80 mV and E_K = −77 mV; τ_d, τ_f and τ_x scaled by α, β and γ.
- Fig. 4: γ = 2.5, chaos at intermediate PCL, LE 0.38 s⁻¹ at PCL 0.907 s.
- The LR1 equations come from the Physiome Model Repository CellML `luo_rudy_1991` (the
  original *Circ Res* paper is blocked). Implemented in `morph/lr1.py`.

**Reason for dropping.** The Fig. 4 values of Ḡ_si, α and β are not stated in the
accessible text. The full texts of the companion papers (Sato 2009 PNAS, PMC2651322;
Sato 2010 Biophys J, PMC2913181) are not available through PMC full-text services.

**Readings tried** (`results/verification/lr1_scan_*.jsonl`, `lr1_traces.png`; 2-trajectory
LE, 150 beats after 120):
- **α = β = 1, γ = 2.5:**
  - Ḡ_si 0.10–0.14 × PCL 0.6–1.5 s;
  - Ḡ_si 0.12 × PCL 1.4–3.0 s, refined to 10 ms steps over 2.20–2.40 s.
- **Reading R2:** all K reversal potentials −77 mV.

**Outcome.**
- LE < 0 everywhere, except values ≥ −0.005 s⁻¹ inside a period-adding cascade at
  Ḡ_si = 0.12, PCL 2.32–2.37 s.
- Under α = β = 1, EAD onset coincides with the AP outlasting the cycle.
- R2 produces no EADs.

The published chaos (0.38 s⁻¹ at 0.907 s) was not reproduced, so the family is dropped. As a
result, **only one morphology family exists**. It is placed in TEST, and **development had no
morphology family**.

## 4. Recording conditions (A3)

**Sampling and amplitude.**
- 360 Hz sampling.
- MIT-BIH amplitude resolution: 200 adu/mV, 11-bit (`noise.quantize_mitbih`).

**nstdb noise** (Moody GB, Muldrow WE, Mark RG, *Computers in Cardiology* 1984; PhysioNet
nstdb 1.0.0).
- Records bw, ma and em. SHA-256 values were checked against PhysioNet's SHA256SUMS.txt.
- SNR uses the definitions of the WFDB program `nst` (physionet.org/physiotools/wag/
  nst-1.htm):
  - S = (QRS peak-to-peak)² / 8, with the trimmed mean of the first 300 normal QRS;
  - N = (RMS of 1-s chunks)², trimmed mean.
- **Split:** DEV uses signal 0 and TEST uses signal 1.
- **Mixture:** bw + ma + em at equal nst power.

**Beat times.**
- **True** beats: generator R events on the 360 Hz grid.
- **Detected** beats: `fp.detect_r_peaks` (the project's Pan–Tompkins), matched within
  150 ms.
- **Jitter:** at detected V beats, the 7,761 Phase 7 measured detected − annotated offsets of
  MIT-BIH V/E/F beats, with the Phase 8 0.15 s guard.
  - On the synthetic PVCs the detector's own offset is ~0 samples (median), so the real late
    PVC fiducial is not reproduced by the synthetic morphology. The measured offsets are
    therefore added explicitly.

**Realism notes** (`realism.py`, `results/realism/realism.md`, `plots/A_*`).

The table gives medians over 5 realism seeds at 512 beats. ECGSYN properties visible under
ectopy, kept because they are properties of the published model:
- **Event size scales with the surrounding interval.** Each event's deflection is
  proportional to the duration of the interval it lies in (paper Fig. 9). A PVC followed by
  a long compensatory pause therefore has enlarged S and inverted-T deflections, up to about
  2–4 mV in couplets.
- **Abrupt RR changes leave a baseline offset.** When ω changes at an R event, the R event's
  forcing integrates to zero in phase but not in time. The offset decays at the model's
  relaxation rate (1/s) over a few beats. The authors' `ecgsyn.c` behaves the same way with
  its piecewise-constant RR.
- **Coupled vdP is extreme.** Its chaotic and several non-chaotic regimes keep CV 0.3–0.5
  after rescaling to 0.8 s, with minimum RR down to about 0.2 s. These are the most extreme
  waveforms, and beat detection there has the most misses and extra detections.

## 5–8. Parts B–D methods

*(added as they are implemented and verified)*

## 9. Dropped methods and models

| item | reason |
|---|---|
| Modified LR1 EAD family (Tran et al. 2009) | not verifiable (Section 3.2) |
