# Phase 8 cardiac model library: sources, equations, parameters, verification

Every model is implemented from an accessible published source. Accessibility was checked
on 2026-10-01; see `HANDOFF.md`. In brief:

- **Reachable:** arXiv, PMC / Europe PMC full text (MathML / LaTeX), and PLOS
  supplements.
- **Blocked from this environment:**
  - APS, AIP, Elsevier / ScienceDirect;
  - physiology.org, MDPI;
  - Springer (bot challenge).

A source is **primary** when it is the paper that defines the formulation, and
**restatement** when it restates another paper's equations. Each model was verified by
reproducing a published result before use (Section "Verification"). Unverifiable
models were dropped.

Open-access source files used here are in `sources/`: the PLOS supplements S1/S2 Text
and S1-S3 Tables. The arXiv and PMC papers are cited by identifier.

## 1. Phase-resetting map: periodically stimulated cardiac oscillator (Glass-Shrier lineage). Native beat-to-beat. TEST

**Source (primary, open).** Diagne K, Bury TM, Pettebone ME, ..., Glass L et al.,
"Phase resetting in human stem cell derived cardiomyocytes explains complex cardiac
arrhythmias", *PLoS Comput Biol* 2026, doi:10.1371/journal.pcbi.1013935, PMC12900431.

**Equations.**
- Main text Eq. (1): φ_{i+1} = φ_i + τ − g(φ_i) (mod 1), with τ = T_stim / T.
- S1 Text Eqs. (1)-(3): the piecewise phase-response curve

| range | g(φ) |
|---|---|
| 0 ≤ φ < φ1 | 1 |
| φ1 ≤ φ < φr | 1 + A(φ−φ1)^4 |
| φr ≤ φ < φ3 | B(φ−φ2)^2 + C |
| φ3 ≤ φ < 1 | 1 + S(φ−1) |

  with φ1 = φr − 0.25, φ2 = φr + 0.04, φ3 = φ2 + S/(2B), and
  C = 1 − S(1−φ3) − B(φ3−φ2)^2.

**Parameters (S1 Table).** Aggregate A is the PRC of the main-text Fig 1B: A = 60,
φr = 0.71, B = 50, S = 0.80, intrinsic cycle 2.07–2.41 s. Aggregates B–F are also
implemented.

**Beat events.** The model's own action potentials. The last AP before stimulus i is at
s_i − φ_i T. The perturbed cycle ends with an AP g(φ_i)T after it, and the oscillator
then fires freely every T until the next stimulus. RR = AP-to-AP intervals.

**Units.** The physiological unit is the intrinsic cycle length T. RR series are
generated with T = 1 and rescaled (Section 6).

**Verification** (published results in the main text, Fig 3B and S8 Fig):

| τ | published | this implementation |
|---|---|---|
| 0.66 | stable period-3 orbit | 3 distinct stimulus phases; LE −0.853 per stimulus |
| 0.38, 0.53, 1.74 | (noisy) periodic | 5, 2, 8 distinct phases; LE −0.045, −0.526, −0.192 |
| 1.2 | chaotic rhythm | aperiodic (999 distinct phases in 1,000); LE +0.135 ± 0.002 |

**Verified.**

**Dropped part.** The modulated-parasystole model (Eq. 2; S2 Text Eqs. 4, 7) needs the
ectopic PRC with parameters m, p, r and N (S2 Table). Their functional form is not given
in the paper or its supplements; it refers to Bury et al., *Chaos* 2023, which is not
accessible. Implementing it would require inventing the form, so it was dropped.

## 2. AV nodal conduction (Sun, Amellal, Glass, Billette 1995). Native beat-to-beat. TEST, non-chaotic only

**Source.** The primary paper (*J Theor Biol* 173:79, 1995) is blocked. Its equations
and parameters are restated in Zhao X, Schaeffer DG, "Alternate pacing of border-collision
period-doubling bifurcations", arXiv:math/0609106, Eq. (51).

**Equations.** A = atrial-His interval, R = drift, H = His-to-atrial interval:
- R_{n+1} = R_n e^{−(A_n+H)/τfat} + γ e^{−H/τfat}
- A_{n+1} = Amin + R_{n+1} + (201 − 0.7 A_n) e^{−H/τrec} if A_n ≤ 130, and
  A_{n+1} = Amin + R_{n+1} + (500 − 3.0 A_n) e^{−H/τrec} if A_n ≥ 130

**Parameters.** Amin = 33 ms, τrec = 70 ms, τfat = 30,000 ms, γ = 0.3 ms.

**Beat events.** His activations. The His-His interval is H + A_{n+1}, in physiological
units (seconds).

**Verification.** Published: border-collision period doubling at H_bif = 56.9078 ms, where
A = 130. This implementation gives the A* = 130 crossing on the A ≤ 130 branch at
H = 56.90784 ms, with period 2 below it and period 1 above it. **Verified.**

**Ground truth.** The largest LE is about −0.006 per beat (the slow fatigue direction) at
every valid H. Below H ≈ 35 ms the map leaves its physiological range. No chaotic regime
exists in this model, so the family supplies non-chaotic regimes only.

## 3. Mackey-Glass delayed feedback. DEV

**Source.**
- **Equation.** The original Science 1977 paper is blocked. The original authors' open
  article is used instead: Glass L, Mackey MC, "Mackey-Glass equation", Scholarpedia
  5(3):6908, Eq. (1):

  dx/dt = β x_τ/(1 + x_τ^n) − γ x

- **Parameters.** β = 0.2, γ = 0.1, n = 10, τ varied: the standard set (Farmer 1982), as
  in Phase 5 G3.

**Verification** against three independent published estimates of the largest LE at
τ = 50 (Table 1 of arXiv:1810.01016):

| source | λ1 |
|---|---|
| Galerkin, N = 50 | 5.456e-3 |
| Breda & Van Vleck | 5.76e-3 |
| Sigeti | 5.83e-3 |
| **this implementation** (RK4 with tangent DDE, dt 0.1 and 0.05, 2 seeds) | **5.65–5.90e-3** |

**Verified.**
- At τ = 17 this implementation gives 5.0–5.3e-3. Farmer's often-quoted 0.0086 (via
  arXiv:1509.06057 Table 3) uses a coarser discretization with unstated units; other open
  sources give about 0.006.
- The period-doubling route is reproduced: period 1 for τ ≤ 13, period 2 at τ = 14,
  higher periods at τ = 15–16, chaos from τ ≈ 17.

**Beat events** (predeclared; the model has no heartbeat): intervals between successive
local maxima of x(t), parabolically interpolated. These are rescaled (Section 6).

## 4. Three coupled modified van der Pol oscillators (SA, AV, His-Purkinje) with delayed coupling. TEST

**Source (primary for this formulation, open).** da Silva Lima G, Savi MA, Bessa WM,
"Adaptive control of cardiac rhythms", *Sci Rep* 14:23446, 2024, PMC11458860.
- Eqs. (1)-(5) are given in the article's LaTeX.
- The parameter table is in Fig. 1b, read from the figure image:

| node | α | ν1 | ν2 | d | e |
|---|---:|---:|---:|---:|---:|
| SA | 3 | 1 | −1.9 | 1.9 | 0.55 |
| AV | 3 | 0.5 | −0.5 | 4 | 0.67 |
| HP | 7 | 1.65 | −2 | 7 | 0.67 |

- Couplings: k_SA-AV = k^τ_SA-AV = 3 and k_AV-HP = k^τ_AV-HP = 55, with τ_SA-AV = 0.8 and
  τ_AV-HP = 0.1. Null parameters are omitted, as the paper states: the unidirectional
  SA→AV→HP coupling is used "for all simulations".
- External stimulus: F_SA = ρ sin(ω t). Initial conditions u0 = (−0.1, −0.6, −3.3),
  u0' = (0.025, 0.1, 2/3). RK4 at 1 kHz. Time scaling β_t = 0.1048.

**Verification.**
- **Normal rhythm.** Published: "normal ECG with heart rate in approximately 90 bpm" with
  β_t = 0.1048. This implementation: the unforced period of 6.403 model units × 0.1048 =
  0.671 s, i.e. **89.4 bpm**. This also fixes the direction of the time scaling
  (real = model × β_t).
- **Pathological rhythms.** Published: "high-frequency rate and highly dispersed" rhythms
  at (ρ, ω) = (5.45, 5.6), (8.625, 2.1), (8, 3.3) and (9.6, 2.1). This implementation
  gives 137, 174, 187 and 256 bpm with RR coefficients of variation 0.44, 0.34, 0.43
  and 0.34.
- **Verified.**

**Beat events** (predeclared): upward crossings of u_HP through the midpoint of its
5th–95th percentile range, i.e. ventricular activations.

**An earlier candidate in this family was dropped.** PMC9938421 (Heliyon 2023) restates
the Gois & Savi 2009 model. Its Table 1, even after correcting three evident typos, makes
the SA oscillator decay to the stable equilibrium at x = −e1 for the published initial
conditions, with no oscillation. It therefore cannot reproduce the paper's own
bifurcation tables (chaos for 0 < k53 ≤ 5.76).

## 5. Dropped models

| model | reason |
|---|---|
| Seidel–Herzel baroreflex | **Not verifiable.** Primary paper (*Physica D* 115:145, 1998) blocked. Equations (1)-(14) and parameters from Dudkowska & Makowiec (arXiv:q-bio/0603016), cross-checked with PMC13305918 (fixed-duration systole); implemented in `models/seidel_herzel.py` (kept for the record). Their Figs 4–5 report a regular rhythm at T ≈ 0.76 s, heart-rate oscillations for θcNa > 2 s at θvNa = 1.65 s, and for θcNa > 0.6 s at θvNa = 3 s. Six readings of the ambiguous details were tried: ĉvNa 1 vs 10; respiration average 2/π, 1 or 0; with or without the systolic dp/dt in the baroreceptor signal. All give a regular T of 0.93–0.98 s and no oscillation at θvNa = 3 s; only the θvNa = 1.65 s transition (between 2.1 and 2.5 s) is reproduced |
| DeBoer beat-to-beat model, incl. the chaotic-respiration variant | AJP 253:H680 (1987) not accessible; the only open restatement (arXiv:physics/0503053) omits the parameter values; no accessible source for the chaotic-respiration variant |
| AV-node model of the ventricular response in AF | no accessible source with equations |
| Gois–Savi via PMC9938421 | see Section 4 |
| Modulated parasystole (Diagne et al. Eq. 2) | PRC form not given; see Section 1 |

## 6. Units and rescaling

- **AV node:** native seconds; no rescaling.
- **Phase-resetting map:** RR in units of the intrinsic cycle T. Each window is rescaled
  linearly to mean 0.8 s with the window's own coefficient of variation kept, i.e.
  T = 0.8 s / mean(RR). The model is in normalized phase units; 0.8 s is a typical adult
  sinus cycle.
- **Mackey–Glass:** no physiological units. Each window is rescaled linearly to mean 0.8 s
  and SD 0.05 s (the Phase 5 RR scale).
- **Coupled vdP:** native seconds after the published time scaling β_t.

The Part A realism table reports the resulting statistics.
