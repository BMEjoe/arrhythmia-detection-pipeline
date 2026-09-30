# Phase 5: RR-interval stress test of the frozen combined detector

Evaluation, not tuning. The detector under test (`detector.py`) is
`fp.analyze_segment(rr, fp.phase4_upo_config(replace(fp.CFG, lle_chaos_test=True)))`
on raw RR (`use_corrected_rr_for_dynamics = False`). The components are
`out["lle_chaos_test"]["detected"]` (Phase 3 winner) and
`out["upo"]["instability_gate_detected"]` (Phase 4 winner C2). The primary
decision is their AND. No MIT-BIH record is loaded in this phase.

Files:
- `config.py`: every generator parameter and the seed scheme;
- `systems.py`: the generators;
- `realism.py`: the realism check (development seeds, no detector);
- `equivalence_check.py`: step 0.

`HANDOFF.md` tracks the status.

## Step 0: detector path

On 70 windows (14 Phase 4 conditions × test seeds 2000–2004, 256 samples),
`analyze_segment` with the config above gives the same C2 UPO decisions and the
same `lle_chaos_test` decisions and p-values as the stored Phase 4 results
(`results/equivalence_check.md`): **0 differences**.

## Signals

Every window is **256 RR intervals in seconds**. Unless a condition is marked
unquantized, the series is rounded to multiples of **1/360 s** (the MIT-BIH
sampling interval) after all other transformations.

**Common RR scale.** Each window draws:
- a mean RR μ ~ U(0.65, 0.95) s;
- an SD σ ~ U(0.035, 0.060) s.

The base process is standardized and set to μ + σ·z. μ is kept inside the
0.6–1.0 s target with a margin for trends, steps and ectopics. The σ range
brackets the short-term SDNN norm, 50 ± 16 ms [Nunan 2010]. It also agrees with
the Task Force 5-min LF + HF power: 1170 + 975 ≈ 2145 ms², i.e. about 46 ms
[Task Force 1996].

**Beat time.** Beat k is taken to occur at k·μ. A frequency f in Hz is
therefore f·μ cycles per beat.

### Null conditions in the PASS rule: not chaotic, the detector should not fire

**N1 `linear_rr`**
- Bimodal RR spectrum of the ECGSYN RR generator [McSharry 2003]: Gaussian bumps
  at 0.1 Hz (LF) and 0.25 Hz (HF), each with SD 0.01 Hz.
- Amplitudes are sqrt(S(f)) and phases are uniform on (0, 2π).
- LF/HF power ratio ~ U(1.5, 2.0), the Task Force 5-min norm. ECGSYN's default
  of 0.5 is not used.
- 1024 beats are synthesized and a random 256-beat segment is kept, so the
  window is not exactly periodic.

**N2 `power_law`**
- Gaussian 1/f^β noise with β ~ U(0.9, 1.1): complex Gaussian Fourier
  coefficients with |a(f)|² ∝ f^−β.
- 4096 beats are synthesized and a random 256-beat segment is kept.
- HRV spectra approximate 1/f [Kobayashi & Musha 1982].

**N3 `linear_rr_trend`**
- N1 plus a linear trend in mean RR.
- Total change across the window is U(5 %, 10 %) of μ, with random sign.
- The trend is centred, so the window mean is unchanged.

**N4 `linear_rr_step`**
- N1 plus an abrupt step in mean RR: U(5 %, 10 %) of μ, random sign.
- The step falls at a beat uniform in [0.3N, 0.7N] and is centred.

**N5 `linear_rr_warped`**
- N1 passed through z → exp(0.5 z), then rescaled to (μ, σ).
- Monotonic static nonlinearity giving right skew (median skewness 1.5).

**N6 `noisy_rsa`**
- A respiratory sinusoid (0.2–0.3 Hz, random phase) carrying 50 % of the
  variance, plus independent white Gaussian noise.
- Periodic plus noise, no chaos.

### Secondary nulls: reported, not in the PASS rule

**S1 `setar`**
- SETAR(2;1,1) [Tong & Lim 1980]:
  - x_t = 0.7 x_{t−1} + e_t if x_{t−1} ≤ 0;
  - x_t = −0.5 x_{t−1} + e_t otherwise;
  - e_t ~ N(0, 1).
- The deterministic skeleton has a stable fixed point at 0. The process is
  noise-driven, nonlinear and not chaotic.

**S2 `ectopic_{2,5,10}pct`**
- N1 with isolated premature beats: round(rate·256) = 5, 13 or 26 ectopics.
- Ectopic indices are at least 3 beats apart.
- Each ectopic interval is RR_i ← c·RR_i with c ~ U(0.6, 0.8).
- A **full compensatory pause** follows: RR_{i+1} ← RR_{i+1} + (1 − c)·RR_i, so
  the pair sums to the two sinus intervals and later beats keep their times.
- 10 % equals the MIT-BIH abnormal-window threshold used by the pipeline:
  `CLASSIFIER_ABNORMAL_FRACTION = 0.10` in `final_pipeline.py`.

**S3 `*_unq`**
- N1, N3 and N4 **without quantization**.
- They use the same base realization and the same modifier draws as their
  quantized twins (checked bit-for-bit), so each pair differs only in
  quantization.

### Positive conditions: chaotic, the detector should fire (power; reported)

- **P1 `henon_rr`**: Hénon x (a = 1.4, b = 0.3) [Hénon 1976], set to (μ, σ),
  then quantized. Phase 2E generator: random initial condition, 1000
  iterations discarded.
- **P2 `logistic_rr`**: logistic map r = 4 [May 1976], likewise.
- **P3**: P1 and P2 plus the N3 trend (the same trend distribution).
- **P4**: P1 and P2 plus white Gaussian measurement noise at 30 and 20 dB
  (relative to σ²), added before quantization.

### Generalization conditions for the fixed m = 2

All are set to (μ, σ) and quantized. RK4 integration throughout.

**G1 `lorenz_maxima`**
- Successive maxima of Lorenz z [Lorenz 1963]: σ = 10, ρ = 28, β = 8/3.
- dt = 0.005; 100 time units discarded.
- Each maximum is refined by a parabola through 3 samples.
- A beat-to-beat-like 1-D map.

**G2 `rossler_flow`**
- Rössler x [Rössler 1976]: a = b = 0.2, c = 5.7.
- dt = 0.01; 500 time units discarded.
- Sampled every 1.0 time unit, about 6 samples per mean orbital period of
  about 6.1.
- A flow, which likely needs m ≥ 3.

**G3 `mackey_glass`**
- Mackey–Glass [Mackey & Glass 1977]: dx/dt = 0.2 x(t−17)/(1 + x(t−17)^10) − 0.1 x.
- τ = 17 is in the chaotic regime [Farmer 1982].
- dt = 0.1, with the delayed term linearly interpolated at half steps.
- Constant initial history ~ U(0.5, 1.3); 1000 time units discarded.
- Sampled every 6.0 time units. This step was chosen for this study: about
  8 samples per quasi-period of about 50 time units.

### Seeds

Phase 2E scheme: `SeedSequence([20260926, code, 256, seed])`.
- The **base** stream (code by base process) gives μ, σ and the base process.
- The **modifier** stream (code by condition; S3 uses its twin's code) gives
  trend, step, ectopic and noise draws.
- Conditions on the same base with the same seed share the realization, which
  gives a paired design.
- Development seeds use the Phase 2E seed numbers: 0–99 for N* and S*, 0–29 for
  P* and G*.
- **Test seeds start at 3000**, with the count fixed in `PREREGISTRATION.md`.

## Realism targets and check

| Measure | Target | Source |
|---|---|---|
| Mean RR | 0.6–1.0 s for every window | task specification (resting rate 60–100 bpm) |
| Frequency bands | LF 0.04–0.15 Hz, HF 0.15–0.40 Hz | [Task Force 1996] |
| LF/HF | 1.5–2.0 (5-min supine) | [Task Force 1996] |
| LF / HF power, 5-min supine | 1170 ± 416 / 975 ± 203 ms²; total power 3466 ± 1018 ms²; LF n.u. 54 ± 4, HF n.u. 29 ± 3 | [Task Force 1996] |
| RMSSD | 27 ± 12 ms (24 h) | [Task Force 1996] |
| Short-term SDNN, RMSSD, LF/HF | 50 ± 16 ms; 42 ± 15 ms; 2.8 ± 2.6 (mean RR 926 ms) | [Nunan 2010] |

`realism.py` computes each measure per window:
- 4 Hz cubic-spline resampling;
- linear detrending;
- Welch with 64 s Hann segments.

The criteria are applied to the median over development windows:
- SDNN within 34–66 ms (Nunan ± 1 SD);
- RMSSD within 27–57 ms;
- LF/HF within 1.5–2.0 (Task Force) or 0.2–5.4 (Nunan ± 1 SD).

The full table is `results/realism/realism_dev.md`; the example plots are in
`plots/`:
- `examples_timeseries.png`;
- `examples_psd.png`, the median PSD over 10 seeds;
- `examples_return_maps.png`.

Summary (medians over development windows):

| Condition | Mean RR 0.6–1.0 s (all windows) | SDNN (ms) | RMSSD (ms) | LF/HF | Targets met |
|---|---|---|---|---|---|
| N1 linear_rr | yes | 47.7 | 37.1 | 1.68 | all |
| N2 power_law | yes | 45.8 | 35.8 | 1.20 | all except the Task Force LF/HF band (within Nunan) |
| N3 trend / N4 step / N5 warped | yes | 50.4 / 55.6 / 47.6 | 37–39 | 1.64–1.70 | all |
| N6 noisy_rsa | yes | 47.0 | 60.2 | 0.11 | mean and SDNN only: white noise and a respiratory tone make it HF-dominated **by design** |
| S1 setar | yes | 49.6 | 54.9 | 0.85 | all except the Task Force LF/HF band |
| S2 ectopic 2 / 5 / 10 % | yes | 68 / 91 / 119 | 92 / 138 / 191 | 0.99 / 0.57 / 0.36 | mean only: ectopics inflate SDNN and RMSSD, as in real recordings, which is why HRV guidelines require NN editing [Task Force 1996] |
| S3 unquantized | yes | as N1 / N3 / N4 | as N1 | as N1 | all (256 distinct values; about 70 when quantized) |
| P1–P4 Hénon / logistic | yes | 45–51 | 67–74 | 0.24–0.37 | mean and SDNN only |
| G1 Lorenz maxima | yes | 47.3 | 61.9 | 0.58 | mean and SDNN only |
| G2 Rössler | yes | 48.0 | 50.5 | 0.03 | mean, SDNN and RMSSD |
| G3 Mackey–Glass | yes | 48.4 | 39.6 | 0.34 | mean, SDNN and RMSSD |

**Realism verdict:**
- **The null conditions N1–N5 meet every target.** N6 meets the mean and SD
  targets. Its spectrum is intentionally a respiratory tone plus white noise.
- **The chaotic conditions match only the mean and SD.** Rescaling a map or
  flow preserves its dynamics and therefore its spectrum:
  - Hénon and logistic are anticorrelated beat to beat (HF-heavy, LF/HF about
    0.3);
  - Rössler at this sampling is LF-dominated.

  Matching the spectrum would require changing the dynamics, so these
  conditions test detection of known chaos at realistic RR scale and
  quantization, not realistic HRV.
- **Quantization at 1/360 s** leaves about 50–110 distinct values per window.

### physionet.org

**Not reachable** from this environment. `curl https://physionet.org/...`
returns `CONNECT tunnel failed, response 403` (a proxy policy denial; checked
2026-09-29). The optional comparison with NSRDB summary statistics was
therefore not done. The MIT-BIH phase needs physionet.org to be allowed, or the
data supplied another way. pypi.org is reachable; `wfdb` is not installed.

The literature values above were checked against web-search results, since the
full texts, on publisher and mirror sites, were blocked by the proxy:
- Task Force 5-min norms: total power 3466 ± 1018, LF 1170 ± 416,
  HF 975 ± 203 ms², LF/HF 1.5–2.0, LF n.u. 54 ± 4, HF n.u. 29 ± 3;
- Task Force 24-h norms: SDNN 141 ± 39, RMSSD 27 ± 12 ms;
- Nunan 2010: SDNN 50 ± 16, RMSSD 42 ± 15 ms, LF/HF 2.8, mean RR 926 ms.

## References

- [Task Force 1996] Task Force of the European Society of Cardiology and the
  North American Society of Pacing and Electrophysiology. Heart rate
  variability: standards of measurement, physiological interpretation and
  clinical use. *Circulation* 93(5):1043–1065, 1996; *Eur Heart J*
  17:354–381, 1996.
- [Nunan 2010] Nunan D, Sandercock GRH, Brodie DA. A quantitative systematic
  review of normal values for short-term heart rate variability in healthy
  adults. *Pacing Clin Electrophysiol* 33(11):1407–1417, 2010.
- [McSharry 2003] McSharry PE, Clifford GD, Tarassenko L, Smith LA. A dynamical
  model for generating synthetic electrocardiogram signals. *IEEE Trans Biomed
  Eng* 50(3):289–294, 2003.
- [Kobayashi & Musha 1982] Kobayashi M, Musha T. 1/f fluctuation of heartbeat
  period. *IEEE Trans Biomed Eng* 29(6):456–457, 1982.
- [Tong & Lim 1980] Tong H, Lim KS. Threshold autoregression, limit cycles and
  cyclical data. *J R Stat Soc B* 42(3):245–292, 1980.
- [Hénon 1976] Hénon M. A two-dimensional mapping with a strange attractor.
  *Commun Math Phys* 50:69–77, 1976.
- [May 1976] May RM. Simple mathematical models with very complicated dynamics.
  *Nature* 261:459–467, 1976.
- [Lorenz 1963] Lorenz EN. Deterministic nonperiodic flow. *J Atmos Sci*
  20:130–141, 1963.
- [Rössler 1976] Rössler OE. An equation for continuous chaos. *Phys Lett A*
  57:397–398, 1976.
- [Mackey & Glass 1977] Mackey MC, Glass L. Oscillation and chaos in
  physiological control systems. *Science* 197:287–289, 1977.
- [Farmer 1982] Farmer JD. Chaotic attractors of an infinite-dimensional
  dynamical system. *Physica D* 4:366–393, 1982.
