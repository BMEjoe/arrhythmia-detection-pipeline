"""
Phase 8 Part A: RR-interval windows of every model regime and variant.

Variants (predeclared; chain for (iv)-(vi) follows the PASS/WINNER rule, where (vi) is
"noise, quantization, ectopy, jitter"):
  i_clean        the model's RR (Section 6 rescaling in MODELS.md), unquantized
  ii_dyn_lo/hi   dynamical noise inside the model, two levels (DYN_NOISE)
  iii_meas30/20  i_clean + white measurement noise at 30 / 20 dB (Phase 5 add_noise; SNR
                 relative to the window's own variance)
  iv_q           quantize(i_clean + 30 dB noise) to 1/360 s
  v_S2_5 / v_E1 / v_E3_10
                 the same 30 dB noisy series with Phase 6 ectopy inserted before quantization:
                 S2 5 % isolated premature beats (Phase 5 add_ectopics), E1 bigeminy, E3 10 %
                 couplets (Phase 6 _apply_group / _place_groups), then quantized
  vi_S2_5 / vi_E1 / vi_E3_10
                 (v) plus the Phase 7 measured R-peak timing jitter at every premature
                 (ventricular-type) beat: offset (samples at 360 Hz) drawn from the 7,761
                 measured detected-minus-annotated offsets of MIT-BIH V/E/F beats
                 (experiments/phase7_mitbih/results/spike_in/v_offsets_samples.npy); the
                 interval ending at the beat gets +offset, the next one -offset; an offset
                 that would make either interval < 0.15 s is redrawn (up to 100 times, else 0).
Variants iii-vi share the clean realization; iv-vi share the 30 dB noise draw; v and vi
share the ectopic draw (paired design, as Phases 5-6).
"""
from __future__ import annotations

import json
import pathlib

import numpy as np

from experiments.phase5_rr import systems as S5
from experiments.phase6_robust import config as C6
from experiments.phase6_robust import systems as S6

HERE = pathlib.Path(__file__).resolve().parent
ENTROPY = 20261003
FAMILY_CODE = {"mackey_glass": 1, "phase_reset": 2, "coupled_vdp": 3, "av_node": 4}
VARIANTS = ("i_clean", "ii_dyn_lo", "ii_dyn_hi", "iii_meas30", "iii_meas20", "iv_q",
            "v_S2_5", "v_E1", "v_E3_10", "vi_S2_5", "vi_E1", "vi_E3_10")
VARIANT_CODE = {v: i + 1 for i, v in enumerate(VARIANTS)}
# dynamical-noise levels (predeclared, MODELS.md): (low, high)
DYN_NOISE = {"phase_reset": (0.002, 0.01),      # SD of Gaussian phase noise per stimulus (cycle units)
             "av_node": (0.5, 2.0),             # SD of H per beat (ms)
             "mackey_glass": (0.002, 0.01),     # additive Euler-Maruyama noise intensity
             "coupled_vdp": (0.02, 0.1)}        # additive noise on each oscillator's u''
JITTER_FILE = HERE.parent / "phase7_mitbih" / "results" / "spike_in" / "v_offsets_samples.npy"
_OFFS = None


def _offsets():
    global _OFFS
    if _OFFS is None:
        _OFFS = np.load(JITTER_FILE)
    return _OFFS


def ss(*key):
    return np.random.SeedSequence([ENTROPY, *[int(k) for k in key]])


def model_rr(family, params, n, rng, noise_sd=0.0):
    """Clean (or dynamically noisy) model RR, rescaled per MODELS.md Section 6."""
    if family == "phase_reset":
        from experiments.phase8_cardiac.models import phase_reset as M
        p = M.prc_params(params.get("aggregate", "A"))
        rr, _ = M.simulate_rr(params["tau"], n, p, phi0=rng.uniform(), noise_sd=noise_sd, rng=rng)
        return rr * (0.8 / rr.mean())
    if family == "av_node":
        from experiments.phase8_cardiac.models import av_node as M
        return M.simulate_rr(params["H"], n, noise_sd=noise_sd, rng=rng, A0=rng.uniform(100.0, 140.0))
    if family == "mackey_glass":
        from experiments.phase8_cardiac.models import mackey_glass as M
        return M.simulate_rr(params["tau"], n, rng, noise_sd=noise_sd)
    if family == "coupled_vdp":
        from experiments.phase8_cardiac.models import coupled_vdp as M
        return M.simulate_rr({"rho": params["rho"], "omega": params["omega"]}, n, rng, noise_sd=noise_sd)
    raise ValueError(family)


def _ectopy(x, kind, rng):
    """Phase 5/6 ectopy on an unquantized series; returns (series, premature interval indices)."""
    x = np.array(x, dtype=float)
    n = len(x)
    labels = np.zeros(n, dtype=bool)
    if kind == "S2_5":
        x, info = S5.add_ectopics(x, rng, 0.05)
        return x, list(info["ectopic_indices"])
    if kind == "E1":
        s = int(rng.integers(1, 3))
        starts = list(range(s, n - 1, 2))
        for i in starts:
            S6._apply_group(x, labels, i, 1, rng)
        return x, starts
    if kind == "E3_10":
        m = int(round(0.10 * n / 2))
        starts = S6._place_groups(rng, [3] * m, n)
        for i in starts:
            S6._apply_group(x, labels, i, 2, rng)
        return x, sorted([i for i in starts] + [i + 1 for i in starts])
    raise ValueError(kind)


JITTER_MIN_RR = 0.15   # s; an offset that would make either affected interval shorter is redrawn


def _jitter(x, premature, rng):
    """Guard (predeclared): R-peak timing error cannot reverse beat order; an offset that
    would make either affected interval < JITTER_MIN_RR is redrawn (up to 100 times, else 0).
    Matters only for the fast rhythms (AV node native ~0.18 s; tachycardic coupled-vdP regimes)."""
    x = np.array(x, dtype=float)
    offs = _offsets()
    for j in premature:
        d = 0.0
        for _ in range(100):
            c = int(rng.choice(offs)) / 360.0
            ok = x[j] + c >= JITTER_MIN_RR and (j + 1 >= len(x) or x[j + 1] - c >= JITTER_MIN_RR)
            if ok:
                d = c
                break
        x[j] += d
        if j + 1 < len(x):
            x[j + 1] -= d
    return x


def window(family, regime_idx, params, n, variant, seed):
    """One window (length n) of a model regime and variant.  Deterministic in its key."""
    fc = FAMILY_CODE[family]
    if variant.startswith("ii_dyn"):
        lvl = DYN_NOISE[family][0 if variant == "ii_dyn_lo" else 1]
        rng = np.random.default_rng(ss(fc, regime_idx, n, VARIANT_CODE[variant], seed))
        return model_rr(family, params, n, rng, noise_sd=lvl)
    clean = model_rr(family, params, n, np.random.default_rng(ss(fc, regime_idx, n, 0, seed)))
    if variant == "i_clean":
        return clean
    var = float(np.var(clean))
    if variant in ("iii_meas30", "iii_meas20"):
        snr = 30.0 if variant == "iii_meas30" else 20.0
        return S5.add_noise(clean, np.random.default_rng(ss(fc, regime_idx, n, VARIANT_CODE[variant], seed)),
                            snr, var)[0]
    noisy30 = S5.add_noise(clean, np.random.default_rng(ss(fc, regime_idx, n, 100, seed)), 30.0, var)[0]
    if variant == "iv_q":
        return S5.quantize(noisy30)
    kind = variant.split("_", 1)[1]
    ect, prem = _ectopy(noisy30, kind, np.random.default_rng(ss(fc, regime_idx, n, 200 + VARIANT_CODE["v_" + kind], seed)))
    q = S5.quantize(ect)
    if variant.startswith("v_"):
        return q
    return _jitter(q, prem, np.random.default_rng(ss(fc, regime_idx, n, 300 + VARIANT_CODE["v_" + kind], seed)))


def load_regimes():
    return json.load(open(HERE / "results" / "ground_truth" / "regimes.json"))
