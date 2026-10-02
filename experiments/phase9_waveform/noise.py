"""
Phase 9 Part A3: recording conditions.

- MIT-BIH amplitude resolution: 200 adu/mV, 11-bit range (-1024 .. 1023 adu), i.e. samples
  quantized to 0.005 mV and clipped to [-5.12, 5.115] mV (MIT-BIH Arrhythmia Database
  directory: "11-bit resolution over a 10 mV range", 360 Hz).
- Real noise from the MIT-BIH Noise Stress Test Database (nstdb 1.0.0, records bw, ma, em,
  360 Hz, 2 signals, 200 adu/mV; files in nstdb_data/, SHA-256 checked against PhysioNet's
  SHA256SUMS.txt), added at a prescribed SNR with the definitions of the WFDB program nst
  (physionet.org/physiotools/wag/nst-1.htm, "Signal-to-noise ratios"):
    S = (peak-to-peak QRS amplitude)^2 / 8, the amplitude being the mean of the middle 90 % of the
        ranges within +/-50 ms of each of the first 300 normal QRS complexes;
    N = (RMS noise amplitude)^2, the amplitude being the mean of the middle 90 % of the RMS
        deviations from the mean in consecutive 1-s chunks of the first 300 s of noise;
    gain a solves SNR = 10 log10(S / (N a^2)).
  Here S is measured on the clean synthetic ECG of the window (normal beats only; all of them if
  fewer than 300) and N on the noise segment actually added (all of its 1-s chunks if < 300 s).
- Noise split (predeclared): DEVELOPMENT windows use nstdb signal 0 ('noise1'), TEST windows
  signal 1 ('noise2'), so no noise sample is shared between development and test.
- Mixture 'mix' (predeclared): bw, ma and em segments are each normalized to unit nst noise
  power N, summed, and the sum is scaled to the target SNR.
"""
from __future__ import annotations

import pathlib

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
NSTDB = HERE / "nstdb_data"
FS = 360
ADU_PER_MV = 200.0
ADU_MIN, ADU_MAX = -1024, 1023
_CACHE = {}


def load(kind):
    if kind not in _CACHE:
        import wfdb
        rec = wfdb.rdrecord(str(NSTDB / kind))
        assert rec.fs == FS
        _CACHE[kind] = np.asarray(rec.p_signal, dtype=float)
    return _CACHE[kind]


def quantize_mitbih(x_mV):
    adu = np.clip(np.round(np.asarray(x_mV) * ADU_PER_MV), ADU_MIN, ADU_MAX)
    return adu / ADU_PER_MV


def _trimmed_mean(v):
    v = np.sort(np.asarray(v, float))
    k = int(round(0.05 * len(v)))
    return float(v[k:len(v) - k].mean()) if len(v) - 2 * k > 0 else float(v.mean())


def signal_power(ecg_mV, r_idx, fs=FS):
    """nst S from the clean ECG and the sample indices of its normal R events."""
    w = int(round(0.05 * fs))
    amps = []
    for i in r_idx[:300]:
        if i - w >= 0 and i + w < len(ecg_mV):
            seg = ecg_mV[i - w:i + w + 1]
            amps.append(seg.max() - seg.min())
    return _trimmed_mean(amps) ** 2 / 8.0


def noise_power(noise_mV, fs=FS):
    n = min(len(noise_mV) // fs, 300)
    chunks = noise_mV[:n * fs].reshape(n, fs)
    rms = np.sqrt(np.mean((chunks - chunks.mean(1, keepdims=True)) ** 2, axis=1))
    return _trimmed_mean(rms) ** 2


def segment(kind, channel, n_samples, rng):
    """A random contiguous segment of nstdb record `kind` (signal `channel`)."""
    sig = load(kind)[:, channel]
    if n_samples > len(sig):
        raise ValueError("segment longer than the noise record")
    start = int(rng.integers(0, len(sig) - n_samples + 1))
    return sig[start:start + n_samples].copy(), start


def make_noise(spec, n_samples, rng, channel):
    """spec in ('bw', 'ma', 'em', 'mix').  Returns (noise_mV with nst power N = 1, info)."""
    if spec == "mix":
        parts, info = [], {}
        for kind in ("bw", "ma", "em"):
            s, st = segment(kind, channel, n_samples, rng)
            parts.append(s / np.sqrt(noise_power(s)))
            info[f"start_{kind}"] = st
        tot = np.sum(parts, axis=0)
        return tot / np.sqrt(noise_power(tot)), info
    s, st = segment(spec, channel, n_samples, rng)
    return s / np.sqrt(noise_power(s)), {f"start_{spec}": st}


def add_noise(ecg_mV, r_idx_normal, spec, snr_db, rng, channel):
    """ecg + a * noise with nst SNR = snr_db; returns (noisy, info)."""
    S = signal_power(ecg_mV, r_idx_normal)
    nz, info = make_noise(spec, len(ecg_mV), rng, channel)      # N = 1 by construction
    a = np.sqrt(S / 10.0 ** (snr_db / 10.0))
    info.update(spec=spec, snr_db=float(snr_db), S=float(S), gain=float(a), channel=int(channel))
    return ecg_mV + a * nz, info
