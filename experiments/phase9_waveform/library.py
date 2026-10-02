"""
Phase 9 Part A: one synthetic recording per window key, from which every input level is derived
(RR, per-beat features, waveform segments).

Window key: (source, name, ectopy, n, seed, noise, split)
  source  'null' (name = Phase 5-6 null), 'model' (name = Phase 8 regime name), 'ktz'
          (name = 'ktz:P=<P_nom>'), 'periodic' (name = 'periodic:<RR s>', constant RR; Part B
          waveform null)
  ectopy  'none' | 'S2_5' | 'E1' | 'E3_10' (models and KTz; nulls carry their own ectopy)
  n       beats in the beat-level window (intervals); waveform windows are cut from its start
  noise   'none' or '<spec><snr>' with spec in mix, bw, ma, em and snr in 24, 12, 6 (dB)
  split   'DEV' or 'TEST' (selects the nstdb signal: DEV 0, TEST 1)

Construction (predeclared):
 1. RR source (families.py) with n_gen = n + MARGIN intervals; the window is the first n.
 2. Lead-in: LEAD = 12 beats at the window's mean RR (type N) precede the first window beat so
    ECGSYN's z transient (relaxation rate 1/s) has decayed; ECGSYN starts at the first lead beat
    with z = 0.04 (ecgsyn.c) and a respiratory baseline phase ~ U(0, 2 pi) per window.
 3. ECGSYN (ecgsyn.py) at fs = 360 Hz (internal 720 Hz, RK4), per-beat morphology
    (morphology.py), fixed mV map.
 4. Optional nstdb noise at nst SNR (noise.py), S measured on the clean record's normal beats.
 5. MIT-BIH amplitude quantization (0.005 mV, 11-bit range).
 6. Beat times:
      'true'      R events of the generator on the sample grid (round(t * fs)), as annotations;
      'detected'  fp.detect_r_peaks (the project's Pan-Tompkins detector) on the recorded ECG;
                  detections are matched to true R events within 150 ms (nearest-unused, as
                  Phase 7) to identify ventricular ectopic beats.
    Jitter (detected beats only, 'jit'): at every detection matched to a true 'V' beat, the
    detection is moved by an offset drawn from the 7,761 Phase 7 measured detected-minus-
    annotated offsets of MIT-BIH V/E/F beats (phase7_mitbih/results/spike_in/
    v_offsets_samples.npy); an offset that would make either adjacent interval < 0.15 s is
    redrawn (up to 100 times, else 0) - the Phase 8 guard.
 7. Window RR: the first n intervals starting at the first beat at or after the first true window
    beat (minus 75 ms for detected beats).  If fewer than n intervals are available the window
    has an error (counted as not detected by every method).
"""
from __future__ import annotations

import hashlib

import numpy as np

from experiments.phase9_waveform import ecgsyn as E
from experiments.phase9_waveform import families as F
from experiments.phase9_waveform import morphology as MO
from experiments.phase9_waveform import noise as NZ

FS = 360
MARGIN = 32
LEAD = 12
MATCH_TOL_S = 0.150
JITTER_MIN_RR = 0.15
NOISE_CODE = {"none": 0}
for _i, _s in enumerate(("mix", "bw", "ma", "em")):
    for _j, _d in enumerate((24, 12, 6)):
        NOISE_CODE[f"{_s}{_d}"] = 1 + 3 * _i + _j
SPLIT_CHANNEL = {"DEV": 0, "TEST": 1}
_KTZ = None


def ktz_labels():
    global _KTZ
    if _KTZ is None:
        import json
        _KTZ = {(r["P_nom"], r["input"]): r for r in json.load(open(F.HERE / "results" / "ground_truth" / "ktz_labels.json"))}
    return _KTZ


def _source(source, name, ectopy, n_gen, seed):
    info = {}
    apd = apd_ref = None
    if source == "periodic":
        # Part B waveform null: strictly periodic beats (constant RR = name's value, e.g. 'periodic:0.8')
        rr = np.full(n_gen, float(name.split(":")[1]))
        types = np.array(["N"] * (n_gen + 1))
    elif source == "null":
        rr, types = F.null_window(name, seed, n_gen)
    elif source == "model":
        reg = {r["name"]: r for r in F.load_regimes()}[name]
        rr, types, info = F.model_window(reg["family"], reg["regime_idx"], reg["params"], n_gen, ectopy, seed)
    elif source == "devmap":
        # DEVELOPMENT-only circle-map control (devmaps.py)
        from experiments.phase9_waveform import devmaps as DM
        r = DM.BY_NAME[name]
        clean = DM.circle_rr(r["K"], r["Om"], n_gen, np.random.default_rng(F.ss9(F.FAMILY9["devmap"], hash_name(name), n_gen, 0, seed)))
        rr, types = F._insert_ectopy(clean, ectopy, np.random.default_rng(F.ss9(F.FAMILY9["devmap"], hash_name(name), n_gen,
                                                                                 F.ECT_CODE[ectopy], seed)))
    elif source == "ktz":
        P = int(name.split("=")[1])
        lab = ktz_labels()[(P, ectopy)]
        rr, types, apd, info = F.ktz_window(P, P, ectopy, n_gen, seed, apd_ref_ts=lab["apd_ref_ts"])
        apd_ref = info["apd_ref_s"]
    else:
        raise ValueError(source)
    return np.asarray(rr, float), np.asarray(types), apd, apd_ref, info


def hash_name(name):
    return int(hashlib.sha256(name.encode()).hexdigest()[:8], 16)


def _offsets():
    from experiments.phase8_cardiac import series as SR8
    return SR8._offsets()


def match(det, true, tol):
    """Nearest-unused matching of detections to true events (indices); returns det -> true map."""
    out = -np.ones(len(det), dtype=int)
    used = np.zeros(len(true), dtype=bool)
    j0 = 0
    for i, d in enumerate(det):
        while j0 < len(true) and true[j0] < d - tol:
            j0 += 1
        best, bd = -1, tol + 1
        for j in range(j0, min(j0 + 3, len(true))):
            if not used[j] and abs(true[j] - d) <= tol and abs(true[j] - d) < bd:
                best, bd = j, abs(true[j] - d)
        if best >= 0:
            out[i] = best
            used[best] = True
    return out


def window_seed(source, name, ectopy, n, seed, noise):
    h = int(hashlib.sha256(f"{source}|{name}".encode()).hexdigest()[:8], 16)
    return F.ss9(7, h, F.ECT_CODE[ectopy], n, NOISE_CODE[noise], seed)


def build(source, name, ectopy, n, seed, noise="none", split="DEV", beats="true", jitter=False,
          keep_clean=False, keep_truth=False):
    """Returns a dict with the recorded ECG (mV, fs 360), true and used beat sample indices, the
    window RR (n intervals), per-beat types of the used beats, and diagnostics."""
    n_gen = n + MARGIN
    rr, types, apd, apd_ref, info = _source(source, name, ectopy, n_gen, seed)
    rng = np.random.default_rng(window_seed(source, name, ectopy, n, seed, noise))
    m = float(np.mean(rr[:n]))
    rr_full = np.concatenate([np.full(LEAD, m), rr])
    types_full = np.concatenate([np.array(["N"] * LEAD), types])
    apd_full = None
    if apd is not None:
        apd_full = np.concatenate([np.full(LEAD, apd_ref), apd])
    bt = np.concatenate([[0.0], np.cumsum(rr_full)])                 # K = LEAD + n_gen + 1 beats
    th, a, b = MO.beat_params(rr_full, types_full, apd_full, apd_ref)
    t_end = bt[-1] + 0.5 * m
    bt_ext = np.concatenate([bt, [bt[-1] + m]])
    th = np.vstack([th, th[-1:]]); a = np.vstack([a, a[-1:]]); b = np.vstack([b, b[-1:]])
    ph = float(rng.uniform(0.0, 2.0 * np.pi))
    z = E.synthesize(bt_ext, th, a, b, fs=FS, q=2, t0=0.0, t1=t_end, resp_phase=ph)
    ecg = E.to_mV(z)
    true_idx = np.round(bt * FS).astype(int)
    true_idx = true_idx[true_idx < len(ecg)]
    out = {"fs": FS, "n": n, "true_idx": true_idx, "types_true": types_full[:len(true_idx)],
           "first_window_beat": LEAD, "info": info, "resp_phase": ph}
    if keep_clean:
        out["ecg_clean"] = ecg.copy()
    if keep_truth:
        out["truth"] = true_fiducials(bt, th, b)
    if noise != "none":
        spec = "".join(c for c in noise if not c.isdigit())
        snr = int("".join(c for c in noise if c.isdigit()))
        normal = true_idx[(types_full[:len(true_idx)] == "N")]
        normal = normal[normal >= LEAD]
        ecg, ninfo = NZ.add_noise(ecg, normal, spec, snr, rng, SPLIT_CHANNEL[split])
        out["noise"] = ninfo
    ecg = NZ.quantize_mitbih(ecg)
    out["ecg"] = ecg
    t0_idx = true_idx[LEAD]
    if beats == "true":
        used = true_idx
        used_types = out["types_true"]
    else:
        import final_pipeline as fp
        det = np.asarray(fp.detect_r_peaks(ecg, FS), dtype=int)
        mp = match(det, true_idx, int(round(MATCH_TOL_S * FS)))
        used_types = np.where(mp >= 0, out["types_true"][np.maximum(mp, 0)], "X")   # X = unmatched
        out["det_offsets_V"] = [int(det[i] - true_idx[mp[i]]) for i in range(len(det)) if mp[i] >= 0
                                and out["types_true"][mp[i]] == "V"]
        out["det_offsets_N"] = [int(det[i] - true_idx[mp[i]]) for i in range(len(det)) if mp[i] >= 0
                                and out["types_true"][mp[i]] == "N"]
        w = (true_idx >= t0_idx) & (true_idx <= true_idx[min(LEAD + n, len(true_idx) - 1)])
        inw = (det >= t0_idx - int(0.075 * FS))
        out["n_missed"] = int(np.sum(w) - np.sum(np.isin(np.nonzero(w)[0], mp[inw & (mp >= 0)])))
        out["n_extra"] = int(np.sum((mp < 0) & inw & (det <= true_idx[min(LEAD + n, len(true_idx) - 1)])))
        used = det.copy()
        if jitter:
            offs = _offsets()
            used = used.astype(float)
            for i in np.nonzero(used_types == "V")[0]:
                d = 0.0
                for _ in range(100):
                    c = float(rng.choice(offs))
                    okp = i == 0 or (used[i] + c - used[i - 1]) / FS >= JITTER_MIN_RR
                    okn = i == len(used) - 1 or (used[i + 1] - used[i] - c) / FS >= JITTER_MIN_RR
                    if okp and okn:
                        d = c
                        break
                used[i] += d
            used = np.round(used).astype(int)
        out["det_idx"] = det
    start = int(np.searchsorted(used, t0_idx - (0 if beats == "true" else int(0.075 * FS))))
    if start + n >= len(used):
        out["error"] = "fewer than n intervals"
        out["rr"] = np.diff(used[start:]) / FS
        out["beat_idx"] = used[start:]
        out["beat_types"] = used_types[start:]
        return out
    sel = slice(start, start + n + 1)
    out["beat_idx"] = used[sel]
    out["beat_types"] = used_types[sel]
    out["rr"] = np.diff(used[sel]) / FS
    out["rr_true"] = np.diff(true_idx[LEAD:LEAD + n + 1]) / FS
    if apd_full is not None:
        out["apd_true_s"] = apd_full[LEAD:LEAD + n + 1]
    return out


def true_fiducials(bt, th, b):
    """Analytic fiducials (seconds) of every beat from the generator parameters (multivariate.py
    docstring): QRS onset = t_Q - 2 sigma_Q, QRS offset = t_S + 2 sigma_S, T peak = t_T,
    T end = t_T + 2 sigma_T, with t_i = t_R + theta_i / omega and sigma_i = b_i / omega, omega = 2 pi / RR
    of the half-cycle containing the event (before R for Q, after R for S and T)."""
    K = len(bt)
    out = {k: np.full(K, np.nan) for k in ("qrs_on", "qrs_off", "t_peak", "t_end", "r")}
    for k in range(1, K - 1):
        w_pre = 2 * np.pi / (bt[k] - bt[k - 1])
        w_post = 2 * np.pi / (bt[k + 1] - bt[k])
        tq, sq = bt[k] + th[k, 1] / w_pre, b[k, 1] / w_pre
        ts, ss = bt[k] + th[k, 3] / w_post, b[k, 3] / w_post
        tt, st = bt[k] + th[k, 4] / w_post, b[k, 4] / w_post
        out["qrs_on"][k] = tq - 2 * sq
        out["qrs_off"][k] = ts + 2 * ss
        out["t_peak"][k] = tt
        out["t_end"][k] = tt + 2 * st
        out["r"][k] = bt[k]
    return out


def waveform(win, minutes):
    """Waveform window: `minutes` of the recorded ECG starting 0.5 s before the first window beat."""
    s = max(0, int(win["true_idx"][win["first_window_beat"]] - 0.5 * FS))
    e = s + int(round(minutes * 60 * FS))
    if e > len(win["ecg"]):
        return None
    return win["ecg"][s:e]
