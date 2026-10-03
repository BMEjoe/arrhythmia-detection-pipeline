"""
Phase 10 shared harness: window sampling, detector wrappers, ectopy masking / editing, spike-in
construction (Question 1), 12-minute titration segments (Question 2), synthetic titration checks
(Question 2e).  Everything here is deterministic in its inputs (RNG streams keyed by task identity or by
a CRC of the data, as the frozen detectors do).

Detectors (frozen; never retuned)
  K1   experiments.phase7_mitbih.detector.evaluate (Phase 6 combined AND detector,
       keep_upo_on_short_lle_embedding = True, D2 detrend) -> and / lle / upo decisions.
  K3   final_pipeline.masked_growth_chaos_test(rr, labels, CFG).
  K4   Phase 9 candidate k4_mnlp_growth_rr (RR-rule mask, its frozen thresholds), Question 4d.
  TIT  titration10 (verified noise titration), arms raw / masked (K3 rule) / edited (Phase 6 edit_nn).
"""
from __future__ import annotations

import zlib
from dataclasses import replace

import numpy as np

ENTROPY = 20261010             # Phase 10 RNG entropy (distinct from every earlier phase's)
LO, HI = 0.25, 3.0             # pipeline hard RR limits (fp.CFG.rr_min_seconds / rr_max_seconds)
N_PER = 12                     # windows per subject
L_PRIMARY = 512
NN_EDIT_MIN = 0.80             # Phase 6/7 editing eligibility (NN fraction)
ROBUST_J = (1, 4, 7, 10)       # Question 4a/4b subset of window indices
SEG_SECONDS = 720.0            # Question 2c: 12-min segments (Wu et al. 2009)


def crc(x):
    return zlib.crc32(np.ascontiguousarray(np.asarray(x, dtype=float)).tobytes())


# ------------------------------------------------------------------------------------ windows
def _clean_start(t, i, L, i_limit):
    """First start >= i whose L intervals all lie in [LO, HI] and that ends before i_limit."""
    n = len(t)
    while i + L < n and i + L <= i_limit:
        rr = np.diff(t[i:i + L + 1])
        bad = np.nonzero((rr < LO) | (rr > HI))[0]
        if not len(bad):
            return int(i)
        i = i + int(bad[-1]) + 1
    return None


def window_starts(t, n_per=N_PER, L=L_PRIMARY):
    """Predeclared rule: window j (1..n_per) nominally starts at the first beat at or after
    T0 + (j - 0.5) / n_per * (T_end - T0); if its L intervals contain one outside [LO, HI] s the start
    moves to the beat after the last such interval (repeatedly); the window must end before the next
    window's nominal start (else it is missing)."""
    T0, T1 = float(t[0]), float(t[-1])
    noms = [int(np.searchsorted(t, T0 + (j - 0.5) / n_per * (T1 - T0))) for j in range(1, n_per + 1)]
    out = {}
    for j in range(1, n_per + 1):
        limit = noms[j] if j < n_per else len(t) - 1
        s = _clean_start(t, noms[j - 1], L, limit)
        if s is not None:
            out[j] = s
    return out


def clock_hour(base_time, t_rel):
    """Clock hour (0-24) of a time t_rel seconds after the record's base time 'HH:MM:SS'."""
    if base_time is None:
        return None
    h, m, s = (float(v) for v in str(base_time).split(":"))
    return ((h * 3600 + m * 60 + s + t_rel) / 3600.0) % 24.0


# ------------------------------------------------------------------------------------ masks
def k3_mask(labels):
    bad = np.asarray(labels).astype(str) != "N"
    return bad[:-1] | bad[1:]


def edit(rr, labels=None, imask=None):
    """Phase 6 edit_nn with the Phase 7 / K3 mask (an interval touching a non-normal beat); the
    interval mask may be given directly (imask)."""
    from experiments.phase6_robust.systems import edit_nn
    m = k3_mask(labels) if imask is None else np.asarray(imask, bool)
    if not m.any():
        return np.asarray(rr, float).copy(), 1.0
    if (~m).sum() < 2:
        return None, float(1 - m.mean())
    x, info = edit_nn(rr, m)
    return x, info["nn_fraction"]


def rr_rule_labels(rr):
    """Labels from the pipeline's RR outlier rule (Phase 8 methods.outlier_mask): the beat ending a
    flagged interval is marked 'X' (then K3's rule masks that interval and the next)."""
    from experiments.phase8_cardiac import methods as MM
    f = MM.outlier_mask(np.asarray(rr, float))
    lab = np.array(["N"] * (len(rr) + 1))
    lab[1:][f] = "X"
    return lab


# ------------------------------------------------------------------------------------ detectors
K1_KEYS = ("and_detected", "lle_detected", "upo_detected", "lle_p", "lle_z", "lle_stat", "upo_score",
           "upo_status", "error", "detrend_applied")


def k1(rr, n_sur=None):
    from experiments.phase7_mitbih import detector as DT
    cfg = DT.DETECTOR if n_sur is None else replace(DT.DETECTOR, lle_chaos_test_surrogates=n_sur,
                                                    so_surrogate_count=n_sur)
    r = DT.evaluate(np.asarray(rr, float), cfg)
    return {k: r.get(k) for k in K1_KEYS}


def k3(rr, labels, n_sur=None):
    import final_pipeline as fp
    cfg = fp.CFG if n_sur is None else replace(fp.CFG, masked_growth_surrogates=n_sur)
    try:
        r = fp.masked_growth_chaos_test(np.asarray(rr, float), np.asarray(labels), cfg)
    except Exception as exc:                                        # noqa: BLE001
        return {"detected": False, "analysable": False, "error": f"{type(exc).__name__}: {exc}"}
    return {k: r.get(k) for k in ("detected", "analysable", "n_masked", "zmax", "G")}


def k4(rr):
    from experiments.phase9_waveform import candidates9 as C9
    r = C9.k4_mnlp_growth_rr(np.asarray(rr, float), None, {})
    return {k: r.get(k) for k in ("detected", "analysable", "n_masked", "zmax", "G")}


def titration(rr, labels=None, arm="raw", salt=0, imask=None):
    """arm: raw (all intervals), masked (K3 rule rows), edited (Phase 6 edit; NN fraction >= 0.80).
    The ectopy-related interval mask is k3_mask(labels), or imask if given."""
    from experiments.phase10_final import titration10 as T
    x = np.asarray(rr, float)
    if arm != "raw":
        m = k3_mask(labels) if imask is None else np.asarray(imask, bool)
    rng = np.random.default_rng(np.random.SeedSequence([ENTROPY, crc(x), 77, salt]))
    if arm == "raw":
        o = T.titrate(x, rng)
    elif arm == "masked":
        o = T.titrate(x, rng, valid_rows=T.rows_from_mask(m))
    elif arm == "edited":
        xe, nnf = edit(x, imask=m)
        if xe is None or nnf < NN_EDIT_MIN:
            return {"eligible": False, "nn_fraction": nnf, "positive": False, "NL": 0.0}
        o = T.titrate(xe, rng)
        o["nn_fraction"] = nnf
    else:
        raise ValueError(arm)
    b = o["base"]
    return {"eligible": True, "analysable": bool(b.get("analysable", False)), "positive": bool(o["positive"]),
            "NL": float(o["NL"]), "p": b.get("p"), "r_lin": b.get("r_lin"), "r_nl": b.get("r_nl"),
            "n_rows": b.get("n_rows"), "nn_fraction": o.get("nn_fraction")}


# ------------------------------------------------------------------------------------ spike-ins (Q1)
FAMILIES = {
    # family: list of (param label, params, lambda, lambda unit)
    "henon": [("a1.08", {"a": 1.08, "b": 0.3}, 0.1364, "per beat"), ("a1.14", {"a": 1.14, "b": 0.3}, 0.2442, "per beat"),
              ("a1.22", {"a": 1.22, "b": 0.3}, 0.3029, "per beat"), ("a1.40", {"a": 1.40, "b": 0.3}, 0.4193, "per beat")],
    "logistic": [("r3.58", {"r": 3.58}, 0.1052, "per beat"), ("r3.65", {"r": 3.65}, 0.2554, "per beat"),
                 ("r3.88", {"r": 3.88}, 0.4642, "per beat"), ("r4.00", {"r": 4.0}, 0.6931, "per beat")],
    "coupled_vdp": [("rho2_om2.7", 22, 0.0850, "per model time unit"), ("rho6_om3.3", 23, 0.0788, "per model time unit"),
                    ("rho6_om4.0", 24, 0.1012, "per model time unit"), ("rho8_om3.3", 25, 0.0751, "per model time unit")],
    "phase_reset": [("tau1.14", 12, 0.0354, "per stimulus"), ("tau0.58", 10, 0.0810, "per stimulus"),
                    ("tau1.20", 14, 0.1379, "per stimulus"), ("tau1.16", 13, 0.1848, "per stimulus")],
}
FAM_CODE = {"henon": 1, "logistic": 2, "coupled_vdp": 3, "phase_reset": 4}
F_GRID = (0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9)


def chaos_component(family, pidx, n, key):
    """Standardised (mean 0, SD 1) chaotic sequence of length n, keyed by `key` (ints)."""
    rng = np.random.default_rng(np.random.SeedSequence([ENTROPY, 1, FAM_CODE[family], pidx, *key]))
    lab, params, _, _ = FAMILIES[family][pidx]
    if family == "henon":
        x, y = rng.uniform(-0.1, 0.1), rng.uniform(-0.1, 0.1)
        out = np.empty(n)
        for i in range(2000 + n):
            x, y = 1 - params["a"] * x * x + y, params["b"] * x
            if i >= 2000:
                out[i - 2000] = x
    elif family == "logistic":
        x = rng.uniform(0.2, 0.8)
        out = np.empty(n)
        for i in range(2000 + n):
            x = params["r"] * x * (1 - x)
            if i >= 2000:
                out[i - 2000] = x
    else:
        from experiments.phase8_cardiac import series as SR8
        reg = next(r for r in SR8.load_regimes() if r["regime_idx"] == params)
        out = np.asarray(SR8.model_rr(reg["family"], reg["params"], n, rng), float)
    sd = out.std()
    return (out - out.mean()) / (sd if sd > 0 else 1.0)


def quantize_times(t0, rr, fs):
    """Beat times t0 + cumsum(rr) rounded to the 1/fs grid; returns the quantized intervals."""
    T = np.round((t0 + np.concatenate([[0.0], np.cumsum(rr)])) * fs) / fs
    return np.diff(T)


def spike_additive(rr, labels, c, f, fs):
    """x' = x + s c, var(s c) / (var(s c) + var_NN(x)) = f; var_NN over the window's non-masked (K3)
    intervals; then quantized to the database resolution 1/fs.  Labels unchanged."""
    rr = np.asarray(rr, float)
    m = k3_mask(labels)
    v = float(np.var(rr[~m])) if (~m).sum() >= 2 else float(np.var(rr))
    s = np.sqrt(f / (1.0 - f) * v)
    x = rr + s * c
    q = quantize_times(0.0, x, fs)
    return q, {"scale_s": float(s), "var_nn": v, "n_out_of_range": int(np.sum((q < LO) | (q > HI))),
               "min_rr": float(q.min())}


def spike_replacement(rr, labels, symbols, c, fs, rng):
    """Phase 7 replacement design: chaos rescaled to the window's NN mean and SD, the window's real
    ectopic beats inserted at their real positions with the Phase 6 Part B model (Phase 7
    spike_in.insert_ectopy / groups), quantized to 1/fs."""
    from experiments.phase7_mitbih import spike_in as SP7
    rr = np.asarray(rr, float)
    m = k3_mask(labels)
    nn = rr[~m] if (~m).sum() >= 2 else rr
    mean, sd = float(nn.mean()), float(nn.std())
    capped = False
    if mean + sd * c.min() < SP7.RR_FLOOR:
        sd, capped = (mean - SP7.RR_FLOOR) / (-c.min()), True
    x = mean + sd * c
    abn = (np.asarray(labels)[1:] != "N").tolist()
    grp = SP7.groups(abn, list(np.asarray(symbols)[1:]))
    xe, n_floor = SP7.insert_ectopy(x, grp, rng)
    q = quantize_times(0.0, xe, fs)
    return q, {"sd_used": sd, "sd_capped": capped, "n_groups": len(grp), "n_comp_floored": n_floor}


# ------------------------------------------------------------------------------------ 12-min segments (Q2c)
SEG_MIN_COVER = 0.90           # in-range intervals must cover >= 90 % of the 12 min
SEG_MIN_N = 400                # and >= 400 in-range intervals must remain


def segments(t, lab):
    """Consecutive 12-min segments from the first beat; per segment the indices of the beats whose
    interval to the next beat STARTS inside the segment.  Intervals outside [LO, HI] s (missing /
    spurious annotations, data gaps) are removed in every arm; a segment is analysable iff the
    in-range intervals cover >= SEG_MIN_COVER of the 12 min and >= SEG_MIN_N remain ('gap' = not)."""
    T0 = float(t[0])
    k = np.floor((t[:-1] - T0) / SEG_SECONDS).astype(int)
    out = []
    nseg = int(np.floor((t[-1] - T0) / SEG_SECONDS))
    for s in range(nseg):
        idx = np.nonzero(k == s)[0]
        if len(idx) < 2:
            continue
        i0, i1 = int(idx[0]), int(idx[-1]) + 1       # beats i0..i1, intervals i0..i1-1
        rr = np.diff(t[i0:i1 + 1])
        keep = (rr >= LO) & (rr <= HI)
        cover = float(rr[keep].sum()) / SEG_SECONDS
        out.append({"seg": s, "i0": i0, "i1": i1, "n": int(keep.sum()), "n_removed": int((~keep).sum()),
                    "cover": cover, "gap": bool(cover < SEG_MIN_COVER or keep.sum() < SEG_MIN_N)})
    return out


def wu_nn(rr, imask):
    """[Wu09] preprocessing: elimination (without interpolation) of every interval that involves a
    non-normal beat; the remaining NN intervals are concatenated."""
    return np.asarray(rr, float)[~np.asarray(imask, bool)]


# ------------------------------------------------------------------------------------ synthetic (Q2e)
VDP_NONCHAOTIC = (26, 27, 28, 29, 30, 31)
ECT = ("none", "S2_5", "E1", "E3_10")


def synthetic_window(kind, cond, seed, n):
    """kind 'null': Phase 5-6 PASS null (families.null_window, 1/360 s as generated);
    kind 'vdp': non-chaotic coupled-vdP regime (families.model_window) with ectopy, quantized 1/360 s."""
    from experiments.phase9_waveform import families as F9
    from experiments.phase5_rr import systems as S5
    if kind == "null":
        rr, types = F9.null_window(cond, seed, n)
        return np.asarray(rr, float), types
    ridx, ect = cond
    reg = next(r for r in F9.load_regimes() if r["regime_idx"] == ridx)
    rr, types, _ = F9.model_window(reg["family"], ridx, reg["params"], n, ect, seed)
    return np.asarray(S5.quantize(rr), float), types
