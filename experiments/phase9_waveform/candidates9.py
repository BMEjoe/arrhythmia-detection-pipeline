"""
Phase 9 Part E candidates (designed on DEVELOPMENT material only; dev_diag.py / DIAGNOSIS.md).

Every candidate maps a 512-interval RR window (plus, for the annotation-masked ones, the beat
labels of its 513 beats) to a record with a boolean "detected".  Errors count as not detected.

Building blocks
  mask_ann(types)      interval k (beat k -> beat k+1) is ECTOPY-RELATED iff beat k or beat k+1 is
                       not a normal beat ('N'); labels are the beat annotations ('V', 'A', or 'X' =
                       a detection not matched to any annotated beat).  Peltola, Front Physiol 3:148
                       (2012), PMC3358711 (read via NCBI efetch): "Normally, both the ectopic beat and
                       the following compensatory pause are edited."
  mask_rr(x)           intervals flagged by the pipeline's own RR outlier rule (Phase 8
                       methods.outlier_mask = fp.correct_rr_intervals), no annotations.
  masked NLP           the Phase 8 robust normalized nonlinear-prediction error (k = 5 nearest
                       neighbours, Theiler 5, unit delay), computed only on delay vectors whose span
                       from the first coordinate to the predicted value contains no masked interval
                       (neighbour search restricted to the same valid vectors) -- the deletion
                       editing reviewed by Peltola (2012) ("the abnormal R-R intervals are removed"),
                       applied to delay vectors instead of concatenating the series (which, per
                       Peltola, introduces step-like artefacts); no interval is invented.  Normalized
                       by the MAD of the unmasked intervals.
                       Horizons h = 1..5; one valid-vector set for all h (span (m-1) + 5).
  masked IAAFT         surrogate = IAAFT (fp.iaaft_surrogate) of the series with masked intervals
                       linearly interpolated (as Phase 8 ep_iaaft); the statistic is computed with
                       the SAME mask (identical processing of data and surrogates).
                       NULL: the unmasked intervals are a monotone static transform of a stationary
                       linear Gaussian process; the masked (ectopy-related) intervals are ignored.
  growth G             G = E(h = 5) - E(h = 2) of the masked NLP at m = 3: the growth of the forecast
                       error with prediction time, the chaos signature of nonparametric forecasting
                       (Sugihara & May, Nature 344:734 (1990) -- primary not accessible here; restated
                       by Mühlbauer et al., Ecol Evol 2022, PMC9760897, who estimate Lyapunov
                       exponents from forecasts 1-10 steps ahead).  Periodic and quasi-periodic
                       dynamics (and their noisy versions) keep the error flat after the first steps
                       (development diagnosis, DIAGNOSIS.md); the threshold is calibrated on
                       development non-chaotic regimes, not taken from the paper.
  FSLE gate            FSLE (measures.fsle, m 3, delta0 0.01 SD, r sqrt 2) at delta_10 = 0.32 SD of the
                       series with masked intervals linearly interpolated (Part D, D2).
Minimum number of valid delay vectors: 100 (else not analysable -> not detected).
"""
from __future__ import annotations

import traceback
import zlib

import numpy as np

HMAX = 5
MIN_VALID = 100
N_SUR = 39
MS = (2, 3, 4, 5)

# Calibrated on DEVELOPMENT data only (calibrate9.py, results/dev/calibration9.json); set before the
# preregistration.  None until calibrated.
Z_DET = None       # masked-NLP determinism threshold (annotation mask)
Z_DET_RR = None    # same, RR-rule mask
G_MIN = None       # growth gate (annotation mask)
G_MIN_RR = None    # growth gate (RR-rule mask)
F_MIN = None       # FSLE gate


def rng_for(x, salt):
    h = zlib.crc32(np.ascontiguousarray(np.asarray(x, dtype=float)).tobytes())
    return np.random.default_rng(np.random.SeedSequence([20261006, int(h), int(salt)]))


def mask_ann(types):
    t = np.asarray(types)
    bad = t != "N"
    return bad[:-1] | bad[1:]


def mask_rr(x):
    from experiments.phase8_cardiac import methods as MM
    return MM.outlier_mask(x)


def interp_masked(x, mask):
    x = np.asarray(x, float).copy()
    if mask.any() and (~mask).sum() >= 2:
        idx = np.arange(len(x))
        x[mask] = np.interp(idx[mask], idx[~mask], x[~mask])
    return x


def valid_vectors(mask, m, hmax=HMAX):
    """Start indices t whose span t .. t + m - 1 + hmax contains no masked interval."""
    n = len(mask)
    span = m + hmax
    c = np.concatenate([[0], np.cumsum(mask.astype(int))])
    t = np.arange(0, n - span + 1)
    return t[(c[t + span] - c[t]) == 0]


def masked_nlp(x, mask, m, hs=(1, 2, 3, 4, 5), k=5, theiler=5):
    """Masked robust NLP errors {h: E}; None if fewer than MIN_VALID valid vectors."""
    x = np.asarray(x, float)
    t = valid_vectors(mask, m)
    if len(t) < MIN_VALID:
        return None
    A = np.stack([x[t + j] for j in range(m)], axis=1)
    d = ((A[:, None, :] - A[None, :, :]) ** 2).sum(-1)
    d[np.abs(t[:, None] - t[None, :]) <= theiler] = np.inf
    nn = np.argpartition(d, k, axis=1)[:, :k]
    xs = x[~mask]
    mad = max(np.median(np.abs(xs - np.median(xs))), 1e-12)
    out = {}
    for h in hs:
        tg = x[t + m - 1 + h]
        out[h] = float(np.median(np.abs(tg - tg[nn].mean(1))) / mad)
    return out


def masked_scan(x, mask, salt):
    """z_m (h = 1) vs masked IAAFT for m = 2..5, the data's growth G (m = 3), and the surrogates'."""
    import final_pipeline as fp
    x = np.asarray(x, float)
    xi = interp_masked(x, mask)
    rng = rng_for(x, salt)
    obs = {m: masked_nlp(x, mask, m) for m in MS}
    if obs[2] is None:
        return {"analysable": False, "n_masked": int(mask.sum())}
    sur = {m: [] for m in MS}
    for _ in range(N_SUR):
        s = fp.iaaft_surrogate(xi, rng)
        for m in MS:
            e = masked_nlp(s, mask, m, hs=(1,))
            sur[m].append(e[1] if e is not None else np.nan)
    out = {"analysable": all(obs[m] is not None for m in MS), "n_masked": int(mask.sum()),
           "n_valid_m3": int(len(valid_vectors(mask, 3)))}
    for m in MS:
        if obs[m] is None:
            out[f"z_m{m}"] = float("nan")
            continue
        sv = np.array(sur[m], float)
        out[f"z_m{m}"] = float((np.nanmean(sv) - obs[m][1]) / max(np.nanstd(sv, ddof=1), 1e-12))
        out[f"E_m{m}"] = [obs[m][h] for h in range(1, HMAX + 1)]
    out["zmax"] = float(np.nanmax([out[f"z_m{m}"] for m in MS]))
    out["G"] = float(obs[3][5] - obs[3][2]) if obs[3] is not None else float("nan")
    return out


def fsle_gate_value(x, mask):
    from experiments.phase9_waveform import measures as M
    dl, lam, cnt, _ = M.fsle(interp_masked(x, mask), m=3, L=1, delta0=0.01, n_levels=14, maxt=50)
    return float(lam[10]) if cnt[10] >= 5 else float("nan")


# ------------------------------------------------------------------------------- candidates
def k1_frozen512(x, types=None):
    from experiments.phase8_cardiac import methods as MM
    return MM._frozen(x)


def _ann(x, types, cache):
    if "ann" not in cache:
        cache["ann"] = masked_scan(x, mask_ann(types), 1)
    return cache["ann"]


def _rrm(x, cache):
    if "rr" not in cache:
        cache["rr"] = masked_scan(x, mask_rr(x), 2)
    return cache["rr"]


def k2_mnlp_ann(x, types, cache):
    r = dict(_ann(x, types, cache))
    r["detected"] = bool(r.get("analysable") and r["zmax"] >= Z_DET)
    return r


def k3_mnlp_growth_ann(x, types, cache):
    r = dict(_ann(x, types, cache))
    r["detected"] = bool(r.get("analysable") and r["zmax"] >= Z_DET and r["G"] >= G_MIN)
    return r


def k4_mnlp_growth_rr(x, types, cache):
    r = dict(_rrm(x, cache))
    r["detected"] = bool(r.get("analysable") and r["zmax"] >= Z_DET_RR and r["G"] >= G_MIN_RR)
    return r


def k5_mnlp_fsle_ann(x, types, cache):
    r = dict(_ann(x, types, cache))
    if "fsle" not in cache:
        cache["fsle"] = fsle_gate_value(x, mask_ann(types))
    r["fsle10"] = cache["fsle"]
    r["detected"] = bool(r.get("analysable") and r["zmax"] >= Z_DET and np.isfinite(r["fsle10"])
                         and r["fsle10"] >= F_MIN)
    return r


CANDIDATES = ("k1_frozen512", "k2_mnlp_ann", "k3_mnlp_growth_ann", "k4_mnlp_growth_rr", "k5_mnlp_fsle_ann")
FUNCS = {"k1_frozen512": lambda x, t, c: k1_frozen512(x), "k2_mnlp_ann": k2_mnlp_ann,
         "k3_mnlp_growth_ann": k3_mnlp_growth_ann, "k4_mnlp_growth_rr": k4_mnlp_growth_rr,
         "k5_mnlp_fsle_ann": k5_mnlp_fsle_ann}


def run_all(x, types, names=CANDIDATES):
    cache, out = {}, {}
    for nm in names:
        try:
            out[nm] = FUNCS[nm](np.asarray(x, float), types, cache)
        except Exception as exc:                        # noqa: BLE001
            out[nm] = {"detected": False, "error": f"{type(exc).__name__}: {exc}",
                       "traceback": traceback.format_exc(limit=2)}
    return out


def raw_stats(x, types):
    """Calibration statistics (no thresholds): both masks' scans and the FSLE gate value."""
    cache = {}
    a = _ann(x, types, cache)
    r = _rrm(x, cache)
    f = fsle_gate_value(x, mask_ann(types)) if a.get("analysable") else float("nan")
    return {"ann": a, "rr": r, "fsle10": f}
