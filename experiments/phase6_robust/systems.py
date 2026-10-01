"""
Phase 6 signals: Phase 5 conditions (experiments/phase5_rr/systems.py, unchanged)
plus the Part B ectopy conditions, all on the N1 linear_rr base.

Ventricular-type ectopic group at interval index i with k premature beats
(k = 1: isolated / bigeminy / trigeminy; 2: couplet; 3-6: run):
  premature intervals  i .. i+k-1: the first is c * x_i (c ~ U(0.6, 0.8)); later
                       ones are c' * x_i (couplet, c' ~ U(0.6, 0.8)) or an
                       absolute run cycle U(0.40, 0.55) s (runs);
  compensatory interval i+k: the sinus node is NOT reset (full compensatory
                       pause), so the next conducted sinus beat keeps its
                       scheduled time: sum(x_i .. x_i+k) - sum(premature).
  The interval count is preserved (k premature beats replace k scheduled
  sinus beats); labels mark intervals i .. i+k as non-NN.
Atrial-type (E5): premature interval c * x_i, then the sinus node is reset: the
  post-ectopic interval is d * x_{i+1}, d ~ U(1.0, 1.1) (less than
  compensatory); later beats shift earlier; count preserved.
Annotation-edited variants: every non-NN interval (premature, couplet/run and
  compensatory) is replaced by linear interpolation (by index) between the
  nearest NN intervals before and after its block; a block at an edge takes
  the nearest NN value.  Editing is applied to the quantized raw series;
  interpolated values are not re-quantized.  In bigeminy only the edge
  intervals are NN, so the edited window is (by construction of the rule) a
  near-straight line between them; this is reported, not patched.
"""
from __future__ import annotations

import numpy as np

from experiments.phase5_rr import config as C5
from experiments.phase5_rr import systems as S5
from experiments.phase6_robust import config as C


def _rng(condition, seed):
    return np.random.default_rng(np.random.SeedSequence([C.ENTROPY, C.PARTB_CODE[condition], C.N, int(seed)]))


def _apply_group(x, labels, i, k, rng, run_cycle=None):
    """Ventricular group of k premature beats starting at interval i (full compensatory)."""
    base = x[i:i + k + 1].copy()
    prem = [rng.uniform(*C.PREMATURE_COUPLING_RANGE) * base[0]]
    for _ in range(k - 1):
        prem.append(run_cycle[0] if run_cycle is not None else rng.uniform(*C.PREMATURE_COUPLING_RANGE) * base[0])
        if run_cycle is not None:
            run_cycle = run_cycle[1:]
    x[i:i + k] = prem
    x[i + k] = base.sum() - np.sum(prem)
    labels[i:i + k + 1] = True


def _place_groups(rng, sizes, n):
    """Start indices for groups occupying sizes[j] intervals, placed one at a time (as
    Phase 5 S2): a uniformly drawn start in 1 .. n-1-size is accepted if the group lies
    inside intervals 1 .. n-2 and leaves >= MIN_NORMAL_GAP normal intervals to every
    group already placed."""
    spans = []
    for size in sizes:
        for _ in range(100000):
            s0 = int(rng.integers(1, n - size))
            s1 = s0 + size - 1
            if s1 <= n - 2 and all(s0 - b1 - 1 >= C.MIN_NORMAL_GAP or a0 - s1 - 1 >= C.MIN_NORMAL_GAP
                                   for a0, b1 in spans):
                spans.append((s0, s1))
                break
        else:
            raise RuntimeError("could not place ectopic groups")
    return [s0 for s0, _ in spans]


def ectopy(kind, params, seed, condition):
    """(unquantized RR, labels, info) for a Part B raw condition."""
    p5 = S5.parts("N1_linear_rr", seed)
    x = p5["unquantized"].copy()
    n = len(x)
    labels = np.zeros(n, dtype=bool)
    rng = _rng(condition, seed)
    info = {}
    if kind in ("bigeminy", "trigeminy"):
        period = 2 if kind == "bigeminy" else 3
        s = int(rng.integers(1, period + 1))
        starts = list(range(s, n - 1, period))
        for i in starts:
            _apply_group(x, labels, i, 1, rng)
        info["n_ectopic"] = len(starts)
    elif kind == "couplets":
        m = int(round(params["rate"] * n / 2))
        starts = _place_groups(rng, [3] * m, n)
        for i in starts:
            _apply_group(x, labels, i, 2, rng)
        info["n_ectopic"] = 2 * m
    elif kind == "runs":
        n_runs = int(rng.integers(C.RUN_COUNT_RANGE[0], C.RUN_COUNT_RANGE[1] + 1))
        lengths = [int(rng.integers(C.RUN_LENGTH_RANGE[0], C.RUN_LENGTH_RANGE[1] + 1)) for _ in range(n_runs)]
        starts = _place_groups(rng, [k + 1 for k in lengths], n)
        for i, k in zip(starts, lengths):
            cyc = list(rng.uniform(*C.RUN_CYCLE_RANGE_S, size=k - 1))
            _apply_group(x, labels, i, k, rng, run_cycle=cyc)
        info.update(n_ectopic=int(sum(lengths)), run_lengths=lengths)
    elif kind == "atrial":
        m = int(round(params["rate"] * n))
        starts = _place_groups(rng, [2] * m, n)
        for i in starts:
            c = rng.uniform(*C.PREMATURE_COUPLING_RANGE)
            d = rng.uniform(*C.ATRIAL_POST_FACTOR_RANGE)
            x[i] = c * x[i]
            x[i + 1] = d * x[i + 1]
            labels[i:i + 2] = True
        info["n_ectopic"] = m
    elif kind == "s2_trend":
        # exactly N3's trend draw followed by S2 10 %'s ectopic draw (Phase 5 modifier streams)
        mr_t = np.random.default_rng(S5._ss(C5.COND_CODE["N3_linear_rr_trend"], seed))
        x, pt = S5.add_trend(x, mr_t, p5["mean_rr"])
        mr_e = np.random.default_rng(S5._ss(C5.COND_CODE["S2_ectopic_10pct"], seed))
        x, pe = S5.add_ectopics(x, mr_e, 0.10)
        for i in pe["ectopic_indices"]:
            labels[i:i + 2] = True
        info.update(pt, n_ectopic=len(pe["ectopic_indices"]))
    else:
        raise ValueError(kind)
    if not np.all(np.isfinite(x)) or np.min(x) <= 0:
        raise FloatingPointError(f"{condition} seed {seed}: invalid RR series")
    return x, labels, info


def _s2_labels(condition, seed):
    parts = S5.parts(condition, seed)
    labels = np.zeros(C.N, dtype=bool)
    for i in parts["params"]["ectopic_indices"]:
        labels[i:i + 2] = True
    return parts["rr"], labels


def edit_nn(rr, labels):
    """Annotation-edited NN series (module docstring).  Returns (series, info) with the
    number and fraction of genuine NN intervals kept."""
    rr = np.asarray(rr, dtype=float).copy()
    good = ~np.asarray(labels, dtype=bool)
    if not good.any():
        raise ValueError("no NN interval in the window")
    idx = np.arange(len(rr))
    rr[~good] = np.interp(idx[~good], idx[good], rr[good])   # linear inside, nearest value at the edges
    return rr, {"n_nn": int(good.sum()), "nn_fraction": float(good.mean())}


def parts(condition, seed):
    if condition in C5.CONDITIONS:
        p = S5.parts(condition, seed)
        return {"rr": p["rr"], "labels": None, "info": p["params"]}
    grp, kind, params = C.PARTB[condition]
    if kind == "edited":
        raw = params["raw"]
        if raw in C5.CONDITIONS:
            rr_raw, labels = _s2_labels(raw, seed)
        else:
            x, labels, _ = ectopy(C.PARTB[raw][1], C.PARTB[raw][2], seed, raw)
            rr_raw = S5.quantize(x)
        rr, info = edit_nn(rr_raw, labels)
        return {"rr": rr, "labels": labels, "info": info, "raw": rr_raw}
    x, labels, info = ectopy(kind, params, seed, condition)
    return {"rr": S5.quantize(x), "labels": labels, "info": info}


def generate(condition, seed):
    return parts(condition, seed)["rr"]


def data_seed(condition, seed):
    if condition in C5.CONDITIONS:
        return S5.data_seed(condition, seed)
    return S5.data_seed("N1_linear_rr", seed)
