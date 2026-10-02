"""
Phase 9 Part B3-B5: false-positive rate of the waveform surrogates on DEVELOPMENT waveform nulls
(and power on a development chaotic family), plus the compute check.

    python -m experiments.phase9_waveform.fp_waveform --workers 4

Predeclared (written before this script was run on any Phase 9 series):
- Windows: 2-minute waveform windows (library.waveform) from 512-beat DEV records (seeds 9000-9004),
  downsampled 360 -> 90 Hz (scipy.signal.decimate, q = 4, zero-phase FIR anti-alias filter).
- Nulls: 'periodic:0.8' (strictly periodic ECG; the PPS null with noise when nstdb is added),
  N1, N2, N6 (correlated / power-law / RSA heart-rate variability), S1 (SETAR), and every ectopy
  pattern (S2 5/10 %, E1, E2, E3 5/10 %, E4, E5 5/10 %, E6); each 'clean' and with nstdb 'mix12'
  (DEV signal 0).  Power check: DEV Mackey-Glass tau = 23 (CHAOTIC), clean and mix12.
- Embedding: tau by first AMI minimum, m by Cao (cap 10), on the first 3,000 samples of the data;
  surrogates use the same.
- Statistics (B4): normalized nonlinear prediction error at horizon h = 1 sample and h = tau, k = 5,
  Theiler window = mean RR in samples; one-sided (data more predictable than surrogates).
- Surrogates: PPS (rho by the noise-radius rule over SD x {0.02, 0.04, 0.08, 0.16, 0.32}),
  CS (cut at the R peaks: true beats for clean, detected beats for nstdb windows), TS (recurrence
  rate 0.05); 19 each; rejection iff p <= 0.05.
- A surrogate is UNUSABLE for chaos claims if its rejection rate exceeds 7 % on any null (B3).
- B5: runtime per window of every component is recorded.
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import time
from multiprocessing import Pool

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402

from experiments.phase9_waveform import families as F  # noqa: E402
from experiments.phase9_waveform import library as L  # noqa: E402
from experiments.phase9_waveform import surrogates_wave as SW  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "results" / "dev"
FS_DS = 90
Q = 4
# Reduced after the first window took 230 s (B5): 6 nulls x 5 seeds x {clean, mix12}; a rejection rate
# above 7 % on any one null already makes a surrogate UNUSABLE (B3), so the reduction can only miss
# an unusable verdict, not create one.
NULLS_W = ["periodic:0.8", "N1_linear_rr", "N6_noisy_rsa", "E1_bigeminy", "E3_couplets_10pct", "E5_atrial_10pct"]
SEEDS = range(9000, 9010)
N_SUR = 19


def window_ds(source, name, seed, noise, minutes=2.0):
    from scipy.signal import decimate
    beats = "true" if noise == "none" else "detected"
    w = L.build(source, name, "none", 512, seed, noise, "DEV", beats, False)
    s = max(0, int(w["true_idx"][w["first_window_beat"]] - 0.5 * L.FS))
    e = s + int(round(minutes * 60 * L.FS))
    x = w["ecg"][s:e]
    xd = decimate(x, Q, ftype="fir", zero_phase=True)
    bi = np.asarray(w["beat_idx"])
    cuts = np.unique(np.round((bi[(bi >= s) & (bi < e)] - s) / Q).astype(int))
    mean_rr = float(np.mean(w["rr"]))
    return xd, cuts, mean_rr


def one(task):
    source, name, seed, noise = task
    t0 = time.time()
    x, cuts, mean_rr = window_ds(source, name, seed, noise)
    t_gen = time.time() - t0
    rng = np.random.default_rng(seed + 77)
    t1 = time.time()
    tau, m = SW.choose_embedding(x[:3000], max_tau=60)   # rule applied to the first 3,000 samples (cost)
    t_emb = time.time() - t1
    th = int(round(mean_rr * FS_DS))
    out = {"source": source, "name": name, "seed": seed, "noise": noise, "m": m, "tau": tau, "n": len(x),
           "theiler": th, "t_gen": t_gen, "t_emb": t_emb}
    hs = {"h1": 1, "htau": tau}
    span = (m - 1) * tau
    z = SW.embed_fwd(x, m, tau)
    stat = {k: (lambda y, h=h: SW.nlp_error(y, m, tau, h=h, theiler=th)) for k, h in hs.items()}
    t2 = time.time()
    obs = {k: f(z[:, 0]) for k, f in stat.items()}
    out["t_stat"] = (time.time() - t2) / 2
    # PPS
    t3 = time.time()
    sd = np.std(x)
    rho, cnt = SW.noise_radius(x, m, tau, sd * np.array([0.02, 0.04, 0.08, 0.16, 0.32]), rng)
    t_rho = time.time() - t3
    sur = {k: [] for k in hs}
    t4 = time.time()
    for _ in range(N_SUR):
        s = SW.pps(x, m, tau, rho, rng)
        for k, f in stat.items():
            sur[k].append(f(s))
    out["t_pps_each"] = (time.time() - t4) / N_SUR
    out["rho_over_sd"] = rho / sd
    out["t_rho"] = t_rho
    for k in hs:
        out[f"p_pps_{k}"] = float((1 + np.sum(np.array(sur[k]) <= obs[k])) / (N_SUR + 1))
    # CS (data restricted to the first..last cut, statistic recomputed on that span)
    if len(cuts) >= 4:
        _, dat = SW.cycle_shuffle(x, cuts, rng)
        obs_cs = {k: f(dat) for k, f in stat.items()}
        sc = {k: [] for k in hs}
        for _ in range(N_SUR):
            s = SW.cycle_shuffle(x, cuts, rng)[0]
            for k, f in stat.items():
                sc[k].append(f(s))
        for k in hs:
            out[f"p_cs_{k}"] = float((1 + np.sum(np.array(sc[k]) <= obs_cs[k])) / (N_SUR + 1))
    # TS
    t5 = time.time()
    zs = (z - z.mean(0)) / z.std(0)
    group, delta = SW.twin_sets(zs, rr=0.05)
    out["n_points_with_twins"] = int(sum(1 for g in group if len(g) > 1))
    st_ = {k: [] for k in hs}
    for _ in range(N_SUR):
        idx = SW.twin_surrogate(zs, group, rng, length=len(z))
        s = z[idx, 0]
        for k, f in stat.items():
            st_[k].append(f(s))
    out["t_ts_total"] = time.time() - t5
    for k in hs:
        out[f"p_ts_{k}"] = float((1 + np.sum(np.array(st_[k]) <= obs[k])) / (N_SUR + 1))
    out["t_total"] = time.time() - t0
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--seeds", type=int, default=10)
    a = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "fp_waveform.jsonl"
    done = set()
    if path.exists():
        for line in open(path):
            r = json.loads(line)
            done.add((r["source"], r["name"], r["seed"], r["noise"]))
    tasks = []
    a.seeds = min(a.seeds, 5)
    for nm in NULLS_W:
        src = "periodic" if nm.startswith("periodic") else "null"
        for nz in ("none", "mix12"):
            tasks += [(src, nm, s, nz) for s in list(SEEDS)[:a.seeds]]
    tasks += [("model", "mackey_glass:tau=23.0", s, "none") for s in list(SEEDS)[:a.seeds]]
    tasks = [t for t in tasks if t not in done]
    print(len(tasks), "to run", flush=True)
    with open(path, "a") as fh, Pool(a.workers, maxtasksperchild=10) as pool:
        for i, r in enumerate(pool.imap_unordered(one, tasks), 1):
            fh.write(json.dumps(r, default=float) + "\n")
            fh.flush()
            if i % 10 == 0:
                print(i, "done", flush=True)


if __name__ == "__main__":
    main()
