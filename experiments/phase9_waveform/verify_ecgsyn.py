"""
Phase 9 Part A1 verification of the externally driven ECGSYN (ecgsyn.py).

    python -m experiments.phase9_waveform.verify_ecgsyn [--cref DIR]

Checks (results/verification/ecgsyn.json, plots/A1_ecgsyn_verification.png):
 V1  Against the authors' C implementation (PhysioNet ecgsyn.c, compiled with stand-in FFT
     and uniform RNG for the two missing Numerical Recipes files, which only feed the RR
     process): constant 60 and 75 bpm, and variable RR (60 +/- 5 and 90 +/- 3 bpm).  The C beat
     times are reconstructed from its piecewise-constant RR (rrpc.dat) by integrating the
     angular frequency as its RK4 step does; ecgsyn.py is driven by those beat times.  Both
     signals are min/max scaled to [-0.4, 1.2] mV on the common span (as ecgsyn.c does).
 V2  Paper Table I: times of the P, Q, R, S, T extrema of one beat at 60 bpm.
 V3  Paper Fig. 9 / Fig. 10 (60 bpm, SD 5 bpm): RS amplitude increases with RR; QT (Q event
     to T event) is linear in RR.
 V4  Beat times: R-peak times (parabolic peak of z near each input beat) reproduce the input
     RR series (constant, variable, and a Phase 8 chaotic RR window with ectopy).
"""
from __future__ import annotations

import argparse
import json
import math
import pathlib

import numpy as np

from experiments.phase9_waveform import ecgsyn as E

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "results" / "verification"
PLOTS = HERE / "plots"
FS = 360
Q = 2


def c_beat_times(rrpc, h):
    """Beat times of ecgsyn.c's limit-cycle phase: theta_{n+1} = theta_n + h (5 w_n + w_{n+1}) / 6
    (RK4 evaluates angfreq at t, t + h/2, t + h/2, t + h; floor((t+h/2)/h) = n), theta_0 = 0."""
    w = 2.0 * np.pi / rrpc
    dth = h * (5.0 * w[:-1] + w[1:]) / 6.0
    theta = np.concatenate([[0.0], np.cumsum(dth)])
    k = np.floor(theta / (2.0 * np.pi))
    idx = np.nonzero(np.diff(k) > 0)[0]
    t = []
    for i in idx:
        target = (k[i] + 1) * 2.0 * np.pi
        frac = (target - theta[i]) / (theta[i + 1] - theta[i])
        t.append((i + frac) * h)
    return np.concatenate([[0.0], np.array(t)])


def mm(v):
    return (v - v.min()) * 1.6 / (v.max() - v.min()) - 0.4


def r_peak_times(z, beat_t, fs=FS, t0=0.0, win=0.05):
    out = []
    for tb in beat_t:
        i0 = int(round((tb - t0 - win) * fs))
        i1 = int(round((tb - t0 + win) * fs)) + 1
        if i0 < 1 or i1 >= len(z) - 1:
            out.append(np.nan)
            continue
        i = i0 + int(np.argmax(z[i0:i1]))
        y0, y1, y2 = z[i - 1], z[i], z[i + 1]
        den = y0 - 2 * y1 + y2
        d = 0.5 * (y0 - y2) / den if den != 0 else 0.0
        out.append(t0 + (i + d) / fs)
    return np.array(out)


def compare_c(cdir, hr, label):
    d = np.loadtxt(cdir / "ecg.dat")
    zc = d[:, 1]
    rrpc = np.loadtxt(cdir / "rrpc.dat")
    h = 1.0 / (Q * FS)
    bt = c_beat_times(rrpc, h)
    th, a, b = E.hr_adjusted(hr)
    T = min(len(zc) / FS, bt[-1]) - 1.0
    z = E.synthesize(bt, th, a, b, fs=FS, q=Q, t0=0.0, t1=T, A=E.A_RESP, resp_phase=0.0, zinit=E.Z_INIT)
    n = len(z)
    zc_s, z_s = mm(zc[:n]), mm(z)
    i0 = int(5 * FS)
    return {"label": label, "hr_mean": hr, "n_samples": n, "n_beats": int(np.sum(bt < T)),
            "max_abs_diff_mV": float(np.abs(zc_s[i0:] - z_s[i0:]).max()),
            "rms_diff_mV": float(np.sqrt(np.mean((zc_s[i0:] - z_s[i0:]) ** 2))),
            "corr": float(np.corrcoef(zc_s[i0:], z_s[i0:])[0, 1])}, (zc_s, z_s, bt)


def extrema_times():
    th, a, b = E.hr_adjusted(60.0)
    bt = np.arange(0.0, 31.0)
    z = E.to_mV(E.synthesize(bt, th, a, b, fs=FS * 10, q=1, A=0.0))
    fs = FS * 10
    tR = 20.0
    seg = lambda lo, hi: slice(int((tR + lo) * fs), int((tR + hi) * fs))
    res = {}
    for name, (lo, hi, fn) in {"P": (-0.35, -0.08, np.argmax), "Q": (-0.08, 0.0, np.argmin),
                               "R": (-0.03, 0.03, np.argmax), "S": (0.0, 0.1, np.argmin),
                               "T": (0.1, 0.6, np.argmax)}.items():
        s = seg(lo, hi)
        i = s.start + int(fn(z[s]))
        res[name] = {"time_s": round(i / fs - tR, 4), "amp_mV": round(float(z[i]), 4)}
    return res, (np.arange(len(z)) / fs - tR, z)


def paper_fig9_fig10(cdir):
    """60 +/- 5 bpm RR from the C RR process; per beat RS amplitude and Q-to-T time."""
    rrpc = np.loadtxt(cdir / "rrpc.dat")
    bt = c_beat_times(rrpc, 1.0 / (Q * FS))
    th, a, b = E.hr_adjusted(60.0)
    fs = FS * 4
    z = E.to_mV(E.synthesize(bt, th, a, b, fs=fs, q=1, A=0.0))
    rows = []
    for k in range(3, len(bt) - 3):
        tk = bt[k]
        if (tk + 0.7) * fs >= len(z):
            break
        sl = lambda lo, hi: slice(int((tk + lo) * fs), int((tk + hi) * fs))
        sR, sS, sQ = sl(-0.03, 0.03), sl(0.0, 0.12), sl(-0.12, 0.0)
        iR = sR.start + int(np.argmax(z[sR]))
        iS = sS.start + int(np.argmin(z[sS]))
        iQ = sQ.start + int(np.argmin(z[sQ]))
        rr_next = bt[k + 1] - bt[k]
        sT = sl(0.1, 0.8 * rr_next)
        iT = sT.start + int(np.argmax(z[sT]))
        rows.append((rr_next, bt[k] - bt[k - 1], z[iR] - z[iS], (iT - iQ) / fs))
    r = np.array(rows)
    rs_corr_next = float(np.corrcoef(r[:, 0], r[:, 2])[0, 1])
    rs_corr_prev = float(np.corrcoef(r[:, 1], r[:, 2])[0, 1])
    slope, icpt = np.polyfit(r[:, 0], r[:, 3], 1)
    pred = slope * r[:, 0] + icpt
    r2 = 1 - np.sum((r[:, 3] - pred) ** 2) / np.sum((r[:, 3] - r[:, 3].mean()) ** 2)
    return {"n_beats": len(r), "rr_mean": float(r[:, 0].mean()), "rr_sd": float(r[:, 0].std()),
            "corr_RS_amp_vs_RR_in_which_beat_sits": rs_corr_next, "corr_RS_amp_vs_preceding_RR": rs_corr_prev,
            "QT_vs_RR_slope": float(slope), "QT_vs_RR_intercept_s": float(icpt), "QT_vs_RR_R2": float(r2)}, r


def beat_time_check(bt, label, params=None):
    if params is None:
        th, a, b = E.hr_adjusted(60.0 / np.mean(np.diff(bt)))
    else:
        th, a, b = params
    z = E.synthesize(bt, th, a, b, fs=FS, q=Q, A=E.A_RESP)
    tp = r_peak_times(z, bt[2:-2], t0=bt[0])
    off = tp - bt[2:-2]
    rr_in = np.diff(bt[2:-2])
    rr_out = np.diff(tp)
    return {"label": label, "n_beats": len(tp), "offset_ms_median": float(np.nanmedian(off) * 1e3),
            "offset_ms_range": [float(np.nanmin(off) * 1e3), float(np.nanmax(off) * 1e3)],
            "rr_error_ms_max_abs": float(np.nanmax(np.abs(rr_out - rr_in)) * 1e3),
            "rr_error_ms_sd": float(np.nanstd(rr_out - rr_in) * 1e3),
            "corr_rr_in_out": float(np.corrcoef(rr_in, rr_out)[0, 1]), "sample_ms": 1e3 / FS}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--cref", default="/tmp/claude-0/ecgsyn_build")
    a = ap.parse_args(argv)
    cref = pathlib.Path(a.cref)
    OUT.mkdir(parents=True, exist_ok=True)
    PLOTS.mkdir(parents=True, exist_ok=True)
    res = {"V1_vs_C": [], "notes": {}}
    traces = {}
    # constant HR: C run with -n 30 -h 60 -H 0 (ref60.dat) and -h 75 -H 0 (r75)
    for sub, hr, lab in (("r60c", 60.0, "60 bpm constant"), ("r75", 75.0, "75 bpm constant"),
                         ("r60v", 60.0, "60 +/- 5 bpm"), ("r90v", 90.0, "90 +/- 3 bpm")):
        if (cref / sub / "ecg.dat").exists():
            r, tr = compare_c(cref / sub, hr, lab)
            res["V1_vs_C"].append(r)
            traces[lab] = tr
    res["V2_table1_60bpm"], ex = extrema_times()
    res["V2_table1_published_times_s"] = {"P": -0.2, "Q": -0.05, "R": 0.0, "S": 0.05, "T": 0.3}
    if (cref / "r60v" / "rrpc.dat").exists():
        res["V3_fig9_fig10"], f910 = paper_fig9_fig10(cref / "r60v")
    rng = np.random.default_rng(9001)
    rr_const = np.full(200, 0.8)
    rr_var = 0.8 + 0.05 * rng.standard_normal(200)
    from experiments.phase8_cardiac import series as P8
    rr_pr = P8.window("phase_reset", 14, {"tau": 1.2}, 300, "vi_E3_10", 9001)
    res["V4_beat_times"] = [beat_time_check(np.concatenate([[0], np.cumsum(r)]), lab)
                            for r, lab in ((rr_const, "constant 0.8 s"), (rr_var, "white 0.8 +/- 0.05 s"),
                                           (rr_pr, "phase_reset tau=1.2 (vi_E3_10), dev seed 9001"))]
    json.dump(res, open(OUT / "ecgsyn.json", "w"), indent=1)
    print(json.dumps(res, indent=1))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(11, 7))
    if "60 +/- 5 bpm" in traces:
        zc, zm, bt = traces["60 +/- 5 bpm"]
        t = np.arange(len(zc)) / FS
        s = slice(int(20 * FS), int(26 * FS))
        ax[0, 0].plot(t[s], zc[s], lw=2.5, color="0.7", label="ecgsyn.c (authors)")
        ax[0, 0].plot(t[s], zm[s], lw=1, color="C3", label="ecgsyn.py (beat-time driven)")
        ax[0, 0].set(title="V1: 60 +/- 5 bpm, same beat times", xlabel="s", ylabel="mV")
        ax[0, 0].legend(fontsize=8)
    tt, zz = ex
    s = (tt > -0.4) & (tt < 0.7)
    ax[0, 1].plot(tt[s], zz[s], color="k")
    for nm, v in res["V2_table1_published_times_s"].items():
        ax[0, 1].axvline(v, color="C0", ls=":", lw=1)
        ax[0, 1].text(v, 1.25, nm, ha="center", color="C0")
    ax[0, 1].set(title="V2: one beat at 60 bpm (dotted: Table I times)", xlabel="s from R", ylabel="mV")
    if "V3_fig9_fig10" in res:
        ax[1, 0].plot(f910[:, 0], f910[:, 2], ".", ms=3)
        ax[1, 0].set(title=f"V3 (Fig. 9): RS amplitude vs RR, r = {res['V3_fig9_fig10']['corr_RS_amp_vs_RR_in_which_beat_sits']:.2f}",
                     xlabel="RR (s)", ylabel="RS amplitude (mV)")
        ax[1, 1].plot(f910[:, 0], f910[:, 3], ".", ms=3)
        ax[1, 1].set(title=f"V3 (Fig. 10): QT vs RR, R2 = {res['V3_fig9_fig10']['QT_vs_RR_R2']:.3f}",
                     xlabel="RR (s)", ylabel="Q-to-T time (s)")
    fig.tight_layout()
    fig.savefig(PLOTS / "A1_ecgsyn_verification.png", dpi=110)


if __name__ == "__main__":
    main()
