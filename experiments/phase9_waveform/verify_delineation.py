"""
Phase 9 Part C1: verification of the delineator (NeuroKit2 0.2.13 ecg_delineate(method='dwt'),
the wavelet delineator of Martinez et al., IEEE TBME 51:570 (2004)) against the published
QT-database (QTDB) accuracy before it is used on any Phase 9 series.

    python -m experiments.phase9_waveform.verify_delineation

Data: PhysioNet QTDB 1.0.0, all 105 records with the first annotator's manual annotations (.q1c),
250 Hz, 2 leads (qtdb_data/, SHA-256 checked against SHA256SUMS.txt; not committed).
Published behaviour (Martinez et al. 2004, QTDB, as restated in Table 6 of Di Marco & Chiari,
BioMed Eng OnLine 10:23 (2011), PMC3076264; mean +/- SD of automatic - manual, ms):
  QRS onset 4.6 +/- 7.7, QRS offset 0.8 +/- 8.7, T peak 0.2 +/- 13.9, T offset (end) -1.6 +/- 18.1.
Procedure (predeclared): the delineator is given the manual QRS peak positions ('N' marks) of each
record (so that only delineation is evaluated); for every manual fiducial (QRS onset '(' before
N, QRS offset ')' after N, T peak 't', T end ')' after t) the delineator output of the same beat is
compared; an output within 150 ms counts as detected (sensitivity); errors are computed over
detected fiducials.  Reported for lead 0 and for the better of the two leads per fiducial (the
QTDB convention of choosing, per annotation, the lead with the smaller error).
Pass criterion (predeclared): best-of-two-leads |mean| <= published SD and SD <= 2 x published SD
for QRS onset, QRS offset and T end, and sensitivity >= 95 % for each.
"""
from __future__ import annotations

import json
import pathlib
import warnings

import numpy as np

warnings.filterwarnings("ignore")
HERE = pathlib.Path(__file__).resolve().parent
QT = HERE / "qtdb_data"
OUT = HERE / "results" / "verification"
PUBLISHED = {"qrs_on": (4.6, 7.7), "qrs_off": (0.8, 8.7), "t_peak": (0.2, 13.9), "t_end": (-1.6, 18.1)}
NK_KEY = {"qrs_on": "ECG_R_Onsets", "qrs_off": "ECG_R_Offsets", "t_peak": "ECG_T_Peaks", "t_end": "ECG_T_Offsets"}
TOL_MS = 150.0


def manual_fiducials(ann):
    """List of beats: dict with 'r' (N sample) and the manual fiducials present."""
    sym, smp = list(ann.symbol), np.asarray(ann.sample)
    beats = []
    for i, s in enumerate(sym):
        if s != "N":
            continue
        b = {"r": int(smp[i])}
        if i >= 1 and sym[i - 1] == "(":
            b["qrs_on"] = int(smp[i - 1])
        if i + 1 < len(sym) and sym[i + 1] == ")":
            b["qrs_off"] = int(smp[i + 1])
        # T wave of this beat: the next 't' before the next 'N'
        for k in range(i + 1, min(i + 8, len(sym))):
            if sym[k] == "N":
                break
            if sym[k] == "t":
                b["t_peak"] = int(smp[k])
                if k + 1 < len(sym) and sym[k + 1] == ")":
                    b["t_end"] = int(smp[k + 1])
                break
        beats.append(b)
    return beats


def delineate(sig, rpk, fs):
    import neurokit2 as nk
    _, w = nk.ecg_delineate(sig, rpk, sampling_rate=fs, method="dwt")
    return {k: np.asarray(w[NK_KEY[k]], dtype=float) for k in NK_KEY}


def main():
    import wfdb
    recs = [r.strip() for r in open(QT / "RECORDS") if r.strip()]
    err = {lead: {k: [] for k in NK_KEY} for lead in ("lead0", "best")}
    n_man = {k: 0 for k in NK_KEY}
    n_det = {lead: {k: 0 for k in NK_KEY} for lead in ("lead0", "best")}
    per_rec = []
    for rec in recs:
        try:
            ann = wfdb.rdann(str(QT / rec), "q1c")
        except Exception:
            continue
        r = wfdb.rdrecord(str(QT / rec))
        fs = r.fs
        beats = manual_fiducials(ann)
        if not beats:
            continue
        # context: delineate a segment around the annotated beats (all beats of the record are given
        # as R peaks: manual 'N' of q1c only cover the annotated part; use the reference 'atr'-free
        # approach: NeuroKit's own peaks over the whole record plus the manual ones)
        import neurokit2 as nk
        outs = []
        for ch in range(r.p_signal.shape[1]):
            sig = np.nan_to_num(r.p_signal[:, ch])
            _, info = nk.ecg_peaks(sig, sampling_rate=fs)
            auto = np.asarray(info["ECG_R_Peaks"], dtype=int)
            man = np.array([b["r"] for b in beats])
            # replace automatic peaks within 150 ms of a manual N by the manual N (delineation only)
            keep = [a for a in auto if np.min(np.abs(man - a)) > TOL_MS * fs / 1000]
            rp = np.sort(np.concatenate([np.array(keep, dtype=int), man]))
            d = delineate(sig, rp, fs)
            idx = {int(v): i for i, v in enumerate(rp)}
            outs.append((d, idx))
        rec_err = {k: [] for k in NK_KEY}
        for b in beats:
            for k in NK_KEY:
                if k not in b:
                    continue
                n_man[k] += 1
                cand = []
                for d, idx in outs:
                    v = d[k][idx[b["r"]]]
                    cand.append((v - b[k]) * 1000.0 / fs if np.isfinite(v) else np.nan)
                c0 = cand[0]
                if np.isfinite(c0) and abs(c0) <= TOL_MS:
                    err["lead0"][k].append(c0)
                    n_det["lead0"][k] += 1
                fin = [c for c in cand if np.isfinite(c) and abs(c) <= TOL_MS]
                if fin:
                    cb = min(fin, key=abs)
                    err["best"][k].append(cb)
                    n_det["best"][k] += 1
                    rec_err[k].append(cb)
        per_rec.append({"record": rec, **{k: [float(np.mean(v)), float(np.std(v)), len(v)] if v else None
                                          for k, v in rec_err.items()}})
    res = {"n_records": len(per_rec), "published_martinez2004": PUBLISHED, "results": {}}
    ok = True
    for lead in ("lead0", "best"):
        res["results"][lead] = {}
        for k in NK_KEY:
            e = np.array(err[lead][k])
            se = n_det[lead][k] / max(n_man[k], 1)
            res["results"][lead][k] = {"n_manual": n_man[k], "sensitivity": se, "mean_ms": float(e.mean()),
                                       "sd_ms": float(e.std(ddof=1))}
            if lead == "best" and k in ("qrs_on", "qrs_off", "t_end"):
                m_pub, s_pub = PUBLISHED[k]
                ok &= abs(e.mean()) <= s_pub and e.std(ddof=1) <= 2 * s_pub and se >= 0.95
    res["criterion"] = ("best-of-two-leads |mean| <= published SD and SD <= 2 x published SD for QRS onset, QRS offset "
                        "and T end; sensitivity >= 95 %")
    res["VERIFIED"] = bool(ok)
    res["per_record"] = per_rec
    OUT.mkdir(parents=True, exist_ok=True)
    json.dump(res, open(OUT / "delineation_qtdb.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "per_record"}, indent=1))


if __name__ == "__main__":
    main()
