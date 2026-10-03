"""
Phase 10 data layer: download (SHA-256 checked against PhysioNet SHA256SUMS.txt), beat series, labels,
window sampling.  DEVELOPMENT databases: mitdb, nsrdb, chfdb.  CONFIRMATORY: nsr2db, chf2db (download is
refused until PREREGISTRATION.md is committed and pushed; see run10.py).

    python -m experiments.phase10_final.data10 --download dev
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import urllib.request

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
DEV_DIR = HERE / "dev_data"
CONF_DIR = HERE / "conf_data"
DEV_DBS = ("mitdb", "nsrdb", "chfdb")
CONF_DBS = ("nsr2db", "chf2db")
ANN_EXT = {"mitdb": "atr", "nsrdb": "atr", "chfdb": "ecg", "nsr2db": "ecg", "chf2db": "ecg"}
BASE = "https://physionet.org/files/{db}/1.0.0/"

# Beat symbols (wfdb beat codes, as Phase 9 Part F); 'N' label iff normal / bundle-branch / escape of
# supraventricular origin with normal conduction sequence (Phase 9 Part F rule).
BEATS = set("NLRBAaJSVrFejnE/fQ?")
NORMAL = set("NLRBej")
RR_MIN, RR_MAX = 0.25, 2.5     # an interval outside marks its ending beat non-normal (Phase 9 Part F)


def _get(url, tries=6):
    import time
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read()
        except Exception:                                   # noqa: BLE001
            if k == tries - 1:
                raise
            time.sleep(2 ** (k + 1))


def download(db, root):
    d = root / db
    d.mkdir(parents=True, exist_ok=True)
    sums = _get(BASE.format(db=db) + "SHA256SUMS.txt").decode()
    (d / "SHA256SUMS.txt").write_text(sums)
    table = {}
    for line in sums.splitlines():
        if line.strip():
            h, name = line.split(None, 1)
            table[name.strip()] = h
    recs = _get(BASE.format(db=db) + "RECORDS").decode().split()
    (d / "RECORDS").write_text("\n".join(recs) + "\n")
    want = ["RECORDS"] + [f"{r}.hea" for r in recs] + [f"{r}.{ANN_EXT[db]}" for r in recs]
    manifest = {}
    for name in want:
        p = d / name
        if not p.exists():
            p.write_bytes(_get(BASE.format(db=db) + name))
        h = hashlib.sha256(p.read_bytes()).hexdigest()
        if table.get(name) != h:
            raise SystemExit(f"SHA-256 mismatch for {db}/{name}")
        manifest[name] = h
    (d / "MANIFEST.json").write_text(json.dumps(manifest, indent=1))
    return recs


def records(db, root):
    return (root / db / "RECORDS").read_text().split()


def beat_series(db, rec, root):
    """Annotation beat times (s), labels ('N' / 'X'), the raw symbols, and fs of the annotation file."""
    import wfdb
    a = wfdb.rdann(str(root / db / rec), ANN_EXT[db])
    fs = float(a.fs) if a.fs else float(wfdb.rdheader(str(root / db / rec)).fs)
    sym = np.array(a.symbol)
    smp = np.asarray(a.sample)
    keep = np.array([s in BEATS for s in sym], bool)
    sym, smp = sym[keep], smp[keep]
    t = smp / fs
    lab = np.array(["N" if s in NORMAL else "X" for s in sym])
    rr = np.diff(t)
    lab[1:][(rr < RR_MIN) | (rr > RR_MAX)] = "X"
    return t, lab, sym, fs


def header_info(db, rec, root):
    import wfdb
    h = wfdb.rdheader(str(root / db / rec))
    return {"fs": float(h.fs), "sig_len": int(h.sig_len), "base_time": str(h.base_time) if h.base_time else None,
            "base_date": str(h.base_date) if h.base_date else None, "comments": list(h.comments or [])}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--download", choices=["dev"], required=True)
    a = ap.parse_args(argv)
    for db in DEV_DBS:
        recs = download(db, DEV_DIR)
        print(db, len(recs), "records; SHA-256 OK")


if __name__ == "__main__":
    main()
