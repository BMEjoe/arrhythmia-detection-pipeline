"""
Download every PhysioNet database used in Phases 7-10, verify each file's SHA-256, and place it where each
phase's code reads it.

    python -m experiments.final_refinement.download_data --list            # what would be downloaded, and where
    python -m experiments.final_refinement.download_data                   # all databases
    python -m experiments.final_refinement.download_data --db nsr2db chf2db
    python -m experiments.final_refinement.download_data --verify-only     # re-check local files
    python -m experiments.final_refinement.download_data --pin             # maintainer: refresh data_manifests/

Every file is checked twice: against PhysioNet's SHA256SUMS.txt for the database version, and against the
pinned manifest committed in experiments/final_refinement/data_manifests/<db>.sha256 (the hashes of the files
this project used; for mitdb they equal experiments/phase7_mitbih/DATA_MANIFEST.json). A mismatch stops the
script. Files go to physionet_data/<db>/ (gitignored), and symbolic links make them visible where each phase
reads them (table below). No file of any database is committed to the repository.

Licences: all seven databases are distributed by PhysioNet under the Open Data Commons Attribution License v1.0
(stated on each database's PhysioNet page, checked 2026-10-03). Attribution is required: see README.md
"Data and licences" for the dataset DOIs and the publications PhysioNet asks users to cite.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import pathlib
import sys
import time
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[2]
CACHE = ROOT / "physionet_data"
PINS = pathlib.Path(__file__).resolve().parent / "data_manifests"
BASE = "https://physionet.org/files/{db}/{ver}/"

# database -> version, records (None = the database's RECORDS file), file extensions used, link locations
DBS = {
    "mitdb": {"ver": "1.0.0", "records": None, "ext": ["hea", "dat", "atr"],
              "links": ["mitdb_data", "experiments/phase10_final/dev_data/mitdb"],
              "used_by": "Phase 7 (MIT-BIH), Phase 9 exploratory, Phase 10 development"},
    "nsrdb": {"ver": "1.0.0", "records": None, "ext": ["hea", "atr"],
              "links": ["experiments/phase9_waveform/partf_data/nsrdb", "experiments/phase10_final/dev_data/nsrdb"],
              "used_by": "Phase 9 Part F, Phase 10 development"},
    "chfdb": {"ver": "1.0.0", "records": None, "ext": ["hea", "ecg"],
              "links": ["experiments/phase9_waveform/partf_data/chfdb", "experiments/phase10_final/dev_data/chfdb"],
              "used_by": "Phase 9 Part F, Phase 10 development"},
    "nsr2db": {"ver": "1.0.0", "records": None, "ext": ["hea", "ecg"],
               "links": ["experiments/phase10_final/conf_data/nsr2db"], "used_by": "Phase 10 confirmatory"},
    "chf2db": {"ver": "1.0.0", "records": None, "ext": ["hea", "ecg"],
               "links": ["experiments/phase10_final/conf_data/chf2db"], "used_by": "Phase 10 confirmatory"},
    "nstdb": {"ver": "1.0.0", "records": ["bw", "ma", "em"], "ext": ["hea", "dat"],
              "links": ["experiments/phase9_waveform/nstdb_data"], "used_by": "Phase 9 (recording noise)"},
    "qtdb": {"ver": "1.0.0", "records": None, "ext": ["hea", "dat", "q1c"],
             "links": ["experiments/phase9_waveform/qtdb_data"], "used_by": "Phase 9 (delineator check)"},
}


def get(url, tries=6):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read()
        except Exception:                                   # noqa: BLE001
            if k == tries - 1:
                raise
            time.sleep(2 ** (k + 1))


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def physionet_sums(db):
    txt = get(BASE.format(db=db, ver=DBS[db]["ver"]) + "SHA256SUMS.txt").decode()
    out = {}
    for line in txt.splitlines():
        if line.strip():
            h, name = line.split(None, 1)
            out[name.strip()] = h
    return out


def wanted(db, sums):
    spec = DBS[db]
    recs = spec["records"]
    if recs is None:
        recs = get(BASE.format(db=db, ver=spec["ver"]) + "RECORDS").decode().split()
    files = ["RECORDS"] if spec["records"] is None else []
    files += [f"{r}.{e}" for r in recs for e in spec["ext"]]
    missing = [f for f in files if f not in sums]
    if missing:
        raise SystemExit(f"{db}: not in SHA256SUMS.txt: {missing[:5]}")
    return files


def pinned(db):
    p = PINS / f"{db}.sha256"
    if not p.exists():
        return None
    return {name: h for h, name in (l.split(None, 1) for l in p.read_text().splitlines() if l.strip())}


def pin(dbs):
    PINS.mkdir(exist_ok=True)
    for db in dbs:
        sums = physionet_sums(db)
        files = wanted(db, sums)
        lines = [f"{sums[f]}  {f}" for f in files]
        (PINS / f"{db}.sha256").write_text("\n".join(lines) + "\n")
        print(f"{db}: pinned {len(files)} files from PhysioNet SHA256SUMS.txt ({DBS[db]['ver']})")


def link(db):
    for rel in DBS[db]["links"]:
        dst = ROOT / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        target = os.path.relpath(CACHE / db, dst.parent)
        if dst.is_symlink():
            dst.unlink()
        elif dst.exists():
            print(f"  {rel} exists and is not a link; left unchanged (files are verified in place by --verify-only)")
            continue
        dst.symlink_to(target)
        print(f"  linked {rel} -> {target}")


def download(dbs, dry=False):
    for db in dbs:
        spec = DBS[db]
        pins = pinned(db)
        if pins is None:
            raise SystemExit(f"{db}: no pinned manifest (run --pin once)")
        sums = physionet_sums(db)
        files = wanted(db, sums)
        print(f"{db} {spec['ver']}: {len(files)} files -> physionet_data/{db}/ ({spec['used_by']})")
        if dry:
            continue
        d = CACHE / db
        d.mkdir(parents=True, exist_ok=True)
        (d / "SHA256SUMS.txt").write_text("".join(f"{sums[f]}  {f}\n" for f in files))
        for i, f in enumerate(files):
            p = d / f
            if not p.exists() or sha256(p) != sums[f]:
                p.write_bytes(get(BASE.format(db=db, ver=spec["ver"]) + f))
            h = sha256(p)
            if h != sums[f]:
                raise SystemExit(f"{db}/{f}: SHA-256 differs from PhysioNet SHA256SUMS.txt")
            if pins.get(f) != h:
                raise SystemExit(f"{db}/{f}: SHA-256 differs from the pinned manifest data_manifests/{db}.sha256")
            if (i + 1) % 50 == 0:
                print(f"  {i + 1}/{len(files)}")
        print(f"  {len(files)} files verified (PhysioNet SHA256SUMS.txt and pinned manifest)")
        link(db)


def verify(dbs):
    bad = 0
    for db in dbs:
        pins = pinned(db) or {}
        d = CACHE / db
        if not d.exists():
            alt = [ROOT / l for l in DBS[db]["links"] if (ROOT / l).exists()]
            d = alt[0] if alt else d
        n = 0
        for f, h in pins.items():
            p = d / f
            if not p.exists():
                print(f"{db}/{f}: missing")
                bad += 1
            elif sha256(p) != h:
                print(f"{db}/{f}: SHA-256 MISMATCH")
                bad += 1
            else:
                n += 1
        print(f"{db}: {n}/{len(pins)} files verified in {d.relative_to(ROOT) if d.is_relative_to(ROOT) else d}")
    return bad


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", nargs="+", choices=list(DBS), default=list(DBS))
    ap.add_argument("--list", action="store_true", help="list databases, files and link locations; no download")
    ap.add_argument("--verify-only", action="store_true")
    ap.add_argument("--pin", action="store_true", help="maintainer: write data_manifests/ from PhysioNet SHA256SUMS")
    a = ap.parse_args(argv)
    if a.pin:
        pin(a.db)
    elif a.list:
        for db in a.db:
            s = DBS[db]
            p = pinned(db)
            print(f"{db} {s['ver']}: {len(p) if p else '?'} files ({', '.join(s['ext'])}); used by {s['used_by']}; "
                  f"linked at {', '.join(s['links'])}")
    elif a.verify_only:
        sys.exit(1 if verify(a.db) else 0)
    else:
        download(a.db)


if __name__ == "__main__":
    main()
