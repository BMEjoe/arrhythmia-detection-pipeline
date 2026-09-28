"""
Phase 4 analysis: evaluate every method (methods.METHODS) on a results file.

    python -m experiments.phase4_upo.analysis --phase dev  --file explore
    python -m experiments.phase4_upo.analysis --phase test --file test

Per method x window x condition:
    det        windows with >= 1 detected peak (the method's rule), k/N + Wilson CI
    levelC     windows with >= 1 detected peak that is Level-C verified_unstable
    and / or   combined binary detector with fp.lle_chaos_test on the same window
    lle        fp.lle_chaos_test detections alone
    loc        median |nearest detected peak - analytic on-attractor fixed point|
               over windows with a detection (chaotic systems only)
    status     run statuses (logistic_p4: embedding_not_saturated etc.)
Writes results/<phase>/tables/<file>_summary.{md,csv}.
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import pathlib
import statistics as st

from experiments.phase3_lle.analysis import frac
from experiments.phase4_upo import config as C
from experiments.phase4_upo import detector as D
from experiments.phase4_upo import methods as MM
from experiments.phase4_upo import systems as S

HERE = pathlib.Path(__file__).resolve().parent


def load(phase, file):
    return [json.loads(line) for line in open(HERE / "results" / phase / f"{file}.jsonl")]


def evaluate(recs, method):
    spec = MM.METHODS[method]
    key = D.run_key(*spec["run"])
    rows = {}
    for r in recs:
        t = r["task"]
        run_out = r["runs"].get(key)
        if run_out is None:
            continue
        lab = C.condition_label(t["system"], t["snr_db"])
        g = rows.setdefault((t["window_length"], lab), collections.defaultdict(list))
        det = spec["detect"](run_out)
        lc = bool(r["lle_chaos_test"]["detected"])
        g["det"].append(bool(det))
        g["levelC"].append(any(p["verified_unstable"] for p in det))
        g["and"].append(bool(det) and lc)
        g["or"].append(bool(det) or lc)
        g["lle"].append(lc)
        g["status"].append(run_out.get("status"))
        g["runtime"].append(r["runtime_s"].get(key))
        ref = S.fixed_point_reference(t["system"])
        if det and ref is not None:
            g["loc"].append(D.loc_error(det, ref))
    return rows


def summarize(recs, methods):
    out = []
    order = [C.condition_label(s, snr) for s, snr in C.conditions()]
    for m in methods:
        rows = evaluate(recs, m)
        for (n, lab) in sorted(rows, key=lambda k: (k[0], order.index(k[1]))):
            g = rows[(n, lab)]
            N = len(g["det"])
            loc = [v for v in g["loc"] if v is not None]
            out.append({
                "method": m, "window": n, "condition": lab, "N": N,
                "det_k": sum(g["det"]), "det": frac(sum(g["det"]), N),
                "levelC_k": sum(g["levelC"]), "and_k": sum(g["and"]), "or_k": sum(g["or"]),
                "lle_k": sum(g["lle"]),
                "loc_median": st.median(loc) if loc else None, "loc_n": len(loc),
                "status": "; ".join(f"{k}: {v}" for k, v in collections.Counter(g["status"]).most_common()),
                "runtime_median": st.median(v for v in g["runtime"] if v is not None),
            })
    return out


def write(phase, file, methods):
    recs = load(phase, file)
    rows = summarize(recs, methods)
    tdir = HERE / "results" / phase / "tables"
    tdir.mkdir(parents=True, exist_ok=True)
    with open(tdir / f"{file}_summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    md = [f"# Phase 4 summary -- {phase} seeds ({file})", "",
          "det = windows with >= 1 detected peak (method rule), k/N [Wilson 95 %]; levelC = detected and "
          "Level-C verified_unstable; AND / OR = combined with fp.lle_chaos_test on the same window; "
          "lle = lle_chaos_test alone; loc = median error of the nearest detected peak to the analytic "
          "on-attractor fixed point.", ""]
    for m in methods:
        md += [f"## {m}", "", "| N | condition | det | levelC | AND | OR | lle | loc (n) | status | s/run |",
               "|---:|---|---|---:|---:|---:|---:|---|---|---:|"]
        for r in rows:
            if r["method"] != m:
                continue
            loc = "" if r["loc_median"] is None else f"{r['loc_median']:.4f} ({r['loc_n']})"
            md.append(f"| {r['window']} | {r['condition']} | {r['det']} | {r['levelC_k']} | {r['and_k']} | "
                      f"{r['or_k']} | {r['lle_k']} | {loc} | {r['status']} | {r['runtime_median']:.1f} |")
        md.append("")
    (tdir / f"{file}_summary.md").write_text("\n".join(md))
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", choices=("dev", "test"), required=True)
    ap.add_argument("--file", required=True)
    ap.add_argument("--methods", default=None)
    a = ap.parse_args(argv)
    recs = load(a.phase, a.file)
    avail = set(recs[0]["runs"])
    methods = a.methods.split(",") if a.methods else [
        m for m, s in MM.METHODS.items() if D.run_key(*s["run"]) in avail]
    write(a.phase, a.file, methods)
    print((HERE / "results" / a.phase / "tables" / f"{a.file}_summary.md").read_text())


if __name__ == "__main__":
    main()
