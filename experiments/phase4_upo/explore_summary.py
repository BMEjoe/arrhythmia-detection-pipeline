"""
Phase 4 development exploration summary (DEVELOPMENT seeds only).

For every run configuration in results/dev/explore.jsonl and every gate rule,
count windows with >= 1 Level-B peak surviving the gate, per condition; plus
median localization error on clean Henon and runtime per run.

    python -m experiments.phase4_upo.explore_summary [--file explore]
"""
from __future__ import annotations

import argparse
import json
import pathlib
import statistics as st

from experiments.phase4_upo import config as C
from experiments.phase4_upo import detector as D
from experiments.phase4_upo import systems as S

HERE = pathlib.Path(__file__).resolve().parent

GATES = [("none", None, None)]
for stat in ("ext_lead", "src_median", "src_trim10"):
    for d in (0.1, 0.2, 0.3, 0.5):
        GATES.append((f"{stat}>={1 + d:.1f}", stat, d))
GATES.append(("levelC_verified", "levelC", None))


def apply(peaks, stat, d):
    if stat is None:
        return peaks
    if stat == "levelC":
        return [p for p in peaks if p["verified_unstable"]]
    return D.gate(peaks, stat, d)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", default="explore")
    ap.add_argument("--window", type=int, default=256)
    a = ap.parse_args(argv)
    R = [json.loads(line) for line in open(HERE / "results" / "dev" / f"{a.file}.jsonl")]
    R = [r for r in R if r["task"]["window_length"] == a.window]
    labs = [C.condition_label(s, snr) for s, snr in C.conditions()]
    present = [lab for lab in labs if any(C.condition_label(r["task"]["system"], r["task"]["snr_db"]) == lab for r in R)]
    runs = list(R[0]["runs"].keys())
    short = {"white_noise": "WN", "logistic_p4": "p4", "skewed_henon": "skH", "two_tone": "2tone"}
    head = [short.get(x, x).replace("logistic", "log").replace("henon", "hen") for x in present]
    print(f"detections (windows with >= 1 gated Level-B peak), dev seeds, N = {a.window}")
    print("N per condition: " + ", ".join(
        f"{h}={sum(1 for r in R if C.condition_label(r['task']['system'], r['task']['snr_db']) == lab)}"
        for h, lab in zip(head, present)))
    for run in runs:
        rt = [r["runtime_s"][run] for r in R]
        print(f"\n== run {run}  (median runtime {st.median(rt):.1f} s)")
        print("gate | " + " | ".join(head) + " | locErr hen (med)")
        for name, stat, d in GATES:
            cells, locs = [], []
            for lab in present:
                g = [r for r in R if C.condition_label(r["task"]["system"], r["task"]["snr_db"]) == lab]
                k = 0
                for r in g:
                    kept = apply(D.peaks_level_b(r["runs"][run]), stat, d)
                    k += bool(kept)
                    if lab == "henon" and kept:
                        locs.append(D.loc_error(kept, S.fixed_point_reference("henon")))
                cells.append(str(k))
            le = f"{st.median(locs):.3f}" if locs else "NA"
            print(f"{name} | " + " | ".join(cells) + f" | {le}")
    lc = {lab: sum(r["lle_chaos_test"]["detected"] for r in R
                   if C.condition_label(r["task"]["system"], r["task"]["snr_db"]) == lab) for lab in present}
    print("\nlle_chaos_test detections: " + ", ".join(f"{h}={lc[lab]}" for h, lab in zip(head, present)))


if __name__ == "__main__":
    main()
