"""
Phase 5 step 0: the analyze_segment path gives the same decisions as the stored
Phase 4 results (C2 UPO detection, lle_chaos_test) on 5 windows per Phase 4
condition (test seeds 2000-2004, 256 samples; experiments/phase4_upo/results/test/test.jsonl).

    python -m experiments.phase5_rr.equivalence_check
"""
from __future__ import annotations

import json
import pathlib
import sys
import warnings

ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.phase4_upo import config as C4  # noqa: E402
from experiments.phase4_upo import methods as MM  # noqa: E402
from experiments.phase4_upo import systems as S4  # noqa: E402
from experiments.phase5_rr import detector as D5  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
STORED = ROOT / "experiments/phase4_upo/results/test/test.jsonl"
SEEDS = range(2000, 2005)
N = 256


def main():
    stored = {}
    for line in open(STORED):
        r = json.loads(line)
        t = r["task"]
        if t["window_length"] == N and t["seed"] in SEEDS:
            stored[(t["system"], t["snr_db"], t["seed"])] = r
    rows, n_diff = [], 0
    for system, snr in C4.conditions():
        for seed in SEEDS:
            r = stored[(system, snr, seed)]
            x = S4.generate(system, N, seed, snr)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                d = D5.decide(x)
            c2_ref = bool(MM.METHODS["c2_m2M15_mediangate"]["detect"](r["runs"]["m2_M15_S50"]))
            lle_ref = bool(r["lle_chaos_test"]["detected"])
            p_ref = r["lle_chaos_test"]["p"]
            same = (d["upo"] == c2_ref and d["lle"] == lle_ref and d["error"] is None
                    and (p_ref is None or abs(d["lle_p"] - p_ref) < 1e-12))
            n_diff += not same
            rows.append({"condition": C4.condition_label(system, snr), "seed": seed,
                         "upo_stored": c2_ref, "upo_analyze_segment": d["upo"],
                         "lle_stored": lle_ref, "lle_analyze_segment": d["lle"],
                         "lle_p_stored": p_ref, "lle_p_analyze_segment": d["lle_p"],
                         "error": d["error"], "same": same})
    out = HERE / "results" / "equivalence_check.json"
    out.write_text(json.dumps({"n_windows": len(rows), "n_differences": n_diff, "rows": rows}, indent=1))
    md = ["# Step 0: analyze_segment path vs stored Phase 4 decisions", "",
          f"{len(rows)} windows (14 Phase 4 conditions x test seeds 2000-2004, 256 samples). "
          f"Differences: **{n_diff}**.", "",
          "| condition | UPO C2 detections (stored / analyze_segment) | lle_chaos_test detections (stored / analyze_segment) | identical windows |",
          "|---|---|---|---|"]
    for lab in dict.fromkeys(r["condition"] for r in rows):
        rs = [r for r in rows if r["condition"] == lab]
        md.append(f"| {lab} | {sum(r['upo_stored'] for r in rs)} / {sum(r['upo_analyze_segment'] for r in rs)} | "
                  f"{sum(r['lle_stored'] for r in rs)} / {sum(r['lle_analyze_segment'] for r in rs)} | "
                  f"{sum(r['same'] for r in rs)}/{len(rs)} |")
    (HERE / "results" / "equivalence_check.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))
    return n_diff


if __name__ == "__main__":
    sys.exit(1 if main() else 0)
