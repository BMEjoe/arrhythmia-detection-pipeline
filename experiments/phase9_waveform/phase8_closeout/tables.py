"""
Phase 9 Part 0: descriptive tables for docs/PHASE8_CARDIAC_CHAOS.md, computed ONLY from the
committed Phase 8 TEST results (experiments/phase8_cardiac/results/test/test_{256,512}.jsonl).
No detector, generator or model is run.  Phase 8 files are read, never modified.

    python -m experiments.phase9_waveform.phase8_closeout.tables
"""
from __future__ import annotations

import json
import pathlib

import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[3]
P8 = ROOT / "experiments" / "phase8_cardiac" / "results" / "test"
OUT = pathlib.Path(__file__).resolve().parent
LENGTH = {"c1_frozen512": 512, "c2_nlp_iaaft512": 512, "c3_nlp_ep256": 256, "c4_titration256": 256,
          "baseline_frozen256": 256}
ORDER = list(LENGTH)


def load():
    rows = []
    for n in (256, 512):
        seen = set()
        for line in open(P8 / f"test_{n}.jsonl"):
            r = json.loads(line)
            if r["task_id"] in seen:
                continue
            seen.add(r["task_id"])
            for m, res in r["results"].items():
                if LENGTH[m] != n:
                    continue
                rows.append({k: r.get(k) for k in ("group", "regime", "label", "family", "split", "variant", "seed", "n")}
                            | {"method": m, "detected": bool(res.get("detected", False)),
                               "error": res.get("error"), "lle": res.get("lle_detected"),
                               "upo": res.get("upo_detected"), "runtime_s": res.get("runtime_s")})
    return pd.DataFrame(rows)


def rate(d):
    return f"{int(d.detected.sum())}/{len(d)}"


def main():
    df = load()
    md = ["# Phase 8 close-out tables (from committed TEST results only)", ""]
    md += [f"Windows: {df.groupby('n').size().to_dict()} method-window rows; "
           f"errors: {int(df.error.notna().sum())}", ""]
    # 1. TEST chaotic by family x variant group
    prim = df[df.group == "chaos_primary"].copy()
    prim["vgrp"] = prim.variant.map(lambda v: "iv" if v == "iv_q" else "vi")
    t = prim.groupby(["family", "vgrp", "method"]).apply(rate, include_groups=False).unstack("method")[ORDER]
    md += ["## TEST CHAOTIC regimes by family and variant (chaos_primary)", "", t.to_markdown(), ""]
    t = prim[prim.vgrp == "vi"].groupby(["variant", "method"]).apply(rate, include_groups=False).unstack("method")[ORDER]
    md += ["## TEST CHAOTIC at (vi) by ectopy pattern", "", t.to_markdown(), ""]
    t = prim.groupby(["regime", "method"]).apply(rate, include_groups=False).unstack("method")[ORDER]
    md += ["## TEST CHAOTIC by regime (iv + vi pooled, 120 windows each)", "", t.to_markdown(), ""]
    # 2. secondary variants
    sec = df[df.group == "chaos_second"]
    t = sec.groupby(["variant", "method"]).apply(rate, include_groups=False).unstack("method")[ORDER]
    md += ["## TEST CHAOTIC, secondary variants (chaos_second, 11 regimes x 10 seeds)", "", t.to_markdown(), ""]
    # 3. DEV chaos at test seeds
    dv = df[df.group == "dev_chaos"]
    t = dv.groupby(["variant", "method"]).apply(rate, include_groups=False).unstack("method")[ORDER]
    md += ["## DEV Mackey-Glass CHAOTIC at TEST seeds (dev_chaos)", "", t.to_markdown(), ""]
    # 4. non-chaotic under ectopy
    ne = df[df.group == "noncha_ect"]
    t = ne.groupby(["family", "variant", "method"]).apply(rate, include_groups=False).unstack("method")[ORDER]
    md += ["## NON-CHAOTIC regimes under ectopy + jitter (noncha_ect; secondary, not in the PASS rule)", "",
           t.to_markdown(), ""]
    # 5. flows
    fl = df[df.group == "flows"]
    t = fl.groupby(["regime", "method"]).apply(rate, include_groups=False).unstack("method")[ORDER]
    md += ["## Phase 5 flows G2 / G3 (not out-of-sample)", "", t.to_markdown(), ""]
    # 6. pooled null detections
    nul = df[df.group.isin(["null", "noncha"])]
    t = nul.groupby(["group", "method"]).apply(rate, include_groups=False).unstack("method")[ORDER]
    md += ["## Pooled detections on PASS conditions", "", t.to_markdown(), ""]
    # 7. components for frozen detectors
    fr = df[df.method.isin(["c1_frozen512", "baseline_frozen256"])].copy()
    fr["grp"] = fr.group.where(fr.group != "chaos_primary", fr.group + ":" + fr.variant.map(lambda v: "iv" if v == "iv_q" else "vi"))
    comp = fr.groupby(["grp", "method"]).agg(N=("detected", "size"), AND=("detected", "sum"),
                                              LLE=("lle", "sum"), UPO=("upo", "sum")).reset_index()
    md += ["## Components of the frozen detector (C1 at 512, baseline at 256)", "", comp.to_markdown(index=False), ""]
    # 8. runtime
    rt = df.groupby("method").runtime_s.median().round(2)
    md += ["## Median runtime per window (s, 4-worker load)", "", rt.to_frame().T[ORDER].to_markdown(index=False), ""]
    (OUT / "phase8_closeout_tables.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
