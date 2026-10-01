"""
Phase 8 preregistered decision rule and tables (PREREGISTRATION.md Sections 4-5).

    python -m experiments.phase8_cardiac.decide            # results/test -> results/test/tables

PASS (per candidate, at its own window length): detections <= floor(0.07 N) on EACH of
  - the 17 Phase 5-6 null conditions (group "null"), and
  - EACH NON-CHAOTIC model regime, TEST and DEV (group "noncha", variant iv_q).
An analysis error counts as not detected.
WINNER: among passing candidates, the highest pooled detection count on the CHAOTIC
regimes of the TEST models at variant (vi) (group "chaos_primary", variants vi_S2_5,
vi_E1, vi_E3_10).  Tie-break: higher pooled detection on TEST chaotic regimes at (iv)
(variant iv_q), then candidate index (C1 < C2 < C3 < C4).
The baseline (frozen detector at 256) is reported, never eligible.
"""
from __future__ import annotations

import json
import math
import pathlib

import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
CANDIDATES = ("c1_frozen512", "c2_nlp_iaaft512", "c3_nlp_ep256", "c4_titration256")
LENGTH = {"c1_frozen512": 512, "c2_nlp_iaaft512": 512, "c3_nlp_ep256": 256, "c4_titration256": 256,
          "baseline_frozen256": 256}
VI = ("vi_S2_5", "vi_E1", "vi_E3_10")


def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return (math.nan, math.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def load(n):
    rows = []
    path = HERE / "results" / "test" / f"test_{n}.jsonl"
    seen = set()
    for line in open(path):
        r = json.loads(line)
        if r["task_id"] in seen:
            continue
        seen.add(r["task_id"])
        for m, res in r["results"].items():
            rows.append({k: r.get(k) for k in ("group", "regime", "label", "family", "split", "variant", "seed", "n")}
                        | {"method": m, "detected": bool(res.get("detected", False)), "error": res.get("error"),
                           "lle": res.get("lle_detected"), "upo": res.get("upo_detected"), "zmax": res.get("zmax"),
                           "NL": res.get("NL")})
    return pd.DataFrame(rows)


def main():
    out = HERE / "results" / "test" / "tables"
    out.mkdir(parents=True, exist_ok=True)
    df = pd.concat([load(256), load(512)], ignore_index=True)
    df = df[df.apply(lambda r: LENGTH[r["method"]] == r["n"], axis=1)]
    summary, passes, md = {}, {}, ["# Phase 8 TEST results (decide.py)", ""]
    for m in CANDIDATES + ("baseline_frozen256",):
        d = df[df.method == m]
        spec = d[d.group.isin(["null", "noncha"])].groupby("regime").agg(N=("detected", "size"), k=("detected", "sum"),
                                                                        errors=("error", lambda s: s.notna().sum()))
        spec["limit"] = (spec.N * 0.07).apply(math.floor)
        spec["within"] = spec.k <= spec.limit
        ok = bool(spec.within.all())
        prim = d[(d.group == "chaos_primary")]
        vi = prim[prim.variant.isin(VI)]
        iv = prim[prim.variant == "iv_q"]
        summary[m] = {"PASS": ok, "n_conditions_failed": int((~spec.within).sum()),
                      "failed": spec[~spec.within].index.tolist(),
                      "pooled_vi": int(vi.detected.sum()), "N_vi": int(len(vi)),
                      "pooled_iv": int(iv.detected.sum()), "N_iv": int(len(iv)), "length": LENGTH[m]}
        passes[m] = ok
        spec.to_csv(out / f"spec_{m}.csv")
        md += [f"## {m} (window {LENGTH[m]})", "",
               f"PASS: **{ok}**; failed conditions: {summary[m]['failed']}", "",
               f"Pooled TEST chaotic at (vi): {summary[m]['pooled_vi']}/{summary[m]['N_vi']}; at (iv): "
               f"{summary[m]['pooled_iv']}/{summary[m]['N_iv']}", "", spec.to_markdown(), ""]
    elig = [m for m in CANDIDATES if passes[m]]
    winner = None
    if elig:
        winner = sorted(elig, key=lambda m: (-summary[m]["pooled_vi"], -summary[m]["pooled_iv"], CANDIDATES.index(m)))[0]
    md = ["# Phase 8 primary decision", "", f"Passing candidates: {elig}", f"WINNER: **{winner}**", ""] + md
    # secondary: by group/regime/variant
    sec = df.groupby(["method", "group", "regime", "variant"]).agg(N=("detected", "size"), k=("detected", "sum")).reset_index()
    sec["rate"] = sec.k / sec.N
    sec.to_csv(out / "by_regime_variant.csv", index=False)
    json.dump({"summary": summary, "passing": elig, "winner": winner}, open(out / "decision.json", "w"), indent=1)
    open(out / "decision.md", "w").write("\n".join(md) + "\n")
    print("\n".join(md[:6]))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
