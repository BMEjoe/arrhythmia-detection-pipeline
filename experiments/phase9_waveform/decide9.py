"""
Phase 9 Part E: the preregistered decision (PREREGISTRATION.md Section 4) from results/test/test9.jsonl.

    python -m experiments.phase9_waveform.decide9

PASS iff detections <= floor(0.07 N) in EACH condition of the groups null_input, null_real,
noncha_input, noncha_real (condition = group + regime; N = 100 each); an error / unanalysable window
counts as not detected.
WINNER: among passing candidates, the highest pooled detection count in chaos_primary (1,320
windows); ties: fewer pooled detections over all PASS conditions, then candidate index.
Outputs results/test/tables/decision9.md, decision9.json, by_condition9.csv, secondary9.csv.
"""
from __future__ import annotations

import json
import math
import pathlib

import pandas as pd

from experiments.phase9_waveform import candidates9 as C

HERE = pathlib.Path(__file__).resolve().parent
PASS_GROUPS = ("null_input", "null_real", "noncha_input", "noncha_real")


def load():
    rows = []
    for line in open(HERE / "results" / "test" / "test9.jsonl"):
        r = json.loads(line)
        base = {k: r.get(k) for k in ("group", "name", "label", "family", "ectopy", "level", "seed", "task_id")}
        base["error"] = r.get("error")
        for c in C.CANDIDATES:
            res = (r.get("results") or {}).get(c, {})
            rows.append({**base, "cand": c, "det": bool(res.get("detected", False)),
                         "cerr": res.get("error"), "analysable": res.get("analysable")})
    return pd.DataFrame(rows)


def main():
    D = load()
    D["cond"] = D["group"] + "|" + D["name"]
    out = {"candidates": {}}
    lines = ["# Phase 9 primary decision (decide9.py)\n"]
    spec_rows = []
    for c in C.CANDIDATES:
        d = D[D.cand == c]
        p = d[d.group.isin(PASS_GROUPS)].groupby("cond")["det"].agg(["sum", "count"])
        p["limit"] = p["count"].map(lambda n: math.floor(0.07 * n))
        p["within"] = p["sum"] <= p["limit"]
        failed = list(p.index[~p["within"]])
        prim = d[d.group == "chaos_primary"]
        out["candidates"][c] = {"PASS": len(failed) == 0, "failed": failed,
                                "primary_detected": int(prim.det.sum()), "primary_N": int(len(prim)),
                                "pass_pool_detected": int(p["sum"].sum()), "pass_pool_N": int(p["count"].sum()),
                                "window_errors": int(d.error.notna().sum()), "candidate_errors": int(d.cerr.notna().sum())}
        for cond, row in p.iterrows():
            spec_rows.append({"cand": c, "cond": cond, "k": int(row["sum"]), "N": int(row["count"]),
                              "limit": int(row["limit"]), "within": bool(row["within"])})
    passing = [c for c in C.CANDIDATES if out["candidates"][c]["PASS"]]
    if passing:
        win = sorted(passing, key=lambda c: (-out["candidates"][c]["primary_detected"],
                                             out["candidates"][c]["pass_pool_detected"], C.CANDIDATES.index(c)))[0]
    else:
        win = None
    out["passing"], out["winner"] = passing, win
    lines.append(f"Passing candidates: {passing}\n\nWINNER: **{win}**\n")
    for c, v in out["candidates"].items():
        lines.append(f"## {c}\n\nPASS: **{v['PASS']}**; failed conditions ({len(v['failed'])}): {v['failed']}\n\n"
                     f"Pooled TEST chaotic (primary): {v['primary_detected']}/{v['primary_N']}; pooled PASS-condition "
                     f"detections: {v['pass_pool_detected']}/{v['pass_pool_N']}; window errors {v['window_errors']}, "
                     f"candidate errors {v['candidate_errors']}\n")
    S = pd.DataFrame(spec_rows)
    sec = (D.groupby(["cand", "group", "family", "name", "ectopy", "level"], dropna=False)["det"]
           .agg(["sum", "count"]).reset_index())
    tdir = HERE / "results" / "test" / "tables"
    tdir.mkdir(parents=True, exist_ok=True)
    S.to_csv(tdir / "by_condition9.csv", index=False)
    sec.to_csv(tdir / "secondary9.csv", index=False)
    json.dump(out, open(tdir / "decision9.json", "w"), indent=1)
    (tdir / "decision9.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
