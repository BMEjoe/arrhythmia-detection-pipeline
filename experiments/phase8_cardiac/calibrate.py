"""
Phase 8 threshold calibration for C2 and C3 (DEVELOPMENT data only).

Z = smallest value such that each of the 17 PASS null conditions has <= 2 % of its
development windows (seeds 0-99 at the candidate's length) at or above Z, plus 0.5.

    python -m experiments.phase8_cardiac.calibrate
"""
from __future__ import annotations

import json
import pathlib

import numpy as np
import pandas as pd

from experiments.phase8_cardiac import nulls as NL

HERE = pathlib.Path(__file__).resolve().parent
DEV = HERE / "results" / "dev"
SPEC = {"c2_nlp_iaaft512": ("nlpred_scan", 512, ["nulls512", "calib_ia512"], "mg512"),
        "c3_nlp_ep256": ("nlpred_ep_scan", 256, ["nulls256", "calib_ep256"], "mg256")}


def load(tag):
    return pd.DataFrame([json.loads(l) for l in open(DEV / f"{tag}.jsonl")])


def main():
    out = {}
    for cand, (meth, n, tags, mgtag) in SPEC.items():
        d = pd.concat([load(t) for t in tags], ignore_index=True)
        d = d[(d.method == meth) & d.regime.isin(NL.PASS_NULLS) & (d.n == n)].drop_duplicates("task_id")
        d["zmax"] = d[[f"np_m{m}_z" for m in (2, 3, 4, 5)]].max(axis=1)
        per = {}
        for c, g in d.groupby("regime"):
            z = np.sort(g.zmax.to_numpy())
            k_allowed = int(np.floor(0.02 * len(z)))
            per[c] = {"N": int(len(z)), "z_needed": float(z[len(z) - 1 - k_allowed] + 1e-9),
                      "max": float(z.max()), "q95": float(np.quantile(z, 0.95))}
        Z = max(v["z_needed"] for v in per.values()) + 0.5
        mg = load(mgtag)
        mg = mg[mg.method == meth]
        mg["zmax"] = mg[[f"np_m{m}_z" for m in (2, 3, 4, 5)]].max(axis=1)
        nonc = mg[(mg.label == "NON-CHAOTIC") & (mg.variant == "iv_q")].zmax
        ch = mg[mg.label == "CHAOTIC"]
        out[cand] = {"Z": round(float(Z), 2), "per_condition": per, "binding_condition":
                     max(per, key=lambda c: per[c]["z_needed"]),
                     "dev_mg_noncha_iv_max": float(nonc.max()),
                     "dev_mg_power_by_variant": ch.groupby("variant").zmax.apply(lambda s: float((s >= Z).mean())).to_dict()}
        print(cand, "Z =", out[cand]["Z"], "binding", out[cand]["binding_condition"],
              "MG noncha iv max", round(out[cand]["dev_mg_noncha_iv_max"], 2))
        print("  dev MG power", {k: round(v, 2) for k, v in out[cand]["dev_mg_power_by_variant"].items()})
    json.dump(out, open(DEV / "calibration.json", "w"), indent=1)


if __name__ == "__main__":
    main()
