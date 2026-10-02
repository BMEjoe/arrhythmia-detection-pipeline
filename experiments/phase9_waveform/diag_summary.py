"""
Summaries of the development diagnosis bank (results/dev/diag.jsonl) -> results/dev/diag_summary.md.

    python -m experiments.phase9_waveform.diag_summary

Groups: CHAOTIC (Mackey-Glass chaotic regimes, Phase 5 positives), MG-NONCHA (Mackey-Glass
non-chaotic), QP (circle map quasi-periodic), LOCKED (circle map locked), NULL-in (nulls at the
input level), NULL-ecg (nulls at the realistic ECG level).
For each statistic: z = (data - surrogate mean) / surrogate SD (19 ectopy-preserving IAAFT), and
the raw data value; group medians, and the fraction of windows with z >= 3 / z <= -3.
"""
from __future__ import annotations

import json
import pathlib

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent


def load():
    R = pd.DataFrame([json.loads(line) for line in open(HERE / "results" / "dev" / "diag.jsonl")])
    R = R[R.get("error").isna()] if "error" in R else R

    def grp(r):
        if r["label"] == "NULL":
            return "NULL-in" if r["kind"] == "rr_null" else "NULL-ecg"
        if r["source"] == "devmap":
            return "QP" if "Om=0.382" in r["name"] or "Om=0.618" in r["name"] else "LOCKED"
        if r["label"] == "CHAOTIC":
            return "CHAOTIC-pos" if r["kind"] == "rr_null" else ("CHAOTIC-mg-ecg" if r["kind"] == "ecg" else "CHAOTIC-mg-rr")
        return "MG-NONCHA-ecg" if r["kind"] == "ecg" else "MG-NONCHA-rr"
    R["group"] = R.apply(grp, axis=1)
    return R


def stat_names(R):
    return sorted(c[:-4] for c in R.columns if c.endswith("__sm"))


def zframe(R):
    Z = {}
    for k in stat_names(R):
        Z[k] = (R[k] - R[k + "__sm"]) / R[k + "__ss"].replace(0, np.nan)
    return pd.DataFrame(Z)


def main():
    R = load()
    Z = zframe(R)
    Z["group"] = R["group"].values
    groups = ["CHAOTIC-mg-rr", "CHAOTIC-mg-ecg", "CHAOTIC-pos", "MG-NONCHA-rr", "MG-NONCHA-ecg", "QP", "LOCKED",
              "NULL-in", "NULL-ecg"]
    med = Z.groupby("group").median(numeric_only=True).T.reindex(columns=groups)
    raw = R.groupby("group")[stat_names(R) + ["c2_zmax", "c3_zmax"]].median(numeric_only=True).T.reindex(columns=groups)
    cnt = R.groupby("group").size().reindex(groups)
    txt = "# Development diagnosis bank summary\n\nWindows per group:\n\n" + cnt.to_frame("n").to_markdown() + "\n\n"
    txt += "## Median z (data vs 19 EP-IAAFT surrogates)\n\n" + med.round(2).to_markdown() + "\n\n"
    txt += "## Median raw data values\n\n" + raw.round(3).to_markdown() + "\n"
    (HERE / "results" / "dev" / "diag_summary.md").write_text(txt)
    print(txt)


if __name__ == "__main__":
    main()
