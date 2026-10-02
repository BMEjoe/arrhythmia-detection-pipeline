"""Summary of the B3 waveform-surrogate false-positive run (results/dev/fp_waveform.jsonl) ->
results/dev/fp_waveform_summary.md.  Rejection iff p <= 0.05 (19 surrogates, one-sided)."""
import json
import pathlib

import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
R = pd.DataFrame([json.loads(line) for line in open(HERE / "results" / "dev" / "fp_waveform.jsonl")])
rows = []
for (nm, nz), g in R.groupby(["name", "noise"]):
    row = {"null": nm, "noise": nz, "windows": len(g)}
    for s in ("pps", "cs", "ts"):
        for h in ("h1", "htau"):
            row[f"{s.upper()} {h}"] = f"{int((g[f'p_{s}_{h}'] <= 0.05).sum())}/{len(g)}"
    row["median s/window"] = round(float(g["t_total"].median()))
    rows.append(row)
t = pd.DataFrame(rows)
txt = ("# B3 waveform surrogate false positives (DEV, 2-min windows at 90 Hz)\n\n"
       "The run was interrupted by a container restart after 8 of the planned 65 windows and was NOT resumed: every surrogate already exceeded the 7 % "
       "UNUSABLE limit on the strictly periodic null (B3 rule: > 7 % on ANY null => UNUSABLE), so the "
       "remaining windows could not change any verdict; the per-window cost (B5) is also recorded.\n\n"
       + t.to_markdown(index=False) + "\n")
(HERE / "results" / "dev" / "fp_waveform_summary.md").write_text(txt)
print(txt)
