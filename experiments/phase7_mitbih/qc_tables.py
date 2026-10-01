"""Markdown tables for QC.md from results/qc/*.csv (no detector output).
    python -m experiments.phase7_mitbih.qc_tables > experiments/phase7_mitbih/results/qc/tables.md"""
import os

import numpy as np
import pandas as pd

Q = os.path.join(os.path.dirname(__file__), "results", "qc")


def main():
    r = pd.read_csv(os.path.join(Q, "records.csv"), dtype={"record": str})
    w = pd.read_csv(os.path.join(Q, "windows.csv"), dtype={"record": str})
    a = pd.read_csv(os.path.join(Q, "windows_annotation_times.csv"), dtype={"record": str})
    c = pd.read_csv(os.path.join(Q, "candidates.csv"), dtype={"record": str})
    c150 = c[c.tol_ms == 150].groupby("record")
    print("## Per-record R-peak quality and windows\n")
    print("| record | subject | ref beats | detected | Se 75 ms | PPV 75 ms | Se 150 ms | PPV 150 ms | candidates | "
          "with extra det. (150) | with missed beat (150) | retained (150): N / A | retained 75 ms | annotation-time: N / A | NN>=0.80: N / A | flags |")
    print("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---|---|---|")
    for _, x in r.iterrows():
        ww = w[w.record == x.record]
        aa = a[a.record == x.record]
        cc = c150.get_group(x.record)
        flags = ", ".join(f for f, b in (("paced", x.paced), ("AF/AFL", x.af_record)) if b)
        print(f"| {x.record} | {x.subject} | {x.n_ref_beats} | {x.n_detected} | {x.sensitivity_75:.4f} | {x.PPV_75:.4f} | "
              f"{x.sensitivity_150:.4f} | {x.PPV_150:.4f} | {x.n_candidate_windows} | {(cc.n_extra_detections > 0).sum()} | "
              f"{(cc.n_missed_beats > 0).sum()} | {(ww.label == 0).sum()} / {(ww.label == 1).sum()} | {ww.retained_tol75.sum()} | "
              f"{(aa.label == 0).sum()} / {(aa.label == 1).sum()} | "
              f"{(ww.edited_eligible & (ww.label == 0)).sum()} / {(ww.edited_eligible & (ww.label == 1)).sum()} | {flags} |")
    for tag in ("75", "150"):
        tp, fpn, fn = r[f"TP_{tag}"].sum(), r[f"FP_{tag}"].sum(), r[f"FN_{tag}"].sum()
        print(f"\nTotal ({tag} ms): TP {tp}, FP {fpn}, FN {fn}, Se {tp/(tp+fn):.4f}, PPV {tp/(tp+fpn):.4f}")
    nonpaced = r[~r.paced]
    for tag in ("75", "150"):
        tp, fpn, fn = nonpaced[f"TP_{tag}"].sum(), nonpaced[f"FP_{tag}"].sum(), nonpaced[f"FN_{tag}"].sum()
        print(f"Non-paced 44 records ({tag} ms): Se {tp/(tp+fn):.4f}, PPV {tp/(tp+fpn):.4f}")
    print("\n## NN fraction (primary windows)\n")
    bins = [0, .2, .4, .6, .7, .8, .9, .95, 1.0001]
    labels = ["[0, 0.2)", "[0.2, 0.4)", "[0.4, 0.6)", "[0.6, 0.7)", "[0.7, 0.8)", "[0.8, 0.9)", "[0.9, 0.95)", "[0.95, 1]"]
    print("| NN fraction | " + " | ".join(labels) + " |")
    print("|---|" + "---:|" * len(labels))
    for lab, name in ((0, "normal"), (1, "abnormal")):
        h, _ = np.histogram(w[w.label == lab].nn_fraction, bins=bins)
        print(f"| {name} | " + " | ".join(str(v) for v in h) + " |")
    print("\n## HRV baseline by label (median [IQR]; primary raw-RR windows)\n")
    print("| label | n | SDNN ms | RMSSD ms | pNN50 % | mean RR s | abnormal fraction |")
    print("|---|---:|---|---|---|---|---|")
    for lab, name in ((0, "normal"), (1, "abnormal")):
        x = w[w.label == lab]
        f = lambda s: f"{s.median():.1f} [{s.quantile(.25):.1f}, {s.quantile(.75):.1f}]"
        g = lambda s: f"{s.median():.3f} [{s.quantile(.25):.3f}, {s.quantile(.75):.3f}]"
        print(f"| {name} | {len(x)} | {f(x.sdnn_ms)} | {f(x.rmssd_ms)} | {f(x.pnn50)} | {g(x.rr_mean)} | {g(x.abnormal_fraction)} |")
    x = w[w.edited_eligible]
    print("\nEdited-NN eligible windows (HRV of the edited series):\n")
    print("| label | n | SDNN ms | RMSSD ms | pNN50 % |")
    print("|---|---:|---|---|---|")
    for lab, name in ((0, "normal"), (1, "abnormal")):
        y = x[x.label == lab]
        f = lambda s: f"{s.median():.1f} [{s.quantile(.25):.1f}, {s.quantile(.75):.1f}]"
        print(f"| {name} | {len(y)} | {f(y.sdnn_ms_edited)} | {f(y.rmssd_ms_edited)} | {f(y.pnn50_edited)} |")


if __name__ == "__main__":
    main()
