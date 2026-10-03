"""
Phase 10 markdown tables for the report, from results/<phase>/analysis*/analysis.json.

    python -m experiments.phase10_final.tables10 --phase conf [--dir analysis_sensitivity_flagged]
"""
from __future__ import annotations

import argparse
import json
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
FAM = {"henon": "Hénon", "logistic": "logistic", "coupled_vdp": "coupled vdP", "phase_reset": "phase-reset"}


def pct(x, d=1):
    return "–" if x is None else f"{100 * x:.{d}f} %"


def ci_s(c, d=1):
    return "–" if not c or c[0] is None else f"[{100 * c[0]:.{d}f}, {100 * c[1]:.{d}f}]"


def fmt_f(v):
    return "none ≤ 0.9" if v is None else f"{v:g}"


def q1_tables(A):
    out = []
    Q = A.get("Q1", {}).get("detectors", {})
    for det in ("K3", "K1", "LLE", "UPO"):
        D = Q.get(det)
        if not D:
            continue
        r = D["real"]
        out.append(f"#### {det}: real detections {r['k']}/{r['n_windows']} windows ({r['n_subjects']} subjects); "
                   f"U = {r['U']:.4f} (cluster bootstrap {r['U_boot95']:.4f}, Clopper–Pearson {r['U_cp95']:.4f})\n")
        fs = sorted(next(iter(D["curves"].values()))["empirical"], key=float)
        out.append(f"| family | parameter | λ | detected / base windows at f = {', '.join(fs)} | "
                   "replacement | π_upper at f = 0.9 | smallest f, π < 0.05 | smallest f, π < 0.20 | "
                   "(secondary: Wilson-lower p) π < 0.20 | (descriptive: Firth fit) π < 0.20 |")
        out.append("|---|---|---|---|---|---|---|---|---|---|")
        keys = sorted(D["curves"], key=lambda k: (list(FAM).index(D["curves"][k]["family"]), D["curves"][k]["lam"]))
        for k in keys:
            c = D["curves"][k]
            e = D["exclusion_L0"][k]
            emp = " · ".join(f"{c['empirical'][f]['k']}/{c['empirical'][f]['n']}" for f in sorted(c["empirical"], key=float))
            rp = c["replacement"]
            out.append(f"| {FAM[c['family']]} | {c['param']} | {c['lam']:.3g} | {emp} | {rp['k']}/{rp['n']} | "
                       f"{e['design']['pi_upper'][-1]:.3f} | {fmt_f(e['design']['f_min']['0.05'])} | "
                       f"{fmt_f(e['design']['f_min']['0.2'])} | {fmt_f(e['design_wilson_lower']['f_min']['0.2'])} | "
                       f"{fmt_f(e['p_hat']['f_min']['0.2'])} |")
        out.append("")
    return "\n".join(out)


def q2_tables(A):
    Q = A.get("Q2", {})
    if not Q:
        return ""
    out = [f"Segments: {Q['n_segments']} ({Q['n_gap']} not analysable: coverage < 90 % or < 400 intervals).\n",
           "| arm | NSR mean subject DR [95 % CI] | CHF mean subject DR [95 % CI] | CHF − NSR [95 % CI] | permutation p |",
           "|---|---|---|---|---|"]
    for key, name in (("P2_raw_by_group", "**raw (P2)**"), ("masked_by_group", "masked (analysable)"),
                      ("edited_by_group", "edited (eligible)"), ("wu_by_group", "Wu NN preprocessing")):
        g = Q.get(key, {})
        if "NSR" not in g:
            continue
        d = g.get("CHF_minus_NSR", {})
        out.append(f"| {name} | {pct(g['NSR']['mean_subject_rate'])} {ci_s(g['NSR']['ci95'])} (n = {g['NSR']['n_subjects']}) | "
                   f"{pct(g.get('CHF', {}).get('mean_subject_rate'))} {ci_s(g.get('CHF', {}).get('ci95'))} "
                   f"(n = {g.get('CHF', {}).get('n_subjects', 0)}) | {pct(d.get('diff'))} {ci_s(d.get('ci95'))} | "
                   f"{d.get('perm_p', float('nan')):.4f} |")
    out += ["", "| change vs raw | group | segments | raw DR | arm DR | mean change [95 % CI] | raw positives removed |",
            "|---|---|---|---|---|---|---|"]
    for key, name in (("P2_change_masked", "**masked (P2)**"), ("P2_change_masked_nonanalysable_negative",
                                                               "masked, non-analysable = negative"),
                      ("change_edited", "edited"), ("change_wu", "Wu NN")):
        for g, v in Q.get(key, {}).items():
            out.append(f"| {name} | {g} | {v['n_units']} | {pct(v['rate_a'])} | {pct(v['rate_b'])} | "
                       f"{pct(v['change'])} {ci_s(v['ci95'])} | {pct(v['a_pos_removed_frac'])} |")
    out += ["", f"Masked arm analysable fraction: {json.dumps({k: round(v, 3) for k, v in Q['masked_analysable_frac'].items()})}; "
                f"edited eligible fraction: {json.dumps({k: round(v, 3) for k, v in Q['edited_eligible_frac'].items()})}.",
            f"Mean NL among positive segments, raw: {json.dumps({k: round(v, 3) for k, v in Q['NL_among_positive_raw'].items()})}; "
            f"Wu: {json.dumps({k: round(v, 3) for k, v in Q['NL_among_positive_wu'].items()})}.",
            f"Raw DR night / day (clock time approximate): {json.dumps(Q['raw_night_day'], default=lambda o: round(o, 3))}."]
    g = Q.get("seg_gee_burden") or {}
    if g.get("estimated"):
        c = g["coef"]
        out.append(f"Segment GEE (raw positive ~ log2(1 + burden) + CHF): burden OR {c['x']['OR']:.2f} "
                   f"[{c['x']['OR_ci95'][0]:.2f}, {c['x']['OR_ci95'][1]:.2f}], p = {c['x']['p']:.2g}; CHF OR {c['chf']['OR']:.2f} "
                   f"[{c['chf']['OR_ci95'][0]:.2f}, {c['chf']['OR_ci95'][1]:.2f}], p = {c['chf']['p']:.2g}.")
    return "\n".join(out)


def gee_row(name, g):
    if not g.get("estimated"):
        return f"| {name} | {g['n']} | {g['n_pos']} | not estimable ({g.get('reason', '')}) | | |"
    c = g["coef"]
    s = f"| {name} | {g['n']} | {g['n_pos']} | {c['x']['OR']:.2f} [{c['x']['OR_ci95'][0]:.2f}, {c['x']['OR_ci95'][1]:.2f}] | {c['x']['p']:.2g} |"
    if "chf" in c:
        s += f" {c['chf']['OR']:.2f} [{c['chf']['OR_ci95'][0]:.2f}, {c['chf']['OR_ci95'][1]:.2f}], p = {c['chf']['p']:.2g} |"
    else:
        s += " – |"
    return s


def q3_tables(A, key="Q3"):
    Q = A.get(key, {})
    if not Q:
        return ""
    out = ["| outcome | windows | positive | OR per doubling of 1 + burden [95 % CI] | p | CHF OR (adjusted for burden) |",
           "|---|---|---|---|---|---|"]
    for c, nm in (("TIT", "**titration (P3)**"), ("K3", "**K3 (P3)**"), ("LLE", "LLE alone"), ("UPO", "UPO alone"), ("K1", "K1")):
        out.append(gee_row(nm, Q["P3_gee"][c]))
    out += ["", "| burden (beats) | windows | subjects | titration | LLE | UPO | K1 | K3 |", "|---|---|---|---|---|---|---|---|"]
    for b in Q["burden_table"]:
        cells = " | ".join(f"{b[c]['k']} ({pct(b[c]['rate'], 0)})" for c in ("TIT", "LLE", "UPO", "K1", "K3"))
        out.append(f"| {b['burden']} | {b['n']} | {b['n_subjects']} | {cells} |")
    out += ["", "| paired comparison (3b) | subset | windows | raw positive | other arm positive | both | raw positives removed | mean subject-level change [95 % CI] |",
            "|---|---|---|---|---|---|---|---|"]
    for name, v in Q["Q3b_paired"].items():
        for sub in ("all", "n_masked_ge1"):
            r = v[sub].get("all")
            if not r:
                continue
            out.append(f"| {name} | {sub} | {r['n_units']} | {r['a_pos']} | {r['b_pos']} | {r['both_pos']} | "
                       f"{pct(r['a_pos_removed_frac'])} | {pct(r['change'])} {ci_s(r['ci95'])} |")
    u = Q.get("Q3c_titration_group_unadjusted") or {}
    if u.get("estimated"):
        c = u["coef"]["chf"]
        out.append(f"\nQ3c: titration CHF OR unadjusted {c['OR']:.2f} [{c['OR_ci95'][0]:.2f}, {c['OR_ci95'][1]:.2f}], p = {c['p']:.2g}.")
    return "\n".join(out)


def q4_tables(A):
    Q = A.get("Q4", {})
    if not Q:
        return ""
    out = []
    if "Q4a_499" in Q:
        out += ["| detector | windows | default + / 499 + | both | default only | 499 only | agreement | κ | McNemar p |",
                "|---|---|---|---|---|---|---|---|---|"]
        for d, v in Q["Q4a_499"].items():
            k = "–" if v["kappa"] is None else f"{v['kappa']:.2f}"
            out.append(f"| {d} | {v['n']} | {v['both'] + v['a_only']} / {v['both'] + v['b_only']} | {v['both']} | "
                       f"{v['a_only']} | {v['b_only']} | {pct(v['agreement'])} | {k} | {v['mcnemar_exact_p']:.3g} |")
        out += ["", "| detector | windows (all 3 lengths) | 256 | 512 | 1024 |", "|---|---|---|---|---|"]
        for d, v in Q["Q4b_lengths"].items():
            out.append(f"| {d} | {v['n_windows']} | {pct(v['rate_256'])} | {pct(v['rate_512'])} | {pct(v['rate_1024'])} |")
    out += ["", "| detector | NSR night | NSR day | CHF night | CHF day |", "|---|---|---|---|---|"]
    for d, v in Q["Q4c_night_day"].items():
        cell = lambda x: f"{x['k']}/{x['n']}"  # noqa: E731
        out.append(f"| {d} | {cell(v['NSR']['night'])} | {cell(v['NSR']['day'])} | {cell(v['CHF']['night'])} | {cell(v['CHF']['day'])} |")
    g = Q.get("Q4c_titration_gee_night", {})
    if g.get("estimated"):
        c = g["coef"]["night_i"]
        out.append(f"\nTitration GEE night effect (adjusted for burden and group): OR {c['OR']:.2f} "
                   f"[{c['OR_ci95'][0]:.2f}, {c['OR_ci95'][1]:.2f}], p = {c['p']:.2g}.")
    L = Q["Q4d_labels"]
    out.append(f"\nQ4d rates: {json.dumps(L['rates'], default=lambda o: round(o, 4))}; K3 vs K4 agreement "
               f"{json.dumps({k: L['K3_vs_K4'][k] for k in ('both', 'a_only', 'b_only', 'neither')})}; K3 vs K3RR "
               f"{json.dumps({k: L['K3_vs_K3RR'][k] for k in ('both', 'a_only', 'b_only', 'neither')})}.")
    return "\n".join(out)


def q2e_table(A):
    Q = A.get("Q2e", {})
    if not Q:
        return ""
    out = ["| condition | n | windows | raw | raw at 1/128 s | masked (analysable) | edited (eligible) |", "|---|---|---|---|---|---|---|"]
    for k in sorted(Q, key=lambda s: (s.split("|")[0].startswith("vdp"), s)):
        c, n = k.split("|")
        v = Q[k]
        out.append(f"| {c} | {n} | {v['n_windows']} | {v['raw']} | {v.get('raw_q128', '–')} | {v['masked']} ({v['masked_analysable']}) | "
                   f"{v['edited']} ({v['edited_eligible']}) |")
    return "\n".join(out)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", choices=["dev", "conf"], required=True)
    ap.add_argument("--dir", default="analysis")
    a = ap.parse_args(argv)
    d = HERE / "results" / a.phase / a.dir
    A = json.load(open(d / "analysis.json"))
    md = [f"# Phase 10 tables ({a.phase}, {a.dir})", "", f"Windows {A['n_windows']}, subjects {A['n_subjects']}.", "",
          "## Q1 detection limits", q1_tables(A), "## Q2 titration (12-min segments)", q2_tables(A),
          "## Q2e synthetic titration checks", q2e_table(A), "## Q3 ectopy dose-response", q3_tables(A),
          "## Q4 robustness", q4_tables(A)]
    if "Q3_mitbih_exploratory" in A:
        md += ["## Q3 MIT-BIH (development, exploratory)", q3_tables(A, "Q3_mitbih_exploratory")]
    (d / "tables.md").write_text("\n\n".join(md) + "\n")
    print("wrote", d / "tables.md")


if __name__ == "__main__":
    main()
