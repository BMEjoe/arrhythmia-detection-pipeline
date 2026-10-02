"""
Phase 9 source list with ground-truth labels and the model-level split (recorded in HANDOFF.md
before any detector or measure was run on a Phase 9 series).

  NULL          the 17 Phase 5-6 PASS nulls (families.NULLS); used in DEV (seeds 9000-9499) and
                TEST (seeds >= 9500) as null conditions.
  Phase 8       all 36 Phase 8 regimes with Phase 8's labels and split (mackey_glass DEV;
                phase_reset, coupled_vdp, av_node TEST).  As Phase 8, the ectopy variants keep
                the model's label (ectopy is added to the beat sequence, it does not drive the
                model).
  KTz (A2)      morphology family, TEST only (the only verified morphology family; the
                modified Luo-Rudy EAD family was dropped, METHODS.md).  A regime is a
                (P_nom, input) pair labelled CHAOTIC / NON-CHAOTIC by ground_truth_ktz.py for
                that exact input; AMBIGUOUS pairs are discarded.
"""
from __future__ import annotations

import json

from experiments.phase9_waveform import families as F


def sources():
    out = []
    for c in F.NULLS:
        out.append({"source": "null", "name": c, "label": "NULL", "split": "BOTH", "family": "null56",
                    "ectopy_fixed": None, "ectopic_type": F.NULL_TYPE[c]})
    for r in F.load_regimes():
        out.append({"source": "model", "name": r["name"], "label": r["label"], "split": r["split"],
                    "family": r["family"], "ectopy_fixed": None, "lam": r["lam"]})
    for r in json.load(open(F.HERE / "results" / "ground_truth" / "ktz_labels.json")):
        if r["label"] == "AMBIGUOUS":
            continue
        out.append({"source": "ktz", "name": f"ktz:P={r['P_nom']}", "label": r["label"], "split": "TEST",
                    "family": "ktz", "ectopy_fixed": r["input"], "lam": r["lam_per_beat"],
                    "apd_ref_ts": r["apd_ref_ts"]})
    return out


def key(s, ectopy=None):
    """Regime key used in tables: name, plus the input for KTz."""
    if s["source"] == "ktz":
        return f"{s['name']},in={s['ectopy_fixed']}"
    return s["name"] if ectopy in (None, "none") else f"{s['name']}|{ectopy}"
