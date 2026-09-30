"""
Checks that the direct component calls used for the m = 3 / 4 arms
(detector.direct) reproduce analyze_segment's component decisions at m = 2 on
development windows (2 seeds per condition).  Writes results/dev/m_path_check.json.

    python -m experiments.phase5_rr.m_path_check
"""
from __future__ import annotations

import json
import pathlib
import sys
import warnings

ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.phase5_rr import config as C  # noqa: E402
from experiments.phase5_rr import detector as D  # noqa: E402
from experiments.phase5_rr import systems as S  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent


def main():
    rows = []
    for cond in C.CONDITIONS:
        for seed in (0, 1):
            rr = S.generate(cond, seed)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                a = D.decide(rr, 2)
                b = D.direct(rr, D.config_m(2))
            same = all(a[k] == b[k] for k in ("upo", "lle", "lle_p", "upo_n_level_b", "upo_n_gated")) \
                if a["error"] is None else None
            rows.append({"condition": cond, "seed": seed, "same": same, "error": a["error"]})
    n_diff = sum(r["same"] is False for r in rows)
    out = HERE / "results" / "dev" / "m_path_check.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"n": len(rows), "n_differences": n_diff, "rows": rows}, indent=1))
    print(f"{len(rows)} windows, {n_diff} differences, {sum(r['error'] is not None for r in rows)} errors")


if __name__ == "__main__":
    main()
