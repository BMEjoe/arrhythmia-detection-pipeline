"""
Checks the runner's exact-sharing shortcut (run_phase6.execute, methods.effective_key)
on DEVELOPMENT windows: every preregistered method is run in full and via the
shortcut; LLE, UPO and AND must agree.  Writes results/dev/shortcut_check.json.

    python -m experiments.phase6_robust.shortcut_check
"""
from __future__ import annotations

import json
import multiprocessing as mp
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.phase6_robust import methods as MM  # noqa: E402
from experiments.phase6_robust import run_phase6 as R  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
WINDOWS = ([("P3_henon_rr_trend", s) for s in range(5)] + [("P3_logistic_rr_trend", s) for s in range(3)]
           + [("G1_lorenz_maxima", s) for s in range(5)] + [("P2_logistic_rr", s) for s in range(5)]
           + [("N2_power_law", s) for s in range(8)] + [("N4_linear_rr_step", s) for s in range(3)]
           + [("S2_ectopic_10pct", s) for s in range(3)] + [("E6_ectopic10_trend", s) for s in range(3)]
           + [("E1_bigeminy", s) for s in range(2)] + [("E1_bigeminy_edited", s) for s in range(2)])


def both(w):
    task = {"task_id": f"check:{w[0]}:{w[1]}", "phase": "dev", "condition": w[0], "seed": w[1]}
    full = R.execute((task, list(MM.PREREGISTERED), False))
    fast = R.execute((task, list(MM.PREREGISTERED), True))
    diffs = [m for m in MM.PREREGISTERED
             if any(full["methods"][m].get(k) != fast["methods"][m].get(k) for k in ("lle", "upo", "and", "lle_p"))]
    return {"condition": w[0], "seed": w[1], "trend_ratio": full["trend_ratio"], "computed_as": fast["computed_as"],
            "differences": diffs}


def main():
    with mp.get_context("fork").Pool(4) as pool:
        rows = pool.map(both, WINDOWS)
    n = sum(bool(r["differences"]) for r in rows)
    shared = sum(v != m for r in rows for m, v in r["computed_as"].items())
    out = HERE / "results" / "dev" / "shortcut_check.json"
    out.write_text(json.dumps({"n_windows": len(rows), "n_windows_with_differences": n,
                               "n_shared_method_runs": shared, "rows": rows}, indent=1))
    print(f"{len(rows)} windows, {n} with differences, {shared} method runs shared")


if __name__ == "__main__":
    main()
