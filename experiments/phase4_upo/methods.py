"""
Phase 4 methods: fixed decision rules applied to saved detector runs.

Each method names the RUN configuration it needs (detector.run arguments) and a
`detect(run_out) -> list of detected peaks` rule.  A window counts as a
detection when that list is non-empty.  Candidates are frozen in
PREREGISTRATION.md; nothing here may change after it is pushed.
"""
from __future__ import annotations

from experiments.phase4_upo import detector as D

PROD_RUN = ("cao", 7, 50)


def _baseline(run_out):
    """Unchanged production detector: >= 1 Level-B (surrogate-significant) peak."""
    return D.peaks_level_b(run_out)


METHODS = {
    "baseline": {"run": PROD_RUN, "detect": _baseline, "eligible": False},
}
