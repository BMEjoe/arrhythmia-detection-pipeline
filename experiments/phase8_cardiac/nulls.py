"""
Phase 5-6 null generators at any window length (Phase 8 PASS rule: N1-N6, S1, S2 5 % and
10 %, E1-E6, all quantized, at the candidate's window length).

The Phase 5 / Phase 6 generators read the window length from their config modules
(experiments.phase5_rr.config.N, experiments.phase6_robust.config.N) at call time and
include N in the seed sequence.  `with length(n):` sets both, so n = 256 reproduces the
Phase 5/6 series exactly and other lengths give independent realizations of the same
conditions.  Nothing in the Phase 5/6 code is modified.
"""
from __future__ import annotations

import contextlib
import threading

from experiments.phase5_rr import config as C5
from experiments.phase5_rr import systems as S5
from experiments.phase6_robust import config as C6
from experiments.phase6_robust import systems as S6

PASS_NULLS = ("N1_linear_rr", "N2_power_law", "N3_linear_rr_trend", "N4_linear_rr_step", "N5_linear_rr_warped",
              "N6_noisy_rsa", "S1_setar", "S2_ectopic_5pct", "S2_ectopic_10pct",
              "E1_bigeminy", "E2_trigeminy", "E3_couplets_5pct", "E3_couplets_10pct", "E4_runs",
              "E5_atrial_5pct", "E5_atrial_10pct", "E6_ectopic10_trend")
SECONDARY = ("G2_rossler_flow", "G3_mackey_glass", "P1_henon_rr", "P2_logistic_rr", "G1_lorenz_maxima")
_LOCK = threading.Lock()


@contextlib.contextmanager
def length(n):
    with _LOCK:
        old5, old6 = C5.N, C6.N
        C5.N = int(n)
        C6.N = int(n)
        try:
            yield
        finally:
            C5.N, C6.N = old5, old6


def generate(condition, seed, n=256):
    with length(n):
        return S6.generate(condition, seed)
