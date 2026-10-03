# Phase 9 primary decision (decide9.py)

Passing candidates: ['k3_mnlp_growth_ann', 'k4_mnlp_growth_rr', 'k5_mnlp_fsle_ann']

WINNER: **k3_mnlp_growth_ann**

## k1_frozen512

PASS: **False**; failed conditions (4): ['noncha_real|coupled_vdp:rho=10.0,omega=3.3', 'noncha_real|coupled_vdp:rho=5.45,omega=5.6', 'noncha_real|coupled_vdp:rho=6.0,omega=5.6', 'noncha_real|coupled_vdp:rho=9.6,omega=2.1']

Pooled TEST chaotic (primary): 74/1320; pooled PASS-condition detections: 97/10100; window errors 63, candidate errors 0

## k2_mnlp_ann

PASS: **False**; failed conditions (12): ['noncha_input|coupled_vdp:rho=2.0,omega=5.6', 'noncha_input|coupled_vdp:rho=4.0,omega=5.6', 'noncha_input|coupled_vdp:rho=5.45,omega=5.6', 'noncha_input|coupled_vdp:rho=9.6,omega=2.1', 'noncha_real|coupled_vdp:rho=10.0,omega=3.3', 'noncha_real|coupled_vdp:rho=2.0,omega=5.6', 'noncha_real|coupled_vdp:rho=4.0,omega=5.6', 'noncha_real|coupled_vdp:rho=5.45,omega=5.6', 'noncha_real|coupled_vdp:rho=6.0,omega=5.6', 'noncha_real|coupled_vdp:rho=9.6,omega=2.1', 'noncha_real|mackey_glass:tau=15.0', 'noncha_real|mackey_glass:tau=16.0']

Pooled TEST chaotic (primary): 304/1320; pooled PASS-condition detections: 607/10100; window errors 63, candidate errors 0

## k3_mnlp_growth_ann

PASS: **True**; failed conditions (0): []

Pooled TEST chaotic (primary): 90/1320; pooled PASS-condition detections: 0/10100; window errors 63, candidate errors 0

## k4_mnlp_growth_rr

PASS: **True**; failed conditions (0): []

Pooled TEST chaotic (primary): 2/1320; pooled PASS-condition detections: 1/10100; window errors 63, candidate errors 0

## k5_mnlp_fsle_ann

PASS: **True**; failed conditions (0): []

Pooled TEST chaotic (primary): 0/1320; pooled PASS-condition detections: 3/10100; window errors 63, candidate errors 0


---
**Erratum E1 (final refinement, 2026-10-03).** This file is the unmodified output of the frozen `decide9.py`,
which grouped conditions by `group|name` and therefore pooled the 31 preregistered KTz non-chaotic regimes
(P × input pattern) into 11 groups (81 conditions instead of the preregistered 101). Under the preregistered
101 conditions, k2_mnlp_ann fails 16 conditions (not 12); every PASS/FAIL and the winner (k3_mnlp_growth_ann)
are unchanged. Corrected tabulation: `experiments/final_refinement/errata/phase9_decision_101.md`; details:
`experiments/final_refinement/ERRATA.md`.
