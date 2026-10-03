# Phase 10 methods and sources

Every method is implemented from a primary published source (or, where the primary is not
accessible, an open restatement, stated as such). Every implementation choice and deviation is
recorded here. Sections are filled in as each method is built.

## 1. Detectors (frozen; no retuning)
- K1 `k1_frozen512`: `experiments.phase7_mitbih.detector.evaluate` (Phase 6 Section 5 combined
  detector, `keep_upo_on_short_lle_embedding = True`, D2 detrend). Components LLE alone, UPO alone.
- K3: `final_pipeline.masked_growth_chaos_test(rr, labels, CFG)`.
