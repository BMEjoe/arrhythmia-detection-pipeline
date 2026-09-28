# Phase 3: LLE as a statistically meaningful, noise-robust chaos detector

Synthetic validation only. The findings are in
[`docs/PHASE3_LLE_VALIDATION.md`](../../docs/PHASE3_LLE_VALIDATION.md); the
frozen design is in [`PREREGISTRATION.md`](PREREGISTRATION.md); the resumable
work log is [`HANDOFF.md`](HANDOFF.md).

## Seed split

- **Development** seeds are the Phase 2E core / noise seeds, used for diagnosis
  and tuning.
- **Test** seeds are 1000–1099 for every condition. `run_phase3.py --phase test`
  refuses to run unless `PREREGISTRATION.md` is committed, unmodified, pushed,
  and lists the method.

## Reproduce

From the repository root, with `requirements.txt` installed (Python 3.13):

```bash
python -m experiments.phase3_lle.run_phase3 --phase dev  --methods all --workers 4
python -m experiments.phase3_lle.diagnose  && python -m experiments.phase3_lle.diagnose --summary
python -m experiments.phase3_lle.tune      && python -m experiments.phase3_lle.tune --summary
python -m experiments.phase3_lle.run_phase3 --phase test --methods all --workers 4
python -m experiments.phase3_lle.analysis --phase test
python -m experiments.phase3_lle.decide   --phase test
```

The runner can be resumed; finished task ids are skipped. BLAS is pinned to one
thread.

## Files

| File | Content |
|---|---|
| `config.py` | Systems, SNR levels, window lengths, LLE references, dev/test seeds |
| `estimators.py` | Baseline (existing pipeline LLE + AAFT test) and candidates C1–C4; vectorized divergence curve, IAAFT, 0–1 test |
| `run_phase3.py` | Checkpointed JSONL runner with the preregistration guard |
| `diagnose.py`, `diagnose_summary.py` | B2 diagnosis (delay, dimension, fit region, noise floor) |
| `tune.py` | B3 tuning grid on development seeds |
| `analysis.py` | Per method × window × condition tables (Wilson CIs, LLE bias) |
| `decide.py` | The preregistered primary decision rule |
| `replicate_check.py`, `groundrule_check.sh` | Checks run after every `final_pipeline.py` change |
| `results/dev/`, `results/test/` | Raw JSONL, manifests, tables |
| `plots/` | B2 divergence curves |

## Pipeline options added in Phase 3 (all default off)

- `keep_upo_on_short_lle_embedding` (A4)
- `upo_report_cao_diagnostics`, and `cao_method(..., return_diagnostics=True)` (A5)
- `lle_chaos_test` (B5): the preregistered winner C1
