# Phase 7 handoff (living document; update after every step)

Read this first when resuming. Branch `claude/amazing-cray-9jbx6p`, from main at
`1b4b7a0` (Phases 2E-6). Never merge into main, never open a pull request.
Commit and push after every step.

## Task (user, Phase 7) in brief
Evaluate the FROZEN combined detector on MIT-BIH (48 records, MITDB_RECORDS).
DETECTOR = combined_chaos_config(replace(CFG, keep_upo_on_short_lle_embedding=True), detrend=True);
decision = combined_chaos_detected(analyze_segment(rr, DETECTOR)).
- Q1 (PRIMARY): AND rate abnormal - normal windows (raw RR), subject-cluster
  bootstrap 95 % CI (10,000); supported iff CI excludes 0. Also sens/spec/PPV/NPV/BA.
- Q2: abnormal AND rate vs Phase 6 ectopy-only FP (max 3.7 %, upper bound 6.4 %);
  "exceeds" iff lower 95 % bound > 6.4 %.
- Q3: LOSO logistic regression (fixed L2 C, standardize in fold): M0 HRV
  (SDNN, RMSSD, pNN50), M1 LLE, M2 UPO, M3 LLE+UPO, M4 all; OOF AUC +
  subject-cluster bootstrap; paired diffs M3-M1, M3-M2, M4-M0; DeLong secondary.
- Secondary: edited NN (NN fraction >= 0.80), per-record table, beat type,
  AF records. Sensitivity: annotation beat times instead of detected peaks.
- Steps: 0 data + manifest; 1 QC (NO detector output on MIT-BIH); 2
  PREREGISTRATION.md pushed alone before any detector run (Phase 3 guard);
  3 run; 4 analysis; 5 docs/PHASE7_MITBIH_RESULTS.md.
- Subject not record for all grouping (201 and 202 = same subject).
- Anything decided after seeing MIT-BIH detector output = EXPLORATORY.

## Ground rules (as Phases 3-6)
- final_pipeline.py: opt-in options only, defaults unchanged; after any change
  `bash experiments/phase3_lle/groundrule_check.sh`. Not expected to change.
- Environment: /root/venv313 (Python 3.13.14, numpy 2.1.3, scipy 1.18.1,
  sklearn 1.9.1, wfdb 4.3.1, pandas 3.0.6); recreate with
  `uv venv -p python3.13 /root/venv313 && uv pip install -p /root/venv313/bin/python -r requirements.txt wfdb`.
  OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1. 4 cores.
- physionet.org IS reachable in this session (HTTP 200).
- Raw data in mitdb_data/ (gitignored). Re-download with the command in
  DATA_MANIFEST.json and compare checksums.

## Status
| Step | State | Commit |
|---|---|---|
| 0 data: 48 records x (hea, dat, atr), mitdb 1.0.0, 144/144 SHA-256 match PhysioNet SHA256SUMS.txt; DATA_MANIFEST.json; mitdb_data/ gitignored | done | (this commit) |
| 1 QC (qc.py, qc_tables.py, data.py, tests/test_phase7_mitbih_helpers.py 43 tests): Se/PPV 75 ms 0.9787/0.9734, 150 ms 0.9965/0.9911; late PVC fiducial (~100 ms) -> decision: label transfer at 150 ms (EC57), beat annotations only; 305 windows (211 N / 94 A), 43 subjects, 19 with abnormal; 75 ms arm 224 (187/37); annotation-time 368 (248/120); NN>=0.80: 211 N / 7 A; Q stays abnormal; paced records 0 windows | done | (this commit) |

## Notes
- NO detector output on MIT-BIH exists yet. Next: PREREGISTRATION.md (commit and
  push alone), then runner with the Phase 3 guard.
- Peak cache: mitdb_data/_phase7_peaks/ (gitignored); committed copy
  results/qc/detected_peaks.npz.
| 2a code before preregistration: detector.py (frozen DETECTOR, lle_z, upo_score), run_phase7.py (guarded; --phase synthetic code check), analysis.py (Q1-Q3, secondary, sensitivity); tests test_phase7_detector.py (5), test_phase7_analysis.py (7); synthetic check results/synthetic_check (72 Phase 5 windows, 53 s) | done | (this commit) |
| 2 PREREGISTRATION.md committed and pushed alone | done | fb212da |

- Step 3 command (resumable): `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /root/venv313/bin/python -m experiments.phase7_mitbih.run_phase7 --phase test --workers 4 > /tmp/claude-0/p7_run.log 2>&1`
  then `python -m experiments.phase7_mitbih.analysis`. Do NOT edit data.py, detector.py,
  run_phase7.py, analysis.py or final_pipeline.py (amendments only).
