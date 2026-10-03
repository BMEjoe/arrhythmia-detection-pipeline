#!/usr/bin/env bash
# Single entry point for reproducing this repository's verification, tables and figures.
#
#   bash reproduce.sh fast            regenerate every verification, table and figure from the STORED results
#                                     (no data needed; ~15 min on 4 CPUs) and check them against the committed files
#   bash reproduce.sh tests           pytest, default and with numpy AVX-512 dispatch disabled (~5 min)
#   bash reproduce.sh data [db ...]   download the PhysioNet databases, SHA-256 verified (see download_data.py)
#   bash reproduce.sh slow <phase>    full rerun of one phase from scratch in a separate clone (hours; needs data)
#   bash reproduce.sh slow list       phases, commands and expected runtimes
#   bash reproduce.sh slow <phase> --dry-run
#
# Environment: Python 3.13 with requirements-lock.txt (README.md, "Environment"). Set PY to the interpreter
# (default: python). Every step runs single-threaded BLAS, as in the original runs.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
PY="${PY:-python}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
WORKERS="${WORKERS:-4}"

step() {  # step "<label>" <python -m module args...>
  local label="$1"; shift
  local t0=$SECONDS
  printf '%-62s ' "$label"
  if "$PY" -m "$@" >> "${LOG:-/dev/null}" 2>&1; then printf 'ok (%ds)\n' $((SECONDS - t0))
  else printf 'FAILED (%ds)\n' $((SECONDS - t0)); FAILED=1; fi
}

fast() {
  FAILED=0
  echo "== Final refinement (independent verification, errata, preregistration integrity, figures)"
  step "verify_all (VERIFICATION.md)"                      experiments.final_refinement.verify_all
  step "errata_phase9_decision (erratum E1 tabulation)"     experiments.final_refinement.errata_phase9_decision
  if [ "$(git rev-parse --is-shallow-repository)" = "true" ]; then
    echo "prereg_integrity: SKIPPED (shallow clone; run 'git fetch --unshallow' first)"
  else
    step "prereg_integrity (PREREG_INTEGRITY.md)"           experiments.final_refinement.prereg_integrity
  fi
  step "make_figures (figures/, tables/, FIGURE_NOTES.md)"  experiments.final_refinement.make_figures
  echo "== Phase tables and decisions from the stored raw results"
  step "Phase 2E analysis (tables, plots)"                  experiments.phase2e.analysis
  step "Phase 3 analysis"                                   experiments.phase3_lle.analysis --phase test
  step "Phase 3 decision"                                   experiments.phase3_lle.decide
  step "Phase 4 analysis"                                   experiments.phase4_upo.analysis --phase test --file test
  step "Phase 4 decision"                                   experiments.phase4_upo.decide --phase test --file test
  step "Phase 5 analysis"                                   experiments.phase5_rr.analysis --phase test --file test
  step "Phase 5 decision"                                   experiments.phase5_rr.decide --phase test --file test
  step "Phase 5 report extras"                              experiments.phase5_rr.report_extras
  step "Phase 6 analysis"                                   experiments.phase6_robust.analysis --phase test --file test
  step "Phase 6 decision"                                   experiments.phase6_robust.decide --phase test --file test
  step "Phase 7 preregistered analysis"                     experiments.phase7_mitbih.analysis
  step "Phase 7 exploratory analysis"                       experiments.phase7_mitbih.exploratory
  step "Phase 7 spike-in summary"                           experiments.phase7_mitbih.spike_in --summary
  step "Phase 8 decision"                                   experiments.phase8_cardiac.decide
  step "Phase 8 close-out tables"                           experiments.phase9_waveform.phase8_closeout.tables
  step "Phase 9 decision (frozen decide9.py)"               experiments.phase9_waveform.decide9
  cat experiments/final_refinement/errata/decision9_erratum_note.md \
      >> experiments/phase9_waveform/results/test/tables/decision9.md   # re-append the erratum E1 note
  step "Phase 9 Part F analysis"                            experiments.phase9_waveform.partf --analyze
  step "Phase 10 analysis, confirmatory (~5 min)"           experiments.phase10_final.analysis10 --phase conf
  step "Phase 10 analysis, sensitivity (~5 min)"            experiments.phase10_final.analysis10 --phase conf --sensitivity-flagged
  step "Phase 10 analysis, development (~2 min)"            experiments.phase10_final.analysis10 --phase dev
  step "Phase 10 tables, confirmatory"                      experiments.phase10_final.tables10 --phase conf
  step "Phase 10 tables, development"                       experiments.phase10_final.tables10 --phase dev
  step "Phase 10 figures, confirmatory"                     experiments.phase10_final.figures10 --phase conf
  step "Phase 10 figures, development"                      experiments.phase10_final.figures10 --phase dev
  echo "== Comparison with the committed files"
  # Phase 2E was produced in an unrecorded pandas / matplotlib environment: its regenerated tables differ only in
  # Markdown column alignment and the last digit of one float column, and its PNGs in bytes (KNOWN_ISSUES.md 7).
  local expected='^experiments/phase2e/(plots/.*\.png|results/tables/(table_E_map_modes\.md|table_H_lle\.csv|table_coverage_vs_count\.csv))$'
  local changed; changed="$(git diff --name-only)"
  local unexpected; unexpected="$(printf '%s\n' "$changed" | grep -Ev "$expected" | grep -v '^$' || true)"
  local known; known="$(printf '%s\n' "$changed" | grep -Ec "$expected" || true)"
  echo "known environment-only differences (Phase 2E formatting / PNG bytes): $known file(s)"
  if [ -n "$unexpected" ]; then
    echo "UNEXPECTED differences from the committed files:"; printf '  %s\n' $unexpected; FAILED=1
  else
    echo "every other regenerated file is byte-identical to the committed version"
  fi
  echo "restore the Phase 2E files with: git checkout -- experiments/phase2e"
  return $FAILED
}

tests() {
  "$PY" -m pytest -q -p no:cacheprovider tests | tail -2 || true
  echo "-- with AVX-512 dispatch disabled (KNOWN_ISSUES.md 1: one test depends on the CPU floating-point path)"
  NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL AVX512_SPR" \
    "$PY" -m pytest -q -p no:cacheprovider tests | tail -1 || true
}

# ------------------------------------------------------------------------------------------------ slow
# phase | raw result files removed in the clone before the rerun | commands | expected wall time (4 workers)
slow_plan() {
  case "$1" in
    2e) RESULTS="experiments/phase2e/results/*.jsonl"
        CMDS=("experiments.phase2e.run_phase2e --experiments all --workers $WORKERS --no-resume"
              "experiments.phase2e.analysis")
        TIME="~31,000 task-seconds (Phase 2E report 1); about 2-3 h on 4 workers";;
    3)  RESULTS="experiments/phase3_lle/results/dev/*.jsonl experiments/phase3_lle/results/test/*.jsonl"
        CMDS=("experiments.phase3_lle.diagnose" "experiments.phase3_lle.tune"
              "experiments.phase3_lle.run_phase3 --phase dev --methods all --workers $WORKERS"
              "experiments.phase3_lle.run_phase3 --phase test --methods all --workers $WORKERS"
              "experiments.phase3_lle.analysis --phase test" "experiments.phase3_lle.decide")
        TIME="test run 14,854 s (4.1 h) plus development runs (Phase 3 report 5.4)";;
    4)  RESULTS="experiments/phase4_upo/results/dev/*.jsonl experiments/phase4_upo/results/test/*.jsonl experiments/phase4_upo/results/m_sensitivity.jsonl"
        CMDS=("experiments.phase4_upo.run_phase4 --phase dev --runs cao:7:50,2:7:50,2:15:50,3:7:50,cao:15:50 --tag explore --windows 256 --workers $WORKERS"
              "experiments.phase4_upo.run_phase4 --phase dev --runs cao:7:50,2:15:50,2:7:50 --tag dev512 --windows 512 --workers $WORKERS"
              "experiments.phase4_upo.run_phase4 --phase test --methods all --workers $WORKERS"
              "experiments.phase4_upo.analysis --phase test --file test" "experiments.phase4_upo.decide --phase test --file test"
              "experiments.phase4_upo.m_sensitivity")
        TIME="test run about 4.9 h (Phase 4 report 3, budget) plus development runs";;
    5)  RESULTS="experiments/phase5_rr/results/dev/*.jsonl experiments/phase5_rr/results/test/*.jsonl experiments/phase5_rr/results/phase4_m_check.jsonl"
        CMDS=("experiments.phase5_rr.equivalence_check" "experiments.phase5_rr.realism"
              "experiments.phase5_rr.run_phase5 --phase dev --seeds 0-9 --tag runtime --workers $WORKERS"
              "experiments.phase5_rr.run_phase5 --phase test --workers $WORKERS"
              "experiments.phase5_rr.analysis --phase test --file test" "experiments.phase5_rr.decide --phase test --file test"
              "experiments.phase5_rr.report_extras" "experiments.phase5_rr.phase4_m_check")
        TIME="test run 11,421 s (3.2 h; Phase 5 report 4)";;
    6)  RESULTS="experiments/phase6_robust/results/test/*.jsonl experiments/phase6_robust/results/partA/*.jsonl"
        CMDS=("experiments.phase6_robust.partA_crash_check"
              "experiments.phase6_robust.run_phase6 --phase test --workers $WORKERS"
              "experiments.phase6_robust.analysis --phase test --file test" "experiments.phase6_robust.decide --phase test --file test")
        TIME="test run about 3.5 h (Phase 6 report 4.3); development tuning commands in the Phase 6 report 7";;
    7)  RESULTS="experiments/phase7_mitbih/results/run/*.jsonl experiments/phase7_mitbih/results/spike_in/spike_in.jsonl"
        NEEDS="mitdb"
        CMDS=("experiments.phase7_mitbih.qc" "experiments.phase7_mitbih.qc_tables"
              "experiments.phase7_mitbih.run_phase7 --phase test --workers $WORKERS"
              "experiments.phase7_mitbih.analysis" "experiments.phase7_mitbih.exploratory" "experiments.phase7_mitbih.audit"
              "experiments.phase7_mitbih.spike_in --workers $WORKERS" "experiments.phase7_mitbih.spike_in --summary")
        TIME="about 1 h (891 detector windows; spike-in 1,331 s; Phase 7 report 12.2)";;
    8)  RESULTS="experiments/phase8_cardiac/results/test/*.jsonl"
        CMDS=("experiments.phase8_cardiac.run_test --workers $WORKERS" "experiments.phase8_cardiac.decide"
              "experiments.phase9_waveform.phase8_closeout.tables")
        TIME="about 8 h (TEST checkpoints 15:23-23:20, Phase 8 report 7); ground-truth scan and calibration in its report 10";;
    9)  RESULTS="experiments/phase9_waveform/results/test/*.jsonl experiments/phase9_waveform/results/partf/*.jsonl"
        NEEDS="nstdb nsrdb chfdb mitdb"
        CMDS=("experiments.phase9_waveform.run_test9 --workers $WORKERS" "experiments.phase9_waveform.decide9"
              "experiments.final_refinement.errata_phase9_decision"
              "experiments.phase9_waveform.partf --run --workers $WORKERS" "experiments.phase9_waveform.partf --analyze"
              "experiments.phase9_waveform.partf --mitbih --workers $WORKERS")
        TIME="TEST 30,052 s (8.3 h; Phase 9 report 6) plus Part F (about 1 h)";;
    10) RESULTS="experiments/phase10_final/results/dev/*.jsonl experiments/phase10_final/results/conf/*.jsonl"
        NEEDS="mitdb nsrdb chfdb nsr2db chf2db"
        CMDS=("experiments.phase10_final.verify_titration"
              "experiments.phase10_final.run10 --phase dev --part all --workers $WORKERS"
              "experiments.phase10_final.overlap10 --amended"
              "experiments.phase10_final.run10 --phase conf --part all --workers $WORKERS"
              "experiments.phase10_final.analysis10 --phase conf" "experiments.phase10_final.analysis10 --phase conf --sensitivity-flagged"
              "experiments.phase10_final.analysis10 --phase dev" "experiments.phase10_final.tables10 --phase conf"
              "experiments.phase10_final.figures10 --phase conf")
        TIME="confirmatory about 15 h projected (PREREGISTRATION.md 8) plus development about 10 h (the original development spike-in and robustness runs were partial)";;
    *) echo "unknown phase '$1' (2e, 3, 4, 5, 6, 7, 8, 9, 10)"; exit 2;;
  esac
}

slow() {
  local phase="${1:-list}"; local dry="${2:-}"
  if [ "$phase" = "list" ]; then
    for p in 2e 3 4 5 6 7 8 9 10; do
      NEEDS=""; slow_plan "$p"; printf 'Phase %-3s %s\n' "$p" "$TIME"
      [ -n "$NEEDS" ] && printf '          data: %s\n' "$NEEDS"
      for c in "${CMDS[@]}"; do printf '          python -m %s\n' "$c"; done
    done
    echo "Total: roughly 55-65 h on 4 workers. Synthetic phases (2E-6, 8) need no data."
    return 0
  fi
  NEEDS=""; slow_plan "$phase"
  local branch; branch="$(git rev-parse --abbrev-ref HEAD)"
  local dir="${RERUN_DIR:-$ROOT/../rerun-phase$phase}"
  echo "Phase $phase full rerun in a separate clone: $dir (branch $branch). Expected: $TIME"
  echo "The clone's runners pass their preregistration guard (the file is committed and identical on 'origin',"
  echo "the local repository). Raw results removed in the clone before the rerun: $RESULTS"
  [ -n "$NEEDS" ] && echo "Data needed: $NEEDS (downloaded into $ROOT/physionet_data and linked into the clone)"
  for c in "${CMDS[@]}"; do echo "  python -m $c"; done
  [ "$dry" = "--dry-run" ] && return 0
  [ -e "$dir" ] && { echo "$dir exists; remove it or set RERUN_DIR"; exit 2; }
  git clone -q --branch "$branch" "$ROOT" "$dir"
  if [ -n "$NEEDS" ]; then
    "$PY" -m experiments.final_refinement.download_data --db $NEEDS
    ln -s "$ROOT/physionet_data" "$dir/physionet_data"
    (cd "$dir" && "$PY" -m experiments.final_refinement.download_data --db $NEEDS)   # links + verification
  fi
  (cd "$dir" && rm -f $RESULTS)
  echo "log: $dir/rerun_phase$phase.log"
  for c in "${CMDS[@]}"; do (cd "$dir" && LOG="$dir/rerun_phase$phase.log" step "python -m $c" $c); done
  echo "Rerun finished. Differences from the committed results (full reruns on another machine reproduce counts"
  echo "approximately, not bit for bit; KNOWN_ISSUES.md 2):"
  (cd "$dir" && git status --short | head -50)
}

case "${1:-}" in
  fast)  fast ;;
  tests) tests ;;
  data)  shift; if [ $# -gt 0 ]; then "$PY" -m experiments.final_refinement.download_data --db "$@";
         else "$PY" -m experiments.final_refinement.download_data; fi ;;
  slow)  shift; slow "$@" ;;
  *) sed -n '2,12p' "$0"; exit 2 ;;
esac
