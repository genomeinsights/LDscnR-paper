#!/bin/bash
## final_analysis/run_rerun_sequence.sh
##
## Sequential rerun of everything downstream of the 2026-09-18 review's
## Stage-2 ordering-unification fix (region_scoring_version 3 -> 4) and the
## bootstrap-pairing / opportunity-control / paired-map-null redesigns that
## went in alongside it. Order matters: 06 depends on 05's fresh
## truth_scores.rds; 11/12 depend on 05's fresh per-combo outputs; 13/14
## depend on 11's freshly regenerated region_details.rds; figures depend on
## 13/14's final TSVs. Stops on the first failure (set -e) rather than
## silently continuing on stale/partial inputs.
set -euo pipefail
cd "$(dirname "$0")"
export PATH="/usr/local/bin:$PATH"
mkdir -p out_final_v1/logs

run_stage() {
  local label="$1"; shift
  local logfile="out_final_v1/logs/rerun_${label}.log"
  echo "=== [$(date '+%H:%M:%S')] starting ${label} ==="
  t0=$(date +%s)
  if "$@" > "$logfile" 2>&1; then
    t1=$(date +%s)
    echo "=== [$(date '+%H:%M:%S')] ${label} OK ($((t1 - t0))s) -- see ${logfile} ==="
  else
    t1=$(date +%s)
    echo "=== [$(date '+%H:%M:%S')] ${label} FAILED ($((t1 - t0))s) -- see ${logfile} ==="
    tail -50 "$logfile"
    exit 1
  fi
}

## step 0: 05_score_truth full-grid rescore -- run separately/already done by
## the time this script is launched (rerun_score_truth_grid.R); NOT repeated
## here to avoid re-running an already-completed multi-hour... actually fast
## ... step from inside this script.

run_stage "06_summarise" Rscript R/06_summarise.R
run_stage "11_stage2_region_details" Rscript R/11_stage2_region_details.R
run_stage "12_floor_decomposition" Rscript R/12_floor_decomposition.R
run_stage "13_region_null_calibration_balanced" Rscript R/13_region_null_calibration.R balanced 200
run_stage "14_summarise_additional_analyses" Rscript R/14_summarise_additional_analyses.R
run_stage "figure_floor_decomposition" Rscript R_figures/figure_simulation_floor_decomposition.R
run_stage "figure_floor_decomposition_by_cell" Rscript R_figures/figure_simulation_floor_decomposition_by_cell.R
run_stage "figure_null_truth_calibration" Rscript R_figures/figure_simulation_null_truth_calibration.R
run_stage "figureS_stage2_size_truth" Rscript R_figures/figureS_simulation_stage2_size_truth.R

echo "RERUN_SEQUENCE_DONE"
