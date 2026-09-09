#!/bin/bash
## final_analysis/run_gate3_grid.sh
##
## Phase 3 gate 3: all 7 cells x both tags, rep=1, env=1 (14 combinations),
## through the full pipeline (parse -> Stage 1/GRM -> EMMAX -> truth score).
set -euo pipefail
cd "$(dirname "$0")"
export PATH="/usr/local/bin:$PATH"

CELLS=(V0.5_c1 V0.5_c1.5 V0.5_c2 V1_c1 V1_c1.5 V2_c1 V2_c1.5)
TAGS=(nobgs bgs)
REP=1
ENVN=1
CONCURRENCY="${1:-7}"

mkdir -p out_final_v1/logs

for tag in "${TAGS[@]}"; do
  for cell in "${CELLS[@]}"; do
    while [ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$CONCURRENCY" ]; do
      sleep 2
    done
    logfile="out_final_v1/logs/${tag}_${cell}_rep${REP}_env${ENVN}.txt"
    ( Rscript R/run_combo.R "$tag" "$cell" "$REP" "$ENVN" > "$logfile" 2>&1 && \
      echo "[$(date +%H:%M:%S)] ${tag}/${cell}: OK" || \
      echo "[$(date +%H:%M:%S)] ${tag}/${cell}: FAILED -- see $logfile" ) &
  done
done
wait
echo "GATE3_GRID_DONE"
