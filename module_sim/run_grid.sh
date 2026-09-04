#!/usr/bin/env bash
## module_sim/run_grid.sh
##
## Run 01_parse_nemo -> 02_bundle -> 03_scan for every (cell, rep) combination
## in CELLS x REPS below, SEQUENTIALLY -- not in parallel. 02_bundle.R's GDS
## file (module_sim/cache/sim.gds) is a single shared path, deleted and rebuilt
## fresh at the start of every run rather than kept per-combination; concurrent
## runs would collide on it mid-computation. The decay/clustering cache IS keyed
## per (tag,cell,rep) as of 2026-09-04, so repeated combinations across separate
## invocations of this script are cheap; a first full pass is not.
##
## Logs per-stage wall time to module_sim/out/grid_timing.csv and each
## invocation's full output to module_sim/out/logs/<cell>_rep<rep>_<stage>.txt.
## TAG and ENV are fixed (SIM_TAG=nobgs, SIM_ENV=1) -- neither has been widened
## in 00_config.R; only CELLS and REPS are looped here.
set -euo pipefail
cd "$(dirname "$0")"
export SIM_TAG=nobgs SIM_ENV=1

CELLS=(V0.5_c1 V0.5_c1.5 V0.5_c2 V1_c1 V1_c1.5 V2_c1 V2_c1.5)
REPS=(1 2 3)

mkdir -p out/logs
LOG=out/grid_timing.csv
NEW_LOG=$([ -f "$LOG" ] && echo 0 || echo 1)
[ "$NEW_LOG" = 1 ] && echo "cell,rep,stage,seconds,status" > "$LOG"

for cell in "${CELLS[@]}"; do
  for rep in "${REPS[@]}"; do
    export SIM_CELL="$cell" SIM_REP="$rep"
    for pair in "R_parsing/01_parse_nemo.R:01_parse_nemo" "R/02_bundle.R:02_bundle" "R/03_scan.R:03_scan"; do
      script="${pair%%:*}"; stage="${pair##*:}"
      logfile="out/logs/${cell}_rep${rep}_${stage}.txt"
      t0=$(date +%s)
      if Rscript "$script" > "$logfile" 2>&1; then status=ok; else status=FAIL; fi
      t1=$(date +%s)
      echo "${cell},${rep},${stage},$((t1 - t0)),${status}" >> "$LOG"
      echo "[$(date +%H:%M:%S)] ${cell} rep${rep} ${stage}: ${status} ($((t1 - t0))s)"
      if [ "$status" = "FAIL" ]; then
        echo "  FAILED -- see $logfile" >&2
        tail -15 "$logfile" >&2
      fi
    done
  done
done
echo "done -> $LOG"
