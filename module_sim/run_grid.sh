#!/usr/bin/env bash
## module_sim/run_grid.sh
##
## Run 01_parse_nemo -> 02_bundle -> 03_scan -> 04_score for every (cell, rep) combination
## in CELLS x REPS below, for one TAG, up to $CONCURRENCY combinations at
## once. Safe to run concurrently across machines too, as long as each
## machine's (tag,rep) combinations are disjoint from the other's --
## 00_config.R's receipt/GDS/untar/edge-list paths are all keyed per
## (tag,cell,rep) as of 2026-09-04 (see R/02_bundle.R, R/03_scan.R,
## R_parsing/01_parse_nemo.R), so combinations never alias each other's
## state; two machines both working through the same (tag,rep) combinations
## against the same NEMO_ROOT is the only thing this does NOT protect
## against. Splitting by TAG (one machine nobgs, the other bgs) is disjoint
## by construction and needs no REPS split on top of it.
##
## Usage:
##   ./run_grid.sh [TAG] [REPS_LIST] [CONCURRENCY]
##   TAG           "nobgs" or "bgs" (see 00_config.R's TAGS comment for why
##                 not "bgs5"). Defaults to nobgs.
##   REPS_LIST     space- or comma-separated rep numbers, e.g. "1 2 3 4 5" or
##                 "6,7,8,9,10". Defaults to 1-10 (all of REPS in 00_config.R).
##   CONCURRENCY   max combinations running at once. Defaults to 4.
##
## SIM_NEMO_ROOT, if already exported (e.g. to mini2's /Volumes/T9/Nemo_sim),
## is passed through untouched -- 00_config.R reads it directly.
##
## Examples:
##   ./run_grid.sh bgs "" 4                                       # local: bgs, all reps, 4-way
##   SIM_NEMO_ROOT=/Volumes/T9/Nemo_sim ./run_grid.sh nobgs "" 4  # mini2: nobgs, all reps, 4-way
##
## Logs per-stage wall time to module_sim/out/grid_timing.csv (merged at the
## end from per-worker temp files, so concurrent appends never interleave) and
## each invocation's full output to
## module_sim/out/logs/<tag>_<cell>_rep<rep>_<stage>.txt. ENV is fixed
## (SIM_ENV=1); TAG, CELLS and REPS vary.
set -euo pipefail
cd "$(dirname "$0")"
export SIM_ENV=1

TAG="${1:-nobgs}"
CELLS=(V0.5_c1 V0.5_c1.5 V0.5_c2 V1_c1 V1_c1.5 V2_c1 V2_c1.5)
REPS_ARG="${2:-1 2 3 4 5 6 7 8 9 10}"
REPS=(${REPS_ARG//,/ })
CONCURRENCY="${3:-4}"
export SIM_TAG="$TAG"

mkdir -p out/logs
LOG=out/grid_timing.csv
NEW_LOG=$([ -f "$LOG" ] && echo 0 || echo 1)
[ "$NEW_LOG" = 1 ] && echo "tag,cell,rep,stage,seconds,status" > "$LOG"

TMPDIR_TIMING="$(mktemp -d out/.grid_timing_XXXXXX)"
trap 'rm -rf "$TMPDIR_TIMING"' EXIT

run_combo() {
  local cell="$1" rep="$2"
  local combo_log="$TMPDIR_TIMING/${TAG}_${cell}_rep${rep}.csv"
  export SIM_CELL="$cell" SIM_REP="$rep"
  for pair in "R_parsing/01_parse_nemo.R:01_parse_nemo" "R/02_bundle.R:02_bundle" "R/03_scan.R:03_scan" "R/04_score.R:04_score"; do
    script="${pair%%:*}"; stage="${pair##*:}"
    logfile="out/logs/${TAG}_${cell}_rep${rep}_${stage}.txt"
    t0=$(date +%s)
    if Rscript "$script" > "$logfile" 2>&1; then status=ok; else status=FAIL; fi
    t1=$(date +%s)
    echo "${TAG},${cell},${rep},${stage},$((t1 - t0)),${status}" >> "$combo_log"
    echo "[$(date +%H:%M:%S)] ${TAG} ${cell} rep${rep} ${stage}: ${status} ($((t1 - t0))s)"
    if [ "$status" = "FAIL" ]; then
      echo "  FAILED -- see $logfile" >&2
      tail -15 "$logfile" >&2
      break
    fi
  done
}

## Polling throttle rather than `wait -n` -- macOS ships bash 3.2, which
## lacks it (added in 4.3), and this needs to run unmodified on whatever bash
## /usr/bin/env finds on either machine.
for cell in "${CELLS[@]}"; do
  for rep in "${REPS[@]}"; do
    while [ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$CONCURRENCY" ]; do
      sleep 2
    done
    run_combo "$cell" "$rep" &
  done
done
wait

cat "$TMPDIR_TIMING"/*.csv >> "$LOG" 2>/dev/null || true
echo "done -> $LOG"
