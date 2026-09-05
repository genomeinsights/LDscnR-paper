#!/usr/bin/env bash
## module_sim/run_grid.sh
##
## Run 01_parse_nemo -> 02_bundle -> 03_scan -> 04_score for every
## (cell, rep, env) combination in CELLS x REPS x ENVS below, for one TAG, up
## to $CONCURRENCY combinations at once. Safe to run concurrently across
## machines too, as long as each machine's (tag,rep,env) combinations are
## disjoint from the other's -- 00_config.R's receipt/GDS/untar/edge-list
## paths are all keyed per (tag,cell,rep,env) as of 2026-09-05 (see
## R/02_bundle.R, R/03_scan.R, R/04_score.R, R_parsing/01_parse_nemo.R), so
## combinations never alias each other's state; two machines both working
## through the same (tag,rep,env) combinations against the same NEMO_ROOT is
## the only thing this does NOT protect against. Splitting by TAG (one
## machine nobgs, the other bgs) is disjoint by construction and needs no
## REPS/ENVS split on top of it.
##
## [!] REP IS NOT THE REPLICATE AXIS -- ENV IS (see 00_config.R's ENVS
## comment: each rep is a different recombination map, each env is an
## independent population draw under the SAME map). Both are looped here
## because PK asked for the full reps x envs cross -- R/05_pool.R averages
## over ENV within each REP, then separately summarizes across REPS.
##
## Usage:
##   ./run_grid.sh [TAG] [REPS_LIST] [ENVS_LIST] [CONCURRENCY]
##   TAG           "nobgs" or "bgs" (see 00_config.R's TAGS comment for why
##                 not "bgs5"). Defaults to nobgs.
##   REPS_LIST     space- or comma-separated rep numbers, e.g. "1 2 3 4 5" or
##                 "6,7,8,9,10". Defaults to 1-10 (all of REPS in 00_config.R).
##   ENVS_LIST     same format, for envs. Defaults to 1-10 (all of ENVS).
##   CONCURRENCY   max combinations running at once. Defaults to 4.
##
## SIM_NEMO_ROOT, if already exported (e.g. to mini2's /Volumes/T9/Nemo_sim),
## is passed through untouched -- 00_config.R reads it directly.
##
## Examples:
##   ./run_grid.sh bgs "" "" 4                                       # local: bgs, all reps x envs, 4-way
##   SIM_NEMO_ROOT=/Volumes/T9/Nemo_sim ./run_grid.sh nobgs "" "" 4  # mini2: nobgs, all reps x envs, 4-way
##
## Logs per-stage wall time to module_sim/out/grid_timing.csv (merged at the
## end from per-worker temp files, so concurrent appends never interleave) and
## each invocation's full output to
## module_sim/out/logs/<tag>_<cell>_rep<rep>_env<env>_<stage>.txt.
set -euo pipefail
cd "$(dirname "$0")"

TAG="${1:-nobgs}"
CELLS=(V0.5_c1 V0.5_c1.5 V0.5_c2 V1_c1 V1_c1.5 V2_c1 V2_c1.5)
REPS_ARG="${2:-1 2 3 4 5 6 7 8 9 10}"
REPS=(${REPS_ARG//,/ })
ENVS_ARG="${3:-1 2 3 4 5 6 7 8 9 10}"
ENVS=(${ENVS_ARG//,/ })
CONCURRENCY="${4:-4}"
export SIM_TAG="$TAG"

mkdir -p out/logs
LOG=out/grid_timing.csv
NEW_LOG=$([ -f "$LOG" ] && echo 0 || echo 1)
[ "$NEW_LOG" = 1 ] && echo "tag,cell,rep,env,stage,seconds,status" > "$LOG"

TMPDIR_TIMING="$(mktemp -d out/.grid_timing_XXXXXX)"
trap 'rm -rf "$TMPDIR_TIMING"' EXIT

## SKIP, not FAIL, a cell bgs5 does not have yet -- 00_config.R's raw_nemo_*
## comment: 3 of the 7 target cells (V0.5_c1.5, V1_c1, V2_c1.5) are still
## being simulated (PK, 2026-09-05). Checking the archive directly rather
## than a hardcoded "known missing" list, so dropping the real archives into
## bgs5/ later makes this script pick them up with no edit here.
NEMO_ROOT_CHECK="${SIM_NEMO_ROOT:-/Volumes/Nemo/Nemo_sim}"
archive_exists() {
  local tag="$1" cell="$2" rep="$3" env="$4"
  [ -f "${NEMO_ROOT_CHECK}/bgs5/adapt_${tag}_chr${rep}_${cell}_env${env}.tgz" ]
}

run_combo() {
  local cell="$1" rep="$2" env="$3"
  if ! archive_exists "$TAG" "$cell" "$rep" "$env"; then
    echo "[$(date +%H:%M:%S)] ${TAG} ${cell} rep${rep} env${env}: SKIP (no archive yet)"
    return 0
  fi
  local combo_log="$TMPDIR_TIMING/${TAG}_${cell}_rep${rep}_env${env}.csv"
  export SIM_CELL="$cell" SIM_REP="$rep" SIM_ENV="$env"
  for pair in "R_parsing/01_parse_nemo.R:01_parse_nemo" "R/02_bundle.R:02_bundle" "R/03_scan.R:03_scan" "R/04_score.R:04_score"; do
    script="${pair%%:*}"; stage="${pair##*:}"
    logfile="out/logs/${TAG}_${cell}_rep${rep}_env${env}_${stage}.txt"
    t0=$(date +%s)
    if Rscript "$script" > "$logfile" 2>&1; then status=ok; else status=FAIL; fi
    t1=$(date +%s)
    echo "${TAG},${cell},${rep},${env},${stage},$((t1 - t0)),${status}" >> "$combo_log"
    echo "[$(date +%H:%M:%S)] ${TAG} ${cell} rep${rep} env${env} ${stage}: ${status} ($((t1 - t0))s)"
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
    for env in "${ENVS[@]}"; do
      while [ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$CONCURRENCY" ]; do
        sleep 2
      done
      run_combo "$cell" "$rep" "$env" &
    done
  done
done
wait

cat "$TMPDIR_TIMING"/*.csv >> "$LOG" 2>/dev/null || true
echo "done -> $LOG"
