#!/usr/bin/env bash
## module_sim_3sp53/run_structured_null_grid.sh
##
## Full-grid version of run_structured_null_sample.sh: runs R/12_structured_
## null.R across ALL (tag, cell, rep, env) combinations, not just the 14
## representative rep=1/env=1 combos. PK, 2026-09-08: "Run the full
## structured null sweep on all 1400 combos overnight."
##
## Resumable: R/12_structured_null.R has its own stage_stale()/receipt
## machinery (like every other stage script here) so a combo that already
## completed is skipped even without the pre-check below -- the pre-check
## just avoids paying the ~1-2s R startup cost per already-done combo on a
## restart. Safe to Ctrl-C and re-launch.
##
## Usage: ./run_structured_null_grid.sh [CONCURRENCY] [REPS_LIST] [ENVS_LIST]
##   CONCURRENCY   default 8 (matches run_grid.sh's main-grid default)
##   REPS_LIST     default "1 2 3 4 5 6 7 8 9 10"
##   ENVS_LIST     default "1 2 3 4 5 6 7 8 9 10"
set -euo pipefail
cd "$(dirname "$0")"
export PATH="/usr/local/bin:$PATH"

TAGS=(nobgs bgs)
CELLS=(V0.5_c1 V0.5_c1.5 V0.5_c2 V1_c1 V1_c1.5 V2_c1 V2_c1.5)
CONCURRENCY="${1:-8}"
REPS_ARG="${2:-1 2 3 4 5 6 7 8 9 10}"
REPS=(${REPS_ARG//,/ })
ENVS_ARG="${3:-1 2 3 4 5 6 7 8 9 10}"
ENVS=(${ENVS_ARG//,/ })

mkdir -p out/logs
LOG=out/structured_null_grid_timing.csv
NEW_LOG=$([ -f "$LOG" ] && echo 0 || echo 1)
[ "$NEW_LOG" = 1 ] && echo "tag,cell,rep,env,seconds,status" > "$LOG"

TMPDIR_TIMING="$(mktemp -d out/.structnull_timing_XXXXXX)"
trap 'rm -rf "$TMPDIR_TIMING"' EXIT

run_combo() {
  local tag="$1" cell="$2" rep="$3" env="$4"
  local bundle="out/02_bundle/bundle_${tag}_rep${rep}_${cell}_env${env}.rds"
  local outfile="out/12_structured_null/structnull_${tag}_rep${rep}_${cell}_env${env}.rds"
  if [ ! -f "$bundle" ]; then
    echo "[$(date +%H:%M:%S)] ${tag} ${cell} rep${rep} env${env}: SKIP (no bundle)"
    return 0
  fi
  if [ -f "$outfile" ]; then
    echo "[$(date +%H:%M:%S)] ${tag} ${cell} rep${rep} env${env}: SKIP (already done)"
    return 0
  fi
  export SIM_TAG="$tag" SIM_CELL="$cell" SIM_REP="$rep" SIM_ENV="$env"
  local logfile="out/logs/structnull_${tag}_${cell}_rep${rep}_env${env}.txt"
  local t0 t1 st
  t0=$(date +%s)
  if Rscript R/12_structured_null.R > "$logfile" 2>&1; then st=ok; else st=FAIL; fi
  t1=$(date +%s)
  echo "${tag},${cell},${rep},${env},$((t1 - t0)),${st}" >> "${TMPDIR_TIMING}/$$_${tag}_${cell}_${rep}_${env}.csv"
  echo "[$(date +%H:%M:%S)] ${tag} ${cell} rep${rep} env${env}: ${st} ($((t1 - t0))s)"
  [ "$st" = FAIL ] && { echo "  FAILED -- see $logfile" >&2; tail -20 "$logfile" >&2; }
}

for tag in "${TAGS[@]}"; do
  for cell in "${CELLS[@]}"; do
    for rep in "${REPS[@]}"; do
      for env in "${ENVS[@]}"; do
        while [ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$CONCURRENCY" ]; do
          sleep 1
        done
        run_combo "$tag" "$cell" "$rep" "$env" &
      done
    done
  done
done
wait

cat "${TMPDIR_TIMING}"/*.csv >> "$LOG" 2>/dev/null || true
echo "done"
