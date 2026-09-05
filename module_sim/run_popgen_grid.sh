#!/usr/bin/env bash
## module_sim/run_popgen_grid.sh
##
## Run R/06_popgen_summary.R across the full TAGS x CELLS x REPS x ENVS grid.
## Independent of run_grid.sh -- only needs R/02_bundle.R's output, which is
## already complete for the whole grid (both tags), so this runs entirely
## locally regardless of which machine originally built which tag.
set -euo pipefail
cd "$(dirname "$0")"

TAGS=(nobgs bgs)
CELLS=(V0.5_c1 V0.5_c1.5 V0.5_c2 V1_c1 V1_c1.5 V2_c1 V2_c1.5)
REPS=(1 2 3 4 5 6 7 8 9 10)
ENVS=(1 2 3 4 5 6 7 8 9 10)
CONCURRENCY="${1:-4}"

mkdir -p out/logs

run_combo() {
  local tag="$1" cell="$2" rep="$3" env="$4"
  export SIM_TAG="$tag" SIM_CELL="$cell" SIM_REP="$rep" SIM_ENV="$env"
  logfile="out/logs/popgen_${tag}_${cell}_rep${rep}_env${env}.txt"
  if Rscript R/06_popgen_summary.R > "$logfile" 2>&1; then st=ok; else st=FAIL; fi
  echo "[$(date +%H:%M:%S)] ${tag} ${cell} rep${rep} env${env}: ${st}"
  [ "$st" = FAIL ] && { echo "  FAILED -- see $logfile" >&2; tail -15 "$logfile" >&2; }
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
echo "done"
