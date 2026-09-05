#!/usr/bin/env bash
## module_sim/run_bgsrecomb_grid.sh
##
## Run R/08_bgs_recomb.R across the full TAGS x CELLS x REPS x ENVS grid.
## Only needs R/02_bundle.R's output (already complete) plus the raw
## rec_map<rep>.rds files -- runs entirely locally.
set -euo pipefail
cd "$(dirname "$0")"

TAGS=(nobgs bgs)
CELLS=(V0.5_c1 V0.5_c1.5 V0.5_c2 V1_c1 V1_c1.5 V2_c1 V2_c1.5)
REPS=(1 2 3 4 5 6 7 8 9 10)
ENVS=(1 2 3 4 5 6 7 8 9 10)
CONCURRENCY="${1:-4}"

mkdir -p out/logs

## SKIP, not FAIL, a (tag,cell,rep,env) whose bundle isn't built yet -- e.g. a
## cell bgs5 doesn't have data for yet. [!] FIXED 2026-09-06: this script was
## missing the same bundle-existence check run_popgen_grid.sh and
## run_bgswindows_grid.sh already have, so the first real bgs5 run (4/7
## cells) FAILED (not skipped) 600/1400 combos -- one Rscript error per
## missing cell -- rather than cleanly skipping them.
run_combo() {
  local tag="$1" cell="$2" rep="$3" env="$4"
  local bundle="out/02_bundle/bundle_${tag}_rep${rep}_${cell}_env${env}.rds"
  if [ ! -f "$bundle" ]; then
    echo "[$(date +%H:%M:%S)] ${tag} ${cell} rep${rep} env${env}: SKIP (no bundle yet)"
    return 0
  fi
  export SIM_TAG="$tag" SIM_CELL="$cell" SIM_REP="$rep" SIM_ENV="$env"
  logfile="out/logs/bgsrecomb_${tag}_${cell}_rep${rep}_env${env}.txt"
  if Rscript R/08_bgs_recomb.R > "$logfile" 2>&1; then st=ok; else st=FAIL; fi
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
