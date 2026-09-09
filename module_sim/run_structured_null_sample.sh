#!/usr/bin/env bash
## module_sim/run_structured_null_sample.sh
##
## Run R/12_structured_null.R across a REPRESENTATIVE 14-combo subset -- all
## 7 cells x 2 tags, rep=1 env=1 held fixed so the comparison across cells
## isn't confounded with which particular rep/env got drawn. Matches
## module_3sp/module_9sp's own scale (one calibration check per dataset),
## not an exhaustive resweep across all 1400 combos.
##
## Usage: ./run_structured_null_sample.sh [CONCURRENCY]
set -euo pipefail
cd "$(dirname "$0")"

TAGS=(nobgs bgs)
CELLS=(V0.5_c1 V0.5_c1.5 V0.5_c2 V1_c1 V1_c1.5 V2_c1 V2_c1.5)
REP=1
ENV=1
CONCURRENCY="${1:-4}"

mkdir -p out/logs

run_combo() {
  local tag="$1" cell="$2"
  local bundle="out/02_bundle/bundle_${tag}_rep${REP}_${cell}_env${ENV}.rds"
  if [ ! -f "$bundle" ]; then
    echo "[$(date +%H:%M:%S)] ${tag} ${cell}: SKIP (no bundle)"
    return 0
  fi
  export SIM_TAG="$tag" SIM_CELL="$cell" SIM_REP="$REP" SIM_ENV="$ENV"
  logfile="out/logs/structnull_${tag}_${cell}.txt"
  t0=$(date +%s)
  if Rscript R/12_structured_null.R > "$logfile" 2>&1; then st=ok; else st=FAIL; fi
  t1=$(date +%s)
  echo "[$(date +%H:%M:%S)] ${tag} ${cell}: ${st} ($((t1 - t0))s)"
  [ "$st" = FAIL ] && { echo "  FAILED -- see $logfile" >&2; tail -20 "$logfile" >&2; }
}

for tag in "${TAGS[@]}"; do
  for cell in "${CELLS[@]}"; do
    while [ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$CONCURRENCY" ]; do
      sleep 1
    done
    run_combo "$tag" "$cell" &
  done
done
wait
echo "done"
