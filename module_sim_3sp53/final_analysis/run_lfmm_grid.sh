#!/bin/bash
## final_analysis/run_lfmm_grid.sh
##
## LFMM secondary-engine grid, all 1,400 combinations -- reuses the already-
## built parse/Stage-1/EMMAX outputs from run_full_grid.sh (missing input
## is a hard failure here too, matching that driver, since by this point
## every combo's ld_units.rds should already exist). "Already done" is
## checked against 04_lfmm/<combo>/lfmm.rds specifically, NOT truth_scores.rds
## (which already exists for every combo from the primary EMMAX grid).
set -euo pipefail
cd "$(dirname "$0")"
export PATH="/usr/local/bin:$PATH"

CELLS=(V0.5_c1 V0.5_c1.5 V0.5_c2 V1_c1 V1_c1.5 V2_c1 V2_c1.5)
TAGS=(nobgs bgs)
REPS=(1 2 3 4 5 6 7 8 9 10)
ENVS=(1 2 3 4 5 6 7 8 9 10)
CONCURRENCY="${1:-7}"

mkdir -p out_final_v1/logs
MANIFEST=out_final_v1/lfmm_grid_manifest.tsv
NEW_MANIFEST=$([ -f "$MANIFEST" ] && echo 0 || echo 1)
[ "$NEW_MANIFEST" = 1 ] && echo -e "tag\tcell\trep\tenv\tstatus\tseconds" > "$MANIFEST"

TMPDIR_STATUS="$(mktemp -d out_final_v1/.lfmm_grid_status_XXXXXX)"
trap 'rm -rf "$TMPDIR_STATUS"' EXIT

run_one() {
  local tag="$1" cell="$2" rep="$3" envn="$4"
  local combo_id="${tag}_${cell}_rep${rep}_env${envn}"
  local logfile="out_final_v1/logs/lfmm_${combo_id}.txt"
  local t0 t1 st

  if [ ! -f "out_final_v1/02_build_ld_units/${combo_id}/ld_units.rds" ]; then
    echo "FAILED: primary-grid input missing: out_final_v1/02_build_ld_units/${combo_id}/ld_units.rds" > "$logfile"
    echo -e "${tag}\t${cell}\t${rep}\t${envn}\tFAILED_MISSING_INPUT\t0" >> "${TMPDIR_STATUS}/$$_${combo_id}.tsv"
    return 0
  fi
  if [ -f "out_final_v1/04_lfmm/${combo_id}/lfmm.rds" ]; then
    echo -e "${tag}\t${cell}\t${rep}\t${envn}\tALREADY_DONE\t0" >> "${TMPDIR_STATUS}/$$_${combo_id}.tsv"
    return 0
  fi

  t0=$(date +%s)
  if Rscript R/run_combo.R "$tag" "$cell" "$rep" "$envn" 1 > "$logfile" 2>&1; then st=OK; else st=FAILED; fi
  t1=$(date +%s)
  echo -e "${tag}\t${cell}\t${rep}\t${envn}\t${st}\t$((t1 - t0))" >> "${TMPDIR_STATUS}/$$_${combo_id}.tsv"
}

for tag in "${TAGS[@]}"; do
  for cell in "${CELLS[@]}"; do
    for rep in "${REPS[@]}"; do
      for envn in "${ENVS[@]}"; do
        while [ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$CONCURRENCY" ]; do
          sleep 1
        done
        run_one "$tag" "$cell" "$rep" "$envn" &
      done
    done
  done
  wait
  cat "${TMPDIR_STATUS}"/*.tsv >> "$MANIFEST" 2>/dev/null || true
  rm -f "${TMPDIR_STATUS}"/*.tsv
done
wait

cat "${TMPDIR_STATUS}"/*.tsv >> "$MANIFEST" 2>/dev/null || true
echo "LFMM_GRID_DONE"
