#!/bin/bash
## final_analysis/run_full_grid.sh
##
## Phase 3 gate 4: the full manifest, 7 cells x 2 tags x 10 reps x 10 envs =
## 1,400 combinations. Missing raw inputs are FAILURES, not silent skips
## (checked explicitly below, before run_combo.R is even invoked) --
## instructions: "For the final run, missing inputs are failures, not
## silent skips." Every combo's status and wall time is recorded, merged
## from per-worker temp files at the end so concurrent appends never
## interleave (same pattern module_sim/run_grid.sh already uses for this).
set -euo pipefail
cd "$(dirname "$0")"
export PATH="/usr/local/bin:$PATH"

CELLS=(V0.5_c1 V0.5_c1.5 V0.5_c2 V1_c1 V1_c1.5 V2_c1 V2_c1.5)
TAGS=(nobgs bgs)
REPS=(1 2 3 4 5 6 7 8 9 10)
ENVS=(1 2 3 4 5 6 7 8 9 10)
CONCURRENCY="${1:-7}"

RAW_ROOT="/Volumes/Large_storage/prod_out"
mkdir -p out_final_v1/logs
MANIFEST=out_final_v1/full_grid_manifest.tsv
NEW_MANIFEST=$([ -f "$MANIFEST" ] && echo 0 || echo 1)
[ "$NEW_MANIFEST" = 1 ] && echo -e "tag\tcell\trep\tenv\tstatus\tseconds" > "$MANIFEST"

TMPDIR_STATUS="$(mktemp -d out_final_v1/.full_grid_status_XXXXXX)"
trap 'rm -rf "$TMPDIR_STATUS"' EXIT

run_one() {
  local tag="$1" cell="$2" rep="$3" envn="$4"
  local combo_id="${tag}_${cell}_rep${rep}_env${envn}"
  local logfile="out_final_v1/logs/${combo_id}.txt"
  local t0 t1 st

  ## missing-input check BEFORE attempting the run -- a failure, never a skip
  local raw_dir="${RAW_ROOT}/adapt_${tag}_chr${rep}_${cell}_env${envn}/GENO"
  if [ ! -d "$raw_dir" ]; then
    echo "FAILED: raw input missing: $raw_dir" > "$logfile"
    echo -e "${tag}\t${cell}\t${rep}\t${envn}\tFAILED_MISSING_INPUT\t0" >> "${TMPDIR_STATUS}/$$_${combo_id}.tsv"
    return 0
  fi

  ## already done (existing truth_scores.rds) -- skip re-running, still recorded
  if [ -f "out_final_v1/05_score_truth/${combo_id}/truth_scores.rds" ]; then
    echo -e "${tag}\t${cell}\t${rep}\t${envn}\tALREADY_DONE\t0" >> "${TMPDIR_STATUS}/$$_${combo_id}.tsv"
    return 0
  fi

  t0=$(date +%s)
  if Rscript R/run_combo.R "$tag" "$cell" "$rep" "$envn" > "$logfile" 2>&1; then st=OK; else st=FAILED; fi
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
  wait   ## drain between tags -- keeps the status-merge below from growing unbounded mid-run
  cat "${TMPDIR_STATUS}"/*.tsv >> "$MANIFEST" 2>/dev/null || true
  rm -f "${TMPDIR_STATUS}"/*.tsv
done
wait

cat "${TMPDIR_STATUS}"/*.tsv >> "$MANIFEST" 2>/dev/null || true
echo "FULL_GRID_DONE"
