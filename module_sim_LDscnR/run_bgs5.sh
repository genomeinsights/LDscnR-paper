#!/bin/bash
# ---------------------------------------------------------------------------
# Regenerate the bgs5 tarballs into bundles.
#
#   ./run_bgs5.sh [nproc]        default 8
#
# bgs5 supersedes bgs4: matched bgs/nobgs at FOUR (V,c) cells -- V0.5_c1 is new --
# across TEN environmental gradients, 10 chromosomes each. 800 files, fully
# crossed with no gaps. The ten gradients are what let the BGS analysis be
# replicate-averaged and paired within gradient, rather than resting on env1.
#
# Parallelises over single chromosomes: one R process per (tag, cell, env, chr).
# Output folders are bgs5-specific so nothing collides with the delivered
# regen_sim_data_nobgs set or with bgs4. Resumable -- finished bundles skip.
#
# COST: ~55 s per chromosome, so ~1.7 h for all 800 on 8 cores.
# ---------------------------------------------------------------------------
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
NPROC=${1:-8}
export SIM_ROOT=${SIM_ROOT:-/Volumes/Nemo/Nemo_sim}
export SIM_RAW=$SIM_ROOT/bgs5
export SIM_OUT=$SIM_ROOT/regen_sim_data_bgs5
export SIM_POPGEN=$SIM_ROOT/popgen_sim_data_bgs5
export SIM_STAGE=${SIM_STAGE:-$HOME/sim_stage}
export SIM_SUBSAMPLE=1
export SIM_CORES=1
# a64356f, not 67bc930: min_r2 is now derived per chromosome at min_r2_rho = 0.5
# instead of being a fixed 0.2 floor. PK standardised every analysis on the
# derived form on 2026-08-29. The older pin is kept for the sets built under the
# fixed floor -- see archive_minr2_0.2/README.md.
export LDSCNR_PATH=${LDSCNR_PATH:-$SIM_ROOT/LDscnR_pinned_a64356f}
LOGS=$SIM_ROOT/logs_bgs5
mkdir -p "$SIM_OUT" "$SIM_POPGEN" "$SIM_STAGE" "$LOGS"

echo "raw    $SIM_RAW"
echo "out    $SIM_OUT"
echo "popgen $SIM_POPGEN"
echo "nproc  $NPROC"
echo "start  $(date)"
if [ ! -d "$SIM_RAW" ]; then echo "!! no raw tarballs at $SIM_RAW"; exit 1; fi

{
  for TAG in bgs nobgs; do
    for cell in "0.5 1" "0.5 2" "1 1.5" "2 1"; do
      set -- $cell
      for e in $(seq 1 10); do
        for ch in $(seq 1 10); do echo "$TAG $1 $2 $e $ch"; done
      done
    done
  done
} | xargs -P "$NPROC" -n 5 sh -c '
  Rscript '"$HERE"'/parse_and_regen_sim_data.R "$2" "$3" "$4" "$5" "$1" \
    > '"$LOGS"'/"$1"_V"$2"_c"$3"_env"$4"_chr"$5".log 2>&1
' sh

echo "BGS5 DONE $(date)"
ls "$SIM_OUT" | wc -l | xargs echo "bundles:"
ls "$SIM_POPGEN"/*.rds 2>/dev/null | wc -l | xargs echo "popgen rds:"
