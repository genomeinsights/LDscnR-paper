#!/bin/bash
# ---------------------------------------------------------------------------
# Regenerate the bgs4 tarballs -- the paired bgs/nobgs design, three (V,c)
# cells, env1 only -- into bundles.
#
#   ./run_bgs4.sh [nproc]        default 8
#
# Unlike run_batch.sh this parallelises over single CHROMOSOMES rather than
# whole cells: there are only five complete (tag, cell) combinations here, so
# cell-level parallelism would leave most cores idle.
#
# Both tags share one output folder on purpose. pool_cell() in run_sim_nulls.R
# globs "^adapt_<TAG>_chr[0-9]+_..." so the tag in the filename keeps them
# apart, and keeping them together means one SIM_DATA for the whole design.
#
# The output folders are bgs4-specific, which matters: SIM_OUT for tag=nobgs
# would otherwise default to regen_sim_data_nobgs (the delivered 100-bundle
# set) and SIM_POPGEN to the shared popgen_sim_data, where this run's
# nobgs_V2_c1_env1 summary would land on top of a DIFFERENT simulation's file
# of the same name.
#
# Resumable: a bundle that already exists is skipped, so re-running after new
# tarballs land does only the missing work. (nobgs V2_c1 was 9/10 until its
# chr7 Nemo run was repeated; all six tag x cell combinations are complete now.)
# ---------------------------------------------------------------------------
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
NPROC=${1:-8}
export SIM_ROOT=${SIM_ROOT:-/Volumes/Nemo/Nemo_sim}
export SIM_RAW=$SIM_ROOT/bgs4
export SIM_OUT=$SIM_ROOT/regen_sim_data_bgs4
export SIM_POPGEN=$SIM_ROOT/popgen_sim_data_bgs4
export SIM_STAGE=${SIM_STAGE:-$HOME/sim_stage}
export SIM_SUBSAMPLE=1
# pin the package to the snapshot on the drive, as the null runs do, so a
# branch switch in a working checkout cannot change what gets regenerated
export LDSCNR_PATH=${LDSCNR_PATH:-$SIM_ROOT/LDscnR_pinned_67bc930}
export SIM_CORES=1
LOGS=$SIM_ROOT/logs_bgs4
mkdir -p "$SIM_OUT" "$SIM_POPGEN" "$SIM_STAGE" "$LOGS"

echo "raw    $SIM_RAW"
echo "out    $SIM_OUT"
echo "popgen $SIM_POPGEN"
echo "stage  $SIM_STAGE"
echo "logs   $LOGS"
echo "nproc  $NPROC"
if [ ! -d "$SIM_RAW" ]; then echo "!! no raw tarballs at $SIM_RAW"; exit 1; fi

{
  for cell in "0.5 2" "1 1.5" "2 1"; do
    set -- $cell
    for ch in $(seq 1 10); do echo "bgs $1 $2 $ch"; done
  done
  for cell in "0.5 2" "1 1.5" "2 1"; do
    set -- $cell
    for ch in $(seq 1 10); do echo "nobgs $1 $2 $ch"; done
  done
} | xargs -P "$NPROC" -n 4 sh -c '
  Rscript '"$HERE"'/parse_and_regen_sim_data.R "$2" "$3" 1 "$4" "$1" \
    > '"$LOGS"'/"$1"_V"$2"_c"$3"_chr"$4"_env1.log 2>&1
' sh

echo "BGS4 DONE $(date)"
ls "$SIM_OUT" | wc -l | xargs echo "bundles in $SIM_OUT:"
ls "$SIM_POPGEN" | wc -l | xargs echo "popgen files:"
