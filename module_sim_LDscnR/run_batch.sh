#!/bin/bash
# ---------------------------------------------------------------------------
# Regenerate one tag x one subsample fraction, in parallel.
#
#   ./run_batch.sh <tag> <subsample> [nproc]
#     tag        bgs | nobgs   (raw folder is picked to match, see below)
#     subsample  1 | 0.75 | 0.5   fraction of analysis individuals per population
#     nproc      concurrent R processes (default 14)
#
#   SIM_CELLS  restrict to a comma-separated list of (V,c) cells, e.g.
#              SIM_CELLS=V0.5_c2,V1_c1.5,V2_c1   (default: all nine)
#   SIM_STAGE  local scratch to write into before moving results to the volume
#              (default $HOME/sim_stage -- keep it on the INTERNAL SSD)
#
# One R process handles one (V, c, env) cell = 10 chromosomes. 90 cells per tag.
# Resumable: a bundle that already exists is skipped, so re-running after an
# interruption picks up where it stopped.
# ---------------------------------------------------------------------------
set -u
TAG=${1:?usage: run_batch.sh <tag: bgs|nobgs> <subsample: 1|0.75|0.5> [nproc]}
SUB=${2:?usage: run_batch.sh <tag> <subsample> [nproc]}
NPROC=${3:-14}

export SIM_ROOT=${SIM_ROOT:-/Volumes/Nemo/Nemo_sim}
HERE=$(cd "$(dirname "$0")" && pwd)

# raw tarballs per tag
if [ "$TAG" = "bgs" ]; then export SIM_RAW=${SIM_RAW:-$SIM_ROOT/bgs2}
else                        export SIM_RAW=${SIM_RAW:-$SIM_ROOT/Nemo_out_nobgs}; fi

# one output folder per tag x fraction; suffix is empty for the full sample
case "$SUB" in
  1|1.0) SUF="" ;;
  0.75)  SUF="_sub75" ;;
  0.5)   SUF="_sub50" ;;
  *)     SUF="_sub$(echo "$SUB" | tr -d '0.')" ;;
esac
export SIM_OUT=${SIM_OUT:-$SIM_ROOT/regen_sim_data_${TAG}${SUF}}
export SIM_SUBSAMPLE=$SUB
export SIM_CORES=1
# stage writes on the internal SSD; each finished file is moved to the volume
export SIM_STAGE=${SIM_STAGE:-$HOME/sim_stage}
mkdir -p "$SIM_STAGE"

# popgen summaries describe the simulation, not the analysis subset: compute them
# once on the full sample, skip them for the subsampled runs
if [ "$SUF" = "" ]; then export SIM_POPGEN=${SIM_POPGEN:-$SIM_ROOT/popgen_sim_data}
else                       export SIM_POPGEN=""; fi

LOGS=$SIM_ROOT/logs_${TAG}${SUF}
mkdir -p "$LOGS" "$SIM_OUT"
echo "tag=$TAG subsample=$SUB nproc=$NPROC"
echo "  raw    $SIM_RAW"
echo "  out    $SIM_OUT"
echo "  popgen ${SIM_POPGEN:-<skipped>}"
echo "  cells  ${SIM_CELLS:-all}"
echo "  stage  $SIM_STAGE"
echo "  logs   $LOGS"

CELLS=${SIM_CELLS:-all}
for V in 0.5 1 2; do for c in 1 1.5 2; do
  if [ "$CELLS" != "all" ] && ! printf ",%s," "$CELLS" | grep -q ",V${V}_c${c},"; then continue; fi
  for e in $(seq 1 10); do echo "$V $c $e"; done
done; done | xargs -P "$NPROC" -n 3 sh -c '
  Rscript '"$HERE"'/parse_and_regen_sim_data.R "$1" "$2" "$3" all '"$TAG"' \
    > '"$LOGS"'/V"$1"_c"$2"_env"$3".log 2>&1
' sh

DONE=$(ls "$SIM_OUT" | wc -l | tr -d ' ')
echo "BATCH DONE $(date) -- $DONE bundles in $SIM_OUT" | tee -a "$LOGS/_batch_status.log"
