#!/bin/bash
# ---------------------------------------------------------------------------
# LFMM p-values for the framework's LFMM bases -- an overnight job, meant to be
# run STRAIGHT FROM THE EXTERNAL DRIVE on a second machine.
#
#   ./run_lfmm_overnight.sh [nproc] [cells]
#     nproc   concurrent surrogate scans   (default 12)
#     cells   comma-separated (V,c) cells  (default V2_c1)
#
# Set LDSCNR_PATH to an LDscnR checkout at 67bc930 or later; omit it only if
# LDscnR is installed as a package there and is that recent. Nothing else needs
# configuring -- SIM_ROOT is derived from this script's own location, so the
# drive can mount anywhere.
#
#   LDSCNR_PATH=~/gitlab/LDscnR ./run_lfmm_overnight.sh 12
#
# Produces, in $SIM_ROOT/analysis_inputs:
#   pvals_<cell>_lfmm_latent_B100.rds     LFMM's home-field null
#   pvals_<cell>_lfmm_env_orth_B100.rds   the method-agnostic arbiter
# The panels are already there from the EMMAX run and are reused untouched.
#
# COST: LFMM is ~12.5 s per chromosome file, so ~2.1 min per pooled surrogate.
# B=100 is ~3.5 h per cell x basis; ten environments x two bases is ~70 core-h,
# roughly 6 h on 12 cores. RESUMABLE -- an existing output file is skipped, so
# if it is interrupted just run it again.
# ---------------------------------------------------------------------------
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
export SIM_ROOT=${SIM_ROOT:-$(cd "$HERE/.." && pwd)}       # the drive, wherever it mounted
NPROC=${1:-12}
CELLS=${2:-V2_c1}

export SIM_DATA=$SIM_ROOT/regen_sim_data_nobgs
export SIM_NULL_OUT=$SIM_ROOT/analysis_inputs
export SIM_NULL_ENGINES=lfmm
export SIM_NULL_TYPES=latent,env_orth
export SIM_NULL_B=${SIM_NULL_B:-100}
export SIM_NULL_CORES=$NPROC
LOGS=$SIM_ROOT/logs_nulls; mkdir -p "$LOGS" "$SIM_NULL_OUT"

echo "SIM_ROOT   $SIM_ROOT"
echo "bundles    $SIM_DATA"
echo "output     $SIM_NULL_OUT"
echo "LDscnR     ${LDSCNR_PATH:-<installed package>}"
echo "cells      $CELLS   B=$SIM_NULL_B   concurrent scans=$NPROC"
if [ ! -d "$SIM_DATA" ]; then echo "!! no bundles at $SIM_DATA"; exit 1; fi

IFS=',' read -ra CELL_ARR <<< "$CELLS"
for cell in "${CELL_ARR[@]}"; do
  V=${cell%%_c*}; V=${V#V}; CC=${cell##*_c}
  for e in $(seq 1 10); do
    done_latent=$SIM_NULL_OUT/pvals_V${V}_c${CC}_env${e}_lfmm_latent_B${SIM_NULL_B}.rds
    done_orth=$SIM_NULL_OUT/pvals_V${V}_c${CC}_env${e}_lfmm_env_orth_B${SIM_NULL_B}.rds
    if [ -f "$done_latent" ] && [ -f "$done_orth" ]; then
      echo "=== V$V c$CC env$e  already done, skipping ==="; continue
    fi
    echo "=== V$V c$CC env$e  ($(date +%H:%M)) ==="
    Rscript "$HERE/run_sim_nulls.R" "$V" "$CC" "$e" nobgs \
      >> "$LOGS/lfmm_V${V}_c${CC}_env${e}.log" 2>&1 || echo "  !! failed V$V c$CC env$e"
  done
done
echo "LFMM DONE $(date)"
ls "$SIM_NULL_OUT"/pvals_*_lfmm_*.rds 2>/dev/null | wc -l | xargs echo "lfmm pvals files:"
