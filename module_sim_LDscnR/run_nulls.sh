#!/bin/bash
# ---------------------------------------------------------------------------
# Surrogate/permutation nulls over a set of cells.
#
#   ./run_nulls.sh <stage> [nproc]
#     stage   1 | 0.75 | 0.5   which subsample stage's bundles to score
#     nproc   concurrent DRAWS within a cell (default 12)
#
# Cells are processed ONE AT A TIME with the draws forked inside, deliberately:
# a cell's 10 pooled bundles are ~2-4 GB once loaded, and mclapply's forked
# workers share that copy-on-write. Running 12 cells side by side would copy it
# 12 times over; this way the memory cost is one cell regardless of nproc.
#
#   SIM_CELLS         (V,c) cells, default V0.5_c2,V1_c1.5,V2_c1
#   SIM_NULL_ENGINES  emmax | emmax,lfmm      (default emmax)
#   SIM_NULL_TYPES    default all five: genetic,latent,global_perm,env_orth,spatial
#   SIM_NULL_B        draws per type (default 100)
#
# Emits RAW P-VALUES only -- an observed vector and one vector per surrogate. No
# C-scores, no thresholds, no regions: those are the analysis's job, so the
# (rho, q*, alpha) grid and tau_C can change without re-running a scan.
#
# Home-field nulls are engine-specific and NOT crossed (genetic -> EMMAX only,
# latent -> LFMM only); env_orth runs on both and is the arbiter. Defaults follow
# the framework: genetic, latent, env_orth. global_perm is dominated by env_orth
# and spatial is supplementary -- request either explicitly via SIM_NULL_TYPES.
#
# Cost (measured: emmax_fast 0.10 s/scan x 10 chromosomes; LFMM 12.5 s/file):
#   EMMAX  ~1 s per draw   -> B=100 ~2 min per cell x null
#   LFMM   ~2.1 min per draw -> B=100 ~3.5 h per cell x null
# Over 30 cells: EMMAX (2 nulls) ~20 min on 12 cores; LFMM (2 nulls) ~17.5 h.
# Output ~800 MB per cell (context 41 MB + ~190 MB per engine x null at B=100).
#
#   ./run_nulls.sh 1    12                          # EMMAX side
#   SIM_NULL_ENGINES=lfmm ./run_nulls.sh 1 12       # LFMM side
# ---------------------------------------------------------------------------
set -u
STAGE=${1:?usage: run_nulls.sh <stage: 1|0.75|0.5> [nproc]}
NPROC=${2:-12}
ROOT=${SIM_ROOT:-/Volumes/Nemo/Nemo_sim}
HERE=$(cd "$(dirname "$0")" && pwd)
CELLS=${SIM_CELLS:-V0.5_c2,V1_c1.5,V2_c1}

case "$STAGE" in
  1|1.0) SUF="" ;;
  0.75)  SUF="_sub75" ;;
  0.5)   SUF="_sub50" ;;
  *) echo "unknown stage: $STAGE"; exit 1 ;;
esac

export SIM_ROOT=$ROOT
export SIM_DATA=$ROOT/regen_sim_data_nobgs${SUF}
export SIM_NULL_OUT=${SIM_NULL_OUT:-$ROOT/nulls_nobgs${SUF}}
export SIM_NULL_CORES=$NPROC
export LDSCNR_PATH=${LDSCNR_PATH:-}
mkdir -p "$SIM_NULL_OUT"
LOGS=$ROOT/logs_nulls${SUF}; mkdir -p "$LOGS"

echo "stage=$STAGE  cells=$CELLS  engines=${SIM_NULL_ENGINES:-emmax}  B=${SIM_NULL_B:-100}  draws-parallel=$NPROC"
echo "  bundles $SIM_DATA"
echo "  nulls   $SIM_NULL_OUT"

IFS=',' read -ra CELL_ARR <<< "$CELLS"
for cell in "${CELL_ARR[@]}"; do
  V=${cell%%_c*}; V=${V#V}
  CC=${cell##*_c}
  for e in $(seq 1 10); do
    echo "=== V$V c$CC env$e  ($(date +%H:%M)) ==="
    Rscript "$HERE/run_sim_nulls.R" "$V" "$CC" "$e" nobgs \
      >> "$LOGS/V${V}_c${CC}_env${e}.log" 2>&1 || echo "  !! failed V$V c$CC env$e"
  done
done
echo "NULLS DONE $(date)"
