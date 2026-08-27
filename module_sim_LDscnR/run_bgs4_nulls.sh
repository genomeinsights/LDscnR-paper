#!/bin/bash
# ---------------------------------------------------------------------------
# Observed + surrogate p-values for the bgs4 paired design: two tags (bgs,
# nobgs) x three (V,c) cells x env1 = six panels.
#
#   ./run_bgs4_nulls.sh [nproc]      default 8
#
# Engines and bases follow the framework, as for the nobgs delivery:
#   emmax genetic     home-field null for the genetic engine
#   lfmm  latent      home-field null for the latent engine
#   both  env_orth    the method-agnostic arbiter
# HOME in run_sim_nulls.R pins each home-field null to its own engine, so the
# cross pairs are never computed. Four p-value files per panel, 24 in all.
#
# OUTPUT DIRECTORY IS PER TAG, and that is not cosmetic. cell_id in
# run_sim_nulls.R is "V<V>_c<c>_env<N>" with NO tag in it, so bgs and nobgs
# would write the same panel_*.rds and pvals_*.rds names -- and the script
# REUSES an existing panel rather than rebuilding it. Sharing a directory
# (with each other, or with the delivered analysis_inputs) would silently scan
# one tag's phenotypes against the other tag's genotypes.
#
# COST: EMMAX is ~0.5 min per basis, LFMM ~35 min per basis at B=100 on 8
# cores, so ~75 min per panel and ~7.5 h for all six. RESUMABLE -- existing
# outputs are skipped.
# ---------------------------------------------------------------------------
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
export SIM_ROOT=${SIM_ROOT:-$(cd "$HERE/.." && pwd)}
NPROC=${1:-8}

export SIM_DATA=$SIM_ROOT/regen_sim_data_bgs4
export SIM_NULL_ENGINES=emmax,lfmm
export SIM_NULL_TYPES=genetic,latent,env_orth
export SIM_NULL_B=${SIM_NULL_B:-100}
export SIM_NULL_CORES=$NPROC
export LDSCNR_PATH=${LDSCNR_PATH:-$SIM_ROOT/LDscnR_pinned_67bc930}
LOGS=$SIM_ROOT/logs_nulls_bgs4; mkdir -p "$LOGS"

echo "bundles  $SIM_DATA"
echo "LDscnR   $LDSCNR_PATH"
echo "engines  $SIM_NULL_ENGINES   bases $SIM_NULL_TYPES   B=$SIM_NULL_B   cores=$NPROC"
if [ ! -d "$SIM_DATA" ]; then echo "!! no bundles at $SIM_DATA"; exit 1; fi

for TAG in bgs nobgs; do
  export SIM_NULL_OUT=$SIM_ROOT/analysis_inputs_bgs4/$TAG
  mkdir -p "$SIM_NULL_OUT"
  for cell in "0.5 2" "1 1.5" "2 1"; do
    set -- $cell
    V=$1; CC=$2
    n=$(ls "$SIM_DATA" | grep -c "^adapt_${TAG}_chr[0-9]*_V${V}_c${CC}_env1\.rds$")
    if [ "$n" -ne 10 ]; then
      echo "=== $TAG V$V c$CC  SKIP: $n/10 chromosome bundles ==="; continue
    fi
    echo "=== $TAG V$V c$CC env1  ($(date +%H:%M)) -> $SIM_NULL_OUT ==="
    Rscript "$HERE/run_sim_nulls.R" "$V" "$CC" 1 "$TAG" \
      >> "$LOGS/${TAG}_V${V}_c${CC}_env1.log" 2>&1 || echo "  !! failed $TAG V$V c$CC"
  done
done
echo "BGS4 NULLS DONE $(date)"
find "$SIM_ROOT/analysis_inputs_bgs4" -name "pvals_*.rds" | wc -l | xargs echo "pvals files:"
