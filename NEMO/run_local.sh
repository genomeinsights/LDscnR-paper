#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_local.sh -- smoke-test one burn-in/adapt pair before submitting to SLURM.
#
#   ./run_local.sh 20        # 20-generation burn-in, then its adapt run
#
# Confirms the maps load, the deleterious trait registers, and the adapt phase
# can chain off the burn-in binary. It does NOT produce usable data -- a real
# burn-in is 10000 generations.
#
# On macOS this needs the binary-loader patch in PATCHES.md; without it the
# adapt step fails with a misleading "Binary file appears corrupted".
# ---------------------------------------------------------------------------
set -euo pipefail

GENS=${1:-20}
NEMO_BIN=${NEMO_BIN:-nemo2.4.2-macARM}
cd "$(dirname "$0")/run"

command -v "$NEMO_BIN" >/dev/null || { echo "not on PATH: $NEMO_BIN"; exit 1; }
[ -d ini ] || { echo "no ini/ -- run Rscript ../R/make_run.R first"; exit 1; }

mkdir -p out logs smoke
burn=$(awk 'NR==2{print $1}' manifest_burnin.tsv)
adapt=$(awk -v b="$burn" '$4==b{print $1; exit}' manifest_adapt.tsv)
echo "burn-in: $burn"
echo "adapt:   $adapt"

sed -e "s/^generations     .*/generations     $GENS/" \
    -e "s/^store_generation  .*/store_generation  $GENS/" \
    -e "s/^stat_log_time       .*/stat_log_time       $GENS/" \
    "ini/burnin_$burn.ini" > "smoke/burnin.ini"
sed -e "s/^generations     .*/generations     $GENS/" \
    -e "s/^store_generation  .*/store_generation  $GENS/" \
    -e "s/^stat_log_time       .*/stat_log_time       $GENS/" \
    -e "s/^files_genotyper_logtime    .*/files_genotyper_logtime    $GENS/" \
    "ini/$adapt.ini" > "smoke/adapt.ini"

echo "== burn-in ($GENS generations)"; "$NEMO_BIN" smoke/burnin.ini | tail -2
echo "== adapt   ($GENS generations)"; "$NEMO_BIN" smoke/adapt.ini   | tail -2

geno="out/$adapt/GENO/${adapt}_${GENS}_1.snp_geno"
if [ -s "$geno" ]; then
  echo "OK: $geno  ($(wc -l < "$geno") lines)"
else
  echo "FAILED: no genotype output at $geno"; exit 1
fi
