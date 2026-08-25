#!/usr/bin/env bash
## kingman2021/R/run_all.sh
## Wait for the 15.3 GB joint VCF download to finish, verify it, then build both cohorts.
## Safe to run while the download is still in flight.
set -uo pipefail

DATA=/Users/petrikem/gitlab/LD-scaling-genome-scans/empirical_data/kingman2021
VCF="$DATA/vcf/227_genomes.final.filtered.vcf.gz"
PAPER=/Users/petrikem/gitlab/LDscnR-paper
HERE="$PAPER/kingman2021/R"
WANT=15288279678
export PATH="/private/tmp/claude-539526166/-Users-petrikem-gitlab-LDscnR/e0f3b9b5-e5d3-4cab-8211-c3f2f61c09ff/scratchpad/ucsc:$PATH"

echo "[$(date +%H:%M:%S)] waiting for $VCF to reach $WANT bytes"
last=0; stall=0
while :; do
  cur=$(stat -f%z "$VCF" 2>/dev/null || echo 0)
  [[ "$cur" -ge "$WANT" ]] && break
  if [[ "$cur" -eq "$last" ]]; then
    stall=$((stall+1))
    if [[ "$stall" -ge 40 ]]; then   # ~20 min with no growth
      echo "[$(date +%H:%M:%S)] ABORT: download stalled at $cur/$WANT bytes"; exit 1
    fi
  else stall=0; fi
  last=$cur
  sleep 30
done
echo "[$(date +%H:%M:%S)] download complete: $(stat -f%z "$VCF") bytes"

## BGZF EOF marker present == not truncated
if ! bcftools view -h "$VCF" >/dev/null 2>bgzf.err; then
  echo "ABORT: bcftools cannot read the VCF header"; cat bgzf.err; exit 1
fi
if bcftools view -h "$VCF" 2>&1 | grep -q "no BGZF EOF marker"; then
  echo "ABORT: BGZF EOF marker missing -- file is truncated"; exit 1
fi
## the .tbi was downloaded before the .gz, so refresh its mtime to stop the
## "index is older than the data file" warning being mistaken for a stale index
touch "$VCF.tbi"
echo "[$(date +%H:%M:%S)] VCF verified"

cd "$PAPER"
for COH in c155_global c150_pacNW; do
  echo "=============================================================="
  echo "[$(date +%H:%M:%S)] $COH : extracting genotypes"
  bash "$HERE/01_extract_gts.sh" "$COH" 0.05 0.20 || { echo "extract failed for $COH"; continue; }
  echo "[$(date +%H:%M:%S)] $COH : building rds"
  Rscript "$HERE/02_build_rds.R" "$COH"          || { echo "build failed for $COH";   continue; }
  echo "[$(date +%H:%M:%S)] $COH : tagging peaks"
  Rscript "$HERE/03_peaks_truth.R" "$COH"        || { echo "truth failed for $COH";   continue; }
done

echo "=============================================================="
echo "[$(date +%H:%M:%S)] DONE"
ls -la "$PAPER"/kingman2021/data/kingman2021_*.rds
