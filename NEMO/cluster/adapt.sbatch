#!/bin/bash
#SBATCH --job-name=nemo-adapt
#SBATCH --output=logs/adapt_%A_%a.out
#SBATCH --error=logs/adapt_%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=14G
#SBATCH --time=01:00:00
# ---------------------------------------------------------------------------
# One array task = one burn-in (10000 generations), read from manifest_adapt.tsv.
#
# Nemo is single-threaded, so this is one core per task. The real resource is
# MEMORY: a 2304 x 8 metapopulation with 1e6 bitstring neutral loci holds ~10 GB
# resident, in this phase and the adapt phase alike. 14G leaves headroom for the
# store step; do not drop it below 12G.
#
# Wall time covers 1000 generations (~4.5 min at 0.265 s/gen) plus loading a
# 4.6 GB source population and writing the genotype output. 1 h is ample.
# ---------------------------------------------------------------------------
set -uo pipefail

NEMO_BIN=${NEMO_BIN:-nemo2.4.2}
cd "${RUN_DIR:?set RUN_DIR to the run/ folder}" || exit 1

read -r id ini expect rest < <(awk -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR==n' manifest_adapt.tsv)
[ -z "${id:-}" ] && { echo "no manifest row for task $SLURM_ARRAY_TASK_ID"; exit 1; }

if [ -e "$expect" ]; then
  echo "SKIP $id (output present)"; exit 0
fi

echo "START $id  $(date -Is)"
"$NEMO_BIN" "$ini" || { echo "FAILED $id"; exit 1; }
[ -e "$expect" ] || { echo "MISSING OUTPUT $id ($expect)"; exit 1; }
echo "DONE $id  $(date -Is)"
