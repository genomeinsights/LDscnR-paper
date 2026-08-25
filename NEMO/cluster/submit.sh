#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# submit.sh -- submit both stages to SLURM with the right dependency.
#
#   ACCOUNT=project_2003847 ./submit.sh                        # the production set
#   ACCOUNT=project_2003847 RUN_DIR=$PWD/../grid ./submit.sh    # the parameter grid
#
# Every adapt run reads its cell's burn-in binary, so the adapt array is held
# until the burn-in array completes. That is coarser than strictly necessary --
# an adapt job only needs its own cell -- but SLURM's element-wise `aftercorr`
# cannot map a 60-task array onto a 600-task one. The cost is that the last
# burn-in gates the first adapt; the gain is that it cannot race.
#
# Re-running is safe: both scripts skip a task whose output already exists, so
# resubmitting after a partial failure only redoes what is missing.
# ---------------------------------------------------------------------------
set -euo pipefail

cd "$(dirname "$0")"
RUN_DIR=${RUN_DIR:-$(cd ../run && pwd)}
ACCOUNT=${ACCOUNT:?set ACCOUNT, e.g. ACCOUNT=project_2003847}
PARTITION=${PARTITION:-small}
THROTTLE=${THROTTLE:-50}          # max concurrent array tasks
NEMO_BIN=${NEMO_BIN:-nemo2.4.2}

[ -f "$RUN_DIR/manifest_burnin.tsv" ] || { echo "no manifests in $RUN_DIR -- run Rscript ../R/make_run.R first"; exit 1; }
mkdir -p "$RUN_DIR/logs" "$RUN_DIR/out"

n_burnin=$(( $(wc -l < "$RUN_DIR/manifest_burnin.tsv") - 1 ))
n_adapt=$((  $(wc -l < "$RUN_DIR/manifest_adapt.tsv")  - 1 ))
echo "submitting $n_burnin burn-ins and $n_adapt adapt runs (throttle $THROTTLE)"

common=(--account="$ACCOUNT" --partition="$PARTITION"
        --chdir="$RUN_DIR"
        --export=ALL,RUN_DIR="$RUN_DIR",NEMO_BIN="$NEMO_BIN")

jid=$(sbatch --parsable "${common[@]}" --array=1-"$n_burnin"%"$THROTTLE" burnin.sbatch)
echo "burn-in array: $jid"

ajid=$(sbatch --parsable "${common[@]}" --dependency=afterok:"$jid" \
              --array=1-"$n_adapt"%"$THROTTLE" adapt.sbatch)
echo "adapt array:   $ajid  (held on $jid)"
