#!/bin/bash
# Full sweep for one or more tags: three subsample fractions each, sequentially.
#
#   ./run_all.sh [nproc] [tag ...]        default: 14 nobgs
#
# Honours SIM_CELLS and SIM_STAGE (see run_batch.sh). The paper's reduced scope is
#   SIM_CELLS=V0.5_c2,V1_c1.5,V2_c1 ./run_all.sh 12
# = 3 (V,c) cells x 10 env x 10 chromosomes x 3 subsample fractions = 900 files.
#
# Defaults to nobgs ONLY. The bgs2 tarballs are the old deleterious settings
# (mutation rate 2e-7, DFE median on the drift barrier) under which background
# selection is unmeasurable -- add "bgs" explicitly once re-simulated data exists:
#   ./run_all.sh 14 nobgs bgs
#
# Each stage is independently resumable; re-running skips finished bundles.
set -u
NPROC=${1:-10}
shift || true
TAGS=${*:-nobgs}
HERE=$(cd "$(dirname "$0")" && pwd)
for TAG in $TAGS; do
  for SUB in 1 0.75 0.5; do
    echo "=== $TAG  subsample=$SUB  $(date) ==="
    "$HERE/run_batch.sh" "$TAG" "$SUB" "$NPROC" || echo "!! stage failed: $TAG $SUB"
  done
done
echo "ALL DONE $(date)"
