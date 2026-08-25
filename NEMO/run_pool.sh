#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_pool.sh -- run a manifest directory locally with a worker pool.
#
#   ./run_pool.sh grid          # the parameter grid
#   ./run_pool.sh run 3         # the production set, 3 workers
#
# For the cluster use cluster/submit.sh instead; this is for local validation.
#
# Burn-ins all finish before the adapt stage starts, because each adapt run reads
# its burn-in binary as source_pop.
#
# MEMORY-BOUND, not core-bound: one Nemo process holds ~10 GB whatever the phase,
# so the worker count comes from RAM, not core count. Over-subscribing does not
# swap gracefully -- macOS kills the process outright (exit 137).
#
# Resumable: a run whose expected output exists is skipped.
# ---------------------------------------------------------------------------
set -uo pipefail

DIR=${1:?usage: $0 <grid|run> [workers]}
cd "$(dirname "$0")/$DIR" || exit 1

NEMO=${NEMO:-nemo2.4.2-macARM}
GB_PER_PROC=${GB_PER_PROC:-11}
ram_gb=$(( $(sysctl -n hw.memsize 2>/dev/null || echo $((16*1073741824))) / 1073741824 ))
# Reserve headroom for the OS, page cache and whatever else is running. Using
# 80% of RAM was too aggressive: 3 x 11 GB on a 36 GB machine drove 18 GB of
# swap and the runs made no progress at all in 2.5 h.
RESERVE_GB=${RESERVE_GB:-12}
workers=${2:-$(( (ram_gb - RESERVE_GB) / GB_PER_PROC ))}
[ "$workers" -lt 1 ] && workers=1

command -v "$NEMO" >/dev/null || { echo "not on PATH: $NEMO"; exit 1; }

# Guard against the failure mode this script exists to avoid: N x 10 GB over
# physical RAM does not swap, it gets the processes killed (exit 137).
if [ $(( workers * GB_PER_PROC + RESERVE_GB )) -gt "$ram_gb" ]; then
  echo "refusing: $workers x ${GB_PER_PROC} GB + ${RESERVE_GB} GB reserve > ${ram_gb} GB RAM"
  echo "Nemo does not degrade gracefully here -- it swaps and stops progressing."; exit 1
fi
mkdir -p logs out
rm -f failures.txt

run_stage() {
  local stage=$1 manifest="manifest_$1.tsv" launched=0 skipped=0
  [ -f "$manifest" ] || { echo "no $manifest"; return 1; }
  echo "== $stage: $(( $(wc -l < "$manifest") - 1 )) runs, $workers workers (${ram_gb} GB RAM)"

  while IFS=$'\t' read -r id ini expect rest; do
    [ "$id" = "id" ] && continue
    if [ -e "$expect" ]; then skipped=$((skipped+1)); continue; fi
    # `wait -n` needs bash 4; macOS ships 3.2, so poll instead
    while [ "$(jobs -rp | wc -l)" -ge "$workers" ]; do sleep 5; done
    (
      s=$(date +%s)
      if "$NEMO" "$ini" > "logs/$id.out" 2>&1 && [ -e "$expect" ]; then
        echo "  done $id ($(( ($(date +%s) - s) / 60 )) min)"
      else
        echo "  FAILED $id" | tee -a failures.txt
      fi
    ) &
    launched=$((launched+1))
  done < "$manifest"
  wait
  echo "== $stage: $launched launched, $skipped skipped, $( [ -f failures.txt ] && wc -l < failures.txt || echo 0) failed"
  [ -f failures.txt ] && return 1 || return 0
}

run_stage burnin && run_stage adapt
