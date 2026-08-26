#!/bin/bash
# ---------------------------------------------------------------------------
# Wait until midnight, then run run_lfmm_overnight.sh to fill in the missing
# LFMM p-value files. Launch it detached, e.g.
#
#   nohup caffeinate -is bash schedule_lfmm_midnight.sh 8 \
#         > /Volumes/Nemo/Nemo_sim/logs_nulls/schedule_midnight.log 2>&1 &
#
#   $1  cores          (default 8)
#   $2  target epoch   (default: the next midnight)
#
# The wait re-reads the wall clock every minute rather than sleeping once for
# the whole interval, so a suspend/resume cannot make it fire hours late; run
# it under caffeinate so the machine does not idle-sleep before midnight.
#
# LDscnR is pinned to the 67bc930 snapshot on the drive -- the same commit that
# produced the env2 files -- so the ten environments stay directly comparable
# and a branch switch in the working checkout cannot change the result.
# ---------------------------------------------------------------------------
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
SIM_ROOT=$(cd "$HERE/.." && pwd)
NPROC=${1:-8}
TARGET=${2:-$(date -v+1d -v0H -v0M -v0S +%s)}

export LDSCNR_PATH=$SIM_ROOT/LDscnR_pinned_67bc930
LOGS=$SIM_ROOT/logs_nulls; mkdir -p "$LOGS"

echo "scheduled  $(date -r "$TARGET")"
echo "now        $(date)"
echo "cores      $NPROC"
echo "LDscnR     $LDSCNR_PATH"
echo "waiting..."

# sleep the whole remaining interval, but never more than a minute at a stretch,
# so the clock is re-read regularly AND the last hop lands on the target second
while :; do
  REM=$(( TARGET - $(date +%s) ))
  [ "$REM" -le 0 ] && break
  [ "$REM" -gt 60 ] && REM=60
  sleep "$REM"
done

echo "=== firing $(date) ==="
if [ ! -d "$SIM_ROOT/regen_sim_data_nobgs" ]; then
  echo "!! bundles not reachable at $SIM_ROOT -- is the drive mounted? aborting"; exit 1
fi
if [ ! -d "$LDSCNR_PATH" ]; then
  echo "!! pinned LDscnR missing at $LDSCNR_PATH -- aborting"; exit 1
fi
exec bash "$HERE/run_lfmm_overnight.sh" "$NPROC"
