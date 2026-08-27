#!/bin/bash
# ---------------------------------------------------------------------------
# Wait until a target time, then exec another script in this folder. Generalises
# schedule_lfmm_midnight.sh, which did the same for one hard-coded runner.
#
#   nohup caffeinate -is bash schedule_run.sh midnight run_bgs4_nulls.sh 8 \
#         > /Volumes/Nemo/Nemo_sim/logs_nulls_bgs4/_schedule.log 2>&1 &
#
#   $1  target: "midnight" (the next one) or an epoch-seconds value
#   $2  script in this folder to run
#   $3+ arguments passed through to it
#
# The wait re-reads the wall clock every minute rather than sleeping once for
# the whole interval, so a suspend/resume cannot make it fire hours late, and
# the final hop is the exact remainder so it lands on the target second. Run it
# under caffeinate so the machine does not idle-sleep before the target.
# ---------------------------------------------------------------------------
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
SIM_ROOT=$(cd "$HERE/.." && pwd)
SPEC=${1:?usage: schedule_run.sh <midnight|epoch> <script> [args...]}; shift
RUNNER=${1:?usage: schedule_run.sh <midnight|epoch> <script> [args...]}; shift

if [ "$SPEC" = "midnight" ]; then TARGET=$(date -v+1d -v0H -v0M -v0S +%s); else TARGET=$SPEC; fi
if [ ! -f "$HERE/$RUNNER" ]; then echo "!! no such runner: $HERE/$RUNNER"; exit 1; fi

echo "scheduled  $(date -r "$TARGET")"
echo "now        $(date)"
echo "runner     $RUNNER $*"
echo "waiting..."

while :; do
  REM=$(( TARGET - $(date +%s) ))
  [ "$REM" -le 0 ] && break
  [ "$REM" -gt 60 ] && REM=60
  sleep "$REM"
done

echo "=== firing $(date) ==="
if [ ! -d "$SIM_ROOT/pipeline" ]; then
  echo "!! $SIM_ROOT not reachable -- is the drive mounted? aborting"; exit 1
fi
exec bash "$HERE/$RUNNER" "$@"
