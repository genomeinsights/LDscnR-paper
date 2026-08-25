#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# setup_machine.sh -- get a fresh machine ready to run these simulations.
#
#   ./setup_machine.sh                 check only, report what is missing
#   ./setup_machine.sh --build-nemo    also clone, patch and build Nemo
#
# Checks, in order: Nemo on PATH and new enough; the macOS binary-loader patch,
# without which the adapt phase cannot read a burn-in; R and data.table; the
# params_V4 inputs, which are NOT in git and have to be copied across; and how
# many workers this machine's RAM allows.
# ---------------------------------------------------------------------------
set -uo pipefail
cd "$(dirname "$0")"

NEMO_BIN=${NEMO_BIN:-nemo2.4.2-macARM}
NEMO_SRC=${NEMO_SRC:-$HOME/Nemo/nemo-release}
GSL_PATH=${GSL_PATH:-/opt/homebrew/opt/gsl/}
GB_PER_PROC=${GB_PER_PROC:-11}
RESERVE_GB=${RESERVE_GB:-12}
fail=0
say() { printf '%-14s %s\n' "$1" "$2"; }

if [ "${1:-}" = "--build-nemo" ]; then
  echo "== building Nemo"
  if [ ! -f "${GSL_PATH}lib/libgsl.a" ]; then
    echo "GSL not found at ${GSL_PATH}lib/libgsl.a"
    echo "  install it:  brew install gsl        (or set GSL_PATH=)"
    echo "  the MAC_ARM build links libgsl.a and libgslcblas.a statically"
    exit 1
  fi
  if [ ! -d "$NEMO_SRC/.git" ]; then
    mkdir -p "$(dirname "$NEMO_SRC")"
    git clone https://bitbucket.org/ecoevo/nemo-release.git "$NEMO_SRC" || exit 1
  fi
  # The patch is required on macOS: a stored population here is 4.6 GB and the
  # loader asks for it in one read(), which macOS rejects with EINVAL. Harmless
  # to apply on Linux, where it is not needed.
  if git -C "$NEMO_SRC" apply --check "$PWD/nemo-macos-bigpop.patch" 2>/dev/null; then
    git -C "$NEMO_SRC" apply "$PWD/nemo-macos-bigpop.patch" && say "patch" "applied"
  else
    say "patch" "already applied (or does not apply -- check manually)"
  fi
  ( cd "$NEMO_SRC" && mkdir -p bin && make MAC_ARM=1 GSL_PATH="$GSL_PATH" -j8 >/tmp/nemo_build.log 2>&1 ) \
    || { say "build" "FAILED -- see /tmp/nemo_build.log"; exit 1; }
  mkdir -p "$HOME/Nemo" && cp "$NEMO_SRC"/bin/nemo* "$HOME/Nemo/" && say "build" "ok -> $HOME/Nemo"
  echo "   add to PATH:  export PATH=\$HOME/Nemo:\$PATH"
  echo
fi

echo "== checks"
if command -v "$NEMO_BIN" >/dev/null; then
  say "nemo" "$(command -v "$NEMO_BIN")"
  # Nemo prints its banner and exits non-zero with no input file; that is fine.
  v=$("$NEMO_BIN" 2>&1 | grep -o 'N E M O [0-9.]*' | head -1)
  say "version" "${v:-unknown}   (need >= 2.4.2: 2.4.0/2.4.1 have the RandBool defect)"
else
  say "nemo" "MISSING -- run with --build-nemo"; fail=1
fi

if grep -q "READ_CHUNK_MAX" "$NEMO_SRC/src/binarydataloader.cc" 2>/dev/null; then
  say "loader patch" "present in $NEMO_SRC"
else
  say "loader patch" "NOT FOUND in $NEMO_SRC -- the adapt phase will fail on macOS"; fail=1
fi

if command -v Rscript >/dev/null; then
  say "R" "$(Rscript -e 'cat(R.version.string)' 2>/dev/null)"
  Rscript -e 'q(status = !requireNamespace("data.table", quietly = TRUE))' 2>/dev/null \
    && say "data.table" "ok" || { say "data.table" "MISSING -- install.packages('data.table')"; fail=1; }
else
  say "R" "MISSING"; fail=1
fi

n=$(ls params_V4/map_ntrl*.txt 2>/dev/null | wc -l | tr -d ' ')
d=$(ls params_V4/disp_mat_*.txt 2>/dev/null | wc -l | tr -d ' ')
e=$(ls params_V4/env_*.txt 2>/dev/null | wc -l | tr -d ' ')
if [ "$n" = 10 ] && [ "$d" = 3 ] && [ "$e" = 10 ]; then
  say "params_V4" "ok (10 maps, 3 dispersal matrices, 10 environments, $(du -sh params_V4 | cut -f1))"
else
  say "params_V4" "INCOMPLETE ($n maps, $d disp, $e env) -- not in git, copy it across"; fail=1
fi

ram=$(( $(sysctl -n hw.memsize 2>/dev/null || echo 0) / 1073741824 ))
w=$(( (ram - RESERVE_GB) / GB_PER_PROC )); [ "$w" -lt 1 ] && w=1
say "RAM" "${ram} GB -> ${w} workers (${GB_PER_PROC} GB each, ${RESERVE_GB} GB reserved)"
say "load" "$(uptime | sed 's/.*averages*: //')"
[ "$ram" -lt 24 ] && { say "WARNING" "under 24 GB cannot hold even one process comfortably"; }

echo
if [ "$fail" = 0 ]; then
  echo "ready. Next:"
  echo "  Rscript R/make_grid.R && ./run_pool.sh grid     # parameter grid first (~2.5 h at 3 workers)"
  echo "  Rscript R/make_run.R  && ./run_pool.sh run      # the production set"
else
  echo "not ready -- fix the items marked above."; exit 1
fi
