#!/usr/bin/env bash
# Regenerate or verify MANIFEST.md for the untracked 3sp_data/ inputs.
#   ./verify_manifest.sh          -- check every file against MANIFEST.md
#   ./verify_manifest.sh --write  -- regenerate MANIFEST.md
# These files are NOT in git (too large). The manifest is, so a fresh clone can
# tell whether a fetched copy is the right one, and which files it is missing.
set -uo pipefail
cd "$(dirname "$0")"
sum() { shasum -a 256 "$1" | awk '{print $1}'; }        # follows symlinks
size() { if [ -L "$1" ]; then stat -f%z "$(readlink "$1")"; else stat -f%z "$1"; fi; }
files=$(ls -1 | grep -vE '^(MANIFEST\.md|verify_manifest\.sh|rec_maps)$' | sort)

if [ "${1:-}" = "--write" ]; then
  { echo "| file | bytes | sha256 | stored at |"
    echo "|---|---|---|---|"
    for f in $files; do
      [ -f "$f" ] || continue
      loc=$([ -L "$f" ] && readlink "$f" | sed "s|$HOME|~|" || echo "here")
      printf "| \`%s\` | %s | \`%s\` | %s |\n" "$f" "$(size "$f")" "$(sum "$f")" "$loc"
    done; } > /tmp/_mtable
  echo "wrote table for $(echo "$files" | wc -l | tr -d ' ') files"; cat /tmp/_mtable
  exit 0
fi

fail=0; miss=0
while IFS='|' read -r _ f b s _; do
  f=$(echo "$f" | tr -d ' `'); b=$(echo "$b" | tr -d ' '); s=$(echo "$s" | tr -d ' `')
  [ -z "$f" ] || [ "$f" = "file" ] && continue
  if [ ! -e "$f" ]; then echo "MISSING  $f"; miss=$((miss+1)); continue; fi
  a=$(sum "$f")
  if [ "$a" = "$s" ]; then echo "ok       $f"
  else echo "MISMATCH $f"; echo "         expected $s"; echo "         got      $a"; fail=$((fail+1)); fi
done < <(grep '^| `' MANIFEST.md)
echo; echo "$fail mismatched, $miss missing"
[ $fail -eq 0 ] && [ $miss -eq 0 ]
