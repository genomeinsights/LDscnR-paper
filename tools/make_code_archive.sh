#!/usr/bin/env bash
## tools/make_code_archive.sh
##
## Builds the code archive for the manuscript deposit from COMMITTED files only
## (git archive), so untracked data, caches, credentials and session files can
## never enter it. Adds the manuscript repository's value/table generators and
## figure scripts under manuscript_scripts/, then checks the result.
##
## Usage: tools/make_code_archive.sh [OUT_ZIP] [COMMIT] [MANUSCRIPT_REPO] [MANUSCRIPT_COMMIT]
##   defaults: ../LDscnR_manuscript/deposit/LDscnR-paper-code.zip, HEAD,
##             ../LDscnR_manuscript, HEAD
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
OUT=${1:-$ROOT/../LDscnR_manuscript/deposit/LDscnR-paper-code.zip}
COMMIT=${2:-HEAD}
MS=${3:-$ROOT/../LDscnR_manuscript}
MS_COMMIT=${4:-HEAD}
NAME=LDscnR-paper-code

## The manuscript's analysis code. Excluded on purpose: superseded/exploratory
## modules (legacy/, module_C2/, module_sim/, module_sim_LDscnR/, module_sim_bgs*/,
## module_sticklebacks_LDscnR/, manuscript/), the kingman2021 side analyses
## (only its lifted EcoPeak BEDs and liftover script are inputs here), the old
## C-score scripts in R/, and working notes (AI-session instructions, audits,
## a superseded write-up with outdated numbers).
PATHS=(LICENSE README.md R/00_fetch_inputs.R tools/make_code_archive.sh
       module_3sp module_9sp module_sim_3sp53 3sp_data/rec_maps
       kingman2021/data/liftover kingman2021/R/08_liftover.sh)
EXCLUDE=(module_sim_3sp53/CLAUDE_REANALYSIS_INSTRUCTIONS.md
         module_sim_3sp53/final_analysis/ADDITIONAL_ANALYSES_AUDIT.md
         module_sim_3sp53/final_analysis/FLOOR_STABILITY_AUDIT_NOTE.md
         module_9sp/PERMUTATION_VS_STRUCTURE.md)
MS_FILES=(generate_empirical_values.R generate_simulation_manuscript.R R_figures)

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP/$NAME/manuscript_scripts"
git -C "$ROOT" archive "$COMMIT" -- "${PATHS[@]}" | tar -x -C "$TMP/$NAME"
for x in "${EXCLUDE[@]}"; do rm -f "$TMP/$NAME/$x"; done
find "$TMP/$NAME" -name .gitignore -delete
git -C "$MS" archive "$MS_COMMIT" -- "${MS_FILES[@]}" | tar -x -C "$TMP/$NAME/manuscript_scripts"
printf 'LDscnR-paper %s\nLDscnR_manuscript %s\nbuilt %s\n' \
  "$(git -C "$ROOT" rev-parse "$COMMIT")" "$(git -C "$MS" rev-parse "$MS_COMMIT")" \
  "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$TMP/$NAME/ARCHIVE_COMMITS.txt"

## ---- checks: nothing that must not be published --------------------------------
fail=0
bad=$(cd "$TMP/$NAME" && find . \( -iname 'key.txt' -o -name '.Renviron' -o -name '.Rhistory' \
        -o -name '.RData' -o -name '.RDataTmp' -o -name '.httr-oauth' -o -name '.netrc' -o -iname '*token*' \
        -o -name 'HANDOFF_*' -o -iname '*claude*' -o -path '*/cache/*' -o -name '.DS_Store' \) -print)
[ -n "$bad" ] && { echo "FORBIDDEN FILES:"; echo "$bad"; fail=1; }
big=$(find "$TMP/$NAME" -type f -size +50M -print)
[ -n "$big" ] && { echo "FILES > 50 MB:"; echo "$big"; fail=1; }
secrets=$(grep -rIlE '(ghp_[A-Za-z0-9]{20,}|glpat-[A-Za-z0-9_-]{20,}|AKIA[0-9A-Z]{16}|BEGIN (RSA|OPENSSH) PRIVATE KEY)' "$TMP/$NAME" || true)
[ -n "$secrets" ] && { echo "POSSIBLE SECRETS:"; echo "$secrets"; fail=1; }
[ "$fail" = 0 ] || { echo "archive NOT written"; exit 1; }

mkdir -p "$(dirname "$OUT")"; rm -f "$OUT"
(cd "$TMP" && zip -qr -X "$OUT" "$NAME")
echo "wrote $OUT ($(du -h "$OUT" | cut -f1), $(cd "$TMP/$NAME" && find . -type f | wc -l | tr -d ' ') files)"
(cd "$TMP/$NAME" && du -sh * | sort -k2)
