## =============================================================================
## module_sim/R_parsing/01_parse_nemo.R
##
## PARSE ONE RAW NEMO OUTPUT FILE INTO A CLEAN, UNSCORED BUNDLE.
##
## Separate from module_sim/R/ on PK's instruction: R_parsing/ turns NEMO's
## native simulator output into the raw substrate this project's analyses
## consume; R/ (02_bundle.R onward, mirroring module_3sp's stage numbers) is
## where decay, ld_w, clustering, kinship and association scans happen. Nothing
## in this file computes a p-value, a kinship, or an LD statistic -- if it does,
## that is a bug in this stage, not a convenience, same contract as
## module_sim/R/00_config.R itself.
##
## WHAT THIS REPLACES. R/Parse_sim_data.R (LDscnR-paper root, 654 lines) reads
## the same NEMO format but is unreproducible in the ways this project spent
## this week fixing elsewhere: a hardcoded "BATCH SELECTOR" line hand-edited
## between runs, relative paths, no seed, no receipt, and (checked directly,
## 2026-09-04) a chr_type assignment of `ifelse(Chr==1,"QTN","ntrl")` that is
## wrong on its face once a chromosome other than 1 carries a QTN -- which every
## checked bundle shows every chromosome does. The NEMO-format-reading logic in
## that script is sound and is reused here (the .map/.snp_geno layout, the
## reference-map index matching, the environment-file format); the surrounding
## practice is not.
##
## SCOPE, DELIBERATELY NARROW. One file: TARGET_CHR/TARGET_CELL/TARGET_ENV/
## TARGET_TAG below, matching module_sim/R/00_config.R's CELLS/TAGS/ENVS (both
## currently pinned to the same single slice). Do not widen this script to loop
## over the grid until this one file is verified end to end through R/'s stages
## too -- the whole point of narrowing scope today.
## =============================================================================
suppressMessages({library(data.table)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim"), "R", "00_config.R"))
STAGE <- "01_parse_nemo"
say("=== %s ===\n\n", STAGE)

TARGET_TAG  <- TAGS[1]     # "nobgs"
TARGET_CELL <- CELLS[1]    # "V2_c1"
TARGET_ENV  <- ENVS[1]     # 1L
TARGET_CHR  <- CHRS[1]     # 1L -- ONE chromosome file; 02_bundle.R pools the ten

raw_dir <- if (TARGET_TAG == "nobgs") PATHS$raw_nemo_nobgs else PATHS$raw_nemo_bgs5
archive <- file.path(raw_dir, sprintf("adapt_%s_chr%d_%s_env%d.tgz",
                                      TARGET_TAG, TARGET_CHR, TARGET_CELL, TARGET_ENV))
recmap_rds <- file.path(PATHS$raw_recmap_dir, sprintf("rec_map%d.rds", TARGET_CHR))
env_txt    <- file.path(PATHS$raw_env_dir, sprintf("env_%d.txt", TARGET_ENV))

INPUTS <- c(archive = archive, recmap = recmap_rds, env = env_txt)
PARAMS <- list(tag = TARGET_TAG, cell = TARGET_CELL, env = TARGET_ENV, chr = TARGET_CHR,
               subsample_step = SUBSAMPLE_STEP, maf_keep = MAF_KEEP)

say("[0] target: tag=%s cell=%s env=%d chr=%d\n", TARGET_TAG, TARGET_CELL, TARGET_ENV, TARGET_CHR)
for (nm in names(INPUTS)) say("    %-8s %s  (%s)\n", nm, INPUTS[[nm]],
                              if (file.exists(INPUTS[[nm]])) "exists" else "MISSING")
## [!] BUG, FOUND 2026-09-04: stage_stale() here was called with the NAMED
## INPUTS vector (kept named for the logging loop above), while write_receipt()
## below uses unname(INPUTS) -- matching module_3sp/01_inputs.R's own
## convention. vapply() over a NAMED vector returns a result named by THOSE
## names ("archive"/"recmap"/"env"); the receipt stores hashes keyed by the
## actual PATH strings. stage_stale()'s lookup (old[names(now)]) then compared
## against names that never existed in the receipt, so it reported every input
## "changed" on every run regardless of whether anything had. Confirmed: a
## rerun immediately after a successful one still said "inputs changed: archive,
## recmap, env". Both calls now use unname(INPUTS), matching write_receipt.
if (!stage_stale(STAGE, unname(INPUTS), PARAMS) && !nzchar(Sys.getenv("FORCE"))) {
  say("\nNothing to do. Set FORCE=1 to rerun anyway.\n"); quit(save = "no")
}
if (any(!file.exists(INPUTS))) stop("missing input(s): ",
    paste(names(INPUTS)[!file.exists(INPUTS)], collapse = ", "))

## ---- 1. unpack the archive ----------------------------------------------------
## Into PATHS$cache/untar, NOT a relative "./tmp_<tag>" that has to be hand-kept
## from colliding with a second terminal's run -- one target file per run of
## this script, so there is nothing to collide with, and the location is fixed
## regardless of cwd.
dir.create(PATHS$untar, recursive = TRUE, showWarnings = FALSE)
say("\n[1] unpack -> %s\n", PATHS$untar)
untar(archive, exdir = PATHS$untar)
files <- list.files(PATHS$untar, recursive = TRUE, full.names = TRUE)
map_file  <- files[grepl("\\.map$", files)]
geno_file <- files[grepl("snp_geno", files, fixed = TRUE)]
stopifnot("expected exactly one .map file" = length(map_file) == 1L,
          "expected exactly one snp_geno file" = length(geno_file) == 1L)
say("    %s\n    %s\n", basename(map_file), basename(geno_file))

## ---- 2. read NEMO's own output --------------------------------------------------
## map_nemo$trait.locus encodes "<type>.<0-based index>" (NEMO's own type/index
## pairing, e.g. "ntrl.412", "quant.3"). GTs' first 5 columns are sample
## metadata (pop, ID, ...); genotype dosages start at column 6. Format verified
## against R/Parse_sim_data.R:315-323, 2026-09-04.
map_nemo <- fread(map_file)
GTs_raw  <- fread(geno_file)
nemo_map <- data.table(marker = map_nemo$trait.locus,
                       do.call(rbind, strsplit(map_nemo$trait.locus, ".", fixed = TRUE)))
setnames(nemo_map, c("V1", "V2"), c("type", "idx"))
nemo_map[, idx := as.numeric(idx) + 1]     # NEMO indexes from 0
sample_info <- GTs_raw[, .(pop, ID)]
GTs_raw <- as.matrix(GTs_raw[, 6:ncol(GTs_raw), with = FALSE])
say("    %d individuals x %d NEMO-native markers ; types: %s\n",
    nrow(GTs_raw), ncol(GTs_raw), paste(sprintf("%s=%d", names(table(nemo_map$type)), table(nemo_map$type)), collapse = ", "))

## ---- 3. join against the reference chromosome map -------------------------------
## rec_map<chr>.rds is the genome MODEL (position, type, allelic_values) shared
## across every replicate of this chromosome -- an input, not a simulation
## output. NEMO's per-type index (ntrl_idx / quanti_idx / delet_idx) selects
## rows of the SAME type within it. "QTN" in the reference map corresponds to
## "quant" in NEMO's own vocabulary.
say("\n[3] joining against reference map: %s\n", basename(recmap_rds))
refmap <- readRDS(recmap_rds)
refmap[, indx := .I]

idx_ntrl  <- refmap[type == "ntrl"][nemo_map[type == "ntrl",  idx], indx]
idx_quant <- refmap[type == "QTN" ][nemo_map[type == "quant", idx], indx]
idx_delet <- refmap[type == "delet"][nemo_map[type == "delet", idx], indx]
refmap[idx_ntrl,  nemo_marker := nemo_map[type == "ntrl",  marker]]
refmap[idx_quant, nemo_marker := nemo_map[type == "quant", marker]]
refmap[idx_delet, nemo_marker := nemo_map[type == "delet", marker]]

map <- refmap[nemo_marker %in% colnames(GTs_raw)]
GTs <- GTs_raw[, map$nemo_marker]
map[, Pos := bp]
map[, marker := paste0("Chr", TARGET_CHR, ":", Pos)]
colnames(GTs) <- map$marker
## [!] SEVERE BUG, FOUND BY CHECKING GROUND TRUTH RATHER THAN TRUSTING A CLEAN
## RUN (2026-09-04). The first version of this stage called `setorder(map, bp)`
## here -- reordering map's ROWS by position without reordering GTs's COLUMNS
## to match. Every map row after that point described a DIFFERENT marker than
## the GTs column at the same position: verified directly, a QTN moved from row
## 705 to row 1498 in map, but GTs column 1498 still held a different marker's
## genotype (MAF 0.028 reported for a marker whose true MAF, looked up by name,
## was 0.075-0.28). This corrupted every marker's map-metadata-to-genotype
## correspondence, not just the QTN it happened to be caught by -- ANY
## downstream decay/clustering/association number from the first run's output
## would have been silently wrong. Fixed by computing the sort order once and
## applying it to both objects together, never one alone.
ord <- order(map$bp)
map <- map[ord]; GTs <- GTs[, ord]
map[, Chr := NULL][, Chr := paste0("Chr", TARGET_CHR)]   # := NULL first: refmap's own
## Chr was numeric, and := on an EXISTING typed column coerces the RHS to match
## rather than replacing the type -- assigning a "ChrN" string silently became
## NA. Deleting the column first before recreating it avoids the coercion.

## duplicate positions (two NEMO loci landing on the same bp -- a real property
## of the simulated recombination map at this resolution, not an error; see
## memory on pos_nemo, not re-litigated here): keep the first, drop the rest,
## consistently in both map and GTs, by POSITION (not by name -- marker strings
## are exactly what collides here, so a name-based operation is the wrong tool).
dup <- duplicated(map$marker)
if (any(dup)) { say("    dropping %d duplicate-position marker(s)\n", sum(dup))
  map <- map[!dup]; GTs <- GTs[, !dup] }
say("    %d markers survive the join (of %d NEMO-native, %d reference)\n",
    nrow(map), ncol(GTs_raw), nrow(refmap))

## Durable guard against this exact failure recurring silently: map's row order
## and GTs's column order must describe the same markers, checked by identity
## (not count) every time this stage runs, not just today.
stopifnot("map/GTs marker correspondence broken" = identical(colnames(GTs), map$marker))

## ---- 4. environment --------------------------------------------------------------
## env_<idx>.txt is a flat 48x48 spatial surface in a custom brace-delimited
## format ("{{v1}}{{v2}}..."), one value per population on a regular grid;
## format verified against R/Parse_sim_data.R:363-372, 2026-09-04.
say("\n[4] environment: %s\n", basename(env_txt))
env_raw  <- scan(env_txt, what = character(), quiet = TRUE)
env_raw  <- strsplit(env_raw, "}{", fixed = TRUE)[[1]]
env_vals <- as.numeric(gsub("}}", "", gsub("{{", "", env_raw, fixed = TRUE), fixed = TRUE))
side <- 48L
env_grid <- data.table(expand.grid(x = 1:side, y = side:1))
env_grid[, pop := .I]
stopifnot("environment surface size does not match the 48x48 grid" = length(env_vals) == nrow(env_grid))
env_grid[, env := env_vals]
env <- env_grid[match(sample_info$pop, env_grid$pop)]
env[, indx := .I]
say("    %d individuals placed on the %dx%d grid ; env range [%.3f, %.3f]\n",
    nrow(env), side, side, min(env$env), max(env$env))

## ---- 5. subsample individuals, MAF filter -----------------------------------------
## keep_inds fixed by SUBSAMPLE_STEP (deterministic sequence, not a random draw
## -- no seed needed for this step).
keep_inds <- seq(1, nrow(GTs), by = SUBSAMPLE_STEP)
GTs <- GTs[keep_inds, ]; env <- env[keep_inds]
say("\n[5] subsampled to %d of %d individuals (step %d)\n",
    length(keep_inds), nrow(GTs_raw), SUBSAMPLE_STEP)

maf <- colSums(GTs) / nrow(GTs) / 2
map[, MAF := pmin(maf, 1 - maf)]
keep_snps <- map$MAF > MAF_KEEP
GTs <- GTs[, keep_snps]; map <- map[keep_snps]
say("    %s of %s markers pass MAF > %.2f\n",
    format(sum(keep_snps), big.mark = ","), format(length(keep_snps), big.mark = ","), MAF_KEEP)

## ---- 6. what this bundle does NOT contain -----------------------------------------
## No ld_w, no emx_p/emx_F, no lfmm_p/lfmm_F, no GRM, no complexity_reduction,
## no chr_type/max_LD_with_QTN/focal_QTN (see 00_config.R's RAW_SCORING note for
## why the last three are excluded rather than carried over from the old
## parser). map$type (ntrl/QTN/delet) and map$allelic_values are the raw truth;
## everything else is 02_bundle.R's job.
map <- map[, .(Chr, Pos, marker, type, allelic_values, MAF)]
map[, true_QTN := type == "QTN"]

## ---- 7. save -----------------------------------------------------------------------
dir.create(PATHS$parsed, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(PATHS$parsed, sprintf("nemo_%s_chr%d_%s_env%d.rds",
                                       TARGET_TAG, TARGET_CHR, TARGET_CELL, TARGET_ENV))
saveRDS(list(GTs = GTs, map = map, env = env,
            source = list(archive = archive, recmap = recmap_rds, env_file = env_txt)), OUT)
write_receipt(STAGE, inputs = unname(INPUTS), params = PARAMS, outputs = OUT)
say("\n[7] wrote %s (%.1f MB)\n    receipt: %s\n", OUT, file.size(OUT) / 1e6, receipt_path(STAGE))
say("\n    Next: R/02_bundle.R picks up from here -- decay, ld_w, stage-1 clustering,\n")
say("    the kinship basis, and both engine scans. Not written yet.\n")
