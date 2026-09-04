## =============================================================================
## module_sim/R_parsing/01_parse_nemo.R
##
## PARSE ONE RAW NEMO SIMULATION REPLICATE INTO A CLEAN, UNSCORED, TWO-CHROMOSOME
## BUNDLE.
##
## Separate from module_sim/R/ on PK's instruction: R_parsing/ turns NEMO's
## native simulator output into the raw substrate this project's analyses
## consume; R/ (02_bundle.R onward, mirroring module_3sp's stage numbers) is
## where decay, ld_w, clustering, kinship and association scans happen. Nothing
## in this file computes a p-value, a kinship, or an LD statistic -- if it does,
## that is a bug in this stage, not a convenience, same contract as
## module_sim/R/00_config.R itself.
##
## [!] CORRECTED 2026-09-04, after the first version ran and PK caught what it
## got wrong. "chr1".."chr10" in NEMO's own filenames are NOT ten genomic
## chromosomes to pool into one panel -- they are ten INDEPENDENT SIMULATION
## REPLICATES. Each replicate simulates TWO chromosomes sharing ONE
## recombination map (PK; confirmed directly: rec_map1.rds contains rows with
## Chr %in% c(1,2), identical bp range for both, Chr 1 carrying 100 potential
## QTN sites in the genome model and Chr 2 carrying 1 -- and rec_map2.rds is a
## genuinely DIFFERENT map, one per replicate, not one per chromosome). The
## first version of this script collapsed both chromosomes into a single
## hardcoded "ChrN" label per replicate, destroying that split entirely --
## every one of the ten replicates LOOKED like it had QTN scattered through a
## single homogeneous chromosome, when the real structure is one QTN-capable
## chromosome and one near-neutral one, per replicate.
##
## A REPLICATE'S BUNDLE THEREFORE ALREADY CONTAINS ITS OWN "NEUTRAL" AND "QTN"
## CHROMOSOME -- there is no cross-replicate pooling to do at this stage.
## Analyses run PER REPLICATE (PK); only the final PR/recall scoring pools
## across all ten (20 chromosomes total). 02_bundle.R therefore reads ONE
## parsed replicate file, not a loop over several.
##
## WHAT THIS REPLACES. R/Parse_sim_data.R (LDscnR-paper root, 654 lines) reads
## the same NEMO format but is unreproducible in the ways this project spent
## this week fixing elsewhere: a hardcoded "BATCH SELECTOR" line hand-edited
## between runs, relative paths, no seed, no receipt, and (checked directly)
## a chr_type assignment of `ifelse(Chr==1,"QTN","ntrl")` that collapses the
## same way the first version of this script did. The NEMO-format-reading logic
## in that script is sound and is reused here (the .map/.snp_geno layout, the
## reference-map index matching, the environment-file format); the surrounding
## practice, including that collapse, is not.
##
## SCOPE, DELIBERATELY NARROW. One replicate: TARGET_REP/TARGET_CELL/
## TARGET_ENV/TARGET_TAG below, matching module_sim/R/00_config.R's
## CELLS/TAGS/ENVS (all currently pinned to a single slice). Only "adapt" runs
## are in scope -- NOT the paired "rxpn" archives also present in
## Nemo_out_nobgs/bgs, which are a range-expansion continuation seeded from an
## adapt replicate's end state (confirmed via its own log: root_dir
## .../range_expansion, source_pop .../nobgs_chr<rep>_..., logfile
## range-xpn-nobgs.log), not a neutral/QTN pairing and not yet in scope (PK).
## =============================================================================
suppressMessages({library(data.table)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim"), "R", "00_config.R"))
STAGE <- "01_parse_nemo"
say("=== %s ===\n\n", STAGE)

TARGET_TAG  <- Sys.getenv("SIM_TAG",  TAGS[1])                    # "nobgs"
TARGET_CELL <- Sys.getenv("SIM_CELL", CELLS[1])                   # "V2_c1"
TARGET_ENV  <- as.integer(Sys.getenv("SIM_ENV",  ENVS[1]))       # 1L
TARGET_REP  <- as.integer(Sys.getenv("SIM_REP",  REPS[1]))       # 1L -- one replicate per run; SIM_REP overrides for looping
combo_id <- sprintf("%s_%s_rep%d_env%d", TARGET_TAG, TARGET_CELL, TARGET_REP, TARGET_ENV)   # keys the receipt subdirectory (see 00_config.R)

raw_dir <- if (TARGET_TAG == "nobgs") PATHS$raw_nemo_nobgs else PATHS$raw_nemo_bgs5
archive <- file.path(raw_dir, sprintf("adapt_%s_chr%d_%s_env%d.tgz",
                                      TARGET_TAG, TARGET_REP, TARGET_CELL, TARGET_ENV))
recmap_rds <- file.path(PATHS$raw_recmap_dir, sprintf("rec_map%d.rds", TARGET_REP))
env_txt    <- file.path(PATHS$raw_env_dir, sprintf("env_%d.txt", TARGET_ENV))

INPUTS <- c(archive = archive, recmap = recmap_rds, env = env_txt)
PARAMS <- list(tag = TARGET_TAG, cell = TARGET_CELL, env = TARGET_ENV, rep = TARGET_REP,
               subsample_step = SUBSAMPLE_STEP, maf_keep = MAF_KEEP)

say("[0] target: tag=%s cell=%s env=%d rep=%d\n", TARGET_TAG, TARGET_CELL, TARGET_ENV, TARGET_REP)
for (nm in names(INPUTS)) say("    %-8s %s  (%s)\n", nm, INPUTS[[nm]],
                              if (file.exists(INPUTS[[nm]])) "exists" else "MISSING")
if (!stage_stale(STAGE, unname(INPUTS), PARAMS, target = combo_id) && !nzchar(Sys.getenv("FORCE"))) {
  say("\nNothing to do. Set FORCE=1 to rerun anyway.\n"); quit(save = "no")
}
if (any(!file.exists(INPUTS))) stop("missing input(s): ",
    paste(names(INPUTS)[!file.exists(INPUTS)], collapse = ", "))

## ---- 1. unpack the archive ----------------------------------------------------
## [!] PER-COMBINATION DIRECTORY, NOT THE SHARED PATHS$untar. First fix
## (2026-09-04, running a second cell for the first time): cleared
## PATHS$untar before every unpack, because a second combination's extraction
## had landed ALONGSIDE the first's leftover files under the one shared path,
## caught immediately by the stopifnot below ("expected exactly one .map
## file", two present) rather than silently reading the wrong one. That fix
## was correct for SEQUENTIAL runs and is NOT SAFE for concurrent ones (PK:
## "4 on mini2 and 4 here") -- two processes racing to unlink+recreate+extract
## into the same directory at once would corrupt each other's extraction
## rather than merely leave stale files. Fixed properly now: a subdirectory
## keyed by (tag, cell, rep), so concurrent combinations never share a path at
## all, regardless of timing.
untar_dir <- file.path(PATHS$untar, combo_id)
unlink(untar_dir, recursive = TRUE)
dir.create(untar_dir, recursive = TRUE, showWarnings = FALSE)
say("\n[1] unpack -> %s\n", untar_dir)
untar(archive, exdir = untar_dir)
files <- list.files(untar_dir, recursive = TRUE, full.names = TRUE)
## [!] Some bgs archives bundle more than one env's output together (e.g.
## adapt_bgs_chr2_V0.5_c1_env1.tgz also contains the full env10 run) --
## found 2026-09-04 launching the grid, 5 of 70 bgs archives affected, 0 of
## 70 nobgs. Filter to TARGET_ENV explicitly rather than assuming "env1.tgz"
## means only-env1 inside; env%d(?!\\d) avoids env1 matching env10's files.
env_pat <- sprintf("env%d(?!\\d)", TARGET_ENV)
files <- files[grepl(env_pat, files, perl = TRUE)]
map_file  <- files[grepl("\\.map$", files)]
geno_file <- files[grepl("snp_geno", files, fixed = TRUE)]
stopifnot("expected exactly one .map file" = length(map_file) == 1L,
          "expected exactly one snp_geno file" = length(geno_file) == 1L)
say("    %s\n    %s\n", basename(map_file), basename(geno_file))

## ---- 2. read NEMO's own output --------------------------------------------------
## map_nemo$trait.locus encodes "<type>.<0-based index>" (NEMO's own type/index
## pairing, e.g. "ntrl.412", "quant.3"). GTs' first 5 columns are sample
## metadata (pop, ID, ...); genotype dosages start at column 6.
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
unlink(untar_dir, recursive = TRUE)   # fully read into memory above; nothing on disk needed past here

## ---- 3. join against the reference map -- BOTH chromosomes, index space shared
## rec_map<rep>.rds is the genome MODEL (position, type, allelic_values) for
## this replicate's two chromosomes, sharing one recombination map -- an input,
## not a simulation output. NEMO's per-type index (ntrl_idx / quanti_idx /
## delet_idx) is a FLAT counter across the whole simulated genome (not reset
## per chromosome), which is why the join below filters refmap by type only,
## never by Chr, and still lands each marker on its correct chromosome: refmap's
## own row order already spans Chr 1 then Chr 2 within each type, and NEMO's
## index respects that same order. "QTN" in the reference map corresponds to
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
## [!] Chr IS PRESERVED FROM refmap, not collapsed to the replicate number. The
## first version of this script did `map[, Chr := paste0("Chr", TARGET_CHR)]`
## here -- overwriting the genuine per-marker Chr (1 or 2) that refmap already
## carries with a single label derived from the REPLICATE index, destroying the
## two-chromosome structure entirely. Fixed: keep refmap's own Chr, just
## reformat it ("1" -> "Chr1") without touching which markers get which value.
map[, Chr := paste0("Chr", Chr)]
map[, marker := paste0(Chr, ":", Pos)]     # per-row Chr, not a replicate-wide constant
colnames(GTs) <- map$marker

## Sort by (Chr, Pos), reordering map and GTs together -- never one alone. This
## is exactly the bug found and fixed in the first version (a QTN silently
## matched to a different marker's genotype after an unsynchronised sort); the
## fix generalises the same way now that there are two chromosomes to keep
## ordered rather than one.
ord <- order(map$Chr, map$bp)
map <- map[ord]; GTs <- GTs[, ord]

## duplicate positions (two NEMO loci landing on the same Chr:bp -- a real
## property of the simulated recombination map at this resolution, not an
## error; see memory on pos_nemo, not re-litigated here): keep the first, drop
## the rest, consistently in both map and GTs, by POSITION.
dup <- duplicated(map$marker)
if (any(dup)) { say("    dropping %d duplicate-position marker(s)\n", sum(dup))
  map <- map[!dup]; GTs <- GTs[, !dup] }
say("    %d markers survive the join (of %d NEMO-native, %d reference)\n",
    nrow(map), ncol(GTs_raw), nrow(refmap))
say("    by chromosome: %s\n", paste(sprintf("%s=%d (QTN=%d)", names(table(map$Chr)),
    table(map$Chr), sapply(split(map$type, map$Chr), function(t) sum(t == "QTN"))), collapse = ", "))

## Durable guard against the exact failure that motivated the fix above.
stopifnot("map/GTs marker correspondence broken" = identical(colnames(GTs), map$marker))

## ---- 4. environment --------------------------------------------------------------
## env_<idx>.txt is a flat 48x48 spatial surface in a custom brace-delimited
## format ("{{v1}}{{v2}}..."), one value per population on a regular grid.
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
keep_inds <- seq(1, nrow(GTs), by = SUBSAMPLE_STEP)
GTs <- GTs[keep_inds, ]; env <- env[keep_inds]
say("\n[5] subsampled to %d of %d individuals (step %d)\n",
    length(keep_inds), nrow(GTs_raw), SUBSAMPLE_STEP)

maf <- colSums(GTs) / nrow(GTs) / 2
map[, MAF := pmin(maf, 1 - maf)]
keep_snps <- map$MAF > MAF_KEEP
GTs <- GTs[, keep_snps]; map <- map[keep_snps]
say("    %s of %s markers pass MAF > %.2f (by chromosome: %s)\n",
    format(sum(keep_snps), big.mark = ","), format(length(keep_snps), big.mark = ","), MAF_KEEP,
    paste(sprintf("%s=%d (QTN=%d)", names(table(map$Chr)), table(map$Chr),
          sapply(split(map$type, map$Chr), function(t) sum(t == "QTN"))), collapse = ", "))

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
OUT <- file.path(PATHS$parsed, sprintf("nemo_%s_rep%d_%s_env%d.rds",
                                       TARGET_TAG, TARGET_REP, TARGET_CELL, TARGET_ENV))
saveRDS(list(GTs = GTs, map = map, env = env,
            source = list(archive = archive, recmap = recmap_rds, env_file = env_txt)), OUT)
write_receipt(STAGE, inputs = unname(INPUTS), params = PARAMS, outputs = OUT, target = combo_id)
say("\n[7] wrote %s (%.1f MB)\n    receipt: %s\n", OUT, file.size(OUT) / 1e6, receipt_path(STAGE, combo_id))
say("\n    Next: R/02_bundle.R reads this ONE replicate's file (both chromosomes\n")
say("    already in it) -- decay, ld_w, stage-1 clustering, the kinship basis.\n")
