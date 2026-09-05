## module_sim/R/10_bgs_windows.R
##
## Per-(tag,cell,rep,env), per-500kb-window SEGREGATING SITE COUNT (n_snp),
## for the paired BGS-vs-recombination-rate estimator PK pointed at:
## /Volumes/Nemo/Nemo_sim/bgs_effect_newmaps/measure_bgs_effect.R --
## B_obs = n_snp(bgs) / n_snp(nobgs), same window, paired at (cell,rep,env),
## binned by local recombination rate.
##
## WHY NOT REUSE R/02_bundle.R's ALREADY-PARSED BUNDLE: 01_parse_nemo.R
## applies MAF_KEEP = 0.1 (00_config.R, flagged there as inherited-not-
## verified) before anything downstream -- every existing bundle already
## excludes exactly the rare variants BGS purges first. measure_bgs_effect.R's
## own comment is explicit about why that matters: "BGS removes rare variants
## first, so survivors skew common ... A count of segregating sites has no
## such bias" -- true only if that count is taken BEFORE any MAF filter.
## R/04_score.R's/R/06_popgen_summary.R's Fst-based BGS estimates (near-zero
## effect) are consistent with this: Fst is a ratio (Hs/Ht) that can stay flat
## even as absolute diversity drops, and was computed on the same
## already-MAF-filtered markers regardless. This stage re-reads the RAW
## archive independently (same untar/join logic as 01_parse_nemo.R, same
## individual subsampling for a fair like-for-like comparison, but NO MAF
## filter) rather than modifying 01_parse_nemo.R's contract for every other
## stage that depends on it.
##
## Windowing/rec_rate: identical scheme to measure_bgs_effect.R -- win =
## floor(bp/5e5)+1 by Chr, rec_rate = mean(rec_rate) over ALL reference
## positions in that window (from rec_map<rep>.rds, not from whichever
## markers happen to survive in this tag's own data) -- both tags at the
## same rep get IDENTICAL window boundaries and rec_rate, since both share
## one recombination map.
suppressMessages({library(data.table)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim"), "R", "00_config.R"))
STAGE <- "10_bgs_windows"
say("=== %s ===\n\n", STAGE)

TARGET_TAG  <- Sys.getenv("SIM_TAG",  TAGS[1])
TARGET_CELL <- Sys.getenv("SIM_CELL", CELLS[1])
TARGET_ENV  <- as.integer(Sys.getenv("SIM_ENV",  ENVS[1]))
TARGET_REP  <- as.integer(Sys.getenv("SIM_REP",  REPS[1]))
combo_id <- sprintf("%s_%s_rep%d_env%d", TARGET_TAG, TARGET_CELL, TARGET_REP, TARGET_ENV)

raw_dir <- if (TARGET_TAG == "nobgs") PATHS$raw_nemo_nobgs else PATHS$raw_nemo_bgs5
archive <- file.path(raw_dir, sprintf("adapt_%s_chr%d_%s_env%d.tgz",
                                      TARGET_TAG, TARGET_REP, TARGET_CELL, TARGET_ENV))
recmap_rds <- file.path(PATHS$raw_recmap_dir, sprintf("rec_map%d.rds", TARGET_REP))

INPUTS <- c(archive = archive, recmap = recmap_rds)
PARAMS <- list(tag = TARGET_TAG, cell = TARGET_CELL, env = TARGET_ENV, rep = TARGET_REP,
               subsample_step = SUBSAMPLE_STEP, maf_filter = "none", win_bp = 5e5)
if (!stage_stale(STAGE, unname(INPUTS), PARAMS, target = combo_id) && !nzchar(Sys.getenv("FORCE"))) {
  say("\nNothing to do. Set FORCE=1 to rerun anyway.\n"); quit(save = "no")
}
if (any(!file.exists(INPUTS))) stop("missing input(s): ", paste(names(INPUTS)[!file.exists(INPUTS)], collapse = ", "))

## ---- 1. unpack (same per-combination scratch dir / env-bundling fix as 01_parse_nemo.R) --
untar_dir <- file.path(PATHS$untar, paste0("bgswin_", combo_id))
unlink(untar_dir, recursive = TRUE)
dir.create(untar_dir, recursive = TRUE, showWarnings = FALSE)
say("[1] unpack -> %s\n", untar_dir)
untar(archive, exdir = untar_dir)
files <- list.files(untar_dir, recursive = TRUE, full.names = TRUE)
env_pat <- sprintf("env%d(?!\\d)", TARGET_ENV)
files <- files[grepl(env_pat, basename(files), perl = TRUE)]
map_file  <- files[grepl("\\.map$", files)]
geno_file <- files[grepl("snp_geno", files, fixed = TRUE)]
stopifnot("expected exactly one .map file" = length(map_file) == 1L,
          "expected exactly one snp_geno file" = length(geno_file) == 1L)

## ---- 2. read + join against the reference map (same logic as 01_parse_nemo.R) --
map_nemo <- fread(map_file)
GTs_raw  <- fread(geno_file)
nemo_map <- data.table(marker = map_nemo$trait.locus,
                       do.call(rbind, strsplit(map_nemo$trait.locus, ".", fixed = TRUE)))
setnames(nemo_map, c("V1", "V2"), c("type", "idx"))
nemo_map[, idx := as.numeric(idx) + 1]
GTs_raw <- as.matrix(GTs_raw[, 6:ncol(GTs_raw), with = FALSE])
unlink(untar_dir, recursive = TRUE)

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
map[, Chr := paste0("Chr", Chr)]
map[, marker := paste0(Chr, ":", Pos)]
colnames(GTs) <- map$marker
ord <- order(map$Chr, map$bp)
map <- map[ord]; GTs <- GTs[, ord]
dup <- duplicated(map$marker)
if (any(dup)) { map <- map[!dup]; GTs <- GTs[, !dup] }
say("[2] %d markers joined (%d NEMO-native)\n", nrow(map), ncol(GTs_raw))

## ---- 3. subsample individuals -- SAME step as 01_parse_nemo.R, NO MAF filter --
keep_inds <- seq(1, nrow(GTs), by = SUBSAMPLE_STEP)
GTs <- GTs[keep_inds, , drop = FALSE]
say("[3] subsampled to %d of %d individuals (step %d) -- no MAF filter\n",
    length(keep_inds), nrow(GTs_raw), SUBSAMPLE_STEP)

## ---- 4. per-window segregating-site count -------------------------------------
WIN_BP <- 5e5
map[, win := floor(Pos / WIN_BP) + 1L]
ac <- colSums(GTs)   # allele count, 0..2N
segregating <- ac > 0 & ac < 2 * nrow(GTs)
say("[4] %d of %d markers segregating in the %d-individual sample\n",
    sum(segregating), length(segregating), nrow(GTs))
n_snp_by_win <- data.table(Chr = map$Chr, win = map$win, segregating = segregating)[
  , .(n_snp = sum(segregating), n_markers = .N), by = .(Chr, win)]

## ---- 5. mean recombination rate per window, from ALL reference positions -----
## Tag-independent: both bgs and nobgs at this rep get this exact table.
refmap[, Chr := paste0("Chr", Chr)][, win := floor(bp / WIN_BP) + 1L]
rate_by_win <- refmap[, .(rate = mean(rec_rate, na.rm = TRUE)), by = .(Chr, win)]
n_snp_by_win <- merge(n_snp_by_win, rate_by_win, by = c("Chr", "win"), all.x = TRUE)
say("[5] %d windows, rate joined for %d\n", nrow(n_snp_by_win), sum(!is.na(n_snp_by_win$rate)))

summary_row <- cbind(data.table(tag = TARGET_TAG, cell = TARGET_CELL, rep = TARGET_REP, env = TARGET_ENV), n_snp_by_win)

OUT <- file.path(stage_dir(STAGE), sprintf("bgswin_%s_rep%d_%s_env%d.rds", TARGET_TAG, TARGET_REP, TARGET_CELL, TARGET_ENV))
dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)
saveRDS(summary_row, OUT)
write_receipt(STAGE, inputs = INPUTS, params = PARAMS, outputs = OUT, target = combo_id)
say("\n[6] wrote %s\n    receipt: %s\n", OUT, receipt_path(STAGE, combo_id))
