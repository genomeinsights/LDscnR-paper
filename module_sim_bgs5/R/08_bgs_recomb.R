## module_sim/R/08_bgs_recomb.R
##
## How BGS affects genetic differentiation AS A FUNCTION OF LOCAL
## RECOMBINATION RATE (PK, following up on fig_bgs_effect.R's genome-wide
## contrast: "some estimates that measure the difference between low and
## high recombination regions to see how bgs affects patterns of genetic
## differentiation as a function of recombination rate"). Background
## selection's classic signature is a STRONGER effect in low-recombination
## regions (linked selection has a bigger footprint there) that fades as
## recombination rises -- this stage computes Fst separately within
## recombination-rate bins so that signature is visible (or not) directly,
## rather than averaged away in one genome-wide number.
##
## rec_rate is NOT in R/02_bundle.R's map (dropped after the reference-map
## join in R_parsing/01_parse_nemo.R) -- rejoined here from the same
## rec_map<rep>.rds both tags share (recombination map is fixed per REP,
## identical for bgs and nobgs at that rep; see 00_config.R's ENVS comment).
## Bin edges are cut from the FULL rec_map (all 1e6 reference positions for
## this rep), not from whichever markers happen to survive this bundle's own
## MAF filter -- keeps bin boundaries identical between bgs and nobgs (and
## across cells/envs) at the same rep, which per-bundle quantiles would not
## guarantee.
##
## rec_rate is markedly zero-inflated (~33% of reference positions are
## exactly 0 -- recombination coldspots) -- binned as its own category
## rather than folded into a quantile that would otherwise span 0 to a
## positive value: {0, then tertiles of the nonzero distribution}.
suppressMessages({library(data.table); library(LDscnR); library(SNPRelate)})
devtools::load_all("~/gitlab/LDscnR")
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim"), "R", "00_config.R"))
STAGE <- "08_bgs_recomb"
say("=== %s ===\n\n", STAGE)
invisible(check_ldscnr())

TARGET_TAG <- Sys.getenv("SIM_TAG", TAGS[1]); TARGET_CELL <- Sys.getenv("SIM_CELL", CELLS[1]); TARGET_ENV <- as.integer(Sys.getenv("SIM_ENV", ENVS[1])); TARGET_REP <- as.integer(Sys.getenv("SIM_REP", REPS[1]))
combo_id <- sprintf("%s_%s_rep%d_env%d", TARGET_TAG, TARGET_CELL, TARGET_REP, TARGET_ENV)

BUNDLE_PATH <- file.path(PATHS$out, "02_bundle",
  sprintf("bundle_%s_rep%d_%s_env%d.rds", TARGET_TAG, TARGET_REP, TARGET_CELL, TARGET_ENV))
RECMAP_PATH <- file.path(PATHS$raw_recmap_dir, sprintf("rec_map%d.rds", TARGET_REP))
if (!file.exists(BUNDLE_PATH)) stop("R/02_bundle.R has not produced: ", basename(BUNDLE_PATH))
if (!file.exists(RECMAP_PATH)) stop("missing recombination map: ", RECMAP_PATH)
b <- readRDS(BUNDLE_PATH)
GTs <- b$GTs; map <- b$map; env <- b$env

INPUTS <- c(BUNDLE_PATH, RECMAP_PATH)
PARAMS <- list(bin_scheme = "zero_plus_tertiles_of_nonzero", fst_method = "W&C84")
if (!stage_stale(STAGE, INPUTS, PARAMS, target = combo_id) && !nzchar(Sys.getenv("FORCE"))) {
  say("\nNothing to do. Set FORCE=1 to rerun anyway.\n"); quit(save = "no")
}

## ---- 1. bin edges from the FULL, shared rec_map (rep-level, tag-independent) --
recmap <- readRDS(RECMAP_PATH)
say("[1] recombination-rate bins from rec_map%d.rds (%d reference positions)\n", TARGET_REP, nrow(recmap))
nz <- recmap$rec_rate[!is.na(recmap$rec_rate) & recmap$rec_rate > 0]
breaks <- c(-Inf, 0, quantile(nz, c(1/3, 2/3), na.rm = TRUE), Inf)
bin_labels <- c("zero", "low", "medium", "high")
say("    edges: zero | (0, %.3f] low | (%.3f, %.3f] medium | (%.3f, Inf) high\n",
    breaks[3], breaks[3], breaks[4], breaks[4])

## ---- 2. join rec_rate onto this bundle's markers, by (Chr, bp) ----------------
recmap[, Chr := paste0("Chr", Chr)]
setkey(recmap, Chr, bp)
map_rt <- copy(map)
lookup_chr <- map_rt$Chr; lookup_pos <- map_rt$Pos
map_rt[, rec_rate := recmap[.(lookup_chr, lookup_pos), rec_rate, on = c("Chr", "bp"), mult = "first"]]
say("[2] rec_rate joined for %d of %d markers (%d unmatched -- dropped from binning)\n",
    sum(!is.na(map_rt$rec_rate)), nrow(map_rt), sum(is.na(map_rt$rec_rate)))
map_rt <- map_rt[!is.na(rec_rate)]
map_rt[, recomb_bin := cut(rec_rate, breaks = breaks, labels = bin_labels, include.lowest = TRUE)]
say("    markers per bin: %s\n", paste(sprintf("%s=%d", names(table(map_rt$recomb_bin)), table(map_rt$recomb_bin)), collapse = ", "))

## ---- 3. Fst (W&C84) within each recombination-rate bin ------------------------
say("[3] Fst per recombination-rate bin, %d populations\n", uniqueN(env$pop))
pop <- factor(env$pop)
fst_by_bin <- rbindlist(lapply(bin_labels, function(bn) {
  markers <- map_rt[recomb_bin == bn, marker]
  if (length(markers) < 10) return(data.table(recomb_bin = bn, n_markers = length(markers), Fst = NA_real_))
  gds_path <- tempfile(fileext = ".gds")
  gds <- create_gds_from_geno(geno = GTs[, markers, drop = FALSE], map = map[marker %in% markers], gds_path)
  on.exit({ try(snpgdsClose(gds), silent = TRUE); unlink(gds_path) }, add = TRUE)
  fst <- snpgdsFst(gds, population = pop, method = "W&C84", autosome.only = FALSE, verbose = FALSE)
  data.table(recomb_bin = bn, n_markers = length(markers), Fst = fst$MeanFst)
}))
fst_by_bin[, recomb_bin := factor(recomb_bin, levels = bin_labels)]
setorder(fst_by_bin, recomb_bin)
say("    %s\n", paste(sprintf("%s: Fst=%.4f (n=%d)", fst_by_bin$recomb_bin, fst_by_bin$Fst, fst_by_bin$n_markers), collapse = " | "))

summary_row <- cbind(data.table(tag = TARGET_TAG, cell = TARGET_CELL, rep = TARGET_REP, env = TARGET_ENV), fst_by_bin)

OUT <- file.path(stage_dir(STAGE), sprintf("bgsrecomb_%s_rep%d_%s_env%d.rds", TARGET_TAG, TARGET_REP, TARGET_CELL, TARGET_ENV))
dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)
saveRDS(summary_row, OUT)
write_receipt(STAGE, inputs = INPUTS, params = PARAMS, outputs = OUT, target = combo_id)
say("\n[4] wrote %s\n    receipt: %s\n", OUT, receipt_path(STAGE, combo_id))
