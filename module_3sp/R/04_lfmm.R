## =============================================================================
## module_3sp/R/04_lfmm.R
##
## THE LFMM ARM. Simes only -- LFMM's p-values are precomputed (LFMM_SOURCE =
## "inherit" in 00_config.R: no script in any repository computes them, they
## come from lfmm_F.rds, dated September 2025), so the consensus test, which
## refits the association on a NEW genotype column per cluster, is not
## available for this engine. That asymmetry is a property of what the engine
## ships, not a choice between statistics, and it is stated on the figure
## rather than hidden -- matching R_3sp/128_joint_manhattan.R's convention.
##
## [!] NO PERMUTATION NULL FOR LFMM, AND THIS IS DELIBERATE, NOT AN OMISSION.
## ld_outlier_perm() needs a fresh p-value per surrogate phenotype, which for
## EMMAX costs one emmax_fast() call (milliseconds once set up). LFMM has no
## such fast-surrogate path in this codebase: computing a genuinely permuted
## LFMM p-value means rerunning the external LFMM program per surrogate,
## which nothing here does. Reporting a permutation p-value anyway -- e.g. by
## reusing EMMAX's null -- would silently claim a calibration LFMM's own
## engine never demonstrated. The rotation null against real EcoPeaks is
## still run: it needs only the assembled regions' genomic positions, not a
## new p-value per draw.
## =============================================================================
suppressMessages({library(data.table); library(LDscnR)})
devtools::load_all("~/gitlab/LDscnR")
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_3sp"), "R", "00_config.R"))
STAGE <- "04_lfmm"
say("=== %s ===\n\n", STAGE)
invisible(check_ldscnr())

BUNDLE_PATH <- file.path(PATHS$out, "02_bundle", "bundle.rds")
b <- readRDS(BUNDLE_PATH)
GTs <- b$GTs; map <- b$map; stage1 <- b$stage1

## ---- align the inherited LFMM F-values to this bundle's marker set --------
## lfmm_F.rds is over the FULL PRE-MAF map (881,786 rows, in map_3sp's own row
## order). 02_bundle.R filtered map_3sp by maf > MAF_KEEP and kept row order
## within the kept set, so recomputing the SAME boolean mask against the raw
## map_3sp and indexing lfmm_F_full by it recovers exact alignment to `map`.
if (!identical(LFMM_SOURCE, "inherit")) stop("LFMM_SOURCE is not \"inherit\" -- this stage has nothing to do.")
e <- new.env(); load(PATHS$raw_3sp, envir = e)
keep <- e$map_3sp$maf > MAF_KEEP
lfmm_F_full <- readRDS(PATHS$lfmm_F)
stopifnot(length(lfmm_F_full) == length(keep))
lfmm_F <- lfmm_F_full[keep]
stopifnot(length(lfmm_F) == nrow(map))
## Length matching proves nothing beyond same COUNT -- the failure class an independent
## audit flagged here (ldscnr-26/fe, 2026-09-04): 02_bundle.R's `map` and this stage's
## re-masked `e$map_3sp` are built from the same source file with the same filter
## expression, so they SHOULD be identical, but that was never actually checked. This does:
## marker-by-marker identity, not just row count.
stopifnot(identical(e$map_3sp$marker[keep], map$marker))
n <- nrow(GTs)
lfmm_p <- stats::pf(lfmm_F, 1, n - 2, lower.tail = FALSE)
say("[0] LFMM F aligned: %s values, min p = %.3g\n", format(length(lfmm_p), big.mark=","), min(lfmm_p))

LIFT <- path.expand("~/gitlab/LDscnR-paper/kingman2021/data/liftover")
ROMAN <- c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV",
           "XV","XVI","XVII","XVIII","XIX","XX","XXI")
.rd <- function(f) { x <- fread(file.path(LIFT, f), header=FALSE,
                                col.names=c("chr","start","end","pv"))
  x[, chr_num := match(sub("^chr","",chr), ROMAN)]
  x[!is.na(chr_num), .(chr = paste0("Chr", chr_num), start, end)] }
KG <- unique(rbindlist(lapply(ECOPEAK_BEDS, .rd)))
LEN <- map[, .(len = max(Pos)), by = Chr]; setnames(LEN, "Chr", "chr")

INPUTS <- c(file.path(PATHS$out, "02_bundle", "_receipt.rds"), PATHS$lfmm_F)
PARAMS <- list(size_floor = SIZE_FLOOR, alpha = ALPHA, region_assembly = REGION_ASSEMBLY,
               n_rotations = N_ROTATIONS, rotation_scheme = ROTATION_SCHEME)
if (!stage_stale(STAGE, INPUTS, PARAMS) && !nzchar(Sys.getenv("FORCE"))) {
  say("\nNothing to do. Set FORCE=1 to rerun anyway.\n"); quit(save = "no")
}
dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)

say("\n[1] LFMM arm -- Simes over the inherited per-marker p-values\n")
t0 <- Sys.time()
test_lfmm <- ld_outlier_test(stage1, map, lfmm_p, statistic = "simes", size_floor = SIZE_FLOOR,
                             alpha = ALPHA, assembly = "stage2_discovered", GTs = GTs,
                             LD_decay = b$LD_decay, score_threshold = REGION_ASSEMBLY$score_threshold,
                             distance_threshold = REGION_ASSEMBLY$distance_threshold)
print(test_lfmm)
rot_lfmm <- ld_region_rotation(test_lfmm$regions, KG, LEN, scheme = ROTATION_SCHEME,
                               n_rotations = N_ROTATIONS, seed = 1L)
print(rot_lfmm)
say("    %.1f min\n", as.numeric(difftime(Sys.time(), t0, units="mins")))

## ---- agreement against EMMAX (both statistics, from 03_EMMAX.R) ------------
scan_path <- file.path(PATHS$out, "03_EMMAX", "scan.rds")
if (file.exists(scan_path)) {
  sc <- readRDS(scan_path)
  sig_lfmm <- test_lfmm$units[significant == TRUE]$unit_id
  for (arm in c("consensus", "simes")) {
    sig_emx <- sc[[arm]]$test$units[significant == TRUE]$unit_id
    jac <- length(intersect(sig_lfmm, sig_emx)) / length(union(sig_lfmm, sig_emx))
    say("[2] LFMM vs EMMAX-%s: LFMM %d, EMMAX %d, shared %d, Jaccard %.3f\n",
        arm, length(sig_lfmm), length(sig_emx), length(intersect(sig_lfmm, sig_emx)), jac)
  }
} else say("[2] 03_EMMAX.R output not found -- skipping the EMMAX comparison\n")

OUT <- file.path(stage_dir(STAGE), "lfmm.rds")
saveRDS(list(test = test_lfmm, rotation = rot_lfmm, lfmm_p = lfmm_p,
            ecopeaks = KG, chrom_lengths = LEN), OUT)
write_receipt(STAGE, inputs = INPUTS, params = PARAMS, outputs = OUT)
say("\n[3] wrote %s\n    receipt: %s\n", OUT, receipt_path(STAGE))
