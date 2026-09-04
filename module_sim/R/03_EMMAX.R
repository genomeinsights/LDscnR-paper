## =============================================================================
## module_sim/R/03_EMMAX.R
##
## THE EMMAX ARM, for one replicate's bundle. Compute p-values, test stage-1
## clusters, assemble regions -- mirrors module_3sp/R/03_EMMAX.R's structure and
## its use of the package's stage-1-cluster outlier API (ld_unit_matrix /
## ld_outlier_test), pinned via LDSCNR_PIN in 00_config.R.
##
## THREE DELIBERATE DEPARTURES FROM module_3sp/R/03_EMMAX.R, each because the
## thing it mirrors does not exist yet on this side rather than because it does
## not apply:
##
##   1. NO EcoPeak ROTATION NULL. There is no external validation set for
##      simulated data -- there is KNOWN TRUTH instead (map$true_QTN). Section
##      3 below is a per-replicate ground-truth sanity check, not a
##      substitute for one.
##   2. NO PERMUTATION NULL. 00_config.R's NULLS section (7) marks this
##      explicitly TBD -- the naive MVN/latent-factor null this project
##      established is anti-conservative, and a structure-aware replacement for
##      sim phenotype permutation has not been designed. Running the observed
##      scan without a null is a real gap, not an oversight; it is stated here
##      so it cannot be mistaken for one later.
##   3. THE RESPONSE IS env$env (continuous), not a binary ecotype. This is a
##      genotype-environment association scan, matching what these simulations
##      were built to test, not a case-control design.
##
## SCOPE: one replicate (REPS[1] in 00_config.R, currently 1L). PK: "all
## analyses run per replicate; only at the end do we pool the ten simulations
## (20 chromosomes) to score PR" -- that pooling is a later, separate stage,
## not here and not per-replicate. This stage's job is the per-replicate scan
## only.
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(digest)})
devtools::load_all("~/gitlab/LDscnR")
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim"), "R", "00_config.R"))
STAGE <- "03_EMMAX"
say("=== %s ===\n\n", STAGE)
invisible(check_ldscnr())

TARGET_TAG <- TAGS[1]; TARGET_CELL <- CELLS[1]; TARGET_ENV <- ENVS[1]; TARGET_REP <- REPS[1]
BUNDLE_PATH <- file.path(PATHS$out, "02_bundle",
  sprintf("bundle_%s_rep%d_%s_env%d.rds", TARGET_TAG, TARGET_REP, TARGET_CELL, TARGET_ENV))
if (!file.exists(BUNDLE_PATH)) stop("R/02_bundle.R has not produced: ", basename(BUNDLE_PATH))
b <- readRDS(BUNDLE_PATH)
GTs <- b$GTs; map <- b$map; env <- b$env; stage1 <- b$stage1; GRM <- b$GRM
y <- env$env
say("[0] bundle: %d individuals x %s markers ; %d true QTN ; %s stage-1 units at floor %d\n",
    nrow(GTs), format(ncol(GTs), big.mark = ","), sum(map$true_QTN),
    format(sum({cl <- as.data.table(stage1$clusters); nl <- if ("n_loci" %in% names(cl)) cl$n_loci else cl$n_snps;
                sum(nl >= SIZE_FLOOR)}), big.mark = ","), SIZE_FLOOR)

INPUTS <- BUNDLE_PATH
PARAMS <- list(size_floor = SIZE_FLOOR, alpha = ALPHA, region_assembly = REGION_ASSEMBLY,
               statistics = STATISTICS, unit_repr = UNIT_REPR)
if (!stage_stale(STAGE, INPUTS, PARAMS) && !nzchar(Sys.getenv("FORCE"))) {
  say("\nNothing to do. Set FORCE=1 to rerun anyway.\n"); quit(save = "no")
}
dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)

RESULTS <- list()

## ---- 1. consensus arm (statistic = "unit") ----------------------------------
say("[1] CONSENSUS arm -- ld_unit_matrix(repr = \"%s\") + emmax_fast\n", UNIT_REPR)
t0 <- Sys.time()
um <- ld_unit_matrix(GTs, stage1, map, size_floor = SIZE_FLOOR, repr = UNIT_REPR)
Pu <- emmax_setup(um, GRM)
pu_obs <- emmax_fast(Pu, y)

test_con <- ld_outlier_test(stage1, map, pu_obs, statistic = "unit", size_floor = SIZE_FLOOR,
                            alpha = ALPHA, assembly = "stage2_discovered", GTs = GTs,
                            LD_decay = b$LD_decay, score_threshold = REGION_ASSEMBLY$score_threshold,
                            distance_threshold = REGION_ASSEMBLY$distance_threshold)
print(test_con)
say("    %.1f min\n\n", as.numeric(difftime(Sys.time(), t0, units = "mins")))
RESULTS$consensus <- list(test = test_con)

## ---- 2. Simes arm (statistic = "simes") -------------------------------------
say("[2] SIMES arm -- direct marker-level scan aggregated per unit\n")
t0 <- Sys.time()
Pm <- emmax_setup(GTs, GRM)
pm_obs <- emmax_fast(Pm, y)

test_sim <- ld_outlier_test(stage1, map, pm_obs, statistic = "simes", size_floor = SIZE_FLOOR,
                            alpha = ALPHA, assembly = "stage2_discovered", GTs = GTs,
                            LD_decay = b$LD_decay, score_threshold = REGION_ASSEMBLY$score_threshold,
                            distance_threshold = REGION_ASSEMBLY$distance_threshold)
print(test_sim)
say("    %.1f min\n\n", as.numeric(difftime(Sys.time(), t0, units = "mins")))
RESULTS$simes <- list(test = test_sim)

## ---- 3. per-replicate ground-truth sanity check -------------------------------
## NOT the PR/recall scoring PK described as a later, cross-replicate stage --
## this is a cheap, immediate check that the pipeline is pointed at something
## real before ten replicates' worth of compute goes into it. A region "hits" a
## true QTN if the QTN's marker falls within the region's [Chr, from, to] span.
.hits_truth <- function(regions, map) {
  truth <- map[true_QTN == TRUE, .(Chr, Pos)]
  if (!nrow(regions) || !nrow(truth)) return(0L)
  sum(vapply(seq_len(nrow(truth)), function(i)
    any(regions$Chr == truth$Chr[i] & regions$from <= truth$Pos[i] & regions$to >= truth$Pos[i]), logical(1)))
}
hit_con <- .hits_truth(test_con$regions, map); hit_sim <- .hits_truth(test_sim$regions, map)
say("[3] ground truth (sanity only -- not PR scoring): %d/%d true QTN recovered (consensus), %d/%d (Simes)\n",
    hit_con, sum(map$true_QTN), hit_sim, sum(map$true_QTN))

## ---- 4. agreement between the two statistics --------------------------------
sig_con <- test_con$units[significant == TRUE]$unit_id
sig_sim <- test_sim$units[significant == TRUE]$unit_id
jac <- if (length(union(sig_con, sig_sim))) length(intersect(sig_con, sig_sim)) / length(union(sig_con, sig_sim)) else NA_real_
say("[4] agreement: consensus %d, Simes %d, shared %d, Jaccard %s\n",
    length(sig_con), length(sig_sim), length(intersect(sig_con, sig_sim)),
    if (is.na(jac)) "NA (no discoveries)" else sprintf("%.3f", jac))

## ---- 5. save ------------------------------------------------------------------
OUT <- file.path(stage_dir(STAGE), sprintf("scan_%s_rep%d_%s_env%d.rds", TARGET_TAG, TARGET_REP, TARGET_CELL, TARGET_ENV))
saveRDS(list(consensus = RESULTS$consensus, simes = RESULTS$simes, jaccard = jac,
            truth_hits = list(consensus = hit_con, simes = hit_sim, n_true_qtn = sum(map$true_QTN))), OUT)
write_receipt(STAGE, inputs = INPUTS, params = PARAMS, outputs = OUT)
say("\n[5] wrote %s\n    receipt: %s\n", OUT, receipt_path(STAGE))
say("\n    Next: LFMM arm and/or widen REPS -- permutation null design still TBD (00_config.R).\n")
