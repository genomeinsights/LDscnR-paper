## =============================================================================
## module_sim/R/03_scan.R
##
## THE SCAN, BOTH ENGINES, for one replicate's bundle. Compute p-values, test
## stage-1 clusters, assemble regions -- uses the package's stage-1-cluster
## outlier API (ld_unit_matrix / ld_outlier_test), pinned via LDSCNR_PIN in
## 00_config.R. Renamed from 03_EMMAX.R (2026-09-04, PK: "analyse LFMM at the
## same time as EMMAX") -- unlike module_3sp, where LFMM is a separate, later
## stage (04_lfmm.R) because its p-values are INHERITED from an external,
## unreproducible source with no consensus arm available, this module's data is
## entirely ours: LFMM is computed fresh, here, from the same GTs and env
## EMMAX uses, with the SAME two-arm structure. There is no structural reason
## to split them.
##
## LFMM CALL PATTERN reused from the superseded module_sim_LDscnR/
## regen_sim_data.R, which already computed it this way for the same kind of
## data (LEA::write.lfmm/write.env -> lfmm2(K = LFMM_K) -> lfmm2.test
## (genomic.control = TRUE, full = TRUE)); not reinvented.
##
## THREE DELIBERATE DEPARTURES FROM module_3sp/R/03_EMMAX.R + 04_lfmm.R, each
## because the thing they mirror does not exist yet on this side rather than
## because it does not apply:
##
##   1. NO EcoPeak ROTATION NULL. There is no external validation set for
##      simulated data -- there is KNOWN TRUTH instead (map$true_QTN). Section
##      4 below is a per-replicate ground-truth sanity check, not a
##      substitute for one.
##   2. NO PERMUTATION NULL, for either engine. 00_config.R's NULLS section (7)
##      marks this explicitly TBD -- the naive MVN/latent-factor null this
##      project established is anti-conservative, and a structure-aware
##      replacement for sim phenotype permutation has not been designed.
##      Running the observed scan without a null is a real gap, not an
##      oversight; it is stated here so it cannot be mistaken for one later.
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
suppressMessages({library(data.table); library(LDscnR); library(digest); library(LEA)})
devtools::load_all("~/gitlab/LDscnR")
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim"), "R", "00_config.R"))
STAGE <- "03_scan"
say("=== %s ===\n\n", STAGE)
invisible(check_ldscnr())
if (!identical(LFMM_SOURCE, "compute")) stop("LFMM_SOURCE is not \"compute\" -- this stage expects to.")

TARGET_TAG <- Sys.getenv("SIM_TAG", TAGS[1]); TARGET_CELL <- Sys.getenv("SIM_CELL", CELLS[1]); TARGET_ENV <- as.integer(Sys.getenv("SIM_ENV", ENVS[1])); TARGET_REP <- as.integer(Sys.getenv("SIM_REP", REPS[1]))
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
               statistics = STATISTICS, unit_repr = UNIT_REPR, lfmm_source = LFMM_SOURCE,
               lfmm_k = LFMM_K)
if (!stage_stale(STAGE, INPUTS, PARAMS) && !nzchar(Sys.getenv("FORCE"))) {
  say("\nNothing to do. Set FORCE=1 to rerun anyway.\n"); quit(save = "no")
}
dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)

## LFMM helper: write.lfmm/write.env expect plain files on disk (LEA's own
## API); genotype matrix can be markers (Simes arm) or stage-1 units
## (consensus arm) -- same shape LEA expects either way, individuals x columns.
## Own tempdir per call, cleaned up immediately rather than accumulating.
##
## [!] pv$gif IS THE PRE-CORRECTION INFLATION FACTOR, NOT A CALIBRATION
## DIAGNOSTIC OF pv$pvalues. Checked directly (2026-09-04): on this replicate
## it reports 3.5-3.9, which reads like a badly miscalibrated scan -- but
## pv$pvalues is ALREADY genomic-control-corrected internally by
## lfmm2.test(genomic.control = TRUE), and re-deriving the calibration of the
## RETURNED p-values (converting back to F, computing the median-based gif
## myself) gives 1.000 exactly, with 6.7% of markers below p < 0.05 against a
## 5% nominal rate -- close to calibrated, not inflated. `gif` below is kept as
## a diagnostic of how much correction the raw scan needed, reported
## explicitly as pre-correction so it is not misread as the state of what is
## actually used downstream.
.lfmm_scan <- function(geno_mat, y, K) {
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  gf <- file.path(tmp, "geno.lfmm"); ef <- file.path(tmp, "grad.env")
  write.lfmm(geno_mat, gf); write.env(y, ef)
  proj <- lfmm2(gf, ef, K = K)
  pv <- suppressWarnings(lfmm2.test(proj, gf, ef, genomic.control = TRUE, full = TRUE))
  list(p = pv$pvalues, F = pv$fscores / pv$gif, gif_precorrection = pv$gif)
}

RESULTS <- list()

## ---- 1. EMMAX, consensus arm (statistic = "unit") -----------------------------
say("[1] EMMAX consensus -- ld_unit_matrix(repr = \"%s\") + emmax_fast\n", UNIT_REPR)
t0 <- Sys.time()
um <- ld_unit_matrix(GTs, stage1, map, size_floor = SIZE_FLOOR, repr = UNIT_REPR)
Pu <- emmax_setup(um, GRM)
pu_obs <- emmax_fast(Pu, y)
test_emx_con <- ld_outlier_test(stage1, map, pu_obs, statistic = "unit", size_floor = SIZE_FLOOR,
                                alpha = ALPHA, assembly = "stage2_discovered", GTs = GTs,
                                LD_decay = b$LD_decay, score_threshold = REGION_ASSEMBLY$score_threshold,
                                distance_threshold = REGION_ASSEMBLY$distance_threshold)
print(test_emx_con)
say("    %.1f min\n\n", as.numeric(difftime(Sys.time(), t0, units = "mins")))
RESULTS$emmax$consensus <- list(test = test_emx_con)

## ---- 2. EMMAX, Simes arm (statistic = "simes") --------------------------------
say("[2] EMMAX Simes -- direct marker-level scan aggregated per unit\n")
t0 <- Sys.time()
Pm <- emmax_setup(GTs, GRM)
pm_obs <- emmax_fast(Pm, y)
test_emx_sim <- ld_outlier_test(stage1, map, pm_obs, statistic = "simes", size_floor = SIZE_FLOOR,
                                alpha = ALPHA, assembly = "stage2_discovered", GTs = GTs,
                                LD_decay = b$LD_decay, score_threshold = REGION_ASSEMBLY$score_threshold,
                                distance_threshold = REGION_ASSEMBLY$distance_threshold)
print(test_emx_sim)
say("    %.1f min\n\n", as.numeric(difftime(Sys.time(), t0, units = "mins")))
RESULTS$emmax$simes <- list(test = test_emx_sim)

## ---- 3. LFMM, consensus arm ----------------------------------------------------
say("[3] LFMM consensus -- K = %d, on the same unit matrix EMMAX used\n", LFMM_K)
t0 <- Sys.time()
lf_con <- .lfmm_scan(um, y, LFMM_K)
say("    pre-correction gif = %.3f (pv$pvalues is already GC-corrected -- see .lfmm_scan)\n", lf_con$gif_precorrection)
test_lf_con <- ld_outlier_test(stage1, map, lf_con$p, statistic = "unit", size_floor = SIZE_FLOOR,
                               alpha = ALPHA, assembly = "stage2_discovered", GTs = GTs,
                               LD_decay = b$LD_decay, score_threshold = REGION_ASSEMBLY$score_threshold,
                               distance_threshold = REGION_ASSEMBLY$distance_threshold)
print(test_lf_con)
say("    %.1f min\n\n", as.numeric(difftime(Sys.time(), t0, units = "mins")))
RESULTS$lfmm$consensus <- list(test = test_lf_con, gif_precorrection = lf_con$gif_precorrection)

## ---- 4. LFMM, Simes arm ---------------------------------------------------------
say("[4] LFMM Simes -- direct marker-level scan, K = %d\n", LFMM_K)
t0 <- Sys.time()
lf_sim <- .lfmm_scan(GTs, y, LFMM_K)
say("    pre-correction gif = %.3f (pv$pvalues is already GC-corrected -- see .lfmm_scan)\n", lf_sim$gif_precorrection)
test_lf_sim <- ld_outlier_test(stage1, map, lf_sim$p, statistic = "simes", size_floor = SIZE_FLOOR,
                               alpha = ALPHA, assembly = "stage2_discovered", GTs = GTs,
                               LD_decay = b$LD_decay, score_threshold = REGION_ASSEMBLY$score_threshold,
                               distance_threshold = REGION_ASSEMBLY$distance_threshold)
print(test_lf_sim)
say("    %.1f min\n\n", as.numeric(difftime(Sys.time(), t0, units = "mins")))
RESULTS$lfmm$simes <- list(test = test_lf_sim, gif_precorrection = lf_sim$gif_precorrection)

## ---- 5. per-replicate ground-truth sanity check, all four arms -----------------
## NOT the PR/recall scoring PK described as a later, cross-replicate stage --
## this is a cheap, immediate check that the pipeline is pointed at something
## real before ten replicates' worth of compute goes into it. A region "hits" a
## true QTN if the QTN's marker falls within the region's [Chr, from, to] span.
##
## [!] VECTORISED, NOT LOOPED. The first version was a per-truth-row vapply
## re-scanning the full regions table on every iteration -- O(n_truth x
## n_regions), harmless only because both are tiny here (single digits: 2-9
## true QTN, 0-9 regions per replicate seen so far). Checked before deciding it
## was fine to leave: the final cross-replicate PR-scoring stage PK described
## will pool ~10x the QTN and ~10x the regions, and this exact shape is the one
## kind of "small enough today" loop that turns into a real cost once copied
## into that stage. Replaced with data.table::foverlaps(), the same tool
## ld_outlier_test()'s own physical-merge branch uses for an identical
## Chr/from/to range join -- one vectorised join instead of a loop, and it
## scales the same way regardless of how large truth or regions gets.
.hits_truth <- function(regions, map) {
  truth <- map[true_QTN == TRUE, .(Chr, Pos)]
  if (!nrow(regions) || !nrow(truth)) return(0L)
  truth[, `:=`(from = Pos, to = Pos)]
  rg <- data.table::copy(regions); data.table::setkey(rg, Chr, from, to)
  ov <- data.table::foverlaps(truth, rg, by.x = c("Chr", "from", "to"), type = "within", nomatch = NULL)
  data.table::uniqueN(ov, by = c("Chr", "Pos"))
}
n_truth <- sum(map$true_QTN)
hits <- list(emmax_consensus = .hits_truth(test_emx_con$regions, map),
            emmax_simes     = .hits_truth(test_emx_sim$regions, map),
            lfmm_consensus  = .hits_truth(test_lf_con$regions, map),
            lfmm_simes      = .hits_truth(test_lf_sim$regions, map))
say("[5] ground truth (sanity only -- not PR scoring), of %d true QTN: %s\n", n_truth,
    paste(sprintf("%s=%d", names(hits), unlist(hits)), collapse = ", "))

## ---- 6. agreement, within and across engines -----------------------------------
sig <- list(emmax_consensus = test_emx_con$units[significant == TRUE]$unit_id,
           emmax_simes     = test_emx_sim$units[significant == TRUE]$unit_id,
           lfmm_consensus  = test_lf_con$units[significant == TRUE]$unit_id,
           lfmm_simes      = test_lf_sim$units[significant == TRUE]$unit_id)
.jac <- function(a, b) if (length(union(a, b))) length(intersect(a, b)) / length(union(a, b)) else NA_real_
.fmt <- function(x) if (is.na(x)) "NA" else sprintf("%.3f", x)
say("[6] discoveries: %s\n", paste(sprintf("%s=%d", names(sig), lengths(sig)), collapse = ", "))
say("    Jaccard, within engine  -- EMMAX consensus/Simes: %s ; LFMM consensus/Simes: %s\n",
    .fmt(.jac(sig$emmax_consensus, sig$emmax_simes)), .fmt(.jac(sig$lfmm_consensus, sig$lfmm_simes)))
say("    Jaccard, across engine  -- consensus EMMAX/LFMM: %s ; Simes EMMAX/LFMM: %s\n",
    .fmt(.jac(sig$emmax_consensus, sig$lfmm_consensus)), .fmt(.jac(sig$emmax_simes, sig$lfmm_simes)))

## ---- 7. save ---------------------------------------------------------------------
OUT <- file.path(stage_dir(STAGE), sprintf("scan_%s_rep%d_%s_env%d.rds", TARGET_TAG, TARGET_REP, TARGET_CELL, TARGET_ENV))
saveRDS(list(results = RESULTS, sig = sig, hits = hits, n_true_qtn = n_truth), OUT)
write_receipt(STAGE, inputs = INPUTS, params = PARAMS, outputs = OUT)
say("\n[7] wrote %s\n    receipt: %s\n", OUT, receipt_path(STAGE))
say("\n    Next: widen REPS -- permutation null design still TBD (00_config.R).\n")
