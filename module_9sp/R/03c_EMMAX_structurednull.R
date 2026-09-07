## =============================================================================
## module_9sp/R/03c_EMMAX_structurednull.R
##
## A STRUCTURE-AWARE NULL, as an alternative to 03_EMMAX.R's within-lineage
## permutation. PK, 2026-09-07: "look into structured_null() -- this could be
## a nice additional worked example for the paper."
##
## NOT A PORT OF structured_null() ITSELF (R/structured_null.R) -- that
## function is built for the C-score/ld_cscore() candidate-selection
## framework, which this module deliberately does NOT use (see
## Supplementary.tex's "Interpretive limits": ld_w-based candidate selection
## misses most of the regions Stage-1 clustering + cluster-level testing
## actually recovers on this data -- 21 of 39 3sp regions have median member
## ld_w < 0.1). What IS reused is the SURROGATE-GENERATION MECHANISM, applied
## to the ld_outlier_perm() framework this module already uses instead:
##
##   s ~ MVN(0, K)                      -- draw a surrogate with the SAME
##                                          genetic-relatedness covariance as
##                                          the real GRM (via K's eigenbasis:
##                                          s = V %*% (sqrt(lambda) * rnorm(n)))
##   s_orth <- resid(lm(s ~ eco_resid)) -- Gram-Schmidt: orthogonalise against
##                                          the ACTUAL tested phenotype, so the
##                                          surrogate carries structure-driven
##                                          false-positive risk but NONE of the
##                                          real signal.
##
## WHY TRY THIS: 03_EMMAX.R's within-LINEAGE, population-LEVEL label
## permutation showed the null typically producing MORE apparent discoveries
## than observed (ratios 122.5%/328.4%), and restricting to the two lineages
## with real ecotype contrast (03b) made this WORSE, not better (1257%/675%)
## -- refuting the hypothesis that small, lopsided permutation strata were the
## problem. The remaining, more fundamental hypothesis (00_config.R, 03b's
## header): 9sp's population structure is simply much stronger than 3sp's
## (background LD 0.301 vs 0.055; GRM off-diagonal SD 0.26-0.32) and a
## DISCRETE label permutation may not represent it well. A CONTINUOUS,
## GRM-covariance-matched surrogate is a direct test of exactly that -- it
## reproduces the relatedness structure the GRM itself estimates, not an
## approximation of it via population/lineage group-swapping.
##
## RUNS ON THE FULL 149-INDIVIDUAL SAMPLE, unlike 03b: the structured null has
## no "monomorphic stratum" problem (WL/WA aren't a blocker here -- a
## continuous MVN(0,K) draw doesn't care that WL/WA have zero ecotype
## variance), so there is no reason to drop them and lose power.
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(digest)})
devtools::load_all("~/gitlab/LDscnR")
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_9sp"), "R", "00_config.R"))
STAGE <- "03c_EMMAX_structurednull"
say("=== %s ===\n\n", STAGE)

invisible(check_ldscnr())
stopifnot(!is.null(EMMAX_COVAR))

BUNDLE_PATH <- file.path(PATHS$out, "02_bundle", "bundle.rds")
SCAN_PATH   <- file.path(PATHS$out, "03_EMMAX", "scan.rds")
b  <- readRDS(BUNDLE_PATH)
sc <- readRDS(SCAN_PATH)
GTs <- b$GTs; map <- b$map; stage1 <- b$stage1; GRM <- b$GRM; pheno <- b$pheno
eco_resid <- sc$eco_resid   ## SAME covariate-adjusted phenotype 03_EMMAX.R tested -- the
                            ## structured null must be orthogonal to what is actually tested,
                            ## not to the raw (pre-covariate) ecotype.
say("[0] bundle: %d individuals ; reusing 03_EMMAX.R's lineage-residualised phenotype\n", nrow(GTs))

## ---- 0b. structure-aware surrogate generator ----------------------------------
## eigendecompose the SAME GRM 03_EMMAX.R corrected for -- reusing structured_null()'s own
## mechanism (R/structured_null.R), not the C-score machinery built around it.
say("[0b] eigendecomposing the GRM for the structured surrogate basis\n")
n <- nrow(GTs)
eK <- eigen(GRM, symmetric = TRUE)
Lv <- pmax(eK$values, 0); Vk <- eK$vectors
gen_structured <- function() {
  s <- as.numeric(Vk %*% (sqrt(Lv) * stats::rnorm(n)))
  as.numeric(stats::resid(stats::lm(s ~ eco_resid)))
}

INPUTS <- c(BUNDLE_PATH, SCAN_PATH)
PARAMS <- list(size_floor = SIZE_FLOOR, alpha = ALPHA, region_assembly = REGION_ASSEMBLY,
               statistics = STATISTICS, unit_repr = UNIT_REPR, nperm_consensus = NPERM_CONSENSUS,
               nperm_simes = NPERM_SIMES, emmax_covar = EMMAX_COVAR, null = "structured_grm")
if (!stage_stale(STAGE, INPUTS, PARAMS) && !nzchar(Sys.getenv("FORCE"))) {
  say("\nNothing to do. Set FORCE=1 to rerun anyway.\n"); quit(save = "no")
}
dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)

## Reuse 03_EMMAX.R's OWN test results (test_con/test_sim) -- the discovered units/regions
## are unchanged by which null we calibrate against; only the null itself is new. Re-deriving
## a fresh ld_outlier_test() from the already-saved observed p-values would be redundant.
test_con <- sc$consensus$test
test_sim <- sc$simes$test

RESULTS <- list()

## ---- 1. consensus arm, structured null -----------------------------------------
say("\n[1] CONSENSUS arm -- structured null (%d surrogates)\n", NPERM_CONSENSUS)
t0 <- Sys.time()
um <- ld_unit_matrix(GTs, stage1, map, size_floor = SIZE_FLOOR, repr = UNIT_REPR)
Pu <- emmax_setup(um, GRM)
p_struct_con <- function(bb) { set.seed(bb); emmax_fast(Pu, gen_structured()) }
null_con_struct <- ld_outlier_perm(test_con, stage1, map, p_struct_con, GTs = GTs,
                                   LD_decay = b$LD_decay, B = NPERM_CONSENSUS,
                                   level = "units", verbose = TRUE)
print(null_con_struct)
say("    %.1f min\n\n", as.numeric(difftime(Sys.time(), t0, units="mins")))
RESULTS$consensus <- list(test = test_con, null = null_con_struct)

## ---- 2. Simes arm, structured null ----------------------------------------------
say("[2] SIMES arm -- structured null (%d surrogates)\n", NPERM_SIMES)
t0 <- Sys.time()
Pm <- emmax_setup(GTs, GRM)
p_struct_sim <- function(bb) { set.seed(bb); emmax_fast(Pm, gen_structured()) }
null_sim_struct <- ld_outlier_perm(test_sim, stage1, map, p_struct_sim, GTs = GTs,
                                   LD_decay = b$LD_decay, B = NPERM_SIMES,
                                   level = "units", verbose = TRUE)
print(null_sim_struct)
say("    %.1f min\n\n", as.numeric(difftime(Sys.time(), t0, units="mins")))
RESULTS$simes <- list(test = test_sim, null = null_sim_struct)

## ---- 3. side-by-side with the label-permutation null ---------------------------
say("[3] label-permutation vs structured-GRM null, side by side\n")
say("    consensus: observed %d | label-perm surrogate mean %.2f (p=%.4f) | structured surrogate mean %.2f (p=%.4f)\n",
    null_con_struct$observed, mean(sc$consensus$null$surrogates), sc$consensus$null$p,
    mean(null_con_struct$surrogates), null_con_struct$p)
say("    simes:     observed %d | label-perm surrogate mean %.2f (p=%.4f) | structured surrogate mean %.2f (p=%.4f)\n",
    null_sim_struct$observed, mean(sc$simes$null$surrogates), sc$simes$null$p,
    mean(null_sim_struct$surrogates), null_sim_struct$p)

## ---- 4. save ---------------------------------------------------------------------
OUT <- file.path(stage_dir(STAGE), "scan.rds")
saveRDS(list(consensus = RESULTS$consensus, simes = RESULTS$simes, null_type = "structured_grm"), OUT)
write_receipt(STAGE, inputs = INPUTS, params = PARAMS, outputs = OUT)
say("\n[4] wrote %s\n    receipt: %s\n", OUT, receipt_path(STAGE))
