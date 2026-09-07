## =============================================================================
## module_9sp/R/03b_EMMAX_sensitivity_ELAdmixed.R
##
## SENSITIVITY CHECK, not the canonical result: restrict to the two lineages
## with real within-lineage ecotype contrast -- EL (54 Freshwater / 6 Marine)
## and Admixed (4 Freshwater / 42 Marine) -- dropping WL and WA (43
## individuals, entirely monomorphic for ecotype, see 00_config.R). PK,
## 2026-09-07: "consider a sensitivity rerun ... to see whether the reported
## signal holds up on the cleaner subset."
##
## SAME DESIGN as 03_EMMAX.R (lineage covariate, within-lineage permutation,
## SIZE_FLOOR, region assembly), same Stage-1 clustering (marker-based, not
## individual-based, so NOT rebuilt for the subset -- same reasoning as any
## subgroup analysis reusing a genome-wide-derived marker partition). The GRM
## IS rebuilt, not sliced from the full-sample one: snpgdsGRM's GCTA estimator
## uses the SAMPLE's own allele frequencies, so a subset GRM is not the same
## object as the corresponding block of the full-sample GRM.
##
## WHAT THIS TESTS: does the null-exceeds-observed pattern (both arms' surrogate
## means > observed in 03_EMMAX.R) persist once WL/WA's zero-variance majority
## is removed from both the test and the permutation, or was it partly an
## artefact of population-level permutation being noisy when the informative
## lineages are this small and lopsided (00_config.R's working hypothesis)?
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(SNPRelate); library(digest)})
devtools::load_all("~/gitlab/LDscnR")
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_9sp"), "R", "00_config.R"))
STAGE <- "03b_EMMAX_sensitivity_ELAdmixed"
say("=== %s ===\n\n", STAGE)

invisible(check_ldscnr())
stopifnot(!is.null(EMMAX_COVAR), !is.null(PERM_SCHEMES))

BUNDLE_PATH <- file.path(PATHS$out, "02_bundle", "bundle.rds")
b <- readRDS(BUNDLE_PATH)
map <- b$map; stage1 <- b$stage1; pheno_full <- b$pheno

## ---- 1. subset individuals ----------------------------------------------------
KEEP_LINEAGES <- c("EL", "Admixed")
keep <- pheno_full[[EMMAX_COVAR]] %in% KEEP_LINEAGES
say("[1] restricting to lineages: %s\n", paste(KEEP_LINEAGES, collapse = ", "))
say("    %d of %d individuals kept ; ecotype counts by lineage:\n", sum(keep), length(keep))
print(table(pheno_full[[EMMAX_COVAR]][keep], b$eco[keep]))

GTs   <- b$GTs[keep, , drop = FALSE]
eco   <- b$eco[keep]
pheno <- pheno_full[keep]
covar_obs <- pheno[[EMMAX_COVAR]]

## ---- 2. rebuild the GRM on JUST this subset (own allele frequencies) ---------
## Same GDS 02_bundle.R already built (module_9sp/cache/9sp.gds), same grm_markers
## (stage-1 representatives) -- only the SAMPLE set differs.
gds_path <- file.path(PATHS$cache, "9sp.gds")
stopifnot(file.exists(gds_path))
gds <- snpgdsOpen(gds_path)
on.exit(try(snpgdsClose(gds), silent = TRUE), add = TRUE)
sample_ids_full <- paste0("ind_", seq_len(nrow(b$GTs)))   # create_gds_from_geno's own convention
sample_ids_keep <- sample_ids_full[keep]
say("\n[2] rebuilding GRM on %d individuals, %s markers (%s)\n",
    sum(keep), format(length(b$grm_markers), big.mark=","), GRM_METHOD)
GRM <- snpgdsGRM(gds, sample.id = sample_ids_keep, snp.id = b$grm_markers,
                 method = GRM_METHOD, verbose = FALSE, autosome.only = FALSE)$grm
stopifnot(nrow(GRM) == sum(keep), ncol(GRM) == sum(keep))
ut <- upper.tri(GRM)
say("    %d x %d ; off-diagonal mean %+.4f sd %.4f\n", nrow(GRM), ncol(GRM), mean(GRM[ut]), sd(GRM[ut]))

## ---- 3. covariate: residualise ecotype on lineage, on the SUBSET -------------
eco_resid <- as.numeric(stats::resid(stats::lm(eco ~ covar_obs)))
say("\n[3] covariate: ecotype ~ %s (subset) ; R^2 = %.3f\n", EMMAX_COVAR,
    summary(lm(eco ~ covar_obs))$r.squared)

POPT <- unique(pheno[, .(pop_ID, ecotype = eco, lineage = get(EMMAX_COVAR))])
perm_lineage <- function() {
  pt <- copy(POPT)[, ep := sample(ecotype), by = lineage]
  y  <- pt$ep[match(pheno$pop_ID, pt$pop_ID)]
  as.numeric(stats::resid(stats::lm(y ~ covar_obs)))
}

INPUTS <- c(BUNDLE_PATH, gds_path)
PARAMS <- list(size_floor = SIZE_FLOOR, alpha = ALPHA, region_assembly = REGION_ASSEMBLY,
               statistics = STATISTICS, unit_repr = UNIT_REPR, nperm_consensus = NPERM_CONSENSUS,
               nperm_simes = NPERM_SIMES, emmax_covar = EMMAX_COVAR, perm_schemes = PERM_SCHEMES,
               keep_lineages = KEEP_LINEAGES)
if (!stage_stale(STAGE, INPUTS, PARAMS) && !nzchar(Sys.getenv("FORCE"))) {
  say("\nNothing to do. Set FORCE=1 to rerun anyway.\n"); quit(save = "no")
}
dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)

RESULTS <- list()

## ---- 4. consensus arm ----------------------------------------------------------
say("\n[4] CONSENSUS arm\n")
t0 <- Sys.time()
um <- ld_unit_matrix(GTs, stage1, map, size_floor = SIZE_FLOOR, repr = UNIT_REPR)
Pu <- emmax_setup(um, GRM)
pu_obs <- emmax_fast(Pu, eco_resid)

test_con <- ld_outlier_test(stage1, map, pu_obs, statistic = "unit", size_floor = SIZE_FLOOR,
                            alpha = ALPHA, assembly = "stage2_discovered", GTs = GTs,
                            LD_decay = b$LD_decay, score_threshold = REGION_ASSEMBLY$score_threshold,
                            distance_threshold = REGION_ASSEMBLY$distance_threshold)
print(test_con)
p_perm_con <- function(bb) { set.seed(bb); y <- perm_lineage(); emmax_fast(Pu, y) }
null_con <- ld_outlier_perm(test_con, stage1, map, p_perm_con, GTs = GTs, LD_decay = b$LD_decay,
                            B = NPERM_CONSENSUS, level = "units", verbose = TRUE)
print(null_con)
say("    %.1f min\n\n", as.numeric(difftime(Sys.time(), t0, units="mins")))
RESULTS$consensus <- list(test = test_con, null = null_con)

## ---- 5. Simes arm ---------------------------------------------------------------
say("[5] SIMES arm\n")
t0 <- Sys.time()
Pm <- emmax_setup(GTs, GRM)
pm_obs <- emmax_fast(Pm, eco_resid)

test_sim <- ld_outlier_test(stage1, map, pm_obs, statistic = "simes", size_floor = SIZE_FLOOR,
                            alpha = ALPHA, assembly = "stage2_discovered", GTs = GTs,
                            LD_decay = b$LD_decay, score_threshold = REGION_ASSEMBLY$score_threshold,
                            distance_threshold = REGION_ASSEMBLY$distance_threshold)
print(test_sim)
p_perm_sim <- function(bb) { set.seed(bb); y <- perm_lineage(); emmax_fast(Pm, y) }
null_sim <- ld_outlier_perm(test_sim, stage1, map, p_perm_sim, GTs = GTs, LD_decay = b$LD_decay,
                            B = NPERM_SIMES, level = "units", verbose = TRUE)
print(null_sim)
say("    %.1f min\n\n", as.numeric(difftime(Sys.time(), t0, units="mins")))
RESULTS$simes <- list(test = test_sim, null = null_sim)

## ---- 6. agreement + comparison to the full-sample result ----------------------
sig_con <- test_con$units[significant == TRUE]$unit_id
sig_sim <- test_sim$units[significant == TRUE]$unit_id
jac <- length(intersect(sig_con, sig_sim)) / length(union(sig_con, sig_sim))
say("[6] agreement: consensus %d, Simes %d, shared %d, Jaccard %.3f\n",
    length(sig_con), length(sig_sim), length(intersect(sig_con, sig_sim)), jac)

FULL_PATH <- file.path(PATHS$out, "03_EMMAX", "scan.rds")
if (file.exists(FULL_PATH)) {
  full <- readRDS(FULL_PATH)
  full_sig_con <- full$consensus$test$units[significant == TRUE]$unit_id
  ov <- length(intersect(sig_con, full_sig_con))
  say("    vs full-sample consensus: %d full-sample units, %d shared with this subset (of %d here)\n",
      length(full_sig_con), ov, length(sig_con))
}

## ---- 7. save ---------------------------------------------------------------------
OUT <- file.path(stage_dir(STAGE), "scan.rds")
saveRDS(list(consensus = RESULTS$consensus, simes = RESULTS$simes, jaccard = jac,
            eco_resid = eco_resid, covar = EMMAX_COVAR, keep_lineages = KEEP_LINEAGES,
            n_individuals = sum(keep)), OUT)
write_receipt(STAGE, inputs = INPUTS, params = PARAMS, outputs = OUT)
say("\n[7] wrote %s\n    receipt: %s\n", OUT, receipt_path(STAGE))
