## =============================================================================
## module_9sp/R/03_EMMAX.R
##
## THE EMMAX ARM. Same structure as module_3sp/R/03_EMMAX.R (consensus + Simes,
## permutation null, size floor from 00_config.R) with two real differences,
## both DECIDED by PK on 2026-09-07 after a direct PCA/crosstab check:
##
##   1. COVARIATE. lineage is regressed out of ecotype BEFORE any test, exactly
##      the mechanism emmax()'s own `Covar` argument uses (Y <- resid(lm(Y ~
##      Covar))) -- replicated manually because emmax_fast(), the 25x-faster
##      scan this module uses throughout, has no Covar parameter (see
##      emmax_fast.R). Applied identically to the observed phenotype and to
##      every permutation draw, so the null is calibrated to the same
##      covariate-adjusted quantity being tested, not a different one.
##
##   2. PERMUTATION. "within lineage", not "within locality" -- shuffles
##      ecotype at the POPULATION level (never individually) inside each
##      lineage stratum. WL and WA are entirely monomorphic for ecotype, so
##      permutation draws within those strata are deterministic (reproduce the
##      observed labels every time) -- expected, not a bug; see 00_config.R.
##
## NO ECOPEAKS / ROTATION STAGE. No external-validation reference exists for
## nine-spined stickleback in this codebase (PK, 2026-09-07: skip it for now).
## This script therefore stops after the permutation null -- no
## ld_region_rotation() call, no KG/LEN, unlike module_3sp's version.
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(digest)})
devtools::load_all("~/gitlab/LDscnR")
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_9sp"), "R", "00_config.R"))
STAGE <- "03_EMMAX"
say("=== %s ===\n\n", STAGE)

## ---- 0. version pin, bundle ---------------------------------------------------
invisible(check_ldscnr())
stopifnot("EMMAX_COVAR must be decided before this stage runs" = !is.null(EMMAX_COVAR))
stopifnot("PERM_SCHEMES must be decided before this stage runs" = !is.null(PERM_SCHEMES))

BUNDLE_PATH <- file.path(PATHS$out, "02_bundle", "bundle.rds")
b <- readRDS(BUNDLE_PATH)
GTs <- b$GTs; map <- b$map; eco <- b$eco; stage1 <- b$stage1
GRM <- b$GRM; pheno <- b$pheno
say("[0] bundle: %d individuals x %s markers ; %s stage-1 units at floor %d\n",
    nrow(GTs), format(ncol(GTs), big.mark=","),
    format(sum({cl<-stage1$clusters; nl<-if("n_loci" %in% names(cl)) cl$n_loci else cl$n_snps;
                sum(nl >= SIZE_FLOOR)}), big.mark=","), SIZE_FLOOR)

## ---- 0b. covariate: residualise ecotype on lineage, ONCE, on the OBSERVED data ----
say("[0b] covariate: ecotype ~ %s (residuals tested in place of raw ecotype)\n", EMMAX_COVAR)
stopifnot(EMMAX_COVAR %in% names(pheno))
covar_obs <- pheno[[EMMAX_COVAR]]
eco_resid <- as.numeric(stats::resid(stats::lm(eco ~ covar_obs)))
say("    lineage counts: %s\n", paste(sprintf("%s=%d", names(table(covar_obs)), table(covar_obs)), collapse = ", "))
say("    ecotype variance explained by %s: R^2 = %.3f\n", EMMAX_COVAR,
    summary(lm(eco ~ covar_obs))$r.squared)

## ---- phenotype permutation, within lineage, at the POPULATION level ----------
## Same granularity as 3sp's perm_regional(): individuals within one population
## always share one permuted label. `sample(ecotype)` within a monomorphic
## lineage (WL, WA) deterministically returns the observed labels -- see
## 00_config.R's PERM_SCHEMES note.
POPT <- unique(pheno[, .(pop_ID, ecotype = eco, lineage = get(EMMAX_COVAR))])
perm_lineage <- function() {
  pt <- copy(POPT)[, ep := sample(ecotype), by = lineage]
  y  <- pt$ep[match(pheno$pop_ID, pt$pop_ID)]
  as.numeric(stats::resid(stats::lm(y ~ covar_obs)))   # SAME covariate adjustment as observed
}

INPUTS <- c(BUNDLE_PATH)
PARAMS <- list(size_floor = SIZE_FLOOR, alpha = ALPHA, region_assembly = REGION_ASSEMBLY,
               statistics = STATISTICS, unit_repr = UNIT_REPR, nperm_consensus = NPERM_CONSENSUS,
               nperm_simes = NPERM_SIMES, emmax_covar = EMMAX_COVAR, perm_schemes = PERM_SCHEMES)
if (!stage_stale(STAGE, INPUTS, PARAMS) && !nzchar(Sys.getenv("FORCE"))) {
  say("\nNothing to do. Set FORCE=1 to rerun anyway.\n"); quit(save = "no")
}
dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)

RESULTS <- list()

## ---- 1. consensus arm (statistic = "unit") ----------------------------------
say("\n[1] CONSENSUS arm -- ld_unit_matrix(repr = \"%s\") + emmax_fast\n", UNIT_REPR)
t0 <- Sys.time()
um <- ld_unit_matrix(GTs, stage1, map, size_floor = SIZE_FLOOR, repr = UNIT_REPR)
Pu <- emmax_setup(um, GRM)
pu_obs <- emmax_fast(Pu, eco_resid)

test_con <- ld_outlier_test(stage1, map, pu_obs, statistic = "unit", size_floor = SIZE_FLOOR,
                            alpha = ALPHA, assembly = "stage2_discovered", GTs = GTs,
                            LD_decay = b$LD_decay, score_threshold = REGION_ASSEMBLY$score_threshold,
                            distance_threshold = REGION_ASSEMBLY$distance_threshold)
print(test_con)

say("    permutation null: %d within-lineage surrogates (consensus is cheap)\n", NPERM_CONSENSUS)
p_perm_con <- function(bb) { set.seed(bb); y <- perm_lineage()
  emmax_fast(Pu, y) }
null_con <- ld_outlier_perm(test_con, stage1, map, p_perm_con, GTs = GTs, LD_decay = b$LD_decay,
                            B = NPERM_CONSENSUS, level = "units", verbose = TRUE)
print(null_con)
say("    %.1f min\n\n", as.numeric(difftime(Sys.time(), t0, units="mins")))
RESULTS$consensus <- list(test = test_con, null = null_con)

## ---- 2. Simes arm (statistic = "simes") -------------------------------------
say("[2] SIMES arm -- direct marker-level scan aggregated per unit\n")
t0 <- Sys.time()
Pm <- emmax_setup(GTs, GRM)
pm_obs <- emmax_fast(Pm, eco_resid)

test_sim <- ld_outlier_test(stage1, map, pm_obs, statistic = "simes", size_floor = SIZE_FLOOR,
                            alpha = ALPHA, assembly = "stage2_discovered", GTs = GTs,
                            LD_decay = b$LD_decay, score_threshold = REGION_ASSEMBLY$score_threshold,
                            distance_threshold = REGION_ASSEMBLY$distance_threshold)
print(test_sim)

say("    permutation null: %d within-lineage surrogates (Simes rescans all %s markers per draw)\n",
    NPERM_SIMES, format(ncol(GTs), big.mark=","))
p_perm_sim <- function(bb) { set.seed(bb); y <- perm_lineage(); emmax_fast(Pm, y) }
null_sim <- ld_outlier_perm(test_sim, stage1, map, p_perm_sim, GTs = GTs, LD_decay = b$LD_decay,
                            B = NPERM_SIMES, level = "units", verbose = TRUE)
print(null_sim)
say("    %.1f min\n\n", as.numeric(difftime(Sys.time(), t0, units="mins")))
RESULTS$simes <- list(test = test_sim, null = null_sim)

## ---- 3. agreement between the two statistics --------------------------------
sig_con <- test_con$units[significant == TRUE]$unit_id
sig_sim <- test_sim$units[significant == TRUE]$unit_id
jac <- length(intersect(sig_con, sig_sim)) / length(union(sig_con, sig_sim))
say("[3] agreement: consensus %d, Simes %d, shared %d, Jaccard %.3f\n",
    length(sig_con), length(sig_sim), length(intersect(sig_con, sig_sim)), jac)

## ---- 4. save ------------------------------------------------------------------
OUT <- file.path(stage_dir(STAGE), "scan.rds")
saveRDS(list(consensus = RESULTS$consensus, simes = RESULTS$simes, jaccard = jac,
            eco_resid = eco_resid, covar = EMMAX_COVAR), OUT)
write_receipt(STAGE, inputs = INPUTS, params = PARAMS, outputs = OUT)
say("\n[4] wrote %s\n    receipt: %s\n", OUT, receipt_path(STAGE))
say("\n    No EcoPeaks-equivalent reference for 9sp -- no rotation-null stage.\n")
say("    Consider a sensitivity rerun restricted to EL+Admixed (the two lineages\n")
say("    with real within-lineage ecotype contrast) once this result is in hand.\n")
