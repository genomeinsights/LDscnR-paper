## =============================================================================
## module_3sp/R/03_scan.R
##
## THE REAL SCAN: compute p-values, test stage-1 clusters, assemble regions,
## validate against a permutation null and against real EcoPeaks. Uses the
## package's new stage-1-cluster outlier API (ld_unit_matrix / ld_outlier_test /
## ld_outlier_perm / ld_region_rotation), pinned via LDSCNR_PIN in 00_config.R.
##
## THERE IS NO SEPARATE "03_clusters.R" STAGE. The original plan had one
## (stage-2 grouping as its own step), but ld_outlier_test(assembly =
## "stage2_discovered") now does that internally -- the region-assembly step
## PK validated (mechanism_over_content DA) is inside the package function, not
## a script here. This picks up directly from the bundle.
##
## TWO ARMS, BOTH REPORTED, NEITHER PRIVILEGED IN THE SCRIPT ITSELF:
##   consensus (statistic = "unit", on ld_unit_matrix(repr = "consensus_dosage"))
##       the arm every number in today's record used -- the more powerful
##       statistic, and the one whose kinship-estimator sensitivity (CR) and
##       false-positive-content divergence from the sims (CY) are already
##       characterised.
##   Simes    (statistic = "simes", directly on the marker-level scan)
##       the comparator. No size penalty vs Simes' n-penalty is the mechanism
##       (block 01 today); reporting both, not picking one, is the point.
##
## COST. consensus is cheap: ~1,356 columns, so NPERM_CONSENSUS = 1000
## surrogates is a few minutes. Simes rescans all 790,578 markers per
## surrogate; NPERM_SIMES = 200 is already the discount rate 00_config.R set
## for exactly this reason.
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(digest)})
devtools::load_all("~/gitlab/LDscnR")
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_3sp"), "R", "00_config.R"))
STAGE <- "03_scan"
say("=== %s ===\n\n", STAGE)

## ---- 0. version pin, bundle, EcoPeaks ---------------------------------------
invisible(check_ldscnr())

BUNDLE_PATH <- file.path(PATHS$out, "02_bundle", "bundle.rds")
b <- readRDS(BUNDLE_PATH)
GTs <- b$GTs; map <- b$map; eco <- b$eco; stage1 <- b$stage1
GRM <- b$GRM; pheno <- b$pheno
say("[0] bundle: %d individuals x %s markers ; %s stage-1 units at floor %d\n",
    nrow(GTs), format(ncol(GTs), big.mark=","),
    format(sum({cl<-stage1$clusters; nl<-if("n_loci" %in% names(cl)) cl$n_loci else cl$n_snps;
                sum(nl >= SIZE_FLOOR)}), big.mark=","), SIZE_FLOOR)

LIFT <- path.expand("~/gitlab/LDscnR-paper/kingman2021/data/liftover")
ROMAN <- c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV",
           "XV","XVI","XVII","XVIII","XIX","XX","XXI")
.rd <- function(f) { x <- fread(file.path(LIFT, f), header=FALSE,
                                col.names=c("chr","start","end","pv"))
  x[, chr_num := match(sub("^chr","",chr), ROMAN)]
  x[!is.na(chr_num), .(chr = paste0("Chr", chr_num), start, end)] }
KG <- unique(rbindlist(lapply(ECOPEAK_BEDS, .rd)))
LEN <- map[, .(len = max(Pos)), by = Chr]; setnames(LEN, "Chr", "chr")
say("[0] EcoPeaks: %d intervals, %.1f%% of assembled sequence\n\n",
    nrow(KG), 100*sum(as.numeric(KG$end-KG$start))/sum(as.numeric(LEN$len)))

## ---- phenotype permutation schemes (regional is primary; population is run too)
POPT <- unique(pheno[, .(pop_ID, ecotype, pop_locality)])
perm_regional <- function() { pt <- copy(POPT)[, ep := sample(ecotype), by = pop_locality]
  as.numeric(pt$ep[match(pheno$pop_ID, pt$pop_ID)] == "Marine") }
perm_population <- function() { pt <- copy(POPT)[, ep := sample(ecotype)]
  as.numeric(pt$ep[match(pheno$pop_ID, pt$pop_ID)] == "Marine") }

INPUTS <- c(BUNDLE_PATH, file.path(LIFT, ECOPEAK_BEDS))
PARAMS <- list(size_floor = SIZE_FLOOR, alpha = ALPHA, region_assembly = REGION_ASSEMBLY,
               statistics = STATISTICS, nperm_consensus = NPERM_CONSENSUS,
               nperm_simes = NPERM_SIMES, n_rotations = N_ROTATIONS,
               rotation_scheme = ROTATION_SCHEME)
if (!stage_stale(STAGE, INPUTS, PARAMS) && !nzchar(Sys.getenv("FORCE"))) {
  say("\nNothing to do. Set FORCE=1 to rerun anyway.\n"); quit(save = "no")
}
dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)

RESULTS <- list()

## ---- 1. consensus arm (statistic = "unit") ----------------------------------
say("[1] CONSENSUS arm -- ld_unit_matrix(repr = \"consensus_dosage\") + emmax_fast\n")
t0 <- Sys.time()
um <- ld_unit_matrix(GTs, stage1, map, size_floor = SIZE_FLOOR, repr = "consensus_dosage")
Pu <- emmax_setup(um, GRM)
pu_obs <- emmax_fast(Pu, eco)

test_con <- ld_outlier_test(stage1, map, pu_obs, statistic = "unit", size_floor = SIZE_FLOOR,
                            alpha = ALPHA, assembly = "stage2_discovered", GTs = GTs,
                            LD_decay = b$LD_decay, score_threshold = REGION_ASSEMBLY$score_threshold,
                            distance_threshold = REGION_ASSEMBLY$distance_threshold)
print(test_con)

say("    permutation null: %d regional surrogates (consensus is cheap)\n", NPERM_CONSENSUS)
p_perm_con <- function(bb) { set.seed(bb); y <- perm_regional()
  emmax_fast(Pu, y) }
null_con <- ld_outlier_perm(test_con, stage1, map, p_perm_con, GTs = GTs, LD_decay = b$LD_decay,
                            B = NPERM_CONSENSUS, level = "units", verbose = TRUE)
print(null_con)

say("    rotation null against real EcoPeaks: %s draws, scheme = %s\n", format(N_ROTATIONS, big.mark=","), ROTATION_SCHEME)
rot_con <- ld_region_rotation(test_con$regions, KG, LEN, scheme = ROTATION_SCHEME,
                              n_rotations = N_ROTATIONS, seed = 1L)
print(rot_con)
say("    %.1f min\n\n", as.numeric(difftime(Sys.time(), t0, units="mins")))
RESULTS$consensus <- list(test = test_con, null = null_con, rotation = rot_con)

## ---- 2. Simes arm (statistic = "simes") -------------------------------------
say("[2] SIMES arm -- direct marker-level scan aggregated per unit\n")
t0 <- Sys.time()
Pm <- emmax_setup(GTs, GRM)
pm_obs <- emmax_fast(Pm, eco)

test_sim <- ld_outlier_test(stage1, map, pm_obs, statistic = "simes", size_floor = SIZE_FLOOR,
                            alpha = ALPHA, assembly = "stage2_discovered", GTs = GTs,
                            LD_decay = b$LD_decay, score_threshold = REGION_ASSEMBLY$score_threshold,
                            distance_threshold = REGION_ASSEMBLY$distance_threshold)
print(test_sim)

say("    permutation null: %d regional surrogates (Simes rescans all %s markers per draw)\n",
    NPERM_SIMES, format(ncol(GTs), big.mark=","))
p_perm_sim <- function(bb) { set.seed(bb); y <- perm_regional(); emmax_fast(Pm, y) }
null_sim <- ld_outlier_perm(test_sim, stage1, map, p_perm_sim, GTs = GTs, LD_decay = b$LD_decay,
                            B = NPERM_SIMES, level = "units", verbose = TRUE)
print(null_sim)

rot_sim <- ld_region_rotation(test_sim$regions, KG, LEN, scheme = ROTATION_SCHEME,
                              n_rotations = N_ROTATIONS, seed = 1L)
print(rot_sim)
say("    %.1f min\n\n", as.numeric(difftime(Sys.time(), t0, units="mins")))
RESULTS$simes <- list(test = test_sim, null = null_sim, rotation = rot_sim)

## ---- 3. agreement between the two statistics --------------------------------
sig_con <- test_con$units[significant == TRUE]$unit_id
sig_sim <- test_sim$units[significant == TRUE]$unit_id
jac <- length(intersect(sig_con, sig_sim)) / length(union(sig_con, sig_sim))
say("[3] agreement: consensus %d, Simes %d, shared %d, Jaccard %.3f\n",
    length(sig_con), length(sig_sim), length(intersect(sig_con, sig_sim)), jac)

## ---- 4. save ------------------------------------------------------------------
OUT <- file.path(stage_dir(STAGE), "scan.rds")
saveRDS(list(consensus = RESULTS$consensus, simes = RESULTS$simes,
            jaccard = jac, ecopeaks = KG, chrom_lengths = LEN), OUT)
write_receipt(STAGE, inputs = INPUTS, params = PARAMS, outputs = OUT)
say("\n[4] wrote %s\n    receipt: %s\n", OUT, receipt_path(STAGE))
