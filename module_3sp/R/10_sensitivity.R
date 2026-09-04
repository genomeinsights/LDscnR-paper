## =============================================================================
## module_3sp/R/10_sensitivity.R
##
## ANALYSIS ONLY. Sweeps the two parameters PK identified as the ones that
## matter: rho and size_floor. rho affects BOTH Stage 1 (ld_complexity_
## reduction's clustering threshold) AND Stage 2 (ld_outlier_test's
## stage2_discovered assembly reads min_r2_rho FROM stage1$params$rho, never
## as an independent argument -- see ld_outlier_test.R's own doc) -- ONE knob
## governs both stages, by the package's own design, so sweeping it here
## covers both without a second, separate Stage-2 sweep. LDW_FLAG,
## SCORE_THRESHOLD and DISTANCE_THRESHOLD are Stage-2-only and, per the
## orthogonality discussion, cannot change which clusters are significant
## (only how already-significant ones get merged into regions) -- explicitly
## OUT of scope here, not an oversight.
##
## GRID (PK): rho in {0.3, 0.5, 0.7}, size_floor in {4, 8, 12, 16}. 0.5/8 is
## the canonical point already computed by 03_EMMAX.R/04_lfmm.R -- reused
## here, not rerun, so this script adds exactly 2 (rho) + 3 (size_floor) = 5
## new grid points, each scored on all three analyses (EMMAX-consensus,
## EMMAX-Simes, LFMM-Simes).
##
## COST ASYMMETRY: a rho value changes which markers cluster together, so it
## needs a real re-cluster (ld_complexity_reduction()) and a rebuilt GRM (the
## kinship basis is stage1's own core_snp representatives -- a different rho
## means a different representative set). size_floor only refilters which
## ALREADY-clustered units count as tested -- same stage1, same GRM, so it
## only needs ld_unit_matrix()/ld_outlier_test() rebuilt, no reclustering.
## Genome-wide reclustering reuses b$LD_decay as-is (no decay refit -- rho
## only changes the threshold applied to the SAME fitted curve), so this is
## "minutes" per 02_bundle.R's own checkpoint comment ("hours for the decay
## fit, minutes for the clustering"), not a second multi-hour run.
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(SNPRelate)})
devtools::load_all("~/gitlab/LDscnR")
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_3sp"), "R", "00_config.R"))
STAGE <- "10_sensitivity"
say("=== %s ===\n\n", STAGE)
invisible(check_ldscnr())

BUNDLE_PATH <- file.path(PATHS$out, "02_bundle", "bundle.rds")
b <- readRDS(BUNDLE_PATH)
GTs <- b$GTs; map <- b$map; eco <- b$eco; pheno <- b$pheno
SCAN_PATH <- file.path(PATHS$out, "03_EMMAX", "scan.rds")
LFMM_PATH <- file.path(PATHS$out, "04_lfmm", "lfmm.rds")
sc <- readRDS(SCAN_PATH); lf <- readRDS(LFMM_PATH)

LIFT <- path.expand("~/gitlab/LDscnR-paper/kingman2021/data/liftover")
ROMAN <- c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV",
           "XV","XVI","XVII","XVIII","XIX","XX","XXI")
.rd <- function(f) { x <- fread(file.path(LIFT, f), header=FALSE,
                                col.names=c("chr","start","end","pv"))
  x[, chr_num := match(sub("^chr","",chr), ROMAN)]
  x[!is.na(chr_num), .(chr = paste0("Chr", chr_num), start, end)] }
KG <- unique(rbindlist(lapply(ECOPEAK_BEDS, .rd)))
LEN <- map[, .(len = max(Pos)), by = Chr]; setnames(LEN, "Chr", "chr")

POPT <- unique(pheno[, .(pop_ID, ecotype, pop_locality)])
perm_regional <- function() { pt <- copy(POPT)[, ep := sample(ecotype), by = pop_locality]
  as.numeric(pt$ep[match(pheno$pop_ID, pt$pop_ID)] == "Marine") }

RHO_GRID <- c(0.3, 0.7)
SF_GRID  <- c(4, 12, 16, 24, 48)

INPUTS <- c(BUNDLE_PATH, file.path(LIFT, ECOPEAK_BEDS), SCAN_PATH, LFMM_PATH)
PARAMS <- list(rho_grid = RHO_GRID, size_floor_grid = SF_GRID, alpha = ALPHA,
               score_threshold = SCORE_THRESHOLD, distance_threshold = DISTANCE_THRESHOLD,
               unit_repr = UNIT_REPR, nperm_consensus = NPERM_CONSENSUS,
               nperm_simes = NPERM_SIMES, n_rotations = N_ROTATIONS, rotation_scheme = ROTATION_SCHEME)
if (!stage_stale(STAGE, INPUTS, PARAMS) && !nzchar(Sys.getenv("FORCE"))) {
  say("\nNothing to do. Set FORCE=1 to rerun anyway.\n"); quit(save = "no")
}
dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)

## ---- reuse already-computed grid points -------------------------------------------------
## Extending SF_GRID (or RHO_GRID) changes PARAMS, so stage_stale() alone can't skip a rerun
## -- but re-running the WHOLE sweep would redo the expensive rho reclustering + GRM rebuild
## for no reason. Grid points already present in a prior sensitivity.rds are reused verbatim;
## only genuinely new (rho, size_floor) combinations get computed.
OLD_RDS <- file.path(stage_dir(STAGE), "sensitivity.rds")
OLD_RESULTS <- if (file.exists(OLD_RDS)) readRDS(OLD_RDS)$results else list()
if (length(OLD_RESULTS)) say("[0] reusing %d already-computed grid point(s) from %s\n",
                             length(OLD_RESULTS), OLD_RDS)

## ---- run all three analyses for one (stage1, GRM, size_floor) grid point ---------------
run_all <- function(stage1, GRM, size_floor, label) {
  t0 <- Sys.time()
  say("  [%s] size_floor = %d\n", label, size_floor)

  ## EMMAX consensus
  um <- ld_unit_matrix(GTs, stage1, map, size_floor = size_floor, repr = UNIT_REPR)
  Pu <- emmax_setup(um, GRM); pu_obs <- emmax_fast(Pu, eco)
  test_con <- ld_outlier_test(stage1, map, pu_obs, statistic = "unit", size_floor = size_floor,
                              alpha = ALPHA, assembly = "stage2_discovered", GTs = GTs,
                              LD_decay = b$LD_decay, score_threshold = SCORE_THRESHOLD,
                              distance_threshold = DISTANCE_THRESHOLD)
  p_perm_con <- function(bb) { set.seed(bb); emmax_fast(Pu, perm_regional()) }
  null_con <- ld_outlier_perm(test_con, stage1, map, p_perm_con, GTs = GTs, LD_decay = b$LD_decay,
                              B = NPERM_CONSENSUS, level = "units", verbose = FALSE)
  rot_con <- ld_region_rotation(test_con$regions, KG, LEN, scheme = ROTATION_SCHEME,
                                n_rotations = N_ROTATIONS, seed = 1L)
  say("    consensus: %d/%d sig -> %d regions (%.1f min)\n", sum(test_con$units$significant),
      nrow(test_con$units), nrow(test_con$regions), as.numeric(difftime(Sys.time(), t0, units="mins")))

  ## EMMAX Simes
  t1 <- Sys.time()
  Pm <- emmax_setup(GTs, GRM); pm_obs <- emmax_fast(Pm, eco)
  test_sim <- ld_outlier_test(stage1, map, pm_obs, statistic = "simes", size_floor = size_floor,
                              alpha = ALPHA, assembly = "stage2_discovered", GTs = GTs,
                              LD_decay = b$LD_decay, score_threshold = SCORE_THRESHOLD,
                              distance_threshold = DISTANCE_THRESHOLD)
  p_perm_sim <- function(bb) { set.seed(bb); emmax_fast(Pm, perm_regional()) }
  null_sim <- ld_outlier_perm(test_sim, stage1, map, p_perm_sim, GTs = GTs, LD_decay = b$LD_decay,
                              B = NPERM_SIMES, level = "units", verbose = FALSE)
  rot_sim <- ld_region_rotation(test_sim$regions, KG, LEN, scheme = ROTATION_SCHEME,
                                n_rotations = N_ROTATIONS, seed = 1L)
  say("    Simes:     %d/%d sig -> %d regions (%.1f min)\n", sum(test_sim$units$significant),
      nrow(test_sim$units), nrow(test_sim$regions), as.numeric(difftime(Sys.time(), t1, units="mins")))

  ## LFMM Simes -- lfmm_p is inherited/precomputed, identical across the whole sweep;
  ## only its aggregation into THIS stage1/size_floor's units changes. No permutation null
  ## (LFMM has none, same reason as 04_lfmm.R).
  t2 <- Sys.time()
  test_lfmm <- ld_outlier_test(stage1, map, lf$lfmm_p, statistic = "simes", size_floor = size_floor,
                               alpha = ALPHA, assembly = "stage2_discovered", GTs = GTs,
                               LD_decay = b$LD_decay, score_threshold = SCORE_THRESHOLD,
                               distance_threshold = DISTANCE_THRESHOLD)
  rot_lfmm <- ld_region_rotation(test_lfmm$regions, KG, LEN, scheme = ROTATION_SCHEME,
                                 n_rotations = N_ROTATIONS, seed = 1L)
  say("    LFMM:      %d/%d sig -> %d regions (%.1f min)\n", sum(test_lfmm$units$significant),
      nrow(test_lfmm$units), nrow(test_lfmm$regions), as.numeric(difftime(Sys.time(), t2, units="mins")))

  list(consensus = list(test = test_con, null = null_con, rotation = rot_con),
       simes     = list(test = test_sim, null = null_sim, rotation = rot_sim),
       lfmm      = list(test = test_lfmm, rotation = rot_lfmm))
}

RESULTS <- list()

## ---- canonical point: REUSED from 03_EMMAX.R/04_lfmm.R, not recomputed -----------------
say("[1] rho=0.50, size_floor=8 (canonical) -- reused from 03_EMMAX.R/04_lfmm.R\n")
RESULTS[["rho=0.50_sf=8"]] <- list(
  consensus = sc$consensus, simes = sc$simes,
  lfmm = list(test = lf$test, rotation = lf$rotation))

## ---- size_floor sweep (rho fixed at canonical 0.5) -- cheap, no reclustering -----------
say("\n[2] size_floor sweep (rho = 0.50 canonical, reusing b$stage1/b$GRM)\n")
for (sf in SF_GRID) {
  key <- sprintf("rho=0.50_sf=%d", sf)
  if (!is.null(OLD_RESULTS[[key]])) {
    say("  [sf=%d] reusing already-computed result\n", sf)
    RESULTS[[key]] <- OLD_RESULTS[[key]]
  } else {
    RESULTS[[key]] <- run_all(b$stage1, b$GRM, sf, sprintf("sf=%d", sf))
  }
}

## ---- rho sweep (size_floor fixed at canonical 8) -- real re-cluster + GRM rebuild ------
say("\n[3] rho sweep (size_floor = 8 canonical)\n")
rho_keys <- sprintf("rho=%.2f_sf=8", RHO_GRID)
if (all(rho_keys %in% names(OLD_RESULTS))) {
  say("  all %d rho grid point(s) already computed -- reusing, no reclustering\n", length(RHO_GRID))
  for (i in seq_along(RHO_GRID)) RESULTS[[rho_keys[i]]] <- OLD_RESULTS[[rho_keys[i]]]
} else {
  gds_sens_path <- file.path(PATHS$cache, "3sp_sensitivity.gds")
  if (file.exists(gds_sens_path)) unlink(gds_sens_path)
  gds_sens <- create_gds_from_geno(geno = GTs, map = map, gds_sens_path)
  on.exit(try(SNPRelate::snpgdsClose(gds_sens), silent = TRUE), add = TRUE)

  for (rho_alt in RHO_GRID) {
    key <- sprintf("rho=%.2f_sf=8", rho_alt)
    if (!is.null(OLD_RESULTS[[key]])) {
      say("  rho = %.1f: reusing already-computed result\n", rho_alt)
      RESULTS[[key]] <- OLD_RESULTS[[key]]
      next
    }
    t0 <- Sys.time()
    say("  rho = %.1f: reclustering (ld_complexity_reduction, gds=genome-wide)\n", rho_alt)
    set.seed(SEEDS[["clusters"]])
    stage1_alt <- ld_complexity_reduction(map = map, LD_decay = b$LD_decay, rho = rho_alt, gds = gds_sens)
    grm_markers_alt <- unique(na.omit(stage1_alt$pruned))
    GRM_alt <- snpgdsGRM(gds_sens, snp.id = grm_markers_alt, method = GRM_METHOD,
                         verbose = FALSE, autosome.only = FALSE)$grm
    cl <- as.data.table(stage1_alt$clusters)
    nl <- if ("n_loci" %in% names(cl)) cl$n_loci else cl$n_snps
    say("    %s clusters ; %s at floor %d ; GRM from %s markers (%.1f min to build)\n",
        format(nrow(cl), big.mark=","), format(sum(nl >= SIZE_FLOOR), big.mark=","), SIZE_FLOOR,
        format(length(grm_markers_alt), big.mark=","),
        as.numeric(difftime(Sys.time(), t0, units="mins")))
    RESULTS[[key]] <- run_all(stage1_alt, GRM_alt, SIZE_FLOOR, sprintf("rho=%.1f", rho_alt))
  }
}

## ---- summary table -----------------------------------------------------------------
say("\n[4] summary\n")
SUMMARY <- rbindlist(lapply(names(RESULTS), function(lab) {
  r <- RESULTS[[lab]]
  rbindlist(lapply(c("consensus","simes","lfmm"), function(arm) {
    a <- r[[arm]]
    data.table(grid = lab, analysis = arm,
              tested = nrow(a$test$units), significant = sum(a$test$units$significant),
              regions = nrow(a$test$regions), on_ecopeak = a$rotation$observed,
              rotation_fold = a$rotation$fold, rotation_p = a$rotation$p,
              perm_p = if (!is.null(a$null)) a$null$p else NA_real_,
              realised_fdr = if (!is.null(a$null)) a$null$realised_fdr else NA_real_)
  }))
}))
print(SUMMARY)

## ---- save + receipt ------------------------------------------------------------------
OUT <- file.path(stage_dir(STAGE), "sensitivity.rds")
saveRDS(list(results = RESULTS, summary = SUMMARY, rho_grid = RHO_GRID, sf_grid = SF_GRID), OUT)
fwrite(SUMMARY, file.path(stage_dir(STAGE), "sensitivity_summary.csv"))
write_receipt(STAGE, inputs = INPUTS, params = PARAMS,
             outputs = c(OUT, file.path(stage_dir(STAGE), "sensitivity_summary.csv")))
say("\n[5] wrote %s\n    receipt: %s\n", OUT, receipt_path(STAGE))
