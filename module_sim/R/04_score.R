## module_sim/R/04_score.R
##
## Per-replicate TP/FP/FN scoring: for one (tag,cell,rep), classify each
## engine/arm's SIGNIFICANT stage-1 clusters (from R/03_scan.R) as true or
## false positive against the simulated QTN, using the package's own
## outlier-region benchmark (flag_true_qtns/qtn_ld_table/score_thresholds/
## evaluate_ors, ~/gitlab/LDscnR/R/ld_benchmark.R) -- the same machinery
## module_sim_LDscnR already uses, reused rather than reinvented.
##
## PK, confirming the design: "each stage 1 cluster gets a TP if it contains
## a SNP with ld_rho>0.75 and d_rho<0.95 with a detectable QTN" -- this is
## exactly score_thresholds(decay_sum, rho_r2=0.75, rho_d=0.95) (its own
## defaults; rho_r2/rho_d are decay-quantiles fed to ld_from_rho()/d_from_rho(),
## NOT literal r^2/bp cutoffs) applied to evaluate_ors() with the SIGNIFICANT
## stage-1 clusters (per engine/arm, at 03_scan.R's fixed alpha=0.05) as the
## discovered regions.
##
## [!] Single alpha, not a threshold sweep. 03_scan.R only tests at alpha=0.05
## (no tau_C/rho sweep exists in this pipeline), so this produces ONE
## (TP,FP,FN,Precision,Recall) operating point per replicate per engine/arm --
## not a PR-AUC curve. R/05_pool.R replicate-averages these across the 10
## reps per (tag,cell), matching every other PR-scoring script in this repo's
## house convention (mean +- SE across replicates -- see module_sim_LDscnR/
## run_sim_LDscnR.R's "ALWAYS replicate-average" rule).
##
## dmax_cap = 5e5: the same cap already standing everywhere else in the repo
## (see memory: "dcap 5e5 already universal").
suppressMessages({library(data.table); library(LDscnR)})
devtools::load_all("~/gitlab/LDscnR")
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim"), "R", "00_config.R"))
STAGE <- "04_score"
say("=== %s ===\n\n", STAGE)
invisible(check_ldscnr())

TARGET_TAG <- Sys.getenv("SIM_TAG", TAGS[1]); TARGET_CELL <- Sys.getenv("SIM_CELL", CELLS[1]); TARGET_ENV <- as.integer(Sys.getenv("SIM_ENV", ENVS[1])); TARGET_REP <- as.integer(Sys.getenv("SIM_REP", REPS[1]))
combo_id <- sprintf("%s_%s_rep%d", TARGET_TAG, TARGET_CELL, TARGET_REP)

BUNDLE_PATH <- file.path(PATHS$out, "02_bundle",
  sprintf("bundle_%s_rep%d_%s_env%d.rds", TARGET_TAG, TARGET_REP, TARGET_CELL, TARGET_ENV))
SCAN_PATH <- file.path(PATHS$out, "03_scan",
  sprintf("scan_%s_rep%d_%s_env%d.rds", TARGET_TAG, TARGET_REP, TARGET_CELL, TARGET_ENV))
if (!file.exists(BUNDLE_PATH)) stop("R/02_bundle.R has not produced: ", basename(BUNDLE_PATH))
if (!file.exists(SCAN_PATH))   stop("R/03_scan.R has not produced: ", basename(SCAN_PATH))
bd <- readRDS(BUNDLE_PATH)
sc <- readRDS(SCAN_PATH)
GTs <- bd$GTs; map <- bd$map; stage1 <- bd$stage1; LD_decay <- bd$LD_decay

DMAX_CAP <- 5e5
RHO_R2 <- 0.75; RHO_D <- 0.95   ## score_thresholds()'s own defaults -- PK confirmed these are the intended values

INPUTS <- c(BUNDLE_PATH, SCAN_PATH)
PARAMS <- list(size_floor = SIZE_FLOOR, dmax_cap = DMAX_CAP, rho_r2 = RHO_R2, rho_d = RHO_D,
               maf_min = 0.1, p_va_min = 0.05, max_bp = 2e6)
if (!stage_stale(STAGE, INPUTS, PARAMS, target = combo_id) && !nzchar(Sys.getenv("FORCE"))) {
  say("\nNothing to do. Set FORCE=1 to rerun anyway.\n"); quit(save = "no")
}

## ---- 1. Va, true_pos_QTN -------------------------------------------------------
## Va = 2*p*(1-p)*allelic_value^2 (standard additive-variance contribution,
## p from the ACTUAL sampled genotypes not the recorded MAF) -- the same
## formula module_sim_LDscnR/parse_and_regen_sim_data.R's get_va() uses.
say("[1] Va, true_pos_QTN (MAF > 0.1, per-chromosome Va share >= 0.05)\n")
qtn_rows <- which(map$type == "QTN")
map[, Va := NA_real_]
map$Va[qtn_rows] <- vapply(qtn_rows, function(i) {
  p <- mean(GTs[, i]) / 2
  2 * p * (1 - p) * map$allelic_values[i]^2
}, 0)
map <- flag_true_qtns(map, va_col = "Va", maf_col = "MAF", maf_min = 0.1, p_va_min = 0.05)
say("    %d of %d QTN pass detectability (true_pos_QTN)\n", sum(map$true_pos_QTN), length(qtn_rows))

## ---- 2. units (same filtered/ordered cluster list 03_scan.R's ld_outlier_test used) --
## .ld_outlier_units() is internal but deterministic given (stage1, map, size_floor) --
## the SAME call ld_outlier_test() makes inside .ld_outlier_tested_units(), so row i
## here is unit_id i in R/03_scan.R's test$units, letting the two align by position
## without reconstructing a CL_id<->unit_id mapping by hand.
units <- .ld_outlier_units(stage1, map, SIZE_FLOOR)
say("[2] %d stage-1 units clear size_floor = %d\n", nrow(units), SIZE_FLOOR)

## ---- 3. QTN-LD lookup + decay-relative match thresholds -----------------------
say("[3] qtn_ld_table + score_thresholds (rho_r2 = %.2f, rho_d = %.2f, dmax_cap = %s)\n",
    RHO_R2, RHO_D, format(DMAX_CAP, big.mark = ","))
candidate_markers <- unique(unlist(units$members, use.names = FALSE))
qtab <- qtn_ld_table(GTs, map, candidate_markers, max_bp = 2e6)
th <- score_thresholds(as.data.table(LD_decay$decay_sum), rho_r2 = RHO_R2, rho_d = RHO_D, dmax_cap = DMAX_CAP)
say("    r2min = %.4f, dmax = %s bp\n", th$r2min, format(round(th$dmax), big.mark = ","))

## ---- 4. evaluate_ors() per engine/arm, significant units only -----------------
ARMS <- list(emmax_consensus = sc$results$emmax$consensus$test,
            emmax_simes     = sc$results$emmax$simes$test,
            lfmm_consensus  = sc$results$lfmm$consensus$test,
            lfmm_simes      = sc$results$lfmm$simes$test)
say("\n[4] evaluate_ors() per engine/arm\n")
scored <- rbindlist(lapply(names(ARMS), function(nm) {
  t <- ARMS[[nm]]
  stopifnot("unit_id must align by position between R/03_scan.R and .ld_outlier_units() here" =
              nrow(t$units) == nrow(units), identical(t$units$unit_id, units$unit_id))
  sig_regions <- units$members[t$units$significant]
  ev <- evaluate_ors(sig_regions, map, qtab, r2_match = th$r2min, d_match = th$dmax)
  say("    %-16s sig=%-4d TP=%-3d FP=%-3d FN=%-3d Precision=%s Recall=%.3f\n",
      nm, length(sig_regions), ev$TP, ev$FP, ev$FN,
      if (is.na(ev$Precision)) "  NA" else sprintf("%.3f", ev$Precision), ev$Recall)
  data.table(tag = TARGET_TAG, cell = TARGET_CELL, rep = TARGET_REP, arm = nm,
             n_sig = length(sig_regions), TP = ev$TP, FP = ev$FP, FN = ev$FN,
             Precision = ev$Precision, Recall = ev$Recall, PR = ev$PR)
}))

OUT <- file.path(stage_dir(STAGE), sprintf("score_%s_rep%d_%s_env%d.rds", TARGET_TAG, TARGET_REP, TARGET_CELL, TARGET_ENV))
dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)
saveRDS(scored, OUT)
invisible(check_ldscnr())
write_receipt(STAGE, inputs = INPUTS, params = PARAMS, outputs = OUT, target = combo_id)
say("\n[5] wrote %s\n    receipt: %s\n", OUT, receipt_path(STAGE, combo_id))
