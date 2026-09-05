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
## (TP,FP,FN,Precision,Recall) operating point per (rep,env) per engine/arm --
## not a PR-AUC curve. R/05_pool.R POOLS these (sums TP/FP/FN, then computes
## one Precision/Recall from the totals) across ENV within each REP, and
## again across REP, per PK's repeated instruction that PR is estimated on
## all chromosomes pooled -- not the mean of each combination's own
## Precision/Recall.
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
combo_id <- sprintf("%s_%s_rep%d_env%d", TARGET_TAG, TARGET_CELL, TARGET_REP, TARGET_ENV)

BUNDLE_PATH <- file.path(PATHS$out, "02_bundle",
  sprintf("bundle_%s_rep%d_%s_env%d.rds", TARGET_TAG, TARGET_REP, TARGET_CELL, TARGET_ENV))
SCAN_PATH <- file.path(PATHS$out, "03_scan",
  sprintf("scan_%s_rep%d_%s_env%d.rds", TARGET_TAG, TARGET_REP, TARGET_CELL, TARGET_ENV))
if (!file.exists(BUNDLE_PATH)) stop("R/02_bundle.R has not produced: ", basename(BUNDLE_PATH))
if (!file.exists(SCAN_PATH))   stop("R/03_scan.R has not produced: ", basename(SCAN_PATH))
bd <- readRDS(BUNDLE_PATH)
sc <- readRDS(SCAN_PATH)
GTs <- bd$GTs; map <- bd$map; stage1 <- bd$stage1; LD_decay <- bd$LD_decay

## [!] WAS 5e5, hardcoded separately from DISTANCE_THRESHOLD. Fixed 2026-09-05:
## a peer session (NEMO simulation work) flagged that DCAP moved to 1e5
## paper-wide (LDscnR-paper 8dbb09a, branch null-rework-two-basis), matching
## the bundles' stage-2 distance_threshold -- module_sim's own
## DISTANCE_THRESHOLD (00_config.R) was already 1e5, but this file's
## qtn_ld_table/score_thresholds cap was a separate, stale constant that
## predated the change. Referencing DISTANCE_THRESHOLD directly instead of a
## second hardcoded value, so the two cannot drift apart again.
DMAX_CAP <- DISTANCE_THRESHOLD
RHO_R2 <- 0.75; RHO_D <- 0.95   ## score_thresholds()'s own defaults -- PK confirmed these are the intended values

INPUTS <- c(BUNDLE_PATH, SCAN_PATH)
PARAMS <- list(size_floor = SIZE_FLOOR, dmax_cap = DMAX_CAP, rho_r2 = RHO_R2, rho_d = RHO_D,
               maf_min = 0.1, p_va_min = 0.05, max_bp = 2e6,
               single_snp_arms = TRUE, single_snp_bh_alpha = ALPHA, lfmm_consensus_arm = FALSE,
               cluster_detail = TRUE,
               ## [!] bumped 2026-09-06: single-SNP arms now score each significant
               ## marker as its own size-1 region (external audit item 1 / PK), a
               ## logic change stage_stale() cannot see any other way -- forces a
               ## rerun instead of silently reusing pre-fix score files.
               single_snp_region_size = "per_marker_size1",
               ## [!] bumped again 2026-09-06: added emmax_snp_clustered/
               ## lfmm_snp_clustered, the singletons-EXCLUDED variant, alongside
               ## emmax_snp/lfmm_snp (singletons-INCLUDED) -- PK: "single SNP
               ## included/excluded in the single SNP analyses, different colors".
               single_snp_clustered_variant = TRUE)
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

## [!] NO QTN WITHIN max_bp OF ANY CANDIDATE MARKER -- a real edge case, not
## a bug: found on a V2 (weak selection) replicate whose one QTN existed but
## fell below the MAF>0.1 detectability filter, leaving true_pos_QTN entirely
## empty. qtn_ld_table() then returns a genuinely EMPTY, COLUMNLESS
## data.table (rbindlist() of an all-NULL list has no r2/dist_bp columns at
## all), and both evaluate_ors() and .diagnose_ors() crash trying to filter
## a column that does not exist -- as soon as ANY arm has >=1 significant
## region (the two zero-significant arms upstream of the crash short-circuit
## before ever touching qtab). No possible TP can exist in this case by
## construction (nothing is in LD with a QTN), so every significant region
## is trivially FP; handled directly here rather than by patching the
## package for a table shape it was never designed to receive.
NO_QTN_POSSIBLE <- nrow(qtab) == 0
if (NO_QTN_POSSIBLE) say("    [!] qtn_ld_table is empty (no QTN within max_bp of any candidate marker) -- every significant region will score as FP\n")
N_TRUE <- sum(map$true_pos_QTN, na.rm = TRUE)

## ---- 4. score one arm: aggregate (evaluate_ors) + per-cluster detail ----------
## PK: also want FP proportion as a function of cluster size -- evaluate_ors()
## alone only returns the aggregate TP/FP/FN, discarding which SIGNIFICANT
## cluster was which size and TP/FP. .diagnose_ors() is the internal function
## evaluate_ors() itself calls (~/gitlab/LDscnR/R/ld_benchmark.R) -- same
## dedup logic, same TP/FP definition, but keeps one row per region with
## n_loci (cluster size) and is_TP, plus dropped_by_dedup/
## candidate_qtn_is_true_positive to reconstruct evaluate_ors()'s
## dedup-neutral "extra" (neither TP nor FP) exactly. Called once per arm
## here instead of building the aggregate a second, differently-scoped way.
.score_arm <- function(sig_regions, nm) {
  if (NO_QTN_POSSIBLE) {
    n_sig <- length(sig_regions)
    ev <- list(TP = 0L, FP = n_sig, FN = N_TRUE,
               Precision = if (n_sig > 0) 0 else NA_real_,
               Recall = if (N_TRUE > 0) 0 else NA_real_, PR = NA_real_)
    detail <- data.table(tag = TARGET_TAG, cell = TARGET_CELL, rep = TARGET_REP, env = TARGET_ENV, arm = nm,
                         n_loci = vapply(sig_regions, length, integer(1)),
                         is_TP = rep(FALSE, n_sig), is_FP = rep(TRUE, n_sig))
  } else {
    ev <- evaluate_ors(sig_regions, map, qtab, r2_match = th$r2min, d_match = th$dmax)
    detail <- if (length(sig_regions)) {
      dg <- .diagnose_ors(sig_regions, map, qtab, r2_match = th$r2min, d_match = th$dmax)
      extra <- dg$dropped_by_dedup == TRUE & dg$candidate_qtn_is_true_positive %in% TRUE
      data.table(tag = TARGET_TAG, cell = TARGET_CELL, rep = TARGET_REP, env = TARGET_ENV, arm = nm,
                 n_loci = dg$n_loci, is_TP = dg$is_TP, is_FP = !dg$is_TP & !extra)
    } else {
      data.table(tag = character(), cell = character(), rep = integer(), env = integer(), arm = character(),
                 n_loci = integer(), is_TP = logical(), is_FP = logical())
    }
  }
  agg <- data.table(tag = TARGET_TAG, cell = TARGET_CELL, rep = TARGET_REP, env = TARGET_ENV, arm = nm,
                    n_sig = length(sig_regions), TP = ev$TP, FP = ev$FP, FN = ev$FN,
                    Precision = ev$Precision, Recall = ev$Recall, PR = ev$PR)
  say("    %-16s sig=%-4d TP=%-3d FP=%-3d FN=%-3d Precision=%s Recall=%s\n",
      nm, length(sig_regions), ev$TP, ev$FP, ev$FN,
      if (is.na(ev$Precision)) "  NA" else sprintf("%.3f", ev$Precision),
      if (is.na(ev$Recall)) "  NA" else sprintf("%.3f", ev$Recall))
  list(agg = agg, detail = detail)
}

## ---- 4a. cluster-based arms, significant units only ---------------------------
ARMS <- list(emmax_consensus = sc$results$emmax$consensus$test,
            emmax_simes     = sc$results$emmax$simes$test,
            lfmm_simes      = sc$results$lfmm$simes$test)
say("\n[4] cluster-based arms\n")
res_cluster <- lapply(names(ARMS), function(nm) {
  t <- ARMS[[nm]]
  stopifnot("unit_id must align by position between R/03_scan.R and .ld_outlier_units() here" =
              nrow(t$units) == nrow(units), identical(t$units$unit_id, units$unit_id))
  .score_arm(units$members[t$units$significant], nm)
})

## ---- 4b. single-SNP arms: each significant marker is its OWN size-1 region ----
## PK: compare against single-SNP analyses using the same machinery. BH-correct
## the single-marker p-values genome-wide (not per-unit -- there is no unit
## combination here, that is the whole point of "single-SNP"), then EACH
## significant marker is its own discovered region of size 1 -- NOT the whole
## enclosing Stage-1 unit it happens to sit inside. .score_arm() -> evaluate_ors()/
## .diagnose_ors() are agnostic to region size (n_loci = length(region) falls
## out naturally as 1), so this needs no package change, only feeding a
## different `regions` list in.
##
## [!] FIXED 2026-09-06 (external audit item 1 + PK direct correction, "Single
## SNPs go all the way to 1 not 2. Unless that one is a QTN it is an FP."):
## the previous version credited the FULL enclosing Stage-1 unit (>= SIZE_FLOOR
## markers, since `units` above was already floor-filtered) as "discovered"
## whenever ANY one of its markers was individually significant, then scored
## that whole multi-marker unit's TP/FP status. That (a) let a single
## significant SNP borrow TP credit from neighbouring markers it never tested,
## inflating single-SNP precision/recall relative to a genuine unrestricted
## single-marker benchmark, (b) meant a truly isolated significant marker --
## one whose own Stage-1 unit does not clear SIZE_FLOOR, so it never appears in
## `units` at all -- was invisible to scoring altogether (neither a TP nor an
## FP), and (c) meant cluster_detail's n_loci for these arms was never 1,
## silently excluding single-SNP arms from the smallest FP-by-size bin. Feeding
## singleton regions through the SAME .score_arm() keeps the TP/FP definition
## (and its dedup-neutral handling) identical across all five arms -- the
## comparison is single-SNP vs complexity-reduced significance calling, not two
## different scoring rules.
say("\n[4b] single-SNP arms (BH genome-wide, each significant marker its own size-1 region)\n")
res_snp <- Map(function(p_vec, nm) {
  q <- stats::p.adjust(p_vec, method = "BH")
  sig_markers <- names(q)[!is.na(q) & q <= ALPHA]
  sig_regions <- as.list(sig_markers)
  .score_arm(sig_regions, nm)
}, list(sc$results$single_snp$emmax_p, sc$results$single_snp$lfmm_p), list("emmax_snp", "lfmm_snp"))

## ---- 4c. single-SNP arms, CLUSTERED-ONLY variant (singletons excluded) --------
## ADDED 2026-09-06 (PK: "redo also some of the other precision, recall
## analyses with single SNP included/excluded in the single SNP analyses,
## different colors"). This is deliberately the PRE-FIX behavior from 4b
## above, kept as an explicit second variant rather than dropped: a
## significant marker only counts if its own Stage-1 unit clears SIZE_FLOOR
## (>=2 markers), and the region scored is that WHOLE unit, not just the one
## marker -- truly isolated (size-1) significant markers are excluded
## entirely, neither TP nor FP. Comparing `emmax_snp`/`lfmm_snp` (singletons
## INCLUDED, 4b) against `emmax_snp_clustered`/`lfmm_snp_clustered` (singletons
## EXCLUDED, here) shows how much of the unrestricted single-marker arm's
## apparent performance comes from truly-isolated calls vs. markers already
## embedded in a real LD cluster.
say("\n[4c] single-SNP arms, clustered-only variant (BH genome-wide, unit = discovered iff it contains a significant SNP)\n")
res_snp_clustered <- Map(function(p_vec, nm) {
  q <- stats::p.adjust(p_vec, method = "BH")
  sig_markers <- names(q)[!is.na(q) & q <= ALPHA]
  sig_regions <- units$members[vapply(units$members, function(mk) any(mk %in% sig_markers), logical(1))]
  .score_arm(sig_regions, nm)
}, list(sc$results$single_snp$emmax_p, sc$results$single_snp$lfmm_p),
   list("emmax_snp_clustered", "lfmm_snp_clustered"))

res_all <- c(res_cluster, res_snp, res_snp_clustered)
scored         <- rbindlist(lapply(res_all, `[[`, "agg"))
cluster_detail <- rbindlist(lapply(res_all, `[[`, "detail"))
say("\n[4d] %d significant clusters/units captured across %d arms for FP-vs-size analysis\n",
    nrow(cluster_detail), length(res_all))

OUT <- file.path(stage_dir(STAGE), sprintf("score_%s_rep%d_%s_env%d.rds", TARGET_TAG, TARGET_REP, TARGET_CELL, TARGET_ENV))
dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(scored = scored, cluster_detail = cluster_detail), OUT)
invisible(check_ldscnr())
write_receipt(STAGE, inputs = INPUTS, params = PARAMS, outputs = OUT, target = combo_id)
say("\n[5] wrote %s\n    receipt: %s\n", OUT, receipt_path(STAGE, combo_id))
