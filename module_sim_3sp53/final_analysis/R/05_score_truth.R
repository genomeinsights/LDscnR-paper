## =============================================================================
## final_analysis/R/05_score_truth.R
##
## Truth definition and scoring, per CLAUDE_REANALYSIS_INSTRUCTIONS.md.
##
## Uses LDscnR's own primitives (flag_true_qtns(), qtn_ld_table(),
## score_thresholds()) but NOT classify_ors()/evaluate_ors() -- those
## implement a DIFFERENT, dedup-neutral region-level metric (one region
## claims at most one QTN; a duplicate claim is dropped from both TP and FP)
## that this project has used extensively elsewhere for exploratory C-score
## work, but the instructions explicitly reject that shape here:
## "Do not silently discard duplicate links to an already recovered QTN from
## this denominator." Precision/recall below are computed directly, at the
## hypothesis level, per spec.
##
## Definitions (verbatim from the instructions):
##  - Va_j = 2 p_j (1-p_j) a_j^2, p_j from the analysed individuals.
##  - Detectable QTN: MAF>0.10 AND >=5% of the chromosome's total QTN Va
##    (flag_true_qtns()'s own definition, called with those exact defaults).
##  - A significant marker is truth-linked when it satisfies rho_r2=0.75/
##    rho_d=0.95-derived (r2min, dmax) to >=1 detectable QTN.
##  - A significant Stage-1 unit is truth-linked when >=1 member marker is.
##  - Precision = truth-linked significant hypotheses / all significant
##    hypotheses (no dedup).
##  - Recall = unique detectable QTN recovered by >=1 significant hypothesis
##    / all detectable QTN.
##  - Conditional recall = same numerator / detectable QTN covered by >=1
##    ELIGIBLE Stage-1 unit (a coverage diagnostic, separate from power loss).
## =============================================================================
suppressMessages({library(data.table); library(LDscnR)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))
STAGE <- "05_score_truth"

score_truth <- function(tag, cell, rep, env, force = FALSE) {
  combo_id <- sprintf("%s_%s_rep%d_env%d", tag, cell, rep, env)
  emmax_file <- file.path(stage_dir("03_emmax", combo_id), "emmax.rds")
  ld_units_file <- file.path(stage_dir("02_build_ld_units", combo_id), "ld_units.rds")
  if (!file.exists(emmax_file)) stop("R/03_emmax.R has not produced: ", emmax_file)

  INPUTS <- c(emmax_file, ld_units_file)
  PARAMS <- list(va_share_detectable = VA_SHARE_DETECTABLE, maf_keep = MAF_KEEP,
                 truth_rho_r2 = TRUTH_RHO_R2, truth_rho_d = TRUTH_RHO_D, truth_dmax_cap = TRUTH_DMAX_CAP,
                 size_floor = SIZE_FLOOR)
  if (!force && !stage_stale(STAGE, INPUTS, PARAMS, target = combo_id)) {
    return(readRDS(file.path(stage_dir(STAGE, combo_id), "truth_scores.rds")))
  }

  say("=== %s: %s/%s/rep%d/env%d ===\n\n", STAGE, tag, cell, rep, env)
  em <- readRDS(emmax_file)
  b  <- readRDS(ld_units_file)
  GTs <- b$GTs; map <- copy(b$map); stage1 <- b$stage1

  ## ---- truth definition ---------------------------------------------------------
  say("[1] Va and detectability (MAF>%.2f, Va share>=%.0f%%)\n", MAF_KEEP, 100 * VA_SHARE_DETECTABLE)
  p <- colSums(GTs) / nrow(GTs) / 2
  map[, p_freq := p[match(marker, colnames(GTs))]]
  map[, Va := 2 * p_freq * (1 - p_freq) * allelic_values^2]
  map <- flag_true_qtns(map, va_col = "Va", maf_col = "MAF", maf_min = MAF_KEEP, p_va_min = VA_SHARE_DETECTABLE)
  detectable_qtn <- map[true_pos_QTN == TRUE, marker]
  say("    %d of %d QTN detectable\n", length(detectable_qtn), sum(map$type == "QTN"))

  say("\n[2] LD to QTN + decay-relative match thresholds (rho_r2=%.2f, rho_d=%.2f, dmax_cap=%.0f)\n",
      TRUTH_RHO_R2, TRUTH_RHO_D, TRUTH_DMAX_CAP)
  thr <- score_thresholds(b$LD_decay$decay_sum, rho_r2 = TRUTH_RHO_R2, rho_d = TRUTH_RHO_D, dmax_cap = TRUTH_DMAX_CAP)
  say("    r2min=%.4f, dmax=%.0f bp\n", thr$r2min, thr$dmax)
  qtn_lut <- qtn_ld_table(GTs, map, candidate_markers = map$marker, cores = 1)
  qtn_lut_match <- qtn_lut[r2 > thr$r2min & dist_bp < thr$dmax & qtn_marker %in% detectable_qtn]
  say("    %s marker-QTN pairs pass the match thresholds (of %s within max_bp)\n",
      format(nrow(qtn_lut_match), big.mark = ","), format(nrow(qtn_lut), big.mark = ","))

  ## ---- eligible-unit coverage (for conditional recall) ---------------------------
  units_base <- LDscnR:::.ld_outlier_units(stage1, map, SIZE_FLOOR)
  eligible_markers <- unique(unlist(units_base$members))
  eligible_marker_fraction <- length(eligible_markers) / nrow(map)
  qtn_covered_by_eligible_unit <- unique(qtn_lut_match[marker %in% eligible_markers, qtn_marker])
  say("    eligible-marker fraction: %.4f ; detectable QTN covered by an eligible unit: %d/%d\n",
      eligible_marker_fraction, length(qtn_covered_by_eligible_unit), length(detectable_qtn))

  ## ---- generic hypothesis-level scorer --------------------------------------------
  ## hyp_members: named list, hypothesis id -> member marker(s) (length 1 for a
  ## plain-marker hypothesis). sig_ids: which hypothesis ids are significant.
  .score <- function(method, hyp_members, sig_ids, n_tested, bh_crit_p) {
    n_sig <- length(sig_ids)
    if (n_sig == 0L) {
      return(data.table(method = method, n_tested = n_tested, n_significant = 0L,
                        TP = 0L, FP = 0L, FN = length(detectable_qtn), bh_crit_p = bh_crit_p,
                        precision = NA_real_, recall = 0,
                        eligible_marker_fraction = eligible_marker_fraction,
                        detectable_qtn_covered = length(qtn_covered_by_eligible_unit),
                        n_detectable_qtn = length(detectable_qtn),
                        conditional_recall = if (length(qtn_covered_by_eligible_unit)) 0 else NA_real_))
    }
    linked_qtn_per_hyp <- lapply(hyp_members[sig_ids], function(mm) unique(qtn_lut_match[marker %in% mm, qtn_marker]))
    truth_linked <- lengths(linked_qtn_per_hyp) > 0
    TP <- sum(truth_linked); FP <- sum(!truth_linked)
    recovered_qtn <- unique(unlist(linked_qtn_per_hyp[truth_linked]))
    FN <- length(detectable_qtn) - length(recovered_qtn)
    precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
    recall <- if (length(detectable_qtn) > 0) length(recovered_qtn) / length(detectable_qtn) else NA_real_
    cond_recall <- if (length(qtn_covered_by_eligible_unit) > 0)
      length(intersect(recovered_qtn, qtn_covered_by_eligible_unit)) / length(qtn_covered_by_eligible_unit) else NA_real_
    data.table(method = method, n_tested = n_tested, n_significant = n_sig,
              TP = TP, FP = FP, FN = FN, bh_crit_p = bh_crit_p,
              precision = precision, recall = recall,
              eligible_marker_fraction = eligible_marker_fraction,
              detectable_qtn_covered = length(qtn_covered_by_eligible_unit),
              n_detectable_qtn = length(detectable_qtn),
              conditional_recall = cond_recall)
  }

  say("\n[3] scoring each method\n")
  ## emmax_snp / emmax_snp_nonsingleton: hypothesis = single marker
  marker_members <- stats::setNames(as.list(em$marker$marker), em$marker$marker)
  res_snp <- .score("emmax_snp", marker_members, em$marker$marker[em$marker$significant_snp],
                    em$marker$n_tested, em$marker$bh_crit_p_snp)
  res_snp_ns <- .score("emmax_snp_nonsingleton", marker_members, em$marker$marker[em$marker$significant_snp_nonsingleton],
                       em$marker$n_tested, em$marker$bh_crit_p_snp)

  ## emmax_simes / emmax_consensus: hypothesis = Stage-1 unit (member markers).
  ## ld_outlier_test()'s returned $units table does NOT carry `members` (just
  ## unit_id/Chr/from/to/n_markers/p/q/significant) -- look membership up from
  ## units_base (already built above for the eligible-marker calc), keyed by
  ## the SAME unit_id (units_base is the fixed candidate table both the
  ## observed test and this lookup share).
  units_base_members <- stats::setNames(units_base$members, as.character(units_base$unit_id))
  .unit_members <- function(units_tbl) units_base_members[as.character(units_tbl$unit_id)]
  sim_members <- .unit_members(em$emmax_simes)
  sig_sim <- as.character(em$emmax_simes$unit_id[em$emmax_simes$significant])
  bh_crit_sim <- { qv <- em$emmax_simes$p[em$emmax_simes$significant]; if (length(qv)) max(qv) else NA_real_ }
  res_sim <- .score("emmax_simes", sim_members, sig_sim, nrow(em$emmax_simes), bh_crit_sim)

  con_members <- .unit_members(em$emmax_consensus)
  sig_con <- as.character(em$emmax_consensus$unit_id[em$emmax_consensus$significant])
  bh_crit_con <- { qv <- em$emmax_consensus$p[em$emmax_consensus$significant]; if (length(qv)) max(qv) else NA_real_ }
  res_con <- .score("emmax_consensus", con_members, sig_con, nrow(em$emmax_consensus), bh_crit_con)

  scores <- rbindlist(list(res_snp, res_snp_ns, res_sim, res_con))
  scores[, `:=`(tag = tag, cell = cell, rep = rep, env = env)]
  print(scores[, .(method, n_tested, n_significant, TP, FP, FN, precision, recall, conditional_recall)])

  OUT_DIR <- stage_dir(STAGE, combo_id)
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  OUT <- file.path(OUT_DIR, "truth_scores.rds")
  saveRDS(list(scores = scores, thresholds = thr, n_detectable_qtn = length(detectable_qtn),
              settings = list(tag = tag, cell = cell, rep = rep, env = env, params = PARAMS)), OUT)
  write_receipt(STAGE, inputs = INPUTS, params = PARAMS, outputs = OUT, target = combo_id)
  say("\n[4] wrote %s\n", OUT)
  scores
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 4) stop("Usage: Rscript R/05_score_truth.R <tag> <cell> <rep> <env>")
  score_truth(args[1], args[2], as.integer(args[3]), as.integer(args[4]))
}
