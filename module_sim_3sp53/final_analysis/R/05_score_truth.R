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
##
## PRIMARY estimand (2026-09-10 update): Stage-2 ASSEMBLED REGIONS.
##  - A reported region is TP when >=1 marker belonging to one of its
##    constituent DISCOVERED Stage-1 clusters (not every marker physically
##    between the region bounds) satisfies rho_r2=0.75/rho_d=0.95-derived
##    (r2min, dmax) to >=1 detectable QTN. Otherwise FP.
##  - When several discovered clusters merge into one region, the region
##    counts ONCE (one TP candidate region, not one TP plus FP calls).
##  - Region precision = TP regions / (TP regions + FP regions); region
##    recall = unique detectable QTN recovered / all detectable QTN.
##  - `emmax_snp_region`/`lfmm_snp_region` make the unrestricted
##    marker-wise comparator's reported-call unit comparable to Simes/
##    consensus's: every phenotype-blind Stage-1 cluster (including
##    singletons) containing >=1 BH-significant marker is "discovered"
##    and fed through the SAME Stage-2 assembly -- this changes only
##    reporting, never marker p-values or marker-wise BH.
##
## DIAGNOSTIC (marker-/Stage-1-unit-level, retained to show why region
## assembly is necessary, NOT the primary estimand):
##  - A significant marker is truth-linked when it satisfies the same
##    (r2min, dmax) criteria to >=1 detectable QTN.
##  - A significant Stage-1 unit is truth-linked when >=1 member marker is.
##  - Diagnostic precision = truth-linked significant hypotheses / all
##    significant hypotheses (no dedup).
##  - Conditional recall = recovered-QTN numerator / detectable QTN covered
##    by >=1 ELIGIBLE Stage-1 unit (a coverage diagnostic, separate from
##    test-power loss).
## =============================================================================
suppressMessages({library(data.table); library(LDscnR)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))

score_truth <- function(tag, cell, rep, env, force = FALSE) {
  STAGE <- "05_score_truth"   ## LOCAL -- see 02_build_ld_units.R's comment on this exact bug
  combo_id <- sprintf("%s_%s_rep%d_env%d", tag, cell, rep, env)
  emmax_file <- file.path(stage_dir("03_emmax", combo_id), "emmax.rds")
  lfmm_file <- file.path(stage_dir("04_lfmm", combo_id), "lfmm.rds")
  ld_units_file <- file.path(stage_dir("02_build_ld_units", combo_id), "ld_units.rds")
  if (!file.exists(emmax_file)) stop("R/03_emmax.R has not produced: ", emmax_file)
  ## LFMM is OPTIONAL -- R/04_lfmm.R is a separate, later stage (instructions:
  ## "so that the primary EMMAX analysis can finish and be audited without
  ## rerunning LFMM"). Score whatever is present; a combo scored before LFMM
  ## was run for it is simply re-scored (force=TRUE, or a fresh receipt once
  ## lfmm.rds first appears -- lfmm_file's presence is part of PARAMS below,
  ## so its mtime/existence changing invalidates the cached truth_scores.rds).
  have_lfmm <- file.exists(lfmm_file)

  INPUTS <- c(emmax_file, ld_units_file, if (have_lfmm) lfmm_file)
  PARAMS <- list(va_share_detectable = VA_SHARE_DETECTABLE, maf_keep = MAF_KEEP,
                 truth_rho_r2 = TRUTH_RHO_R2, truth_rho_d = TRUTH_RHO_D, truth_dmax_cap = TRUTH_DMAX_CAP,
                 size_floor = SIZE_FLOOR, have_lfmm = have_lfmm,
                 region_assembly = REGION_ASSEMBLY, region_scoring_version = 2L)
  if (!force && !stage_stale(STAGE, INPUTS, PARAMS, target = combo_id)) {
    return(readRDS(file.path(stage_dir(STAGE, combo_id), "truth_scores.rds")))
  }

  say("=== %s: %s/%s/rep%d/env%d ===\n\n", STAGE, tag, cell, rep, env)
  em <- readRDS(emmax_file)
  lf <- if (have_lfmm) readRDS(lfmm_file) else NULL
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
  ## [!] GUARDED 2026-09-09 -- found in gate 3 (bgs/V1_c1/rep1/env1): a raw
  ## run can save ZERO QTN loci at all (confirmed separately -- each run
  ## saves only a small, variable subset of the reference map's 101 QTN, see
  ## qc/offset_propagation_check.R's header). When there are no QTN on the
  ## map, qtn_ld_table() has nothing to compute distances/r2 to and returns
  ## a genuinely EMPTY data.table with NO COLUMNS (not just zero rows), so
  ## qtn_lut[r2 > ...] errors ("Object 'r2' not found amongst []") instead
  ## of just returning zero matches. Skip the call entirely in that case.
  if (sum(map$type == "QTN") == 0L) {
    qtn_lut <- data.table(marker = character(), qtn_marker = character(), r2 = numeric(), dist_bp = numeric())
  } else {
    qtn_lut <- qtn_ld_table(GTs, map, candidate_markers = map$marker, cores = 1)
  }
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

  score_list <- list(res_snp, res_snp_ns, res_sim, res_con)

  ## ---- Stage-2 assembled-region scoring -- the PRIMARY estimand (PK 2026-09-10) ---
  ## Hypothesis = one WHOLE Stage-2 region; member markers = the ACTUAL
  ## constituent discovered-cluster members (via ld_prune_and_eMLG()'s own
  ## $groups$members), NOT every marker physically between the region's
  ## [Chr,from,to] bounds. The bounds-sweep version used until now was a
  ## real bug: "an intervening, untested marker could otherwise lend truth
  ## credit to the region" (CLAUDE_REANALYSIS_INSTRUCTIONS.md, "Update
  ## after inspection of the Stage-2 Manhattan figures"). .run_stage2()
  ## below is ld_outlier_test()'s own "stage2_discovered" branch, called
  ## directly (same cl_sig/mk_sig/ms_sig/sub construction, same
  ## ld_prune_and_eMLG() call) so we get $groups$members verbatim instead
  ## of reconstructing membership by physical-position approximation.
  ##
  ## Two seeding routes, per the instructions' "Testing unit versus
  ## reported-call unit" section:
  ##  - .assemble_regions_from_units(): Simes/consensus -- seed Stage 2
  ##    with the BH-significant STAGE-1 UNITS (as ld_outlier_test() does).
  ##  - .assemble_regions_from_markers(): unrestricted marker-wise
  ##    (emmax_snp/lfmm_snp) -- mark EVERY phenotype-blind Stage-1 cluster,
  ##    INCLUDING SINGLETONS, as discovered when it contains >=1
  ##    BH-significant marker, then run the SAME assembly. This changes
  ##    only reporting/scoring, never marker p-values or marker-wise BH --
  ##    it makes the unrestricted comparator's reported-call unit
  ##    genuinely comparable to Simes/consensus's, instead of comparing
  ##    "every significant SNP" against "one assembled region."
  ##
  ## Every region is already a significant call by construction, so ALL
  ## region ids are "significant" hypotheses for the generic .score()
  ## scorer above.
  .run_stage2 <- function(cl_sub) {
    if (!nrow(cl_sub)) return(list())
    mk_sub <- unlist(cl_sub$members, use.names = FALSE)
    ms_sub <- as.data.table(stage1$map_snp)[marker %chin% mk_sub]
    sub <- structure(list(map_snp = ms_sub, clusters = cl_sub, pruned = cl_sub$core_snp),
                     class = "ld_complexity_reduction")
    pr <- ld_prune_and_eMLG(GTs = GTs[, mk_sub, drop = FALSE], stage1 = sub,
                            ld_w_col = "ld_w_095", ld_w_threshold = 0,
                            LD_decay = b$LD_decay, min_r2_rho = stage1$params$rho,
                            score_threshold = REGION_ASSEMBLY$score_threshold,
                            distance_threshold = REGION_ASSEMBLY$distance_threshold,
                            compute_unflagged_eMLG = FALSE, min_n_loci_eMLG = 1,
                            min_n_loci_flag = 1, cores = 1)
    g <- as.data.table(pr$groups)
    stats::setNames(g$members, as.character(seq_len(nrow(g))))
  }
  .assemble_regions_from_units <- function(sig_units_tbl) {
    sig_units <- sig_units_tbl[sig_units_tbl$significant == TRUE]
    if (!nrow(sig_units)) return(list())
    cl <- as.data.table(stage1$clusters)
    nl <- if ("n_loci" %in% names(cl)) cl$n_loci else cl$n_snps
    .run_stage2(cl[nl >= SIZE_FLOOR][sig_units$unit_id])
  }
  .assemble_regions_from_markers <- function(sig_markers) {
    if (!length(sig_markers)) return(list())
    cl <- as.data.table(stage1$clusters)
    disc <- vapply(cl$members, function(mm) any(mm %chin% sig_markers), logical(1))
    .run_stage2(cl[disc])
  }

  snp_region_members <- .assemble_regions_from_markers(em$marker$marker[em$marker$significant_snp])
  res_snp_region <- .score("emmax_snp_region", snp_region_members, names(snp_region_members),
                           em$marker$n_tested, NA_real_)

  sim_region_members <- .assemble_regions_from_units(em$emmax_simes)
  res_sim_region <- .score("emmax_simes_region", sim_region_members, names(sim_region_members),
                           nrow(em$emmax_simes), NA_real_)

  con_region_members <- .assemble_regions_from_units(em$emmax_consensus)
  res_con_region <- .score("emmax_consensus_region", con_region_members, names(con_region_members),
                           nrow(em$emmax_consensus), NA_real_)

  score_list <- c(score_list, list(res_snp_region, res_sim_region, res_con_region))

  ## lfmm_snp / lfmm_simes / lfmm_snp_region / lfmm_simes_region -- only
  ## when R/04_lfmm.R has been run for this combo. NO lfmm_consensus
  ## (instructions: "no LFMM consensus-dosage analysis") and no
  ## lfmm_snp_nonsingleton (not in the instructions' LFMM scope, unlike
  ## emmax_snp_nonsingleton).
  if (have_lfmm) {
    lfmm_marker_members <- stats::setNames(as.list(lf$marker$marker), lf$marker$marker)
    res_lfmm_snp <- .score("lfmm_snp", lfmm_marker_members, lf$marker$marker[lf$marker$significant_snp],
                           lf$marker$n_tested, lf$marker$bh_crit_p_snp)
    lfmm_sim_members <- .unit_members(lf$lfmm_simes)
    sig_lfmm_sim <- as.character(lf$lfmm_simes$unit_id[lf$lfmm_simes$significant])
    bh_crit_lfmm_sim <- { qv <- lf$lfmm_simes$p[lf$lfmm_simes$significant]; if (length(qv)) max(qv) else NA_real_ }
    res_lfmm_sim <- .score("lfmm_simes", lfmm_sim_members, sig_lfmm_sim, nrow(lf$lfmm_simes), bh_crit_lfmm_sim)

    lfmm_snp_region_members <- .assemble_regions_from_markers(lf$marker$marker[lf$marker$significant_snp])
    res_lfmm_snp_region <- .score("lfmm_snp_region", lfmm_snp_region_members, names(lfmm_snp_region_members),
                                  lf$marker$n_tested, NA_real_)

    lfmm_sim_region_members <- .assemble_regions_from_units(lf$lfmm_simes)
    res_lfmm_sim_region <- .score("lfmm_simes_region", lfmm_sim_region_members, names(lfmm_sim_region_members),
                                  nrow(lf$lfmm_simes), NA_real_)

    score_list <- c(score_list, list(res_lfmm_snp, res_lfmm_sim, res_lfmm_snp_region, res_lfmm_sim_region))
  }

  scores <- rbindlist(score_list)
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
