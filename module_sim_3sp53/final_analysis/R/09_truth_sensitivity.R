## =============================================================================
## final_analysis/R/09_truth_sensitivity.R
##
## Sensitivity analysis: rescore truth (never association statistics) across
## the prespecified detectable-QTN additive-variance-share threshold grid
## (00_config.R's VA_SHARE_SENSITIVITY_GRID, includes the primary 5% value),
## per CLAUDE_REANALYSIS_INSTRUCTIONS.md's "Sensitivity analyses worth
## retaining": "This is cheap because association statistics do not change."
##
## Cheap by construction here too: Stage-2 region assembly (the expensive
## part of 05_score_truth.R -- ld_prune_and_eMLG() via .run_stage2()) does
## NOT depend on the truth threshold at all, only on which hypotheses are
## significant -- so region_members for each of the 5 PRIMARY region
## methods (emmax_snp_region/emmax_simes_region/emmax_consensus_region/
## lfmm_snp_region/lfmm_simes_region) is computed ONCE per combo and reused
## across every threshold in the grid. Likewise qtn_ld_table()'s full
## marker-QTN (r2, dist_bp) table is threshold-independent and computed
## once; only which QTN count as "detectable" (flag_true_qtns(p_va_min=...))
## and hence which rows of that table survive varies per threshold.
##
## Reuses 05_score_truth.R's region-assembly/scoring logic verbatim (kept
## as a duplicated, self-contained copy here rather than a shared function,
## matching this project's established convention -- see
## manhattan_example_data.R's own header comment on the same choice).
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(parallel)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))
say("=== truth_sensitivity ===\n\n")

.score_one <- function(hyp_members, sig_ids, qtn_lut_match, detectable_qtn) {
  n_sig <- length(sig_ids)
  if (n_sig == 0L) return(list(TP = 0L, FP = 0L, n_recovered = 0L))
  linked_qtn_per_hyp <- lapply(hyp_members[sig_ids], function(mm) unique(qtn_lut_match[marker %in% mm, qtn_marker]))
  truth_linked <- lengths(linked_qtn_per_hyp) > 0
  recovered_qtn <- unique(unlist(linked_qtn_per_hyp[truth_linked]))
  list(TP = sum(truth_linked), FP = sum(!truth_linked), n_recovered = length(recovered_qtn))
}

score_truth_sensitivity <- function(tag, cell, rep, env) {
  combo_id <- sprintf("%s_%s_rep%d_env%d", tag, cell, rep, env)
  em <- readRDS(file.path(stage_dir("03_emmax", combo_id), "emmax.rds"))
  lfmm_file <- file.path(stage_dir("04_lfmm", combo_id), "lfmm.rds")
  have_lfmm <- file.exists(lfmm_file)
  lf <- if (have_lfmm) readRDS(lfmm_file) else NULL
  b <- readRDS(file.path(stage_dir("02_build_ld_units", combo_id), "ld_units.rds"))
  GTs <- b$GTs; map0 <- copy(b$map); stage1 <- b$stage1

  ## ---- threshold-INDEPENDENT: Va per marker, full qtn_lut, region_members ----
  p <- colSums(GTs) / nrow(GTs) / 2
  map0[, p_freq := p[match(marker, colnames(GTs))]]
  map0[, Va := 2 * p_freq * (1 - p_freq) * allelic_values^2]

  thr <- score_thresholds(b$LD_decay$decay_sum, rho_r2 = TRUTH_RHO_R2, rho_d = TRUTH_RHO_D, dmax_cap = TRUTH_DMAX_CAP)
  qtn_lut_full <- if (sum(map0$type == "QTN") == 0L) {
    data.table(marker = character(), qtn_marker = character(), r2 = numeric(), dist_bp = numeric())
  } else {
    qtn_ld_table(GTs, map0, candidate_markers = map0$marker, cores = 1)
  }
  qtn_lut_pass <- qtn_lut_full[r2 > thr$r2min & dist_bp < thr$dmax]   ## LD/distance pass, truth-independent

  .run_stage2 <- function(cl_sub) {
    if (!nrow(cl_sub)) return(list())
    mk_sub <- unlist(cl_sub$members, use.names = FALSE)
    ms_sub <- as.data.table(stage1$map_snp)[marker %chin% mk_sub]
    sub <- structure(list(map_snp = ms_sub, clusters = cl_sub, pruned = cl_sub$core_snp), class = "ld_complexity_reduction")
    pr <- ld_prune_and_eMLG(GTs = GTs[, mk_sub, drop = FALSE], stage1 = sub, ld_w_col = "ld_w_095", ld_w_threshold = 0,
                            LD_decay = b$LD_decay, min_r2_rho = stage1$params$rho,
                            score_threshold = REGION_ASSEMBLY$score_threshold, distance_threshold = REGION_ASSEMBLY$distance_threshold,
                            compute_unflagged_eMLG = FALSE, min_n_loci_eMLG = 1, min_n_loci_flag = 1, cores = 1)
    g <- as.data.table(pr$groups)
    stats::setNames(g$members, as.character(seq_len(nrow(g))))
  }
  cl_all <- as.data.table(stage1$clusters)
  nl_all <- if ("n_loci" %in% names(cl_all)) cl_all$n_loci else cl_all$n_snps
  .from_units <- function(units_tbl) {
    sig <- units_tbl[units_tbl$significant == TRUE]
    if (!nrow(sig)) return(list())
    .run_stage2(cl_all[nl_all >= SIZE_FLOOR][as.integer(sig$unit_id)])
  }
  .from_markers <- function(sig_markers) {
    if (!length(sig_markers)) return(list())
    disc <- vapply(cl_all$members, function(mm) any(mm %chin% sig_markers), logical(1))
    .run_stage2(cl_all[disc])
  }

  region_members <- list(
    emmax_snp_region = .from_markers(em$marker$marker[em$marker$significant_snp]),
    emmax_simes_region = .from_units(em$emmax_simes),
    emmax_consensus_region = .from_units(em$emmax_consensus)
  )
  if (have_lfmm) {
    region_members$lfmm_snp_region <- .from_markers(lf$marker$marker[lf$marker$significant_snp])
    region_members$lfmm_simes_region <- .from_units(lf$lfmm_simes)
  }

  ## ---- threshold-DEPENDENT: detectable_qtn, qtn_lut_match, scoring -----------
  rows <- list()
  for (va_thr in VA_SHARE_SENSITIVITY_GRID) {
    map_t <- flag_true_qtns(map0, va_col = "Va", maf_col = "MAF", maf_min = MAF_KEEP, p_va_min = va_thr)
    detectable_qtn <- map_t[true_pos_QTN == TRUE, marker]
    qtn_lut_match <- qtn_lut_pass[qtn_marker %in% detectable_qtn]
    for (m in names(region_members)) {
      rm_ <- region_members[[m]]
      sc <- .score_one(rm_, names(rm_), qtn_lut_match, detectable_qtn)
      precision <- if ((sc$TP + sc$FP) > 0) sc$TP / (sc$TP + sc$FP) else NA_real_
      recall <- if (length(detectable_qtn) > 0) sc$n_recovered / length(detectable_qtn) else NA_real_
      rows[[length(rows) + 1]] <- data.table(
        tag = tag, cell = cell, rep = rep, env = env, method = m, va_share_threshold = va_thr,
        n_significant = length(rm_), TP = sc$TP, FP = sc$FP,
        n_detectable_qtn = length(detectable_qtn), n_recovered = sc$n_recovered,
        precision = precision, recall = recall)
    }
  }
  rbindlist(rows)
}

if (sys.nframe() == 0L) {
  combos <- CJ(tag = TAGS_ALL, cell = CELLS_ALL, rep = REPS_ALL, env = ENVS_ALL)
  say("[1] rescoring %d combos x %d Va-share thresholds x 5 region methods\n", nrow(combos), length(VA_SHARE_SENSITIVITY_GRID))
  rows <- mclapply(seq_len(nrow(combos)), function(i) {
    tryCatch(score_truth_sensitivity(combos$tag[i], combos$cell[i], combos$rep[i], combos$env[i]),
             error = function(e) { message(sprintf("[truth_sensitivity] %s_%s_rep%d_env%d: %s",
                                                    combos$tag[i], combos$cell[i], combos$rep[i], combos$env[i], conditionMessage(e))); NULL })
  }, mc.cores = 7)
  dt <- rbindlist(rows)
  say("    %d rows\n", nrow(dt))

  say("\n[2] grand-pooled (ratio of pooled counts) by method x Va-share threshold\n")
  ci <- function(x, probs = c(0.025, 0.975)) stats::quantile(x, probs, na.rm = TRUE, names = FALSE)
  n_rep <- length(REPS_ALL); B <- N_BOOTSTRAP
  set.seed(SEEDS[["bootstrap"]])
  draws <- matrix(sample.int(n_rep, size = n_rep * B, replace = TRUE), nrow = n_rep, ncol = B)
  mult <- apply(draws, 2, tabulate, nbins = n_rep)

  rep_dt <- dt[, .(TP = sum(TP), FP = sum(FP), n_detectable_qtn = sum(n_detectable_qtn), n_recovered = sum(n_recovered)),
              by = .(method, va_share_threshold, rep)]
  point <- dt[, .(TP = sum(TP), FP = sum(FP), n_detectable_qtn = sum(n_detectable_qtn), n_recovered = sum(n_recovered)),
             by = .(method, va_share_threshold)]
  point[, `:=`(precision = TP / pmax(TP + FP, 1), recall = n_recovered / pmax(n_detectable_qtn, 1))]

  pooled_rows <- list()
  for (m in unique(dt$method)) for (va_thr in VA_SHARE_SENSITIVITY_GRID) {
    rd <- rep_dt[method == m & va_share_threshold == va_thr][match(REPS_ALL, rep)]
    TP_b <- as.numeric(crossprod(mult, rd$TP)); FP_b <- as.numeric(crossprod(mult, rd$FP))
    ndq_b <- as.numeric(crossprod(mult, rd$n_detectable_qtn)); nrec_b <- as.numeric(crossprod(mult, rd$n_recovered))
    prec_b <- TP_b / pmax(TP_b + FP_b, 1); rec_b <- nrec_b / pmax(ndq_b, 1)
    pt <- point[method == m & va_share_threshold == va_thr]
    pooled_rows[[length(pooled_rows) + 1]] <- data.table(
      method = m, va_share_threshold = va_thr, TP = pt$TP, FP = pt$FP,
      n_detectable_qtn = pt$n_detectable_qtn, n_recovered = pt$n_recovered,
      precision = pt$precision, precision_ci_lo = ci(prec_b)[1], precision_ci_hi = ci(prec_b)[2],
      recall = pt$recall, recall_ci_lo = ci(rec_b)[1], recall_ci_hi = ci(rec_b)[2])
  }
  pooled <- rbindlist(pooled_rows)
  setorder(pooled, method, va_share_threshold)
  print(pooled)

  dir.create("results", showWarnings = FALSE)
  fwrite(pooled, "results/simulation_truth_sensitivity.tsv", sep = "\t")
  say("\nwrote results/simulation_truth_sensitivity.tsv\n")
}
