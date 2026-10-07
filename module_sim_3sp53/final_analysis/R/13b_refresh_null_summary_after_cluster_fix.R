## =============================================================================
## final_analysis/R/13b_refresh_null_summary_after_cluster_fix.R
##
## One-off refresh, not a new pipeline stage: regenerates results/
## simulation_null_truth_summary.tsv using R/13's FIXED map-cluster bootstrap
## (clustered by `rep` alone, not (cell, rep) -- PK's 2026-09-26 third
## review), without rerunning R/13's own expensive null-generation stages
## 1-5 (surrogate phenotype draws, Stage-2 assembly per surrogate).
##
## Reads `map_level` back from results/simulation_null_truth_calibration.tsv
## (R/13's own stage-5 output, already on disk, untouched by this fix -- the
## fix only changes how the summary's bootstrap CI resamples rows already
## computed, not the rows themselves) and re-runs R/13's own summary section
## [6]-[7] verbatim (copied, not reimplemented) with the corrected cluster
## key. Point estimates (spearman_rho, pooled_R_null_obs, pooled_FDP_truth,
## slope, intercept, etc.) are IDENTICAL to before -- only the four CI
## columns (pooled_R_null_obs_ci_lo/hi, pooled_FDP_truth_ci_lo/hi) change.
## =============================================================================
suppressMessages({library(data.table)})
MODULE_ROOT <- path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis")
source(file.path(MODULE_ROOT, "R", "00_config.R"))
source(file.path(MODULE_ROOT, "R", "helpers_stage2_truth.R"))   ## ci_quantile(), bootstrap_rep_matrix()
say("=== 13b_refresh_null_summary_after_cluster_fix ===\n\n")

NULL_METHODS <- c("emmax_consensus", "emmax_simes")
NULL_SCHEMES <- c("group", "mvn", "spatial")

map_level <- fread(file.path(MODULE_ROOT, "results", "simulation_null_truth_calibration.tsv"))
say("[1] loaded cached map_level (%d rows), unchanged by this fix\n", nrow(map_level))

## ---- [6]-[7] copied verbatim from R/13, only cluster_level's `by=` differs --
summary_rows <- list()
for (m in NULL_METHODS) for (s in NULL_SCHEMES) {
  ml <- map_level[method == m & scheme == s]
  complete <- ml[!is.na(R_null_obs) & !is.na(FDP_truth)]
  rho <- if (nrow(complete) >= 3) stats::cor(complete$R_null_obs, complete$FDP_truth, method = "spearman") else NA_real_
  fit <- if (nrow(complete) >= 3) stats::lm(FDP_truth ~ R_null_obs, data = complete) else NULL
  slope <- if (!is.null(fit)) unname(coef(fit)[2]) else NA_real_
  intercept <- if (!is.null(fit)) unname(coef(fit)[1]) else NA_real_
  mean_signed_diff <- if (nrow(complete)) mean(complete$R_null_obs - complete$FDP_truth) else NA_real_
  median_abs_diff <- if (nrow(complete)) stats::median(abs(complete$R_null_obs - complete$FDP_truth)) else NA_real_

  adj_rho <- NA_real_
  if (nrow(complete) >= 3 && uniqueN(complete$cell) >= 2) {
    rx <- rank(complete$R_null_obs); ry <- rank(complete$FDP_truth)
    rx_resid <- stats::resid(stats::lm(rx ~ factor(cell) + factor(tag), data = complete))
    ry_resid <- stats::resid(stats::lm(ry ~ factor(cell) + factor(tag), data = complete))
    adj_rho <- suppressWarnings(stats::cor(rx_resid, ry_resid))
  }

  ## FIXED: clustered by `rep` alone (was (cell, rep) -- see R/13's own
  ## updated comment at this same point for the full justification).
  cluster_level <- ml[, .(sum_E_null = sum(sum_E_null), sum_n_obs = sum(sum_n_obs),
                          sum_TP = sum(sum_TP), sum_FP = sum(sum_FP)), by = .(rep)]
  n_r <- nrow(cluster_level)
  mat <- as.matrix(cluster_level[, .(sum_E_null, sum_n_obs, sum_TP, sum_FP)])
  bs <- if (n_r >= 2) bootstrap_rep_matrix(mat, N_BOOTSTRAP, SEEDS[["bootstrap"]]) else NULL
  ratio_b <- if (!is.null(bs)) bs[, "sum_E_null"] / bs[, "sum_n_obs"] else NA_real_
  fdp_b <- if (!is.null(bs)) bs[, "sum_FP"] / (bs[, "sum_TP"] + bs[, "sum_FP"]) else NA_real_
  ratio_b <- ratio_b[is.finite(ratio_b)]; fdp_b <- fdp_b[is.finite(fdp_b)]

  null_role <- c(group = "candidate_null", mvn = "negative_control", spatial = "candidate_null")[[s]]

  summary_rows[[length(summary_rows) + 1]] <- data.table(
    method = m, scheme = s, null_role = null_role,
    n_map_groups = nrow(ml), n_complete = nrow(complete),
    frac_zero_obs_regions = mean(ml$sum_n_obs == 0),
    spearman_rho = rho, cell_tag_adjusted_rho = adj_rho, slope = slope, intercept = intercept,
    mean_signed_diff = mean_signed_diff, median_abs_diff = median_abs_diff,
    pooled_E_null = sum(ml$sum_E_null), pooled_n_obs = sum(ml$sum_n_obs),
    pooled_TP = sum(ml$sum_TP), pooled_FP = sum(ml$sum_FP),
    pooled_R_null_obs = if (sum(ml$sum_n_obs) > 0) sum(ml$sum_E_null) / sum(ml$sum_n_obs) else NA_real_,
    pooled_R_null_obs_ci_lo = if (length(ratio_b)) ci_quantile(ratio_b)[1] else NA_real_,
    pooled_R_null_obs_ci_hi = if (length(ratio_b)) ci_quantile(ratio_b)[2] else NA_real_,
    pooled_FDP_truth = if (sum(ml$sum_TP) + sum(ml$sum_FP) > 0) sum(ml$sum_FP) / (sum(ml$sum_TP) + sum(ml$sum_FP)) else NA_real_,
    pooled_FDP_truth_ci_lo = if (length(fdp_b)) ci_quantile(fdp_b)[1] else NA_real_,
    pooled_FDP_truth_ci_hi = if (length(fdp_b)) ci_quantile(fdp_b)[2] else NA_real_)
}
summary_dt <- rbindlist(summary_rows)
fwrite(summary_dt, file.path(MODULE_ROOT, "results", "simulation_null_truth_summary.tsv"), sep = "\t")
say("[2] overwrote results/simulation_null_truth_summary.tsv (rep-clustered bootstrap)\n")
print(summary_dt[, .(method, scheme, null_role, n_map_groups, frac_zero_obs_regions,
                     spearman_rho, cell_tag_adjusted_rho, slope, pooled_R_null_obs,
                     pooled_R_null_obs_ci_lo, pooled_R_null_obs_ci_hi, pooled_FDP_truth,
                     pooled_FDP_truth_ci_lo, pooled_FDP_truth_ci_hi)])
cat("REFRESH_DONE\n")
