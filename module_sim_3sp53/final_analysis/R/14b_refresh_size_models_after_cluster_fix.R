## =============================================================================
## final_analysis/R/14b_refresh_size_models_after_cluster_fix.R
##
## One-off refresh, not a new pipeline stage: regenerates results/
## simulation_stage2_size_models.tsv using R/14's FIXED `.fit_size_model()`
## (clustered by `rep` alone, not (cell, rep) -- PK's 2026-09-26 third
## review), without rerunning R/14's own opportunity-effect relocation
## section, which needs the 02_build_ld_units bundles -- NOT reachable on
## this host right now (external volume unmounted).
##
## Reads `dt` back from results/simulation_stage2_region_details.rds (R/11's
## output, already on disk, bundle-independent, untouched by this fix) and
## re-runs R/14's own `.fit_size_model()` verbatim (copied, not
## reimplemented) with the corrected cluster key. Point estimates are
## unchanged -- only ci_lo/ci_hi/n_boot_converged are affected.
##
## Also fixes Supplementary.tex:614 (LDscnR_manuscript), which quotes these
## intervals directly -- see that edit's own commit/note.
## =============================================================================
suppressMessages({library(data.table)})
MODULE_ROOT <- path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis")
source(file.path(MODULE_ROOT, "R", "00_config.R"))
source(file.path(MODULE_ROOT, "R", "helpers_stage2_truth.R"))   ## ci_quantile()
say("=== 14b_refresh_size_models_after_cluster_fix ===\n\n")

REGION_METHODS <- c("emmax_snp_region", "emmax_simes_region", "emmax_consensus_region",
                    "lfmm_snp_region", "lfmm_simes_region")

region_detail_file <- file.path(MODULE_ROOT, "results", "simulation_stage2_region_details.rds")
dt <- readRDS(region_detail_file)
dt[, FP := !TP]
say("[1] loaded cached region details (%d rows), unchanged by this fix\n", nrow(dt))

## ---- copied verbatim from R/14's (now-fixed) .fit_size_model() -------------
.fit_size_model <- function(size_col) {
  d <- copy(dt)[n_markers > 0]
  d[, log2_size := log2(get(size_col))]
  d[, method := factor(method, levels = REGION_METHODS)]
  d[, cell := factor(cell, levels = CELLS_ALL)]
  d[, tag := factor(tag, levels = TAGS_ALL)]
  fit <- stats::glm(FP ~ log2_size + method + cell + tag, data = d, family = stats::binomial())
  co <- summary(fit)$coefficients
  point <- data.table(term = rownames(co), estimate = co[, "Estimate"])

  clusters <- sort(unique(d$rep))
  n_cl <- length(clusters)
  set.seed(SEEDS[["bootstrap"]])
  boot_coef <- matrix(NA_real_, nrow = N_BOOTSTRAP, ncol = nrow(point), dimnames = list(NULL, point$term))
  cl_rows <- split(seq_len(nrow(d)), d$rep)
  for (bb in seq_len(N_BOOTSTRAP)) {
    draw <- sample(clusters, n_cl, replace = TRUE)
    idx <- unlist(cl_rows[as.character(draw)], use.names = FALSE)
    fit_b <- tryCatch(stats::glm(FP ~ log2_size + method + cell + tag, data = d[idx], family = stats::binomial()),
                      error = function(e) NULL, warning = function(w) NULL)
    if (!is.null(fit_b)) {
      cb <- coef(fit_b)
      boot_coef[bb, names(cb)[names(cb) %in% colnames(boot_coef)]] <- cb[names(cb) %in% colnames(boot_coef)]
    }
  }
  point[, `:=`(ci_lo = apply(boot_coef, 2, ci_quantile)[1, ][term],
              ci_hi = apply(boot_coef, 2, ci_quantile)[2, ][term],
              n_boot_converged = colSums(!is.na(boot_coef))[term])]
  list(point = point, n_obs = nrow(d), n_events = sum(d$FP), aic = AIC(fit))
}

model_markers <- .fit_size_model("n_markers")
model_units <- .fit_size_model("n_units")
models <- rbind(model_markers$point[, size_measure := "n_markers"],
               model_units$point[, size_measure := "n_units"])
fwrite(models, file.path(MODULE_ROOT, "results", "simulation_stage2_size_models.tsv"), sep = "\t")
say("[2] overwrote results/simulation_stage2_size_models.tsv (rep-clustered bootstrap)\n")
say("    n_markers model: n=%d, FP events=%d, AIC=%.1f\n", model_markers$n_obs, model_markers$n_events, model_markers$aic)
say("    n_units   model: n=%d, FP events=%d, AIC=%.1f\n", model_units$n_obs, model_units$n_events, model_units$aic)
print(models[term == "log2_size"])
cat("REFRESH_DONE\n")
