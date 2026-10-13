## =============================================================================
## final_analysis/R/30_truth_sensitivity_by_group.R
##
## Truth-threshold sensitivity pooled separately for the primary high-gene-flow
## scenarios (c = 1, 600 runs) and the lower-gene-flow scenarios (c = 1.5 and
## c = 2, 800 runs), from the per-run rows written by R/09_truth_sensitivity.R
## (R09_OUT_TAG = "c1" and "lowgf"). Same estimator and map-cluster bootstrap as
## R/09: counts pooled within a group, replicate maps resampled.
## Output: results/simulation_truth_sensitivity_by_group.tsv
## =============================================================================
suppressMessages(library(data.table))
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))
say("=== 30_truth_sensitivity_by_group ===\n\n")

files <- c(high_gene_flow = "results/simulation_truth_sensitivity_per_run_c1.tsv",
           lower_gene_flow = "results/simulation_truth_sensitivity_per_run_lowgf.tsv")
files <- files[file.exists(files)]
stopifnot(length(files) > 0)
ci <- function(x, probs = c(0.025, 0.975)) stats::quantile(x, probs, na.rm = TRUE, names = FALSE)
n_rep <- length(REPS_ALL); B <- N_BOOTSTRAP
set.seed(SEEDS[["bootstrap"]])
draws <- matrix(sample.int(n_rep, size = n_rep * B, replace = TRUE), nrow = n_rep, ncol = B)
mult <- apply(draws, 2, tabulate, nbins = n_rep)

out <- rbindlist(lapply(names(files), function(g) {
  dt <- fread(files[[g]])
  say("[%s] %d runs\n", g, uniqueN(dt[, .(tag, cell, rep, env)]))
  rep_dt <- dt[, .(TP = sum(TP), FP = sum(FP), n_detectable_qtn = sum(n_detectable_qtn), n_recovered = sum(n_recovered)),
               by = .(method, va_share_threshold, rep)]
  rbindlist(lapply(split(rep_dt, by = c("method", "va_share_threshold")), function(rd) {
    rd <- rd[match(REPS_ALL, rep)]
    rd[is.na(TP), `:=`(TP = 0L, FP = 0L, n_detectable_qtn = 0L, n_recovered = 0L)]
    TP_b <- as.numeric(crossprod(mult, rd$TP)); FP_b <- as.numeric(crossprod(mult, rd$FP))
    ndq_b <- as.numeric(crossprod(mult, rd$n_detectable_qtn)); nrec_b <- as.numeric(crossprod(mult, rd$n_recovered))
    prec_b <- TP_b / pmax(TP_b + FP_b, 1); rec_b <- nrec_b / pmax(ndq_b, 1)
    TP <- sum(rd$TP); FP <- sum(rd$FP); ndq <- sum(rd$n_detectable_qtn); nrec <- sum(rd$n_recovered)
    data.table(group = g, method = rd$method[!is.na(rd$method)][1], va_share_threshold = rd$va_share_threshold[!is.na(rd$va_share_threshold)][1],
               TP = TP, FP = FP, n_detectable_qtn = ndq, n_recovered = nrec,
               precision = TP / max(TP + FP, 1), precision_ci_lo = ci(prec_b)[1], precision_ci_hi = ci(prec_b)[2],
               recall = nrec / max(ndq, 1), recall_ci_lo = ci(rec_b)[1], recall_ci_hi = ci(rec_b)[2])
  }))
}))
setorder(out, group, method, va_share_threshold)
print(out)
fwrite(out, "results/simulation_truth_sensitivity_by_group.tsv", sep = "\t")
say("wrote results/simulation_truth_sensitivity_by_group.tsv\n")
