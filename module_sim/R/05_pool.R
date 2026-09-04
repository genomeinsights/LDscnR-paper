## module_sim/R/05_pool.R
##
## Replicate-average R/04_score.R's per-replicate TP/FP/FN/Precision/Recall
## into one Precision/Recall (mean +- SE) per (tag, cell, arm), pooling across
## all REPS present for that (tag, cell) -- the house convention used
## everywhere else PR is scored in this repo (module_sim_LDscnR/run_sim_LDscnR.R:
## "ALWAYS replicate-average (mean +- SE) -- env1 alone repeatedly flukes"; same
## pattern in pr_curves.R, score_c2_against_truth.R, bgs_vs_nobgs_prauc.R).
## Macro-averaging (mean of each replicate's own Precision), NOT the legacy
## pooled-count micro-averaging in legacy/R_LDscnR/PR_AUC.R -- that was
## superseded specifically because this repo settled on replicate-averaging
## instead. The tradeoff macro-averaging makes: a replicate with zero
## significant discoveries has Precision = NA (0/0), dropped via na.rm rather
## than treated as 0 or excluded outright -- n_sig_reps below reports how many
## of the REPS_N replicates actually contributed a Precision value, so a mean
## built from few of them is visible rather than silently presented as solid.
##
## No target= subdirectory: unlike 01-04, this stage pools ACROSS reps, so its
## receipt/output are keyed by (tag, cell) only, not (tag, cell, rep).
suppressMessages({library(data.table)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim"), "R", "00_config.R"))
STAGE <- "05_pool"
say("=== %s ===\n\n", STAGE)

CELLS_ALL <- c("V0.5_c1", "V0.5_c1.5", "V0.5_c2", "V1_c1", "V1_c1.5", "V2_c1", "V2_c1.5")
TAGS_ALL  <- c("nobgs", "bgs")
REPS_N    <- 10L

score_files <- data.table(expand.grid(tag = TAGS_ALL, cell = CELLS_ALL, rep = seq_len(REPS_N),
                                       stringsAsFactors = FALSE))
score_files[, path := file.path(PATHS$out, "04_score",
  sprintf("score_%s_rep%d_%s_env1.rds", tag, rep, cell))]
score_files[, exists := file.exists(path)]
say("[0] %d/%d expected (tag,cell,rep) score files present\n", sum(score_files$exists), nrow(score_files))
if (any(!score_files$exists)) {
  missing <- score_files[exists == FALSE]
  say("    missing (run R/04_score.R for these first):\n")
  for (i in seq_len(nrow(missing))) say("      %s %s rep%d\n", missing$tag[i], missing$cell[i], missing$rep[i])
}

INPUTS <- score_files[exists == TRUE, path]
PARAMS <- list(cells = CELLS_ALL, tags = TAGS_ALL, reps_n = REPS_N)
if (!stage_stale(STAGE, INPUTS, PARAMS) && !nzchar(Sys.getenv("FORCE"))) {
  say("\nNothing to do. Set FORCE=1 to rerun anyway.\n"); quit(save = "no")
}

say("\n[1] reading %d score files\n", length(INPUTS))
all_scored <- rbindlist(lapply(INPUTS, readRDS))

say("[2] replicate-averaging (mean +- SE) by (tag, cell, arm)\n")
pooled <- all_scored[, .(
  n_reps        = .N,
  n_sig_reps    = sum(!is.na(Precision)),
  mean_n_sig    = mean(n_sig),
  mean_TP       = mean(TP), mean_FP = mean(FP), mean_FN = mean(FN),
  Precision     = mean(Precision, na.rm = TRUE),
  Precision_SE  = if (sum(!is.na(Precision)) > 1) sd(Precision, na.rm = TRUE) / sqrt(sum(!is.na(Precision))) else NA_real_,
  Recall        = mean(Recall),
  Recall_SE     = if (.N > 1) sd(Recall) / sqrt(.N) else NA_real_,
  PR            = mean(PR, na.rm = TRUE)
), by = .(tag, cell, arm)]
setorder(pooled, tag, cell, arm)

say("\n%-6s %-10s %-16s %5s %6s %6s %8s %8s\n", "tag", "cell", "arm", "nreps", "nsig", "mTP", "Precision", "Recall")
for (i in seq_len(nrow(pooled))) with(pooled[i], say(
  "%-6s %-10s %-16s %5d %6d %6.2f %8s %5.3f+-%.3f\n",
  tag, cell, arm, n_reps, n_sig_reps, mean_TP,
  if (is.na(Precision)) "NA" else sprintf("%.3f", Precision), Recall, Recall_SE))

OUT <- file.path(stage_dir(STAGE), "pooled_pr.rds")
dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(per_replicate = all_scored, pooled = pooled), OUT)
write_receipt(STAGE, inputs = INPUTS, params = PARAMS, outputs = OUT)
say("\n[3] wrote %s\n    receipt: %s\n", OUT, receipt_path(STAGE))
