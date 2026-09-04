## module_sim/R/05_pool.R
##
## Two-level POOLED-COUNT scoring of R/04_score.R's per-(tag,cell,rep,env,arm)
## TP/FP/FN (each already summed over both chromosomes in that combo's bundle,
## from a single evaluate_ors() call over all its stage-1 units -- see
## R/04_score.R). PK, repeated explicitly: "precision recall should be
## estimated for all chromosomes pooled" -- i.e. SUM TP/FP/FN first, then
## compute ONE Precision/Recall from the totals, not the mean of each
## combination's own Precision. Summation is associative, so summing
## combo-level TP/FP/FN across reps/envs is exactly equal to summing the
## underlying 20 (or up to 200, once envs are crossed) chromosomes' own
## TP/FP/FN directly -- nothing per-chromosome needs to change upstream.
##
## Two pooling levels, because ENV and REP are not the same kind of axis
## (00_config.R's ENVS comment): ENV is the true replicate under one fixed
## recombination map; REP is a different map entirely.
##   1. per_rep: pool over ENV within each (tag, cell, rep, arm) -- the
##      chromosomes-pooled PR for ONE genomic architecture.
##   2. pooled:  pool per_rep's TP/FP/FN further, over REP, by (tag, cell,
##      arm) -- the chromosomes-and-maps-pooled PR PK's original "10
##      simulations, 20 chromosomes" description describes directly, once
##      REPS is understood as the map axis rather than the replicate axis.
##
## This SUPERSEDES the previous macro-averaging (mean +- SE of each
## replicate's own Precision) design -- dropped 2026-09-05 on PK's repeated,
## explicit instruction, not merely refined. Pooled counts have no natural SE
## the way a mean-of-independent-estimates does; none is reported here.
##
## No target= subdirectory: unlike 01-04, this stage pools ACROSS reps and
## envs, so its receipt/output are keyed by (tag, cell) only.
suppressMessages({library(data.table)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim"), "R", "00_config.R"))
STAGE <- "05_pool"
say("=== %s ===\n\n", STAGE)

CELLS_ALL <- c("V0.5_c1", "V0.5_c1.5", "V0.5_c2", "V1_c1", "V1_c1.5", "V2_c1", "V2_c1.5")
TAGS_ALL  <- c("nobgs", "bgs")
REPS_N    <- 10L
ENVS_N    <- 10L

score_files <- data.table(expand.grid(tag = TAGS_ALL, cell = CELLS_ALL, rep = seq_len(REPS_N),
                                       env = seq_len(ENVS_N), stringsAsFactors = FALSE))
score_files[, path := file.path(PATHS$out, "04_score",
  sprintf("score_%s_rep%d_%s_env%d.rds", tag, rep, cell, env))]
score_files[, exists := file.exists(path)]
say("[0] %d/%d expected (tag,cell,rep,env) score files present\n", sum(score_files$exists), nrow(score_files))
if (any(!score_files$exists)) {
  say("    %d missing (run R/04_score.R for these first) -- e.g.:\n", sum(!score_files$exists))
  missing <- score_files[exists == FALSE]
  for (i in seq_len(min(10, nrow(missing)))) say("      %s %s rep%d env%d\n",
    missing$tag[i], missing$cell[i], missing$rep[i], missing$env[i])
}

INPUTS <- score_files[exists == TRUE, path]
PARAMS <- list(cells = CELLS_ALL, tags = TAGS_ALL, reps_n = REPS_N, envs_n = ENVS_N, pooling = "counts",
               fp_by_size = TRUE)
if (!stage_stale(STAGE, INPUTS, PARAMS) && !nzchar(Sys.getenv("FORCE"))) {
  say("\nNothing to do. Set FORCE=1 to rerun anyway.\n"); quit(save = "no")
}

say("\n[1] reading %d score files\n", length(INPUTS))
raw <- lapply(INPUTS, readRDS)
all_scored     <- rbindlist(lapply(raw, `[[`, "scored"))
cluster_detail <- rbindlist(lapply(raw, `[[`, "cluster_detail"))

## sum TP/FP/FN, THEN divide -- not mean of each row's own Precision/Recall.
.pool_counts <- function(TP, FP, FN) {
  TP <- sum(TP); FP <- sum(FP); FN <- sum(FN)
  list(TP = TP, FP = FP, FN = FN,
       Precision = if ((TP + FP) > 0) TP / (TP + FP) else NA_real_,
       Recall    = if ((TP + FN) > 0) TP / (TP + FN) else NA_real_)
}

say("[2] per-rep pooling: sum TP/FP/FN over ENV (both chromosomes each), by (tag, cell, rep, arm)\n")
per_rep <- all_scored[, {p <- .pool_counts(TP, FP, FN)
  list(n_envs = .N, n_sig = sum(n_sig), TP = p$TP, FP = p$FP, FN = p$FN,
       Precision = p$Precision, Recall = p$Recall)}, by = .(tag, cell, rep, arm)]
setorder(per_rep, tag, cell, rep, arm)

say("[3] across-rep pooling: sum per_rep's TP/FP/FN further over REP, by (tag, cell, arm)\n")
pooled <- per_rep[, {p <- .pool_counts(TP, FP, FN)
  list(n_reps = .N, n_envs_total = sum(n_envs), n_sig = sum(n_sig),
       TP = p$TP, FP = p$FP, FN = p$FN, Precision = p$Precision, Recall = p$Recall)}, by = .(tag, cell, arm)]
setorder(pooled, tag, cell, arm)

say("\n%-6s %-10s %-16s %6s %5s %5s %5s %9s %9s\n", "tag", "cell", "arm", "nenvs", "TP", "FP", "FN", "Precision", "Recall")
for (i in seq_len(nrow(pooled))) with(pooled[i], say(
  "%-6s %-10s %-16s %6d %5d %5d %5d %9s %9s\n",
  tag, cell, arm, n_envs_total, TP, FP, FN,
  if (is.na(Precision)) "NA" else sprintf("%.3f", Precision),
  if (is.na(Recall)) "NA" else sprintf("%.3f", Recall)))

## ---- 4. FP proportion as a function of cluster size (PK) -----------------------
## Pooled across every (rep,env) -- one row per significant cluster/unit
## already carries its own size (n_loci) and TP/FP status from
## R/04_score.R's .diagnose_ors() call, so this is a straight pooled count by
## size bin, same "sum first, then divide" rule as the PR pooling above. Bins
## fixed rather than data-driven (quantile bins would shift under -- and so
## be incomparable across -- different cells/arms/tags): SIZE_FLOOR=2 is the
## smallest scored cluster; the rest are round-number thresholds wide enough
## to keep bin counts usable while size 2 and 3 (the most common, and where
## detection differs most) get their own bins rather than being folded in.
say("[4] FP proportion by cluster size (pooled across all rep,env), by (tag, arm, size_bin)\n")
SIZE_BREAKS <- c(1, 2, 3, 5, 10, 20, 50, Inf)
SIZE_LABELS <- c("2", "3", "4-5", "6-10", "11-20", "21-50", "50+")
cluster_detail[, size_bin := cut(n_loci, breaks = SIZE_BREAKS, labels = SIZE_LABELS)]
fp_by_size <- cluster_detail[, .(
  n_sig = .N, TP = sum(is_TP), FP = sum(is_FP),
  FP_proportion = sum(is_FP) / .N
), by = .(tag, arm, size_bin)]
setorder(fp_by_size, tag, arm, size_bin)

OUT <- file.path(stage_dir(STAGE), "pooled_pr.rds")
dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(per_replicate = all_scored, per_rep = per_rep, pooled = pooled,
            cluster_detail = cluster_detail, fp_by_size = fp_by_size), OUT)
write_receipt(STAGE, inputs = INPUTS, params = PARAMS, outputs = OUT)
say("\n[5] wrote %s\n    receipt: %s\n", OUT, receipt_path(STAGE))
