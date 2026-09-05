## module_sim_bgs2/R_audit/fp_by_size_fixed.R
##
## Audit response (audit_simulations.tex, item 3): recompute FP-by-cluster-size
## with a MATCHED denominator at both levels, and add a cell-stratified
## version, using the already-archived bgs2 cluster_detail. NOBGS ONLY (PK,
## 2026-09-06: continue bgs2 work with nobgs only while the bgs5 grid runs).
##
## The bug: the archived fp_by_size's point estimate was FP/(TP+FP) (excluding
## dedup-neutral rows -- a duplicate claim on an already-claimed true QTN,
## neither TP nor FP, see R/04_score.R's .score_arm()), but its SE was built
## from a per-env proportion that divided by ALL significant rows (.N),
## INCLUDING those neutral ones -- point and SE were different estimators.
## Same fix as R/05_pool.R (live pipeline, fixed 2026-09-06): both now divide
## by (TP+FP) at every level. Also added: cell-stratified estimates, since
## pooling all 7 cells together (as the archived version did) lets cluster-
## size trends partly reflect changing cell composition across the grid
## rather than a within-cell size effect.
suppressMessages({library(data.table); library(ggplot2)})

POOL_PATH <- path.expand("~/gitlab/LDscnR-paper/module_sim_bgs2/results/pooled_pr.rds")
x <- readRDS(POOL_PATH)
cluster_detail <- as.data.table(x$cluster_detail)[tag == "nobgs"]
cat(sprintf("[0] %d significant-cluster rows, nobgs only\n", nrow(cluster_detail)))

.se <- function(v) { v <- v[!is.na(v)]; if (length(v) > 1) sd(v) / sqrt(length(v)) else NA_real_ }
SIZE_BREAKS <- c(1, 2, 3, 5, 10, 20, 50, Inf)
SIZE_LABELS <- c("2", "3", "4-5", "6-10", "11-20", "21-50", "50+")
cluster_detail[, size_bin := cut(n_loci, breaks = SIZE_BREAKS, labels = SIZE_LABELS)]

## pooled across cells -- matched denominator at both levels
per_env_size <- cluster_detail[, .(
  n_sig = .N, TP = sum(is_TP), FP = sum(is_FP),
  FP_proportion = if ((sum(is_TP) + sum(is_FP)) > 0) sum(is_FP) / (sum(is_TP) + sum(is_FP)) else NA_real_
), by = .(env, arm, size_bin)]
fp_by_size <- per_env_size[, .(
  n_envs = sum(!is.na(FP_proportion)), n_sig = sum(n_sig),
  FP_proportion = sum(FP) / (sum(TP) + sum(FP)),
  FP_proportion_SE = .se(FP_proportion)
), by = .(arm, size_bin)]

## cell-stratified -- same rule
per_env_size_cell <- cluster_detail[, .(
  n_sig = .N, TP = sum(is_TP), FP = sum(is_FP),
  FP_proportion = if ((sum(is_TP) + sum(is_FP)) > 0) sum(is_FP) / (sum(is_TP) + sum(is_FP)) else NA_real_
), by = .(cell, env, arm, size_bin)]
fp_by_size_cell <- per_env_size_cell[, .(
  n_envs = sum(!is.na(FP_proportion)), n_sig = sum(n_sig),
  FP_proportion = sum(FP) / (sum(TP) + sum(FP)),
  FP_proportion_SE = .se(FP_proportion)
), by = .(cell, arm, size_bin)]

cat("[1] pooled (nobgs only, matched denominator):\n")
print(fp_by_size[order(arm, size_bin)])

OUT <- path.expand("~/gitlab/LDscnR-paper/module_sim_bgs2/results/fp_by_size_fixed_nobgs.rds")
saveRDS(list(cluster_detail = cluster_detail, per_env_size = per_env_size, fp_by_size = fp_by_size,
            per_env_size_cell = per_env_size_cell, fp_by_size_cell = fp_by_size_cell), OUT)
cat(sprintf("\n[2] wrote %s\n", OUT))
