## module_sim_bgs2/R_audit/pooled_boot_fixed.R
##
## Audit response (audit_simulations.tex, item 4): cluster bootstrap CI for
## the pooled Precision/Recall/PR, NOBGS ONLY, using bgs2's already-archived
## per_replicate table -- pure recompute, no rerun needed. Mirrors the
## bootstrap just added to the live R/05_pool.R: resamples ENVIRONMENTS (with
## replacement, the replicate axis) and, within each resampled environment,
## resamples its RECOMBINATION MAPS/reps (with replacement), recomputing the
## pooled ratio each time. Percentile (2.5%/97.5%) interval, alongside (not
## replacing) the existing across-env SE in `pooled`.
suppressMessages({library(data.table)})

POOL_PATH <- path.expand("~/gitlab/LDscnR-paper/module_sim_bgs2/results/pooled_pr.rds")
x <- readRDS(POOL_PATH)
all_scored <- as.data.table(x$per_replicate)[tag == "nobgs"]
pooled <- as.data.table(x$pooled)[tag == "nobgs"]
cat(sprintf("[0] %d per-replicate rows, nobgs only\n", nrow(all_scored)))

N_BOOT <- 2000L
set.seed(20260906)
.boot_pooled <- function(sub) {
  env_idx <- split(seq_len(nrow(sub)), sub$env)
  envs <- names(env_idx); n_env <- length(envs)
  TPv <- sub$TP; FPv <- sub$FP; FNv <- sub$FN
  Precision <- Recall <- PR <- numeric(N_BOOT)
  for (b in seq_len(N_BOOT)) {
    boot_envs <- sample(envs, n_env, replace = TRUE)
    TP <- 0; FP <- 0; FN <- 0
    for (e in boot_envs) {
      idx <- env_idx[[e]]; ridx <- sample(idx, length(idx), replace = TRUE)
      TP <- TP + sum(TPv[ridx]); FP <- FP + sum(FPv[ridx]); FN <- FN + sum(FNv[ridx])
    }
    Precision[b] <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
    Recall[b]    <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
    PR[b] <- Precision[b] * Recall[b]
  }
  list(Precision = Precision, Recall = Recall, PR = PR)
}
.ci <- function(v) stats::quantile(v, c(0.025, 0.975), na.rm = TRUE, names = FALSE)

pooled_boot <- all_scored[, {
  bt <- .boot_pooled(.SD)
  pci <- .ci(bt$Precision); rci <- .ci(bt$Recall); prci <- .ci(bt$PR)
  list(n_boot = N_BOOT,
       Precision_lo = pci[1], Precision_hi = pci[2],
       Recall_lo = rci[1], Recall_hi = rci[2],
       PR_lo = prci[1], PR_hi = prci[2])
}, by = .(cell, arm)]

pooled_fixed <- merge(pooled, pooled_boot, by = c("cell", "arm"), all.x = TRUE)
setorder(pooled_fixed, cell, arm)
cat("[1] pooled + bootstrap CI (nobgs only):\n")
print(pooled_fixed[, .(cell, arm, Precision, Precision_SE, Precision_lo, Precision_hi,
                        Recall, Recall_SE, Recall_lo, Recall_hi)])

OUT <- path.expand("~/gitlab/LDscnR-paper/module_sim_bgs2/results/pooled_boot_fixed_nobgs.rds")
saveRDS(pooled_fixed, OUT)
cat(sprintf("\n[2] wrote %s\n", OUT))
