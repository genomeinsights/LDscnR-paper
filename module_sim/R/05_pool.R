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
## ENV is the replicate axis for SE (PK, 2026-09-05: "use the environments
## for SE"); REP is a different recombination map each time, not an
## exchangeable replicate (00_config.R's ENVS comment) -- but env_<N>.txt is
## the SAME file for every rep, so pooling REP into each env's estimate is a
## coherent operation (same environmental predictor, different genomic
## architectures), unlike pooling ACROSS env would be.
##   1. per_env: pool over REP within each (tag, cell, env, arm) -- one
##      chromosomes-and-maps-pooled Precision/Recall per environment.
##   2. pooled:  mean +- SE, ACROSS ENV, of per_env's Precision/Recall, by
##      (tag, cell, arm) -- N=10 (the number of environments), the genuine
##      replicate count.
##
## Point estimate is still the pooled-count Precision/Recall (not a mean of
## means) -- the SE is layered on top of that same pooled quantity, computed
## per env, not a second, different estimator.
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
               se_axis = "env", fp_by_size = TRUE,
               ## [!] bumped 2026-09-06: matched fp_by_size denominator + cell
               ## stratification, size_bin now starts at "1", and the new pooled
               ## bootstrap CI -- logic changes stage_stale() cannot see any other
               ## way -- forces a rerun instead of silently reusing pre-fix output.
               pool_version = "2026-09-06-audit")
if (!stage_stale(STAGE, INPUTS, PARAMS) && !nzchar(Sys.getenv("FORCE"))) {
  say("\nNothing to do. Set FORCE=1 to rerun anyway.\n"); quit(save = "no")
}

say("\n[1] reading %d score files\n", length(INPUTS))
raw <- lapply(INPUTS, readRDS)
all_scored     <- rbindlist(lapply(raw, `[[`, "scored"))
cluster_detail <- rbindlist(lapply(raw, `[[`, "cluster_detail"))

## sum TP/FP/FN, THEN divide -- not mean of each row's own Precision/Recall.
## PR = Precision*Recall; beta = FN/(TP+FN) = 1-Recall, the Type II error
## rate (the natural counterpart to ALPHA, 00_config.R's significance
## threshold) -- PK asked for both alongside Precision/Recall.
.pool_counts <- function(TP, FP, FN) {
  TP <- sum(TP); FP <- sum(FP); FN <- sum(FN)
  Precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  Recall    <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  list(TP = TP, FP = FP, FN = FN, Precision = Precision, Recall = Recall,
       PR = Precision * Recall, beta = if ((TP + FN) > 0) FN / (TP + FN) else NA_real_)
}
.se <- function(x) { x <- x[!is.na(x)]; if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_ }

say("[2] per-env pooling: sum TP/FP/FN over REP (both chromosomes each), by (tag, cell, env, arm)\n")
per_env <- all_scored[, {p <- .pool_counts(TP, FP, FN)
  list(n_reps = .N, n_sig = sum(n_sig), TP = p$TP, FP = p$FP, FN = p$FN,
       Precision = p$Precision, Recall = p$Recall, PR = p$PR, beta = p$beta)}, by = .(tag, cell, env, arm)]
setorder(per_env, tag, cell, env, arm)

say("[3] across-env pooling: mean +- SE, over ENV (N=%d), of per_env's Precision/Recall/PR/beta, by (tag, cell, arm)\n", ENVS_N)
pooled <- per_env[, {pc <- .pool_counts(TP, FP, FN)
  list(n_envs = .N, n_reps_total = sum(n_reps), n_sig = sum(n_sig),
       TP = pc$TP, FP = pc$FP, FN = pc$FN,
       Precision = pc$Precision, Precision_SE = .se(Precision),
       Recall = pc$Recall, Recall_SE = .se(Recall),
       PR = pc$PR, PR_SE = .se(PR),
       beta = pc$beta, beta_SE = .se(beta))}, by = .(tag, cell, arm)]
setorder(pooled, tag, cell, arm)

say("\n%-6s %-10s %-16s %6s %5s %5s %5s %14s %14s\n", "tag", "cell", "arm", "nenvs", "TP", "FP", "FN", "Precision", "Recall")
for (i in seq_len(nrow(pooled))) with(pooled[i], say(
  "%-6s %-10s %-16s %6d %5d %5d %5d %6s+-%-6s %6s+-%-6s\n",
  tag, cell, arm, n_envs, TP, FP, FN,
  if (is.na(Precision)) "NA" else sprintf("%.3f", Precision), if (is.na(Precision_SE)) "NA" else sprintf("%.3f", Precision_SE),
  if (is.na(Recall)) "NA" else sprintf("%.3f", Recall), if (is.na(Recall_SE)) "NA" else sprintf("%.3f", Recall_SE)))

## ---- 3b. cluster bootstrap CI for the pooled ratios (external audit item 4) ----
## ADDED 2026-09-06: `pooled`'s Precision/Recall/PR above is a ratio of counts
## summed over ALL environments and reps, but its SE is the SE of the 10
## per-environment ratios -- a fine descriptive summary, but not a formal CI
## on the pooled ratio itself (the audit's point). Added alongside it, not
## replacing it: a two-stage cluster bootstrap that resamples ENVIRONMENTS
## (with replacement -- the genuine replicate axis, as established above) and,
## within each resampled environment, resamples its RECOMBINATION MAPS/reps
## (with replacement) -- exactly the audit's "resampling environments and
## recombination maps" recommendation -- recomputing the pooled ratio each
## time. Percentile (2.5%/97.5%) interval reported.
say("[3b] cluster bootstrap (envs, then reps within env), N_BOOT = %d\n", N_BOOT <- 2000L)
set.seed(SEEDS[["bundle"]])
.boot_pooled <- function(sub) {
  env_idx <- split(seq_len(nrow(sub)), sub$env)
  envs <- names(env_idx)
  n_env <- length(envs)
  TPv <- sub$TP; FPv <- sub$FP; FNv <- sub$FN
  Precision <- Recall <- PR <- numeric(N_BOOT)
  for (b in seq_len(N_BOOT)) {
    boot_envs <- sample(envs, n_env, replace = TRUE)
    TP <- 0; FP <- 0; FN <- 0
    for (e in boot_envs) {
      idx <- env_idx[[e]]
      ridx <- sample(idx, length(idx), replace = TRUE)
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
}, by = .(tag, cell, arm)]
setorder(pooled_boot, tag, cell, arm)
pooled <- merge(pooled, pooled_boot, by = c("tag", "cell", "arm"), all.x = TRUE)
setorder(pooled, tag, cell, arm)

## ---- 4. FP proportion as a function of cluster size (PK) -----------------------
## Same env-as-replicate-axis logic as the PR pooling above: pool over REP
## within (tag, [cell,] env, arm, size_bin) to get one FP proportion per
## environment, then mean +- SE across the 10 environments for the final
## point. Bins fixed rather than data-driven (quantile bins would shift
## under -- and so be incomparable across -- different cells/arms/tags): the
## rest are round-number thresholds wide enough to keep bin counts usable
## while size 1-3 (the most common, and where detection differs most) each
## get their own bin rather than being folded in.
##
## [!] FIXED 2026-09-06 (PK direct correction + external audit item 1):
## SIZE_FLOOR=2 is the smallest CLUSTER-based unit (R/04_score.R's `units`,
## floor-filtered before clustering even starts), but single-SNP arms (fixed
## alongside this in 04_score.R -- each significant marker now its own size-1
## region, not the enclosing multi-marker unit) really do produce n_loci=1
## rows. The breaks used to start at 1 with cut()'s default right-closed
## intervals, under which x==1 exactly falls in NO bin (silently dropped as
## NA) -- "Single SNPs go all the way to 1 not 2." Breaks now start at 0 with
## an explicit "1" label.
##
## [!] FIXED 2026-09-06 (external audit): the point estimate FP/(TP+FP)
## excludes dedup-neutral rows (neither is_TP nor is_FP -- a duplicate claim
## on an already-claimed true QTN, dropped rather than counted either way,
## see R/04_score.R's .score_arm()), but the previous per-env proportion used
## for the SE divided by .N -- ALL significant rows, INCLUDING those neutral
## ones. Point and SE were therefore different estimators. Both now divide by
## (TP+FP) at every level. Also now computed BOTH pooled-across-cells (as
## before, for a single headline number) and cell-stratified (fp_by_size_cell)
## -- the audit's other point, that pooling all selection/dispersal cells
## together lets cluster-size trends partly reflect changing cell composition
## rather than a within-cell size effect.
say("\n[4] FP proportion by cluster size, mean +- SE over ENV, pooled and by-cell\n")
SIZE_BREAKS <- c(0, 1, 2, 3, 5, 10, 20, 50, Inf)
SIZE_LABELS <- c("1", "2", "3", "4-5", "6-10", "11-20", "21-50", "50+")
cluster_detail[, size_bin := cut(n_loci, breaks = SIZE_BREAKS, labels = SIZE_LABELS)]

## pooled across cells
per_env_size <- cluster_detail[, .(
  n_sig = .N, TP = sum(is_TP), FP = sum(is_FP),
  FP_proportion = if ((sum(is_TP) + sum(is_FP)) > 0) sum(is_FP) / (sum(is_TP) + sum(is_FP)) else NA_real_
), by = .(tag, env, arm, size_bin)]
fp_by_size <- per_env_size[, .(
  n_envs = sum(!is.na(FP_proportion)), n_sig = sum(n_sig),
  FP_proportion = sum(FP) / (sum(TP) + sum(FP)),
  FP_proportion_SE = .se(FP_proportion)
), by = .(tag, arm, size_bin)]

## cell-stratified (same denominator rule)
per_env_size_cell <- cluster_detail[, .(
  n_sig = .N, TP = sum(is_TP), FP = sum(is_FP),
  FP_proportion = if ((sum(is_TP) + sum(is_FP)) > 0) sum(is_FP) / (sum(is_TP) + sum(is_FP)) else NA_real_
), by = .(tag, cell, env, arm, size_bin)]
fp_by_size_cell <- per_env_size_cell[, .(
  n_envs = sum(!is.na(FP_proportion)), n_sig = sum(n_sig),
  FP_proportion = sum(FP) / (sum(TP) + sum(FP)),
  FP_proportion_SE = .se(FP_proportion)
), by = .(tag, cell, arm, size_bin)]
setorder(fp_by_size, tag, arm, size_bin)

OUT <- file.path(stage_dir(STAGE), "pooled_pr.rds")
dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(per_replicate = all_scored, per_env = per_env, pooled = pooled,
            cluster_detail = cluster_detail, per_env_size = per_env_size, fp_by_size = fp_by_size,
            per_env_size_cell = per_env_size_cell, fp_by_size_cell = fp_by_size_cell), OUT)
write_receipt(STAGE, inputs = INPUTS, params = PARAMS, outputs = OUT)
say("\n[5] wrote %s\n    receipt: %s\n", OUT, receipt_path(STAGE))
