## module_sim/R/07_popgen_pool.R
##
## Pool R/06_popgen_summary.R's per-(tag,cell,rep,env) simulation-property
## summary (Fst, detectable QTN, Va, local adaptation) across the 10
## environments (the replicate axis -- 00_config.R's ENVS comment: REP is a
## different recombination map each time, ENV is the true replicate under a
## fixed map) into mean +- SE by (tag, cell). Unlike the TP/FP counts in
## 05_pool.R, these are already per-replicate CONTINUOUS quantities (not
## counts to sum first) -- ordinary mean +- SE is the right operation, no
## pooled-then-divide step needed.
##
## [!] FIXED 2026-09-05: the first version of this file grouped per_rep by
## (tag,cell,rep) [averaging over env] then computed SE across those 10
## REP-level values by (tag,cell) -- backwards from 05_pool.R's established
## env-is-the-replicate-axis rule, the same mistake found and fixed there on
## 2026-09-05. Swapped: per_env averages over REP within (tag,cell,env), and
## the final SE is computed across the 10 ENVIRONMENTS' per_env values.
##
## Also computes the BGS EFFECT (PK: "a measure of the effect of bgs, a
## separate figure") as a PAIRED contrast -- every (cell,rep,env) has both a
## bgs and a nobgs run, so log2(bgs/nobgs) is computed at that finest grain
## first (properly paired, not two independently-pooled tag summaries
## differenced afterward), then pooled the same way: mean over REP within
## (cell,env), then mean +- SE over the 10 ENVIRONMENTS, by cell.
##
## No target= subdirectory: pools across reps and envs, keyed by (tag, cell)
## (or just cell, for bgs_effect) only, same convention as 05_pool.R.
suppressMessages({library(data.table)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim"), "R", "00_config.R"))
STAGE <- "07_popgen_pool"
say("=== %s ===\n\n", STAGE)

CELLS_ALL <- c("V0.5_c1", "V0.5_c1.5", "V0.5_c2", "V1_c1", "V1_c1.5", "V2_c1", "V2_c1.5")
TAGS_ALL  <- c("nobgs", "bgs")
REPS_N    <- 10L
ENVS_N    <- 10L

files <- data.table(expand.grid(tag = TAGS_ALL, cell = CELLS_ALL, rep = seq_len(REPS_N),
                                 env = seq_len(ENVS_N), stringsAsFactors = FALSE))
files[, path := file.path(PATHS$out, "06_popgen_summary",
  sprintf("popgen_%s_rep%d_%s_env%d.rds", tag, rep, cell, env))]
files[, exists := file.exists(path)]
say("[0] %d/%d expected (tag,cell,rep,env) popgen summaries present\n", sum(files$exists), nrow(files))
if (any(!files$exists)) {
  say("    %d missing -- e.g.:\n", sum(!files$exists))
  missing <- files[exists == FALSE]
  for (i in seq_len(min(10, nrow(missing)))) say("      %s %s rep%d env%d\n",
    missing$tag[i], missing$cell[i], missing$rep[i], missing$env[i])
}

INPUTS <- files[exists == TRUE, path]
PARAMS <- list(cells = CELLS_ALL, tags = TAGS_ALL, reps_n = REPS_N, envs_n = ENVS_N, bgs_effect = TRUE)
if (!stage_stale(STAGE, INPUTS, PARAMS) && !nzchar(Sys.getenv("FORCE"))) {
  say("\nNothing to do. Set FORCE=1 to rerun anyway.\n"); quit(save = "no")
}

say("\n[1] reading %d popgen summaries\n", length(INPUTS))
all_summary <- rbindlist(lapply(INPUTS, readRDS))

.se <- function(x) { x <- x[!is.na(x)]; if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_ }
METRICS <- c("n_qtn_total", "n_qtn_detectable", "Va_total", "Va_detectable", "Fst", "local_adapt_r2")

say("[2] per-env mean, over REP within each (tag, cell, env) -- 10 reps each\n")
per_env <- all_summary[, {
  vals <- lapply(METRICS, function(m) mean(get(m), na.rm = TRUE))
  names(vals) <- METRICS
  c(list(n_reps = .N), vals)
}, by = .(tag, cell, env)]

say("[3] across-env pooling: mean +- SE, over the 10 environments, by (tag, cell)\n")
pooled <- per_env[, {
  out <- list(n_envs = .N)
  for (m in METRICS) {
    out[[m]] <- mean(get(m), na.rm = TRUE)
    out[[paste0(m, "_SE")]] <- .se(get(m))
  }
  out
}, by = .(tag, cell)]
setorder(pooled, tag, cell)

say("\n%-6s %-10s %8s %8s %10s %8s %10s\n", "tag", "cell", "n_QTN_d", "Va_tot", "Fst", "LA_R2", "")
for (i in seq_len(nrow(pooled))) with(pooled[i], say(
  "%-6s %-10s %5.2f+-%.2f %6.4f+-%.4f %6.4f+-%.4f %6.3f+-%.3f\n",
  tag, cell, n_qtn_detectable, n_qtn_detectable_SE, Va_total, Va_total_SE,
  Fst, Fst_SE, local_adapt_r2, local_adapt_r2_SE))

## ---- 4. BGS effect: paired log2(bgs/nobgs), matched at (cell,rep,env) ---------
## PAIRED, not two pooled tag summaries subtracted afterward -- every
## (cell,rep,env) has both a bgs and a nobgs run, so the contrast is computed
## at that finest grain first (removing shared cell/rep/env variation the
## way a paired t-test would), then pooled the same way as everything else
## here: mean over REP within (cell,env), then mean +- SE over the 10
## ENVIRONMENTS, by cell. log2 so a null effect is 0 and the sign is
## symmetric (negative = bgs lower than nobgs); n_qtn_detectable can be 0 for
## a single (cell,rep,env) (small counts, 1-3 QTN), making its ratio
## undefined -- those pairs are dropped (both a 0/x and x/0 are NA, never a
## silent 0 or Inf) and n_valid_reps tracks how many of the 10 reps' pairs
## actually contributed, per (cell,env).
say("\n[4] BGS effect: paired log2(bgs/nobgs) at (cell,rep,env), by cell\n")
wide <- dcast(all_summary, cell + rep + env ~ tag, value.var = METRICS)
log2r <- function(bgs, nobgs) { r <- bgs / nobgs; r[!is.finite(r) | r <= 0] <- NA_real_; log2(r) }
for (m in METRICS) wide[, (paste0("log2_", m)) := log2r(get(paste0(m, "_bgs")), get(paste0(m, "_nobgs")))]
LOG2_COLS <- paste0("log2_", METRICS)

per_env_effect <- wide[, {
  out <- list(n_valid_reps = sum(!is.na(log2_Fst)))  # Fst always defined; representative denominator
  for (col in LOG2_COLS) out[[col]] <- mean(get(col), na.rm = TRUE)
  out
}, by = .(cell, env)]

bgs_effect <- per_env_effect[, {
  out <- list(n_envs = .N)
  for (col in LOG2_COLS) {
    out[[col]] <- mean(get(col), na.rm = TRUE)
    out[[paste0(col, "_SE")]] <- .se(get(col))
  }
  out
}, by = .(cell)]
setorder(bgs_effect, cell)

say("\n%-10s %10s %10s %10s\n", "cell", "log2 Fst", "log2 Va", "log2 LA_R2")
for (i in seq_len(nrow(bgs_effect))) with(bgs_effect[i], say(
  "%-10s %6.3f+-%.3f %6.3f+-%.3f %6.3f+-%.3f\n",
  cell, log2_Fst, log2_Fst_SE, log2_Va_total, log2_Va_total_SE, log2_local_adapt_r2, log2_local_adapt_r2_SE))

OUT <- file.path(stage_dir(STAGE), "popgen_summary.rds")
dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(per_replicate = all_summary, per_env = per_env, pooled = pooled,
            per_env_effect = per_env_effect, bgs_effect = bgs_effect), OUT)
write_receipt(STAGE, inputs = INPUTS, params = PARAMS, outputs = OUT)
say("\n[5] wrote %s\n    receipt: %s\n", OUT, receipt_path(STAGE))
