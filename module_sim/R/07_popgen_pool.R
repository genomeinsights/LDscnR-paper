## module_sim/R/07_popgen_pool.R
##
## Pool R/06_popgen_summary.R's per-(tag,cell,rep,env) simulation-property
## summary (Fst, detectable QTN, Va, local adaptation) across the 10
## environments (the replicate axis, same reasoning as R/05_pool.R) into
## mean +- SE by (tag, cell). Unlike the TP/FP counts in 05_pool.R, these are
## already per-replicate CONTINUOUS quantities (not counts to sum first) --
## ordinary mean +- SE across environments is the right operation here, no
## pooled-then-divide step needed.
##
## No target= subdirectory: pools across reps and envs, keyed by (tag, cell)
## only, same convention as 05_pool.R.
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
PARAMS <- list(cells = CELLS_ALL, tags = TAGS_ALL, reps_n = REPS_N, envs_n = ENVS_N)
if (!stage_stale(STAGE, INPUTS, PARAMS) && !nzchar(Sys.getenv("FORCE"))) {
  say("\nNothing to do. Set FORCE=1 to rerun anyway.\n"); quit(save = "no")
}

say("\n[1] reading %d popgen summaries\n", length(INPUTS))
all_summary <- rbindlist(lapply(INPUTS, readRDS))

.se <- function(x) { x <- x[!is.na(x)]; if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_ }
METRICS <- c("n_qtn_total", "n_qtn_detectable", "Va_total", "Va_detectable", "Fst", "local_adapt_r2")

say("[2] per-rep mean, over ENV within each (tag, cell, rep) -- 10 envs each\n")
per_rep <- all_summary[, {
  vals <- lapply(METRICS, function(m) mean(get(m), na.rm = TRUE))
  names(vals) <- METRICS
  c(list(n_envs = .N), vals)
}, by = .(tag, cell, rep)]

say("[3] across-rep pooling: mean +- SE, over the 10 environments' rep-level means, by (tag, cell)\n")
pooled <- per_rep[, {
  out <- list(n_reps = .N)
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

OUT <- file.path(stage_dir(STAGE), "popgen_summary.rds")
dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(per_replicate = all_summary, per_rep = per_rep, pooled = pooled), OUT)
write_receipt(STAGE, inputs = INPUTS, params = PARAMS, outputs = OUT)
say("\n[4] wrote %s\n    receipt: %s\n", OUT, receipt_path(STAGE))
