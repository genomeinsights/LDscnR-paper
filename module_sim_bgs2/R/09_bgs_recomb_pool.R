## module_sim/R/09_bgs_recomb_pool.R
##
## Pool R/08_bgs_recomb.R's per-(tag,cell,rep,env,recomb_bin) Fst into the
## BGS effect (paired log2(bgs/nobgs), matched at cell,rep,env,bin -- same
## logic as R/07_popgen_pool.R's genome-wide bgs_effect) BY RECOMBINATION-RATE
## BIN. PK: "measure the difference between low and high recombination
## regions to see how bgs affects patterns of genetic differentiation as a
## function of recombination rate."
##
## Pooled across ALL CELLS as well as reps/envs for this one, unlike
## 05_pool.R/07_popgen_pool.R's per-cell breakdown -- this is a question
## about the BGS/recombination MECHANISM itself (expected to hold regardless
## of selection intensity or dispersal), and pooling across cells gives it
## much more power (up to 700 (cell,rep,env) triples per bin) to see clearly
## rather than splitting into noisy per-cell estimates.
suppressMessages({library(data.table)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim"), "R", "00_config.R"))
STAGE <- "09_bgs_recomb_pool"
say("=== %s ===\n\n", STAGE)

CELLS_ALL <- c("V0.5_c1", "V0.5_c1.5", "V0.5_c2", "V1_c1", "V1_c1.5", "V2_c1", "V2_c1.5")
TAGS_ALL  <- c("nobgs", "bgs")
REPS_N    <- 10L
ENVS_N    <- 10L

files <- data.table(expand.grid(tag = TAGS_ALL, cell = CELLS_ALL, rep = seq_len(REPS_N),
                                 env = seq_len(ENVS_N), stringsAsFactors = FALSE))
files[, path := file.path(PATHS$out, "08_bgs_recomb",
  sprintf("bgsrecomb_%s_rep%d_%s_env%d.rds", tag, rep, cell, env))]
files[, exists := file.exists(path)]
say("[0] %d/%d expected (tag,cell,rep,env) bgs-recomb summaries present\n", sum(files$exists), nrow(files))
if (any(!files$exists)) {
  say("    %d missing -- e.g.:\n", sum(!files$exists))
  missing <- files[exists == FALSE]
  for (i in seq_len(min(10, nrow(missing)))) say("      %s %s rep%d env%d\n",
    missing$tag[i], missing$cell[i], missing$rep[i], missing$env[i])
}

INPUTS <- files[exists == TRUE, path]
PARAMS <- list(cells = CELLS_ALL, tags = TAGS_ALL, reps_n = REPS_N, envs_n = ENVS_N, pooled_across_cells = TRUE)
if (!stage_stale(STAGE, INPUTS, PARAMS) && !nzchar(Sys.getenv("FORCE"))) {
  say("\nNothing to do. Set FORCE=1 to rerun anyway.\n"); quit(save = "no")
}

say("\n[1] reading %d bgs-recomb summaries\n", length(INPUTS))
all_summary <- rbindlist(lapply(INPUTS, readRDS))
BIN_LEVELS <- c("zero", "low", "medium", "high")
all_summary[, recomb_bin := factor(recomb_bin, levels = BIN_LEVELS)]

.se <- function(x) { x <- x[!is.na(x)]; if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_ }

say("[2] paired log2(bgs/nobgs) at (cell,rep,env,bin)\n")
wide <- dcast(all_summary, cell + rep + env + recomb_bin ~ tag, value.var = "Fst")
log2r <- function(bgs, nobgs) { r <- bgs / nobgs; r[!is.finite(r) | r <= 0] <- NA_real_; log2(r) }
wide[, log2_Fst := log2r(bgs, nobgs)]

say("[3] pool over REP within (cell,env,bin), then mean +- SE over the 10 x 7 (env,cell) units, by bin\n")
per_env_cell <- wide[, .(log2_Fst = mean(log2_Fst, na.rm = TRUE)), by = .(cell, env, recomb_bin)]
bgs_effect_by_bin <- per_env_cell[, .(
  n_units = sum(!is.na(log2_Fst)),
  log2_Fst = mean(log2_Fst, na.rm = TRUE),
  log2_Fst_SE = .se(log2_Fst)
), by = .(recomb_bin)]
setorder(bgs_effect_by_bin, recomb_bin)

say("\n%-8s %6s %10s\n", "bin", "n", "log2 Fst")
for (i in seq_len(nrow(bgs_effect_by_bin))) with(bgs_effect_by_bin[i], say(
  "%-8s %6d %6.3f+-%.3f\n", recomb_bin, n_units, log2_Fst, log2_Fst_SE))

## also keep the raw (unpaired) mean Fst per tag x bin, for a figure showing
## both curves alongside the effect-size figure if useful.
say("\n[4] mean Fst per (tag,bin), for reference\n")
fst_by_tag_bin <- all_summary[, .(Fst = mean(Fst, na.rm = TRUE), Fst_SE = .se(Fst)), by = .(tag, recomb_bin)]
setorder(fst_by_tag_bin, tag, recomb_bin)
for (i in seq_len(nrow(fst_by_tag_bin))) with(fst_by_tag_bin[i], say(
  "%-6s %-8s %6.4f+-%.4f\n", tag, recomb_bin, Fst, Fst_SE))

OUT <- file.path(stage_dir(STAGE), "bgs_recomb_summary.rds")
dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(per_replicate = all_summary, per_env_cell = per_env_cell,
            bgs_effect_by_bin = bgs_effect_by_bin, fst_by_tag_bin = fst_by_tag_bin), OUT)
write_receipt(STAGE, inputs = INPUTS, params = PARAMS, outputs = OUT)
say("\n[5] wrote %s\n    receipt: %s\n", OUT, receipt_path(STAGE))
