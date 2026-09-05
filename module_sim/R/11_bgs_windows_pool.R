## module_sim/R/11_bgs_windows_pool.R
##
## Pool R/10_bgs_windows.R's per-(tag,cell,rep,env,Chr,win) segregating-site
## counts into the paired estimator from
## /Volumes/Nemo/Nemo_sim/bgs_effect_newmaps/measure_bgs_effect.R:
## B_obs = n_snp(bgs)/n_snp(nobgs) for the SAME window (Chr,win,rep,env,cell),
## binned into recombination-rate quintiles PER CELL. Same design as that
## script, applied to this repo's own grid (7 cells x 10 reps x 10 envs,
## rec_map resolution matching its "old maps" 32-44% zero-rate category).
##
## PER CELL, NEVER POOLED ACROSS CELLS -- matching measure_bgs_effect.R's own
## stated reason: BGS strength varies with map length/selection/dispersal.
## (This differs from R/09_bgs_recomb_pool.R's Fst-based analysis, which
## pooled across cells for power on what was already a null result; there is
## no such motivation to override the reference script's own convention here,
## now that the effect is expected to be real and cell-dependent.)
suppressMessages({library(data.table)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim"), "R", "00_config.R"))
STAGE <- "11_bgs_windows_pool"
say("=== %s ===\n\n", STAGE)

CELLS_ALL <- c("V0.5_c1", "V0.5_c1.5", "V0.5_c2", "V1_c1", "V1_c1.5", "V2_c1", "V2_c1.5")
TAGS_ALL  <- c("nobgs", "bgs")
REPS_N    <- 10L
ENVS_N    <- 10L

files <- data.table(expand.grid(tag = TAGS_ALL, cell = CELLS_ALL, rep = seq_len(REPS_N),
                                 env = seq_len(ENVS_N), stringsAsFactors = FALSE))
files[, path := file.path(PATHS$out, "10_bgs_windows",
  sprintf("bgswin_%s_rep%d_%s_env%d.rds", tag, rep, cell, env))]
files[, exists := file.exists(path)]
say("[0] %d/%d expected (tag,cell,rep,env) window summaries present\n", sum(files$exists), nrow(files))
if (any(!files$exists)) {
  say("    %d missing -- e.g.:\n", sum(!files$exists))
  missing <- files[exists == FALSE]
  for (i in seq_len(min(10, nrow(missing)))) say("      %s %s rep%d env%d\n",
    missing$tag[i], missing$cell[i], missing$rep[i], missing$env[i])
}

INPUTS <- files[exists == TRUE, path]
PARAMS <- list(cells = CELLS_ALL, tags = TAGS_ALL, reps_n = REPS_N, envs_n = ENVS_N,
               estimator = "n_snp_ratio_per_window", quintiles = "per_cell")
if (!stage_stale(STAGE, INPUTS, PARAMS) && !nzchar(Sys.getenv("FORCE"))) {
  say("\nNothing to do. Set FORCE=1 to rerun anyway.\n"); quit(save = "no")
}

say("\n[1] reading %d window summaries\n", length(INPUTS))
all_windows <- rbindlist(lapply(INPUTS, readRDS))

say("[2] pairing bgs/nobgs at (cell,rep,env,Chr,win)\n")
wide <- dcast(all_windows, cell + rep + env + Chr + win + rate ~ tag, value.var = "n_snp")
wide <- wide[!is.na(bgs) & !is.na(nobgs) & nobgs > 0 & !is.na(rate)]
wide[, B_obs := bgs / nobgs]
say("    %d paired windows total, across %d cells\n", nrow(wide), uniqueN(wide$cell))

say("[3] per-cell quintile binning (Q1 = lowest recombination rate)\n")
quintiles <- rbindlist(lapply(CELLS_ALL, function(cl) {
  d <- copy(wide[cell == cl])
  if (nrow(d) < 20) return(NULL)
  d[, q := cut(frank(rate, ties.method = "first"), 5, labels = FALSE)]
  qq <- d[, .(n_windows = .N, rate = median(rate), B = median(B_obs)), by = q][order(q)]
  qq[, cell := cl]
  qq
}))
setcolorder(quintiles, c("cell", "q", "n_windows", "rate", "B"))

genome_wide <- wide[, .(n_paired_windows = .N, genome_ratio = median(B_obs)), by = cell]
contrast <- quintiles[, .(Q1 = B[q == 1], Q5 = B[q == 5]), by = cell]
contrast[, contrast_pct := 100 * (Q5 - Q1) / Q5]
summary_by_cell <- merge(genome_wide, contrast, by = "cell")
setorder(summary_by_cell, cell)

say("\n%-10s %8s %8s %8s %8s %8s\n", "cell", "n_win", "genome", "Q1", "Q5", "Q1/Q5 %")
for (i in seq_len(nrow(summary_by_cell))) with(summary_by_cell[i], say(
  "%-10s %8d %8.3f %8.3f %8.3f %7.1f%%\n", cell, n_paired_windows, genome_ratio, Q1, Q5, contrast_pct))

OUT <- file.path(stage_dir(STAGE), "bgs_windows_summary.rds")
dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(paired_windows = wide, quintiles = quintiles, summary_by_cell = summary_by_cell), OUT)
write_receipt(STAGE, inputs = INPUTS, params = PARAMS, outputs = OUT)
say("\n[4] wrote %s\n    receipt: %s\n", OUT, receipt_path(STAGE))
