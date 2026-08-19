## =====================================================================
## module_sim_LDscnR / stage_cache.R
##
## Copy the per-file regen bundles from the external drive to a LOCAL slim cache
## so the analysis scripts read at SSD speed (the runs are I/O-bound on the USB
## drive, not RAM-bound) and can be parallelised across envs. Drops the bulky
## LD_decay edge lists (only decay_sum is used downstream) and keeps everything
## the analysis needs. The cache is regenerable and should not be committed.
##
## Run once from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/stage_cache.R  [V c]        # default V2 c1, env1-5, chr1-10
## Then point the analysis at it:
##   SIM_DATA=module_sim_LDscnR/data/cache Rscript module_sim_LDscnR/grm_comparison.R
## =====================================================================

suppressMessages(library(data.table))
a  <- commandArgs(trailingOnly = TRUE)
V  <- if (length(a) >= 1) a[1] else "2"
CC <- if (length(a) >= 2) a[2] else "1"
TAG  <- "nobgs"; ENVS <- 1:5; CHRS <- 1:10
SRC  <- Sys.getenv("SIM_SRC",   "/Volumes/Nemo/Nemo_sim/regen_sim_data")
DST  <- Sys.getenv("SIM_CACHE", "module_sim_LDscnR/data/cache")
if (!dir.exists(DST)) dir.create(DST, recursive = TRUE)

## fields the analysis reads (GTs, ld_ws, map, env, decay_sum, saved GRM + marker sets)
slim_of <- function(d) list(
  GTs = d$GTs, map = d$map, env = d$env, ld_ws = d$ld_ws,
  LD_decay = list(decay_sum = d$LD_decay$decay_sum),       # drop the edge lists (bulky, unused)
  GRM = d$GRM, pruned_markers = d$pruned_markers,
  grm_markers = d$grm_markers, emx_gif = d$emx_gif)

nfile <- 0L; src_mb <- dst_mb <- 0
for (env in ENVS) for (ch in CHRS) {
  fn  <- sprintf("adapt_%s_chr%d_V%s_c%s_env%d.rds", TAG, ch, V, CC, env)
  src <- file.path(SRC, fn); dst <- file.path(DST, fn)
  if (!file.exists(src)) next
  if (file.exists(dst)) { message("skip (cached): ", fn); next }
  d <- readRDS(src); saveRDS(slim_of(d), dst)
  nfile <- nfile + 1L
  src_mb <- src_mb + file.size(src) / 1e6; dst_mb <- dst_mb + file.size(dst) / 1e6
  message(sprintf("  %s : %.0f -> %.0f MB", fn, file.size(src)/1e6, file.size(dst)/1e6))
}
message(sprintf("\nstaged %d files to %s  (%.0f MB source -> %.0f MB cache)", nfile, DST, src_mb, dst_mb))
