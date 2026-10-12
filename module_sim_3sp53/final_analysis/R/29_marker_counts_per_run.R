## =============================================================================
## final_analysis/R/29_marker_counts_per_run.R
##
## Per simulated data set (one primary c = 1 run): analysed markers, markers
## retained at the floor of two (members of Stage-1 units with >= 2 markers)
## and the number of test units they form. The sums reproduce the pooled
## counts used in the manuscript (9,120,379 / 5,534,487 / 1,796,934).
## Output: results/simulation_marker_counts_per_run_c1.tsv
## =============================================================================
suppressMessages({library(data.table); library(parallel)})
MODULE_ROOT <- path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis")
setwd(MODULE_ROOT)
source(file.path(MODULE_ROOT, "R", "00_config.R"))
say("=== 29_marker_counts_per_run ===\n\n")

runs <- list.files(stage_dir("02_build_ld_units", ""), full.names = TRUE)
runs <- runs[grepl("_c1_rep", basename(runs))]
r <- rbindlist(mclapply(runs, function(x) {
  cl <- as.data.table(readRDS(file.path(x, "ld_units.rds"))$stage1$clusters)
  data.table(run = basename(x), markers = sum(cl$n_snps),
             retained = sum(cl[n_snps >= 2]$n_snps), units = sum(cl$n_snps >= 2))
}, mc.cores = as.integer(Sys.getenv("R29_CORES", "8"))))
stopifnot(nrow(r) == 600L)
say("[1] %d runs; totals %s markers, %s retained, %s units\n", nrow(r),
    format(sum(r$markers), big.mark = ","), format(sum(r$retained), big.mark = ","),
    format(sum(r$units), big.mark = ","))
fwrite(r, "results/simulation_marker_counts_per_run_c1.tsv", sep = "\t")
say("[2] wrote results/simulation_marker_counts_per_run_c1.tsv\n")
