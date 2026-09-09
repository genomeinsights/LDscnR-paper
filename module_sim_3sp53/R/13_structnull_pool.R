## module_sim_3sp53/R/13_structnull_pool.R
##
## Pool R/12_structured_null.R's per-(tag,cell,rep,env) structured-null
## output into one archive file: results/structured_null_summary.rds. Run
## after run_structured_null_grid.sh (or run_structured_null_sample.sh)
## finishes -- picks up whatever is in out/12_structured_null/ at the time,
## so it's safe to run on a partial grid too.

suppressMessages(library(data.table))
files <- list.files("out/12_structured_null", pattern = "^structnull_.*\\.rds$", full.names = TRUE)
cat(length(files), "combo files found\n")
all_summary <- rbindlist(lapply(files, function(f) readRDS(f)$summary))
all_sweep   <- rbindlist(lapply(files, function(f) readRDS(f)$size_sweep))
saveRDS(list(summary = all_summary, size_sweep = all_sweep, n_combos = length(files)),
        "results/structured_null_summary.rds")
cat("wrote results/structured_null_summary.rds:", nrow(all_summary), "summary rows,",
    nrow(all_sweep), "sweep rows\n")
