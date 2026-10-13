## =============================================================================
## module_3sp/R/17_span_vs_simulation.R
##
## Are three-spined candidate regions more extended than regions around
## simulated QTN? Compares EMMAX consensus region spans (three-spined data)
## with true-positive EMMAX consensus regions of the 600 primary c = 1
## simulations (floor 2), in base pairs and relative to each data set's LD
## decay distance d_0.95 (from the size-predicted decay rate a_pred; for the
## simulations, the median over the run's chromosomes).
## Output: out/17_span_vs_simulation/span_comparison.tsv
## =============================================================================
suppressMessages(library(data.table))
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_3sp"), "R", "00_config.R"))
SIM <- path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis")
d95 <- function(a) 0.95 / (a * 0.05) / 1e3   # d_from_rho(a, 0.95), in kb

d <- as.data.table(readRDS(file.path(SIM, "results/simulation_stage2_region_details.rds")))
d <- d[method == "emmax_consensus_region" & grepl("_c1$", cell) & size_floor == 2 & TP == TRUE]
runs <- unique(d[, .(tag, cell, rep, env)])
dd <- rbindlist(lapply(seq_len(nrow(runs)), function(i) { r <- runs[i]
  f <- file.path(SIM, sprintf("out_final_v1/02_build_ld_units/%s_%s_rep%d_env%d/ld_units.rds", r$tag, r$cell, r$rep, r$env))
  r[, d95_kb := median(d95(as.data.table(readRDS(f)$LD_decay$decay_sum)$a_pred))] }))
d <- merge(d, dd, by = c("tag", "cell", "rep", "env"))[, `:=`(span_kb = (to - from) / 1e3)][, rel := span_kb / d95_kb]

b <- readRDS(file.path(stage_dir("02_bundle"), "bundle.rds"))
sc <- readRDS(file.path(stage_dir("03_EMMAX"), "scan.rds"))
e <- copy(sc$consensus$test$regions)[, span_kb := (to - from) / 1e3]
e <- merge(e, as.data.table(b$LD_decay$decay_sum)[, .(Chr = as.character(Chr), d95_kb = d95(a_pred))], by = "Chr")[, rel := span_kb / d95_kb]

out <- rbind(
  data.table(set = "simulated_TP", n = nrow(d), median_span_kb = median(d$span_kb), median_d95_kb = median(dd$d95_kb),
             median_rel = median(d$rel), p99_rel = quantile(d$rel, .99, names = FALSE), max_rel = max(d$rel), n_above_sim_max = NA_integer_),
  data.table(set = "threespined", n = nrow(e), median_span_kb = median(e$span_kb), median_d95_kb = median(e$d95_kb),
             median_rel = median(e$rel), p99_rel = quantile(e$rel, .99, names = FALSE), max_rel = max(e$rel), n_above_sim_max = sum(e$rel > max(d$rel))))
print(out)
o <- stage_dir("17_span_vs_simulation"); dir.create(o, showWarnings = FALSE, recursive = TRUE)
fwrite(out, file.path(o, "span_comparison.tsv"), sep = "\t")
