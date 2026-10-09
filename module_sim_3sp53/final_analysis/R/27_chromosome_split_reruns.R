## =============================================================================
## final_analysis/R/27_chromosome_split_reruns.R
##
## PK (2026-10-09): separate the QTN chromosome (Chr1) from the neutral one
## (Chr2) for the two analyses whose saved outputs keep only per-run counts:
## the floor decomposition (R/12 -> R/16) and the marker-QTN redundancy /
## singleton counts (R/17 -> R/18). Primary c=1 runs only (600).
##
## Nothing else changes: this calls R/12's and R/17's own per-run functions,
## which now also return per-chromosome columns (additive; their existing
## columns are untouched), and writes NEW files only:
##   out_final_v1/27_floor_decomposition_raw_c1_chrsplit.rds
##   results/simulation_marker_qtn_redundancy_raw_c1_chrsplit.tsv
## Validation: the existing per-run columns must reproduce the saved
## results/simulation_marker_qtn_redundancy_raw.tsv (c=1 rows) exactly, and
## the floor-decomposition pooled TP/FP must reproduce
## results/simulation_floor_sweep_c1.tsv. Summaries are in R/26.
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(parallel)})
MODULE_ROOT <- path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis")
setwd(MODULE_ROOT)
source(file.path(MODULE_ROOT, "R", "12_floor_decomposition.R"))   ## floor_decomposition_one_combo()
source(file.path(MODULE_ROOT, "R", "17_marker_qtn_redundancy.R")) ## redundancy_one_combo()
say("=== 27_chromosome_split_reruns ===\n\n")

C1_CELLS <- c("V0.5_c1", "V1_c1", "V2_c1")
combos <- CJ(tag = TAGS_ALL, cell = C1_CELLS, rep = REPS_ALL, env = ENVS_ALL)
CORES <- as.integer(Sys.getenv("R27_CORES", "10"))
say("[1] %d c=1 combos, %d cores\n", nrow(combos), CORES)

t0 <- Sys.time()
res <- mclapply(seq_len(nrow(combos)), function(i) {
  cb <- combos[i]
  tryCatch(list(ok = TRUE,
                floor = floor_decomposition_one_combo(cb$tag, cb$cell, cb$rep, cb$env),
                red = redundancy_one_combo(cb$tag, cb$cell, cb$rep, cb$env)),
           error = function(e) list(ok = FALSE, error = sprintf("%s_%s_rep%d_env%d: %s", cb$tag, cb$cell,
                                                               cb$rep, cb$env, conditionMessage(e))))
}, mc.cores = CORES)
ok <- vapply(res, `[[`, logical(1), "ok")
if (!all(ok)) { message(paste(vapply(res[!ok], `[[`, character(1), "error"), collapse = "\n"))
  stop(sprintf("%d/%d combos failed", sum(!ok), length(ok))) }
fl <- rbindlist(lapply(res, `[[`, "floor"), fill = TRUE)
red <- rbindlist(lapply(res, `[[`, "red"), fill = TRUE)
say("[2] all %d combos OK in %.1f min\n", nrow(combos), as.numeric(difftime(Sys.time(), t0, units = "mins")))

## ---- validation against the saved outputs ------------------------------------
stopifnot(all(fl$TP == fl$TP_qtnchr), all(fl$FP == fl$FP_ntrl + fl$FP_qtnchr), all(fl$n_regions_ntrl == fl$FP_ntrl))
old_red <- fread("results/simulation_marker_qtn_redundancy_raw.tsv")[cell %chin% C1_CELLS]
key <- c("tag", "cell", "rep", "env", "method")
cols <- c("n_significant", "n_qtn_linked_naive", "n_false_positive", "n_unique_qtn_recovered", "n_redundant")
chk <- merge(old_red[, c(key, cols), with = FALSE], red[, c(key, cols), with = FALSE], by = key, suffixes = c(".old", ".new"))
stopifnot(nrow(chk) == nrow(old_red), nrow(chk) == nrow(red))
for (cc in cols) stopifnot(identical(chk[[paste0(cc, ".old")]], chk[[paste0(cc, ".new")]]))
say("[3] redundancy counts reproduce simulation_marker_qtn_redundancy_raw.tsv exactly (%d rows)\n", nrow(chk))
sweep_c1 <- fread("results/simulation_floor_sweep_c1.tsv")
pooled <- fl[, .(TP = sum(TP), FP = sum(FP)), by = .(method, arm, floor)]
cmp <- merge(sweep_c1, pooled, by = c("method", "arm", "floor"), suffixes = c(".saved", ".new"))
say("[4] floor decomposition: %d method x arm x floor cells compared; max |TP diff| = %d, max |FP diff| = %d\n",
    nrow(cmp), max(abs(cmp$TP.saved - cmp$TP.new)), max(abs(cmp$FP.saved - cmp$FP.new)))
stopifnot(nrow(cmp) > 0, all(cmp$TP.saved == cmp$TP.new), all(cmp$FP.saved == cmp$FP.new))

saveRDS(fl, "out_final_v1/27_floor_decomposition_raw_c1_chrsplit.rds", compress = "xz")
fwrite(red, "results/simulation_marker_qtn_redundancy_raw_c1_chrsplit.tsv", sep = "\t")
write_receipt("27_chromosome_split_reruns", inputs = character(),
              params = list(cells = C1_CELLS, neutral_rule = "no QTN with non-zero effect"),
              outputs = c("out_final_v1/27_floor_decomposition_raw_c1_chrsplit.rds",
                          "results/simulation_marker_qtn_redundancy_raw_c1_chrsplit.tsv"))
say("[5] wrote out_final_v1/27_floor_decomposition_raw_c1_chrsplit.rds, results/simulation_marker_qtn_redundancy_raw_c1_chrsplit.tsv\n")
