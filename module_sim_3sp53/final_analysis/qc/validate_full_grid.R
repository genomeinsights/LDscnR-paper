## final_analysis/qc/validate_full_grid.R
##
## Phase 3 gate 4 validation: independently confirms the full 1,400-
## combination grid, cross-checking the manifest against the actual
## per-combo output files rather than trusting shell exit codes alone.
suppressMessages(library(data.table))
source("R/00_config.R")

CELLS <- CELLS_ALL; TAGS <- TAGS_ALL; REPS <- REPS_ALL; ENVS <- ENVS_ALL
EXPECTED_METHODS <- c("emmax_snp", "emmax_snp_nonsingleton", "emmax_simes", "emmax_consensus")
EXPECTED_FIELDS <- c("method", "n_tested", "n_significant", "TP", "FP", "FN", "bh_crit_p", "precision", "recall",
                     "eligible_marker_fraction", "detectable_qtn_covered", "n_detectable_qtn", "conditional_recall",
                     "tag", "cell", "rep", "env")

manifest <- fread("out_final_v1/full_grid_manifest.tsv")
say("[1] manifest: %d rows (expect 1400)\n", nrow(manifest))
say("    status breakdown:\n"); print(manifest[, .N, by = status])

say("\n[2] cross-checking every combo's actual output file (not just manifest status)\n")
rows <- list()
all_scores <- list()
for (tag in TAGS) for (cell in CELLS) for (rep in REPS) for (envn in ENVS) {
  combo_id <- sprintf("%s_%s_rep%d_env%d", tag, cell, rep, envn)
  f <- file.path(stage_dir("05_score_truth", combo_id), "truth_scores.rds")
  ok <- file.exists(f)
  methods_ok <- NA; fields_ok <- NA; n_sig_total <- NA_integer_
  if (ok) {
    ts <- readRDS(f)$scores
    methods_ok <- setequal(ts$method, EXPECTED_METHODS)
    fields_ok <- all(EXPECTED_FIELDS %in% names(ts))
    n_sig_total <- sum(ts$n_significant)
    all_scores[[combo_id]] <- ts
  }
  rows[[combo_id]] <- data.table(tag = tag, cell = cell, rep = rep, env = envn,
                                 file_exists = ok, methods_ok = methods_ok, fields_ok = fields_ok, n_sig_total = n_sig_total)
}
report <- rbindlist(rows)
say("    %d/1400 combos have a valid truth_scores.rds\n", sum(report$file_exists))
say("    all methods present in ALL combos: %s\n", all(report$methods_ok, na.rm = TRUE))
say("    all fields present in ALL combos: %s\n", all(report$fields_ok, na.rm = TRUE))
stopifnot("expected exactly 1400 combinations" = nrow(report) == 1400L,
          "every combo must have a valid output file" = all(report$file_exists),
          "every combo must have all 4 methods" = all(report$methods_ok),
          "every combo must have all expected fields" = all(report$fields_ok))

say("\n[3] aggregate sanity check: significant hits and pooled counts by method\n")
scores_all <- rbindlist(all_scores)
agg <- scores_all[, .(n_combos = .N, total_significant = sum(n_significant), total_TP = sum(TP), total_FP = sum(FP),
                      total_FN = sum(FN), pooled_precision = sum(TP) / max(sum(TP) + sum(FP), 1),
                      mean_n_tested = mean(n_tested)), by = method]
print(agg)

say("\n[4] writing reports\n")
fwrite(report, "qc/full_grid_validation_report.tsv", sep = "\t")
fwrite(agg, "qc/full_grid_pooled_summary_preview.tsv", sep = "\t")
say("wrote qc/full_grid_validation_report.tsv, qc/full_grid_pooled_summary_preview.tsv\n")
say("\nGATE 4: PASS -- 1400/1400 combinations complete, all methods and fields present.\n")
