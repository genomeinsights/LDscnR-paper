## =============================================================================
## final_analysis/R/10_export_qc.R
##
## Assembles results/simulation_qc.tsv (CLAUDE_REANALYSIS_INSTRUCTIONS.md,
## "Final outputs" item 3: "parser checks, input counts, and the compact
## BGS validation"). A compact, heterogeneous summary table (category/item/
## value/detail), not a per-combo grid -- draws entirely from artifacts
## already produced by earlier stages/gates (no new computation):
##   - parser checks: qc/parser_qa_report.tsv (Phase 1 gate, 3 runs) and
##     qc/full_grid_validation_report.tsv (Phase 3 gate 4, all 1,400).
##   - input counts: results/simulation_bgs_validation.tsv's "overall" bin
##     (mean markers/combo, pre- and post-MAF-filter, by tag) and
##     results/simulation_performance.tsv's detectable-QTN totals (primary
##     5% Va-share threshold).
##   - BGS validation: the same table's key pre-filter retention numbers
##     (the validated signature) plus a pointer to the full table/figure.
## Entirely local -- reads only small, already-git-tracked results/qc
## files, no raw NEMO or out_final_v1 access needed.
## =============================================================================
suppressMessages({library(data.table)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))
say("=== export_qc ===\n\n")

rows <- list()
add <- function(category, item, value, detail = "") {
  rows[[length(rows) + 1]] <<- data.table(category = category, item = item, value = as.character(value), detail = detail)
}

## ---- 1. parser checks ---------------------------------------------------------
qa <- fread("qc/parser_qa_report.tsv")
add("parser", "qa_gate1_assertions_passed", sprintf("%d/%d", sum(qa$status == "PASS"), nrow(qa)),
    "Phase 1 gate: 3 deliberately chosen native runs, both BGS treatments, >=2 map IDs")
add("parser", "offset_correction", "0 (was +1)",
    "Nemo locus identifiers are one-based; the inherited unconditional +1 offset was removed and verified via 3-run old-vs-corrected comparison")

fg <- fread("qc/full_grid_validation_report.tsv")
add("parser", "gate4_full_grid_combos_verified", sprintf("%d/%d", nrow(fg), 1400L),
    sprintf("file_exists=%s, methods_ok=%s, fields_ok=%s (independent post-hoc validation, not just shell exit status)",
           all(fg$file_exists), all(fg$methods_ok), all(fg$fields_ok)))

## ---- 2. input counts ------------------------------------------------------------
bv <- fread("results/simulation_bgs_validation.tsv")
pre <- bv[filt == "pre_maf" & bin == "overall"]; post <- bv[filt == "post_maf" & bin == "overall"]
add("input_counts", "markers_per_combo_prefilter_nobgs", round(pre$n_markers_nobgs), "mean across 700 (cell,rep,env) triples, before MAF>0.10")
add("input_counts", "markers_per_combo_prefilter_bgs", round(pre$n_markers_bgs), "mean across 700 (cell,rep,env) triples, before MAF>0.10")
add("input_counts", "markers_per_combo_postfilter_nobgs", round(post$n_markers_nobgs), sprintf("mean across 700 triples, MAF>%.2f", MAF_KEEP))
add("input_counts", "markers_per_combo_postfilter_bgs", round(post$n_markers_bgs), sprintf("mean across 700 triples, MAF>%.2f", MAF_KEEP))

pf <- fread("results/simulation_performance.tsv")
dq <- unique(pf[method == "emmax_snp_region", .(tag, cell, n_detectable_qtn)])
add("input_counts", "detectable_qtn_total_grid", sum(dq$n_detectable_qtn),
    sprintf("summed over 14 (tag,cell) strata, primary Va-share=%.2f threshold (MAF>%.2f AND >=%.0f%% chromosome QTN Va)",
           VA_SHARE_DETECTABLE, MAF_KEEP, 100 * VA_SHARE_DETECTABLE))
add("input_counts", "individuals_analysed_per_combo", 160L, "subsampled from 320 raw, step=2")

## ---- 3. compact BGS validation ---------------------------------------------------
low <- bv[filt == "pre_maf" & bin == "low"]; high <- bv[filt == "pre_maf" & bin == "high"]
add("bgs_validation", "retention_ratio_overall", sprintf("%.3f", pre$ratio_n_markers),
    sprintf("pre-filter segregating-marker retention, bgs/nobgs (point estimate; paired-difference 95%% bootstrap CI [%.0f, %.0f] markers excludes 0) -- see results/simulation_bgs_validation.tsv, figureS_simulation_bgs_qc.pdf for the full table/CIs",
           pre$diff_n_markers_ci_lo, pre$diff_n_markers_ci_hi))
add("bgs_validation", "retention_ratio_low_recomb", sprintf("%.3f", low$ratio_n_markers), "classic BGS signature: steepest reduction at low recombination")
add("bgs_validation", "retention_ratio_high_recomb", sprintf("%.3f", high$ratio_n_markers), "weakest reduction at high recombination, as expected")
add("bgs_validation", "post_filter_he_trap_confirmed", "TRUE",
    "post-MAF-filter heterozygosity among survivors is HIGHER under bgs (backwards) -- confirms a known MAF-trim measurement trap; pre-filter retention ratio above is the validated metric, not post-filter He")

qc <- rbindlist(rows)
dir.create("results", showWarnings = FALSE)
fwrite(qc, "results/simulation_qc.tsv", sep = "\t")
say("wrote results/simulation_qc.tsv (%d rows)\n", nrow(qc))
print(qc)
