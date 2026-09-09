## =============================================================================
## final_analysis/R/01_parse_nemo_qa.R
##
## Phase 3 gate 1 / Phase 1's own validation requirement: "On at least three
## deliberately chosen native runs, compare the old and corrected join. Print
## the mappings for the final few QTN indices so that the correction is
## human-verifiable."
##
## Three runs chosen to span both BGS treatments and >= 2 map IDs (both
## genome chromosomes are present within any single run, by design). A
## preliminary survey (6 runs across cells/tags/reps) found that each raw run
## saves only a small, VARIABLE subset of the 101 reference QTN (2-14 of 101
## in that sample) -- locus 101 itself (the single Chr2 QTN, the sharpest
## boundary case for this fix) appeared in only 1 of those 6. Run A below is
## deliberately the one confirmed to include it, so the QTN-index printout
## actually exercises the Chr2 case rather than relying on the reference-map-
## only check to cover it.
##   A: nobgs / V0.5_c2 / rep1 / env1   (confirmed: includes quant.100 AND quant.101)
##   B: bgs   / V0.5_c2 / rep1 / env1   (same cell+rep as A, other tag)
##   C: bgs   / V1_c1   / rep2 / env1   (different cell AND rep)
##
## For each run: parse under BOTH offset=0 (corrected, production) and
## offset=1 (old, reproduced here ONLY for this comparison -- never used to
## produce a real bundle). Outputs:
##   qc/parser_qa_report.tsv            -- every assertion, corrected parse only
##   qc/three_run_offset_comparison.tsv -- old vs corrected, per run
##   qc/qtn_index_mapping_<run>.tsv     -- full QTN-by-QTN old vs corrected mapping
## =============================================================================
suppressMessages({library(data.table)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "01_parse_nemo.R"))

RUNS <- list(
  A = list(tag = "nobgs", cell = "V0.5_c2", rep = 1L, env = 1L),
  B = list(tag = "bgs",   cell = "V0.5_c2", rep = 1L, env = 1L),
  C = list(tag = "bgs",   cell = "V1_c1",   rep = 2L, env = 1L)
)

all_qa <- list()
all_compare <- list()

for (run_name in names(RUNS)) {
  r <- RUNS[[run_name]]
  say("\n============================================================\n")
  say("=== run %s: %s/%s/rep%d/env%d ===\n", run_name, r$tag, r$cell, r$rep, r$env)
  say("============================================================\n")

  say("\n--- corrected (offset=0) ---\n")
  res0 <- parse_nemo_run(r$tag, r$cell, r$rep, r$env, offset = 0L)
  res0$qa[, run := run_name]
  all_qa[[run_name]] <- res0$qa
  say("  %d assertions PASSED\n", nrow(res0$qa))
  say("  markers pre-filter: %d (QTN=%d) ; post-MAF: %d (QTN=%d, %d QTN lost to MAF)\n",
      sum(res0$counts_by_chr_type$prefilter$N), sum(res0$map$true_QTN) + res0$n_qtn_lost_to_maf,
      nrow(res0$map), sum(res0$map$true_QTN), res0$n_qtn_lost_to_maf)
  print(res0$counts_by_chr_type$prefilter)
  print(res0$counts_by_chr_type$postfilter)

  say("\n--- old (offset=1), QA comparison ONLY -- not a production parse ---\n")
  ## Deliberately NOT calling parse_nemo_run(offset=1) here: for run A
  ## specifically, raw idx 101 + offset 1 = 102, out of range for a 101-row
  ## QTN block -- parse_nemo_run's own "idx in [1, n_loci_of_type]" assertion
  ## correctly stops on that (it's the bug, demonstrated), which would abort
  ## this whole comparison script rather than just this one number. Instead,
  ## reproduce the OLD join directly via the same .join_nemo_to_refmap()
  ## helper the production parser uses, at offset=1, and compute the same
  ## local_adapt_r2 diagnostic from it without the full assertion suite.
  bf <- file.path(PATHS$raw_nemo, sprintf("adapt_%s_chr%d_%s_env%d", r$tag, r$rep, r$cell, r$env), "GENO")
  files <- list.files(bf, full.names = TRUE)
  map_nemo <- fread(files[grepl("\\.map$", files)])
  nemo_map <- data.table(marker = map_nemo$trait.locus,
                         do.call(rbind, strsplit(map_nemo$trait.locus, ".", fixed = TRUE)))
  setnames(nemo_map, c("V1", "V2"), c("type", "idx_raw"))
  nemo_map[, idx_raw := as.numeric(idx_raw)]
  refmap0 <- readRDS(file.path(PATHS$raw_recmap_dir, sprintf("rec_map%d.rds", r$rep)))
  n_qtn_total <- refmap0[type == "QTN", .N]
  quant_idx_here <- nemo_map[type == "quant", idx_raw]
  n_out_of_range_old <- sum((quant_idx_here + 1) > n_qtn_total)
  say("  old offset (+1): %d of %d QTN loci in this run go OUT OF RANGE (idx+1 > %d total QTN)\n",
      n_out_of_range_old, length(quant_idx_here), n_qtn_total)

  j0 <- .join_nemo_to_refmap(nemo_map, refmap0, 0L)$refmap
  j1 <- .join_nemo_to_refmap(nemo_map, refmap0, 1L)$refmap

  bv_r2_from_join <- function(joined_refmap, GTs_raw_local, sample_info_local, env_local) {
    mm <- joined_refmap[!is.na(nemo_marker) & type == "QTN"]
    if (!nrow(mm)) return(list(r2 = NA_real_, n_qtn = 0L, n_bv_rows = 0L))
    bv_rows <- mm[!is.na(allelic_values)]
    if (!nrow(bv_rows)) return(list(r2 = NA_real_, n_qtn = nrow(mm), n_bv_rows = 0L))
    GT <- GTs_raw_local[, bv_rows$nemo_marker, drop = FALSE]
    bv <- as.numeric(GT %*% bv_rows$allelic_values)
    if (stats::sd(bv) == 0) return(list(r2 = NA_real_, n_qtn = nrow(mm), n_bv_rows = nrow(bv_rows)))
    list(r2 = stats::cor(bv, env_local)^2, n_qtn = nrow(mm), n_bv_rows = nrow(bv_rows))
  }
  GTs_raw_full <- fread(files[grepl("snp_geno", files)])
  sample_info_full <- GTs_raw_full[, .(pop, ID)]
  GTs_raw_mat <- as.matrix(GTs_raw_full[, 6:ncol(GTs_raw_full), with = FALSE])
  colnames(GTs_raw_mat) <- nemo_map$marker
  env_txt_full <- file.path(PATHS$raw_env_dir, sprintf("env_%d.txt", r$env))
  ev <- scan(env_txt_full, what = character(), quiet = TRUE)
  ev <- strsplit(ev, "}{", fixed = TRUE)[[1]]
  ev <- as.numeric(gsub("}}", "", gsub("{{", "", ev, fixed = TRUE), fixed = TRUE))
  eg <- data.table(expand.grid(x = 1:GRID_SIDE, y = GRID_SIDE:1)); eg[, pop := .I]; eg[, env := ev]
  env_full <- eg[match(sample_info_full$pop, eg$pop), env]

  bv0 <- bv_r2_from_join(j0, GTs_raw_mat, sample_info_full, env_full)
  bv1 <- bv_r2_from_join(j1, GTs_raw_mat, sample_info_full, env_full)

  cmp <- data.table(run = run_name, tag = r$tag, cell = r$cell, rep = r$rep, env = r$env,
                    n_qtn_in_run = length(quant_idx_here), n_out_of_range_old = n_out_of_range_old,
                    n_qtn_corrected = bv0$n_qtn, n_bv_rows_corrected = bv0$n_bv_rows, local_adapt_r2_corrected = bv0$r2,
                    n_qtn_old = bv1$n_qtn, n_bv_rows_old = bv1$n_bv_rows, local_adapt_r2_old = bv1$r2)
  all_compare[[run_name]] <- cmp
  say("\n  local_adapt_r2 (full genome-wide QTN, no MAF filter): corrected=%.4f (n_bv=%d/%d QTN) vs old=%.4f (n_bv=%d/%d QTN)\n",
      bv0$r2, bv0$n_bv_rows, bv0$n_qtn, bv1$r2, bv1$n_bv_rows, bv1$n_qtn)

  ## ---- full QTN-by-QTN mapping table, both offsets, for human verification ----
  quant_markers <- nemo_map[type == "quant", marker]
  qtn_tab <- data.table(
    quant_marker = quant_markers,
    corrected_Chr = j0[match(quant_markers, nemo_marker), Chr],
    corrected_Pos = j0[match(quant_markers, nemo_marker), bp],
    old_Chr = j1[match(quant_markers, nemo_marker), Chr],
    old_Pos = j1[match(quant_markers, nemo_marker), bp])
  qtn_tab[, quant_num := as.integer(sub("quant\\.", "", quant_marker))]
  setorder(qtn_tab, quant_num)
  say("\n  --- QTN index mapping, last 6 (of %d) ---\n", nrow(qtn_tab))
  print(tail(qtn_tab, 6))
  fwrite(qtn_tab, file.path(PATHS$qc, sprintf("qtn_index_mapping_%s.tsv", run_name)), sep = "\t")
}

qa_report <- rbindlist(all_qa)
compare_report <- rbindlist(all_compare)
fwrite(qa_report, file.path(PATHS$qc, "parser_qa_report.tsv"), sep = "\t")
fwrite(compare_report, file.path(PATHS$qc, "three_run_offset_comparison.tsv"), sep = "\t")

say("\n\n============================================================\n")
say("=== FINAL SUMMARY ===\n")
say("============================================================\n")
say("Assertions: %d PASS, 0 FAIL (any failure would have stopped this script) across %d runs, %d checks/run\n",
    nrow(qa_report), length(RUNS), nrow(qa_report) / length(RUNS))
say("\nWrote:\n  %s\n  %s\n  %s\n",
    file.path(PATHS$qc, "parser_qa_report.tsv"),
    file.path(PATHS$qc, "three_run_offset_comparison.tsv"),
    paste0(file.path(PATHS$qc, "qtn_index_mapping_<A|B|C>.tsv")))
print(compare_report)
