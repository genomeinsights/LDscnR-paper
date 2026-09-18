## =============================================================================
## final_analysis/R/11_stage2_region_details.R
##
## Full detailed record of every Stage-2 reported region at the CANONICAL
## floor (SIZE_FLOOR), for the primary region methods (emmax_snp_region,
## emmax_simes_region, emmax_consensus_region, lfmm_snp_region,
## lfmm_simes_region when LFMM is present). Two purposes at once:
##
##  1. This IS the validation gate for R/helpers_stage2_truth.R's refactor:
##     pooling this file's TP/FP by (tag,cell,method) must reproduce
##     results/simulation_performance*.tsv exactly (verified separately, see
##     ADDITIONAL_ANALYSES_AUDIT.md -- the full 1,400-combo rescore via
##     05_score_truth.R already confirmed this at the pooled level; this
##     script exercises the SAME assemble_stage2()/score_stage2_regions()
##     functions, at the SAME floor, so it is the same code path again, not
##     an independent implementation that could drift).
##  2. It is Analysis 2's (Stage-2 size vs. TP/FP) primary input -- the
##     per-region detail that 05_score_truth.R deliberately does NOT retain
##     (only pooled per-combo counts, to keep truth_scores.rds small).
##
## Per CLAUDE_REANALYSIS_INSTRUCTIONS.md / the additional-analyses brief:
## constituent markers/units come from assemble_stage2()'s actual
## ld_prune_and_eMLG() output, never a [Chr,from,to] bounds sweep; every
## region here is a real, individually reported Stage-2 candidate.
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(parallel)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "helpers_stage2_truth.R"))
say("=== 11_stage2_region_details ===\n\n")

## One combo -> a data.table of detailed region rows across every primary
## region method present for that combo, at SIZE_FLOOR. Mirrors
## 05_score_truth.R's truth-table preamble (Va, detectable QTN, qtn_lut_match,
## units_base) -- duplicated rather than shared, matching this project's
## established convention for the truth-table preamble specifically (see
## 09_truth_sensitivity.R's header comment); only the assembly/scoring core
## (the part that must not drift) is the shared helper.
region_details_one_combo <- function(tag, cell, rep, env) {
  combo_id <- sprintf("%s_%s_rep%d_env%d", tag, cell, rep, env)
  em <- readRDS(file.path(stage_dir("03_emmax", combo_id), "emmax.rds"))
  lfmm_file <- file.path(stage_dir("04_lfmm", combo_id), "lfmm.rds")
  have_lfmm <- file.exists(lfmm_file)
  lf <- if (have_lfmm) readRDS(lfmm_file) else NULL
  b <- readRDS(file.path(stage_dir("02_build_ld_units", combo_id), "ld_units.rds"))
  GTs <- b$GTs; map <- copy(b$map); stage1 <- b$stage1

  p <- colSums(GTs) / nrow(GTs) / 2
  map[, p_freq := p[match(marker, colnames(GTs))]]
  map[, Va := 2 * p_freq * (1 - p_freq) * allelic_values^2]
  map <- flag_true_qtns(map, va_col = "Va", maf_col = "MAF", maf_min = MAF_KEEP, p_va_min = VA_SHARE_DETECTABLE)
  detectable_qtn <- map[true_pos_QTN == TRUE, marker]

  thr <- score_thresholds(b$LD_decay$decay_sum, rho_r2 = TRUTH_RHO_R2, rho_d = TRUTH_RHO_D, dmax_cap = TRUTH_DMAX_CAP)
  qtn_lut <- if (sum(map$type == "QTN") == 0L) {
    data.table(marker = character(), qtn_marker = character(), r2 = numeric(), dist_bp = numeric())
  } else {
    qtn_ld_table(GTs, map, candidate_markers = map$marker, cores = 1)
  }
  qtn_lut_match <- qtn_lut[r2 > thr$r2min & dist_bp < thr$dmax & qtn_marker %in% detectable_qtn]

  units_base <- LDscnR:::.ld_outlier_units(stage1, map, SIZE_FLOOR)

  .detail_for <- function(method, cl_sub) {
    d <- assemble_stage2(stage1, map, GTs, b$LD_decay, cl_sub,
                         REGION_ASSEMBLY$score_threshold, REGION_ASSEMBLY$distance_threshold)
    d <- score_stage2_regions(d, qtn_lut_match)
    if (!nrow(d)) return(NULL)
    d[, `:=`(tag = tag, cell = cell, rep = rep, env = env, method = method, size_floor = SIZE_FLOOR)]
    d[]
  }

  rows <- list()
  rows$emmax_snp_region <- .detail_for("emmax_snp_region",
    stage2_seed_from_markers(stage1, map, em$marker$marker[em$marker$significant_snp]))

  sig_sim <- em$emmax_simes[em$emmax_simes$significant == TRUE]
  sig_sim_core <- units_base$core_snp[match(as.character(sig_sim$unit_id), as.character(units_base$unit_id))]
  rows$emmax_simes_region <- .detail_for("emmax_simes_region", stage2_seed_from_units(stage1, map, sig_sim_core))

  sig_con <- em$emmax_consensus[em$emmax_consensus$significant == TRUE]
  sig_con_core <- units_base$core_snp[match(as.character(sig_con$unit_id), as.character(units_base$unit_id))]
  rows$emmax_consensus_region <- .detail_for("emmax_consensus_region", stage2_seed_from_units(stage1, map, sig_con_core))

  if (have_lfmm) {
    rows$lfmm_snp_region <- .detail_for("lfmm_snp_region",
      stage2_seed_from_markers(stage1, map, lf$marker$marker[lf$marker$significant_snp]))
    sig_lsim <- lf$lfmm_simes[lf$lfmm_simes$significant == TRUE]
    sig_lsim_core <- units_base$core_snp[match(as.character(sig_lsim$unit_id), as.character(units_base$unit_id))]
    rows$lfmm_simes_region <- .detail_for("lfmm_simes_region", stage2_seed_from_units(stage1, map, sig_lsim_core))
  }

  out <- rbindlist(Filter(Negate(is.null), rows), fill = TRUE)
  if (!nrow(out)) return(NULL)
  out[]
}

if (sys.nframe() == 0L) {
  combos <- CJ(tag = TAGS_ALL, cell = CELLS_ALL, rep = REPS_ALL, env = ENVS_ALL)
  say("[1] region detail for %d combos x up to 5 region methods (floor=%d)\n", nrow(combos), SIZE_FLOOR)
  t0 <- Sys.time()
  ## [!] {ok, data, error} worker result (PK, 2026-09-18 second review):
  ## region_details_one_combo() itself legitimately returns NULL when a
  ## combo has zero reported regions across every method -- a real, valid
  ## outcome, not a failure. Plain tryCatch(..., error=function(e) NULL)
  ## made that indistinguishable from an actual worker error (a crash, an
  ## OOM kill under the higher mc.cores, a corrupted read), and either one
  ## would just be silently dropped by Filter(Negate(is.null), ...) with the
  ## script still reporting success on the remainder. Wrapping every outcome
  ## in an explicit {ok, data} tag, and hard-stopping if any ok is FALSE,
  ## makes a genuine failure impossible to miss while still allowing the
  ## legitimate empty-combo case through.
  res <- mclapply(seq_len(nrow(combos)), function(i) {
    tryCatch(list(ok = TRUE, data = region_details_one_combo(combos$tag[i], combos$cell[i], combos$rep[i], combos$env[i])),
             error = function(e) {
               msg <- sprintf("[11_stage2_region_details] %s_%s_rep%d_env%d: %s",
                              combos$tag[i], combos$cell[i], combos$rep[i], combos$env[i], conditionMessage(e))
               message(msg)
               list(ok = FALSE, data = NULL, error = msg)
             })
  }, mc.cores = 12)   ## bumped from 7 (PK, 2026-09-18): mini has 14 physical cores, 2 left for the system
  ok_flags <- vapply(res, `[[`, logical(1), "ok")
  if (!all(ok_flags))
    stop(sprintf("%d/%d combos failed -- see [11_stage2_region_details] messages above; refusing to pool a partial result set", sum(!ok_flags), nrow(combos)))
  dt <- rbindlist(Filter(Negate(is.null), lapply(res, `[[`, "data")), fill = TRUE)
  say("[2] %s region rows from %d combos in %.1f min\n", format(nrow(dt), big.mark = ","), nrow(combos),
      as.numeric(difftime(Sys.time(), t0, units = "mins")))

  dt[, region_uid := sprintf("%s_%s_rep%d_env%d_%s_r%d", tag, cell, rep, env, method, region_id)]
  setcolorder(dt, c("region_uid", "tag", "cell", "rep", "env", "method", "size_floor", "region_id",
                    "Chr", "from", "to", "n_markers", "n_units", "TP"))

  RDS_OUT <- "results/simulation_stage2_region_details.rds"
  saveRDS(dt, RDS_OUT, compress = "xz")
  say("[3] wrote %s (%.1f MB)\n", RDS_OUT, file.size(RDS_OUT) / 1e6)

  ## Compact TSV: list-cols (member_markers/member_units/linked_qtn) joined
  ## with ";" for human inspection; the full list-columns survive intact in
  ## the RDS above -- this file is a readable summary, not the sole record.
  tsv <- copy(dt)
  tsv[, member_markers := vapply(member_markers, paste, character(1), collapse = ";")]
  tsv[, member_units := vapply(member_units, paste, character(1), collapse = ";")]
  tsv[, linked_qtn := vapply(linked_qtn, paste, character(1), collapse = ";")]
  TSV_OUT <- "results/simulation_stage2_region_details.tsv"
  fwrite(tsv, TSV_OUT, sep = "\t")
  say("[4] wrote %s\n", TSV_OUT)

  write_receipt("11_stage2_region_details",
                inputs = character(),  # depends on the whole 03/04/02 grid; receipted at combo level upstream
                params = list(size_floor = SIZE_FLOOR, region_assembly = REGION_ASSEMBLY,
                              truth = list(va_share_detectable = VA_SHARE_DETECTABLE, maf_keep = MAF_KEEP,
                                          truth_rho_r2 = TRUTH_RHO_R2, truth_rho_d = TRUTH_RHO_D)),
                outputs = c(RDS_OUT, TSV_OUT))

  say("\n[5] rows by method:\n")
  print(dt[, .N, by = method])
}
