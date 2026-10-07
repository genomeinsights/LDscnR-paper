## =============================================================================
## final_analysis/R/22_floor_stability_regions.R
##
## EXPLORATORY, additive analysis (PK, 2026-09-27) -- part 2 of 3 (see
## R/21's header for the full question/scope). Does NOT touch
## R/12_floor_decomposition.R or its outputs.
##
## For each of the 600 c=1 combos, at each of its (up to 3, often fewer once
## collapsed floors are deduplicated -- see R/21) selected floors, rerun BH +
## Stage-2 assembly + truth scoring for EMMAX Simes and EMMAX consensus
## (LFMM Simes too, if that combo has an LFMM stage -- secondary/optional per
## PK's own framing, reported separately in R/23, not folded into the primary
## EMMAX tables). Reuses, never recomputes:
##   - raw marker p-values (03_emmax.R's em$marker$p / 04_lfmm.R's lf$marker$p)
##   - the fixed relationship matrix (GRM, from 02_build_ld_units.R)
##   - Stage-1 clusters (stage1) and LD-decay (LD_decay)
##   - a SINGLE floor=1 per-unit Simes/consensus p-value (simes_p1/con_p1),
##     subset and re-BH'd at each floor -- never refit per floor.
## This is R/12_floor_decomposition.R's own arms C/D machinery verbatim (same
## core_snp-keyed stable-identity discipline, same reasons -- see that
## file's header on why unit_id is NOT stable across floors), just evaluated
## at THIS analysis's own floor set instead of FLOOR_GRID, and returning full
## per-region detail (not just summary counts), since R/23 needs individual
## regions to match across floors.
##
## Truth prerequisites (flag_true_qtns(), score_thresholds(), qtn_ld_table())
## are identical to 05_score_truth.R/11_stage2_region_details.R/R/12's own
## copies -- computed once per combo, floor-independent, never recomputed
## per floor.
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(parallel)})
MODULE_ROOT <- path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis")
source(file.path(MODULE_ROOT, "R", "00_config.R"))
source(file.path(MODULE_ROOT, "R", "helpers_stage2_truth.R"))
source(file.path(MODULE_ROOT, "R", "helpers_floor_stability.R"))
say("=== 22_floor_stability_regions ===\n\n")

## ---- one combo -> region-detail rows at every DISTINCT floor it needs -------
floor_stability_regions_one_combo <- function(tag, cell, rep, env, floors_needed) {
  combo_id <- sprintf("%s_%s_rep%d_env%d", tag, cell, rep, env)
  em <- readRDS(file.path(stage_dir("03_emmax", combo_id), "emmax.rds"))
  lfmm_file <- file.path(stage_dir("04_lfmm", combo_id), "lfmm.rds")
  have_lfmm <- file.exists(lfmm_file)
  lf <- if (have_lfmm) readRDS(lfmm_file) else NULL
  b <- readRDS(file.path(stage_dir("02_build_ld_units", combo_id), "ld_units.rds"))
  GTs <- b$GTs; map <- copy(b$map); stage1 <- b$stage1; GRM <- b$GRM
  y <- b$env$env

  ## ---- truth prerequisites (verbatim from R/12/05_score_truth.R/R/11) -----
  p <- colSums(GTs) / nrow(GTs) / 2
  map[, p_freq := p[match(marker, colnames(GTs))]]
  map[, Va := 2 * p_freq * (1 - p_freq) * allelic_values^2]
  map <- flag_true_qtns(map, va_col = "Va", maf_col = "MAF", maf_min = MAF_KEEP, p_va_min = VA_SHARE_DETECTABLE)
  detectable_qtn <- map[true_pos_QTN == TRUE, marker]
  n_dq <- length(detectable_qtn)

  thr <- score_thresholds(b$LD_decay$decay_sum, rho_r2 = TRUTH_RHO_R2, rho_d = TRUTH_RHO_D, dmax_cap = TRUTH_DMAX_CAP)
  qtn_lut <- if (sum(map$type == "QTN") == 0L) {
    data.table(marker = character(), qtn_marker = character(), r2 = numeric(), dist_bp = numeric())
  } else {
    qtn_ld_table(GTs, map, candidate_markers = map$marker, cores = 1)
  }
  qtn_lut_match <- qtn_lut[r2 > thr$r2min & dist_bp < thr$dmax & qtn_marker %in% detectable_qtn]

  ## ---- floor=1 base + floor=1 raw per-unit p (verbatim from R/12) --------
  units1 <- LDscnR:::.ld_outlier_units(stage1, map, 1L)
  pm_obs <- em$marker$p; names(pm_obs) <- em$marker$marker
  simes_p1 <- stats::setNames(vapply(units1$members, function(mm) LDscnR:::.simes(pm_obs[mm]), numeric(1)),
                              units1$core_snp)
  um1 <- ld_unit_matrix(GTs, stage1, map, size_floor = 1L, repr = "consensus_dosage")
  stopifnot("ld_unit_matrix(floor=1) unit_id must match units1$unit_id, in row order" =
              identical(colnames(um1), as.character(units1$unit_id)))
  Pu1 <- emmax_setup(um1, GRM)
  con_p1 <- stats::setNames(emmax_fast(Pu1, y), units1$core_snp)
  if (have_lfmm) {
    lp_obs <- lf$marker$p; names(lp_obs) <- lf$marker$marker
    lsimes_p1 <- stats::setNames(vapply(units1$members, function(mm) LDscnR:::.simes(lp_obs[mm]), numeric(1)),
                                 units1$core_snp)
  }

  ## ---- one method's region set at one floor --------------------------------
  .region_set <- function(method, floor, p_vec) {
    q <- stats::p.adjust(p_vec, "BH")
    sig_core <- names(p_vec)[!is.na(q) & q <= ALPHA]
    cl_sub <- stage2_seed_from_units(stage1, map, sig_core)
    detail <- assemble_stage2(stage1, map, GTs, b$LD_decay, cl_sub,
                              REGION_ASSEMBLY$score_threshold, REGION_ASSEMBLY$distance_threshold)
    detail <- score_stage2_regions(detail, qtn_lut_match)
    if (nrow(detail)) detail[, `:=`(tag = tag, cell = cell, rep = rep, env = env, combo_id = combo_id,
                                    method = method, floor = floor,
                                    n_tests = length(p_vec), n_significant = length(sig_core))]
    detail
  }

  rows <- list()
  for (floor in floors_needed) {
    eligible1 <- units1[n_markers >= floor]
    rows[[length(rows) + 1]] <- .region_set("emmax_simes", floor, simes_p1[eligible1$core_snp])
    rows[[length(rows) + 1]] <- .region_set("emmax_consensus", floor, con_p1[eligible1$core_snp])
    if (have_lfmm) rows[[length(rows) + 1]] <- .region_set("lfmm_simes", floor, lsimes_p1[eligible1$core_snp])
  }
  out <- rbindlist(Filter(function(d) !is.null(d) && ncol(d) > 0, rows), fill = TRUE)
  if (nrow(out)) out[, n_detectable_qtn := n_dq]
  ## n_detectable_qtn is combo-level, not region-level. A combo with ZERO
  ## reported regions at every floor/method would otherwise silently lose
  ## this value entirely (it only ever gets attached to region ROWS above,
  ## and a zero-region combo produces none) -- but R/23 needs it as a fixed
  ## recall denominator for EVERY combo, including zero-call ones. Returned
  ## separately here rather than folded into `out` so it survives regardless
  ## of how many regions this combo produced.
  list(regions = out, combo_meta = data.table(tag = tag, cell = cell, rep = rep, env = env,
                                              combo_id = combo_id, n_detectable_qtn = n_dq,
                                              have_lfmm = have_lfmm))
}

if (sys.nframe() == 0L) {
  sel_wide <- fread(file.path(MODULE_ROOT, "results", fs_name("floor_selection_wide.tsv")))
  say("[1] loaded floor selection for %d combos\n", nrow(sel_wide))

  say("[2] rerunning association + Stage-2 assembly at each combo's distinct floor(s)\n")
  t0 <- Sys.time()
  res <- mclapply(seq_len(nrow(sel_wide)), function(i) {
    s <- sel_wide[i]
    floors_needed <- sort(unique(c(s$floor_995, s$floor_997, s$floor_999)))
    tryCatch(list(ok = TRUE, data = floor_stability_regions_one_combo(s$tag, s$cell, s$rep, s$env, floors_needed), error = NA_character_),
             error = function(e) list(ok = FALSE, data = NULL,
                                      error = sprintf("%s: %s", s$combo_id, conditionMessage(e))))
  }, mc.cores = 12)
  ok_flags <- vapply(res, `[[`, logical(1), "ok")
  if (!all(ok_flags)) {
    errs <- vapply(res[!ok_flags], `[[`, character(1), "error")
    message(paste(errs, collapse = "\n"))
    stop(sprintf("%d/%d combos failed region generation -- see messages above; refusing to pool a partial result set",
                 sum(!ok_flags), nrow(sel_wide)))
  }
  combo_results <- lapply(res, `[[`, "data")
  dt <- rbindlist(lapply(combo_results, `[[`, "regions"), fill = TRUE)
  combo_meta <- rbindlist(lapply(combo_results, `[[`, "combo_meta"))
  say("[3] %s region rows from %d combos in %.1f min\n", format(nrow(dt), big.mark = ","), nrow(sel_wide),
      as.numeric(difftime(Sys.time(), t0, units = "mins")))
  stopifnot("combo_meta must have exactly one row per combo (n_detectable_qtn denominator source)" =
              nrow(combo_meta) == nrow(sel_wide))

  saveRDS(dt, file.path(MODULE_ROOT, "results", fs_name("region_details.rds")), compress = "xz")
  compact <- dt[, .(tag, cell, rep, env, combo_id, method, floor, region_id, Chr, from, to,
                    n_markers, n_units, TP, n_tests, n_significant, n_detectable_qtn)]
  fwrite(compact, file.path(MODULE_ROOT, "results", fs_name("region_details.tsv")), sep = "\t")
  fwrite(combo_meta, file.path(MODULE_ROOT, "results", fs_name("combo_meta.tsv")), sep = "\t")
  say("[3b] wrote results/floor_stability_combo_meta.tsv (%d combos, n_detectable_qtn for every combo including zero-region ones)\n",
      nrow(combo_meta))
  say("[4] wrote results/floor_stability_region_details.rds (full, list-cols) + .tsv (compact)\n")

  ## ---- completeness diagnostic: this is a LONG table -- a (combo, method,
  ## floor) cell with zero reported regions simply contributes no rows here,
  ## it is not "missing". R/23 must compute its own zero-call denominator
  ## from sel_wide x method (the full intended grid), never from this
  ## table's own row count, to avoid silently dropping zero-region cells.
  say("[5] %d distinct (combo, method, floor) cells produced >=1 region (out of an intended grid R/23 will compute from floor selection x method)\n",
      uniqueN(dt[, .(combo_id, method, floor)]))
  say("    methods present: %s\n", paste(sort(unique(dt$method)), collapse = ", "))

  write_receipt(fs_stage("22_floor_stability_regions"), inputs = paste0("results/", fs_name("floor_selection_wide.tsv")),
                params = list(alpha = ALPHA, maf_keep = MAF_KEEP, va_share_detectable = VA_SHARE_DETECTABLE,
                              truth_rho_r2 = TRUTH_RHO_R2, truth_rho_d = TRUTH_RHO_D),
                outputs = c(paste0("results/", fs_name("region_details.rds")), paste0("results/", fs_name("region_details.tsv")),
                           paste0("results/", fs_name("combo_meta.tsv"))))
  cat("REGION_GENERATION_DONE\n")
}
