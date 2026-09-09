## final_analysis/R_figures/manhattan_example_data.R
##
## Shared scaffold for the two illustrative example-combo Manhattan figures
## (figureS_simulation_manhattan_example.R, figureS_simulation_manhattan_tpfp.R).
## Not a standalone figure script -- sourced by both, after 00_config.R.
##
## Concatenates ALL 10 REPS (recombination-map replicates) of ONE (tag,
## cell, env) into a single 20-"chromosome" illustrative genome: rep r's
## own Chr1 (the QTN-bearing map) becomes chromosome 2r-1, its Chr2 (the
## near-neutral map -- 0 or 1 QTN by construction, see
## CLAUDE_REANALYSIS_INSTRUCTIONS.md) becomes chromosome 2r, so odd = real,
## even = near-neutral, alternating across the 10 reps. Each rep is an
## independent simulated genome (own reference map/burn-in) -- LD, truth
## linkage, etc. are always computed WITHIN a rep, never across reps.
##
## env=3 (of nobgs/V0.5_c1's 10 environments) was picked as the illustrative
## example because a one-time scan found it is the ONLY environment where
## every one of the 10 reps has >=1 significant EMMAX marker (972 total
## significant markers, 77 significant emmax_simes units, vs e.g. env=9's
## 4/10 reps and 60 markers). This is a representative example for a
## figure, not a systematic result -- the actual pooled performance numbers
## are in results/simulation_performance.tsv (R/06_summarise.R).
suppressMessages({library(data.table); library(LDscnR)})

EXAMPLE_TAG <- "nobgs"; EXAMPLE_CELL <- "V0.5_c1"; EXAMPLE_ENV <- 3L; EXAMPLE_REPS <- 1:10

## One rep's ld_units.rds/emmax.rds/lfmm.rds, all already produced by the
## primary grid -- this only reads existing stage output, never rebuilds it.
.load_example_rep <- function(tag, cell, r, envn) {
  combo_id <- sprintf("%s_%s_rep%d_env%d", tag, cell, r, envn)
  list(combo_id = combo_id,
      lu = readRDS(file.path(stage_dir("02_build_ld_units", combo_id), "ld_units.rds")),
      em = readRDS(file.path(stage_dir("03_emmax", combo_id), "emmax.rds")),
      lf = readRDS(file.path(stage_dir("04_lfmm", combo_id), "lfmm.rds")))
}

## compute_stage2=FALSE skips the (fairly slow) Stage-2 region-assembly
## re-call for figures that don't need it (e.g. the TP/FP figure, which
## only needs the already-saved Stage-1 emmax_simes/lfmm_simes tables).
build_example_genome <- function(tag = EXAMPLE_TAG, cell = EXAMPLE_CELL, envn = EXAMPLE_ENV,
                                 reps = EXAMPLE_REPS, compute_stage2 = TRUE) {
  rep_data <- list(); mk_list <- list()
  for (r in reps) {
    rd <- .load_example_rep(tag, cell, r, envn)
    map <- as.data.table(rd$lu$map)
    qtns <- map[type == "QTN"]

    test_em <- test_lf <- NULL
    if (compute_stage2) {
      test_em <- ld_outlier_test(rd$lu$stage1, map, rd$em$marker$p, statistic = "simes", size_floor = SIZE_FLOOR,
                                 alpha = ALPHA, assembly = "stage2_discovered", GTs = rd$lu$GTs,
                                 LD_decay = rd$lu$LD_decay, score_threshold = REGION_ASSEMBLY$score_threshold,
                                 distance_threshold = REGION_ASSEMBLY$distance_threshold)
      test_lf <- ld_outlier_test(rd$lu$stage1, map, rd$lf$marker$p, statistic = "simes", size_floor = SIZE_FLOOR,
                                 alpha = ALPHA, assembly = "stage2_discovered", GTs = rd$lu$GTs,
                                 LD_decay = rd$lu$LD_decay, score_threshold = REGION_ASSEMBLY$score_threshold,
                                 distance_threshold = REGION_ASSEMBLY$distance_threshold)
    }

    em_mk <- as.data.table(rd$em$marker)[, .(Chr, Pos, marker, q, engine = "EMMAX")]
    lf_mk <- as.data.table(rd$lf$marker)[, .(Chr, Pos, marker, q, engine = "LFMM")]
    mk <- rbindlist(list(em_mk, lf_mk))
    mk[, `:=`(rep = r, is_qtn = marker %in% qtns$marker)]
    mk[, global_chr := 2L * r - 1L + as.integer(Chr == "Chr2")]

    rep_data[[as.character(r)]] <- c(rd, list(map = map, qtns = qtns, test_em = test_em, test_lf = test_lf))
    mk_list[[as.character(r)]] <- mk
  }

  dt <- rbindlist(mk_list)
  dt[, engine := factor(engine, levels = c("EMMAX", "LFMM"))]

  ## genome-cumulative x position across all 20 chromosomes, shared by both
  ## engine rows (facet_grid) since the marker set/positions are identical
  chr_max <- dt[engine == "EMMAX", .(max_pos = max(Pos)), by = global_chr][order(global_chr)]
  chr_max[, offset := c(0, cumsum(max_pos)[-.N]) + c(0, cumsum(rep(2e5, .N - 1)))]
  dt <- merge(dt, chr_max[, .(global_chr, offset)], by = "global_chr")
  dt[, pos_cum := Pos + offset]
  chr_mid <- dt[engine == "EMMAX", .(mid = mean(range(pos_cum))), by = global_chr][order(global_chr)]
  chr_mid[, neutral := global_chr %% 2 == 0]
  band_rects <- dt[engine == "EMMAX", .(xmin = min(pos_cum), xmax = max(pos_cum)), by = global_chr]
  band_rects <- merge(band_rects, chr_mid[, .(global_chr, neutral)], by = "global_chr")
  band_rects <- band_rects[neutral == TRUE]

  list(rep_data = rep_data, dt = dt, chr_mid = chr_mid, band_rects = band_rects,
       tag = tag, cell = cell, envn = envn, reps = reps)
}

## Physical-position region-membership label, "rep<r>.<i>", for one rep's
## Stage-2 $regions table (from build_example_genome(compute_stage2=TRUE)).
assign_stage2_region <- function(map, regions, r) {
  lab <- rep(NA_character_, nrow(map))
  if (nrow(regions) > 0) {
    for (i in seq_len(nrow(regions))) {
      hit <- map$Chr == regions$Chr[i] & map$Pos >= regions$from[i] & map$Pos <= regions$to[i]
      lab[hit] <- sprintf("rep%d.%d", r, i)
    }
  }
  lab
}

## Shared truth-linkage computation, reproducing R/05_score_truth.R's setup
## EXACTLY (same primitives: flag_true_qtns/qtn_ld_table/score_thresholds/
## LDscnR:::.ld_outlier_units, same PARAMS from 00_config.R) for one rep --
## factored out so both the Stage-1-unit scorer and the unrestricted
## marker-wise scorer below start from the identical truth definition.
.truth_linkage_for_rep <- function(rep_bundle) {
  lu <- rep_bundle$lu; map <- copy(rep_bundle$map)
  GTs <- lu$GTs; stage1 <- lu$stage1

  p <- colSums(GTs) / nrow(GTs) / 2
  map[, p_freq := p[match(marker, colnames(GTs))]]
  map[, Va := 2 * p_freq * (1 - p_freq) * allelic_values^2]
  map <- flag_true_qtns(map, va_col = "Va", maf_col = "MAF", maf_min = MAF_KEEP, p_va_min = VA_SHARE_DETECTABLE)
  detectable_qtn <- map[true_pos_QTN == TRUE, marker]

  thr <- score_thresholds(lu$LD_decay$decay_sum, rho_r2 = TRUTH_RHO_R2, rho_d = TRUTH_RHO_D, dmax_cap = TRUTH_DMAX_CAP)
  if (sum(map$type == "QTN") == 0L) {
    qtn_lut <- data.table(marker = character(), qtn_marker = character(), r2 = numeric(), dist_bp = numeric())
  } else {
    qtn_lut <- qtn_ld_table(GTs, map, candidate_markers = map$marker, cores = 1)
  }
  qtn_lut_match <- qtn_lut[r2 > thr$r2min & dist_bp < thr$dmax & qtn_marker %in% detectable_qtn]
  list(map = map, stage1 = stage1, qtn_lut_match = qtn_lut_match)
}

## Per-marker TP/FP label for ONE rep's Stage-1 emmax_simes/lfmm_simes
## significant UNITS -- i.e. precisely the hypotheses R/06_summarise.R's
## pooled precision/recall are counted from for the RESTRICTED methods, not
## the (unscored) Stage-2 assembled regions. A marker is coloured by its
## enclosing significant unit's truth-linkage status, not its own.
score_units_tpfp <- function(rep_bundle) {
  tl <- .truth_linkage_for_rep(rep_bundle)
  em <- rep_bundle$em; lf <- rep_bundle$lf
  units_base <- LDscnR:::.ld_outlier_units(tl$stage1, tl$map, SIZE_FLOOR)
  units_base_members <- stats::setNames(units_base$members, as.character(units_base$unit_id))

  .tpfp_for <- function(units_tbl) {
    sig_ids <- as.character(units_tbl$unit_id[units_tbl$significant])
    if (!length(sig_ids)) return(data.table(marker = character(), status = character()))
    members <- units_base_members[sig_ids]
    linked_qtn <- lapply(members, function(mm) unique(tl$qtn_lut_match[marker %in% mm, qtn_marker]))
    truth_linked <- lengths(linked_qtn) > 0
    rbindlist(Map(function(mm, tlk) data.table(marker = mm, status = if (tlk) "TP" else "FP"),
                  members, truth_linked))
  }

  list(em = .tpfp_for(em$emmax_simes), lf = .tpfp_for(lf$lfmm_simes))
}

## Per-marker TP/FP label for ONE rep's UNRESTRICTED marker-wise
## significant markers (emmax_snp/lfmm_snp) -- each significant marker is
## its OWN hypothesis (no unit borrowing), matching 05_score_truth.R's
## marker_members <- setNames(as.list(marker), marker) exactly. This is the
## comparison case: does the unrestricted engine call many more markers
## significant with no truth link at all, unlike the Stage-1-restricted
## figure above?
score_markers_tpfp <- function(rep_bundle) {
  tl <- .truth_linkage_for_rep(rep_bundle)
  em <- rep_bundle$em; lf <- rep_bundle$lf
  truth_markers <- unique(tl$qtn_lut_match$marker)

  .tpfp_for <- function(marker_tbl) {
    sig <- marker_tbl$marker[marker_tbl$significant_snp]
    if (!length(sig)) return(data.table(marker = character(), status = character()))
    data.table(marker = sig, status = ifelse(sig %in% truth_markers, "TP", "FP"))
  }

  list(em = .tpfp_for(as.data.table(em$marker)), lf = .tpfp_for(as.data.table(lf$marker)))
}

## Per-marker TP/FP label at the STAGE-2 ASSEMBLED-REGION level: one whole
## region is ONE hypothesis (05_score_truth.R's PRIMARY estimand as of
## 2026-09-10), so every marker belonging to one of the region's ACTUAL
## constituent discovered-cluster members shares that region's status --
## NOT every marker physically inside its [Chr,from,to] bounds, which
## could let an untested intervening marker lend truth credit to the
## region (a real bug, fixed here and in 05_score_truth.R together:
## CLAUDE_REANALYSIS_INSTRUCTIONS.md, "Update after inspection of the
## Stage-2 Manhattan figures"). Reproduces 05_score_truth.R's
## .run_stage2()/.assemble_regions_from_units() exactly (same
## ld_prune_and_eMLG() call, same cl_sig/mk_sig/ms_sig/sub construction)
## rather than calling the ld_outlier_test() wrapper and approximating
## membership from its $regions bounds -- does NOT need
## build_example_genome(compute_stage2 = TRUE).
##
## This is deliberately DIFFERENT from score_units_tpfp(): a Stage-2
## region can merge several Stage-1 units, and scoring at the unit level
## can split one visual "peak" into mixed TP (truth-linked unit) and FP (a
## merged-in neighbouring unit that individually isn't) -- an artefact of
## Stage-1 granularity, not evidence the region call itself is wrong.
## Returns region-level TP/FP counts (n_tp_em/n_fp_em/n_tp_lf/n_fp_lf) in
## addition to the per-marker $em/$lf status tables used for colouring --
## the counts are what the figure's subtitle must report prominently, not
## coloured-marker totals (CLAUDE_REANALYSIS_INSTRUCTIONS.md: "Print TP/FP
## region counts prominently; label coloured-marker totals separately if
## retained").
score_stage2_tpfp <- function(rep_bundle) {
  tl <- .truth_linkage_for_rep(rep_bundle)
  lu <- rep_bundle$lu; em <- rep_bundle$em; lf <- rep_bundle$lf
  GTs <- lu$GTs; stage1 <- tl$stage1

  .run_stage2 <- function(cl_sub) {
    if (!nrow(cl_sub)) return(list())
    mk_sub <- unlist(cl_sub$members, use.names = FALSE)
    ms_sub <- as.data.table(stage1$map_snp)[marker %chin% mk_sub]
    sub <- structure(list(map_snp = ms_sub, clusters = cl_sub, pruned = cl_sub$core_snp),
                     class = "ld_complexity_reduction")
    pr <- ld_prune_and_eMLG(GTs = GTs[, mk_sub, drop = FALSE], stage1 = sub,
                            ld_w_col = "ld_w_095", ld_w_threshold = 0,
                            LD_decay = lu$LD_decay, min_r2_rho = stage1$params$rho,
                            score_threshold = REGION_ASSEMBLY$score_threshold,
                            distance_threshold = REGION_ASSEMBLY$distance_threshold,
                            compute_unflagged_eMLG = FALSE, min_n_loci_eMLG = 1,
                            min_n_loci_flag = 1, cores = 1)
    g <- as.data.table(pr$groups)
    stats::setNames(g$members, as.character(seq_len(nrow(g))))
  }
  .region_members_for <- function(units_tbl) {
    sig <- units_tbl[units_tbl$significant == TRUE]
    if (!nrow(sig)) return(list())
    cl <- as.data.table(stage1$clusters)
    nl <- if ("n_loci" %in% names(cl)) cl$n_loci else cl$n_snps
    .run_stage2(cl[nl >= SIZE_FLOOR][as.integer(sig$unit_id)])
  }
  .tpfp_for <- function(region_members) {
    if (!length(region_members)) return(list(status = data.table(marker = character(), status = character()), n_tp = 0L, n_fp = 0L))
    truth_linked <- vapply(region_members, function(mm) length(unique(tl$qtn_lut_match[marker %in% mm, qtn_marker])) > 0, logical(1))
    status <- rbindlist(Map(function(mm, tlk) data.table(marker = mm, status = if (tlk) "TP" else "FP"),
                            region_members, truth_linked))
    list(status = status, n_tp = sum(truth_linked), n_fp = sum(!truth_linked))
  }

  em_out <- .tpfp_for(.region_members_for(em$emmax_simes))
  lf_out <- .tpfp_for(.region_members_for(lf$lfmm_simes))
  list(em = em_out$status, lf = lf_out$status,
      n_tp_em = em_out$n_tp, n_fp_em = em_out$n_fp, n_tp_lf = lf_out$n_tp, n_fp_lf = lf_out$n_fp)
}
