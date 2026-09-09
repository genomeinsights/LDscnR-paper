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

## Per-marker TP/FP/"ns" (not significant) label, reproducing
## R/05_score_truth.R's hypothesis-level truth-linkage logic EXACTLY (same
## primitives: flag_true_qtns/qtn_ld_table/score_thresholds/
## LDscnR:::.ld_outlier_units, same PARAMS from 00_config.R) for ONE rep's
## Stage-1 emmax_simes/lfmm_simes significant units -- i.e. precisely the
## hypotheses R/06_summarise.R's pooled precision/recall are counted from,
## not the (unscored) Stage-2 assembled regions.
score_units_tpfp <- function(rep_bundle) {
  lu <- rep_bundle$lu; map <- copy(rep_bundle$map); em <- rep_bundle$em; lf <- rep_bundle$lf
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

  units_base <- LDscnR:::.ld_outlier_units(stage1, map, SIZE_FLOOR)
  units_base_members <- stats::setNames(units_base$members, as.character(units_base$unit_id))

  .tpfp_for <- function(units_tbl) {
    sig_ids <- as.character(units_tbl$unit_id[units_tbl$significant])
    if (!length(sig_ids)) return(data.table(marker = character(), status = character()))
    members <- units_base_members[sig_ids]
    linked_qtn <- lapply(members, function(mm) unique(qtn_lut_match[marker %in% mm, qtn_marker]))
    truth_linked <- lengths(linked_qtn) > 0
    rbindlist(Map(function(mm, tl) data.table(marker = mm, status = if (tl) "TP" else "FP"),
                  members, truth_linked))
  }

  list(em = .tpfp_for(em$emmax_simes), lf = .tpfp_for(lf$lfmm_simes))
}
