## =============================================================================
## final_analysis/R/12_floor_decomposition.R
##
## Analysis 3: separate the effect of excluding small Stage-1 units (a
## filtering choice) from the effect of LD-complexity reduction itself
## (combining correlated markers into fewer tests).
##
## Four arms at every floor in FLOOR_GRID, all fed through the SAME
## assemble_stage2()/score_stage2_regions() (R/helpers_stage2_truth.R) the
## canonical pipeline uses -- never a separate implementation:
##   A. unrestricted marker (fixed baseline, identical at every floor)
##   B. floor-filtered marker (BH RECOMPUTED on the eligible subset, never
##      a post-hoc filter of the unrestricted BH result)
##   C. Stage-1 Simes (raw per-unit Simes p computed ONCE at floor=1 from
##      already-computed marker p; re-BH per floor)
##   D. Stage-1 consensus, EMMAX only (raw per-unit consensus p computed
##      ONCE at floor=1 on the FIXED relationship matrix; re-BH per floor)
## LFMM gets arms A/B/C only (no consensus arm, matching 04_lfmm.R).
##
## [!] STABLE UNIT KEY: floor=1 keeps every Stage-1 cluster (singletons
## included), so .ld_outlier_units(stage1, map, 1L)'s own unit_id -- assigned
## ONCE here, over the complete cluster set -- is used as this script's own
## stable within-combo reference; the SAME core_snp-keyed seeding functions
## the canonical pipeline uses are called for Stage-2 assembly at every
## floor (see helpers_stage2_truth.R's header comment on why unit_id itself
## is NOT stable across different size_floor values).
##
## GRM identity across floors (validation check 4): never recomputed here --
## b$GRM (built once in 02_build_ld_units.R from stage1-pruned
## representatives, ALREADY including sub-floor units) is reused unchanged.
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(parallel)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "helpers_stage2_truth.R"))
say("=== 12_floor_decomposition ===\n\n")

FLOOR_GRID <- c(1L, 2L, 3L, 5L, 10L, 20L)

## ---- one combo -> one row per (floor, arm) --------------------------------
floor_decomposition_one_combo <- function(tag, cell, rep, env) {
  combo_id <- sprintf("%s_%s_rep%d_env%d", tag, cell, rep, env)
  em <- readRDS(file.path(stage_dir("03_emmax", combo_id), "emmax.rds"))
  lfmm_file <- file.path(stage_dir("04_lfmm", combo_id), "lfmm.rds")
  have_lfmm <- file.exists(lfmm_file)
  lf <- if (have_lfmm) readRDS(lfmm_file) else NULL
  b <- readRDS(file.path(stage_dir("02_build_ld_units", combo_id), "ld_units.rds"))
  GTs <- b$GTs; map <- copy(b$map); stage1 <- b$stage1; GRM <- b$GRM
  y <- b$env$env

  ## ---- truth prerequisites (as in 05_score_truth.R/11_stage2_region_details.R) --
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

  ## ---- floor=1 base: every Stage-1 cluster, singletons included -----------
  units1 <- LDscnR:::.ld_outlier_units(stage1, map, 1L)   ## stable unit_id 1..N_total for THIS combo

  ## Arm B basis: raw marker p (already computed, floor-independent)
  pm_obs <- em$marker$p; names(pm_obs) <- em$marker$marker

  ## [!] Arms C/D are keyed by core_snp EVERYWHERE below, never by unit_id
  ## VALUE. Found by validating this script: .ld_outlier_units() assigns
  ## unit_id = seq_len(nrow(cl)) BEFORE its own internal setorder(Chr, from),
  ## so units1$unit_id is a PERMUTATION of 1:N, not units1's row position --
  ## vapply(units1$members, ...) produces a vector in units1's ROW order,
  ## and indexing that result BY unit_id VALUE (simes_p1[eligible$unit_id])
  ## silently reads the wrong unit's p-value for every row where value !=
  ## position (i.e. almost always). core_snp has no such ambiguity: it is
  ## assigned once per cluster and never reordered independently of its own
  ## row, so naming every per-unit vector by core_snp and indexing by
  ## core_snp throughout removes the failure mode rather than working around
  ## one instance of it.

  ## Arm C basis: raw per-unit Simes p at floor=1 (grouping already-computed
  ## marker p; no refit -- Simes only combines existing p-values)
  simes_p1 <- stats::setNames(vapply(units1$members, function(mm) LDscnR:::.simes(pm_obs[mm]), numeric(1)),
                              units1$core_snp)

  ## Arm D basis: raw per-unit consensus EMMAX p at floor=1 on the FIXED GRM
  ## -- the one genuinely new EMMAX fit this script needs (em$emmax_consensus
  ## was built at SIZE_FLOOR=2, excluding singletons, so it cannot be reused
  ## for a floor=1 baseline). emmax_fast() returns a plain, unnamed,
  ## POSITION-aligned vector (same convention as pm_obs/03_emmax.R) -- name
  ## it from units1$core_snp using that SAME row order (verified, not
  ## assumed: the stopifnot below checks colnames(um1) against units1$unit_id
  ## in row order, which is what licenses attaching units1$core_snp the same
  ## way).
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

  ## ---- generic per-arm scorer at one floor ---------------------------------
  ## `eligible_ids`/`p_vec`/`seed_fun` differ by arm; returns one summary row.
  .arm_row <- function(arm, floor, eligible_markers, n_units_eligible, p_vec, seed) {
    n_tested <- length(p_vec)
    q <- stats::p.adjust(p_vec, "BH")
    sig_ids <- names(p_vec)[!is.na(q) & q <= ALPHA]
    cl_sub <- seed(sig_ids)
    detail <- assemble_stage2(stage1, map, GTs, b$LD_decay, cl_sub,
                              REGION_ASSEMBLY$score_threshold, REGION_ASSEMBLY$distance_threshold)
    detail <- score_stage2_regions(detail, qtn_lut_match)
    n_regions <- nrow(detail)
    TP <- if (n_regions) sum(detail$TP) else 0L
    FP <- n_regions - TP
    recovered_qtn <- if (n_regions) unique(unlist(detail$linked_qtn[detail$TP])) else character()
    qtn_covered <- unique(qtn_lut_match[marker %in% eligible_markers, qtn_marker])
    data.table(
      tag = tag, cell = cell, rep = rep, env = env, arm = arm, floor = floor,
      n_eligible_markers = length(eligible_markers), n_eligible_units = n_units_eligible,
      n_tests = n_tested, frac_markers_retained = length(eligible_markers) / nrow(map),
      frac_dqtn_covered = if (n_dq) length(qtn_covered) / n_dq else NA_real_,
      n_significant = length(sig_ids), n_regions = n_regions, TP = TP, FP = FP,
      n_detectable_qtn = n_dq, n_recovered = length(recovered_qtn),
      n_qtn_covered_by_eligible = length(qtn_covered))
  }

  rows <- list()
  for (floor in FLOOR_GRID) {
    eligible_units1 <- units1[n_markers >= floor]
    eligible_markers <- unique(unlist(eligible_units1$members, use.names = FALSE))

    ## Arm A: unrestricted marker -- IDENTICAL at every floor (read, not recomputed)
    qA <- em$marker$q; names(qA) <- em$marker$marker
    sig_A <- em$marker$marker[em$marker$significant_snp]
    detail_A <- score_stage2_regions(
      assemble_stage2(stage1, map, GTs, b$LD_decay, stage2_seed_from_markers(stage1, map, sig_A),
                      REGION_ASSEMBLY$score_threshold, REGION_ASSEMBLY$distance_threshold), qtn_lut_match)
    TP_A <- if (nrow(detail_A)) sum(detail_A$TP) else 0L; FP_A <- nrow(detail_A) - TP_A
    rec_A <- if (nrow(detail_A)) unique(unlist(detail_A$linked_qtn[detail_A$TP])) else character()
    rows[[length(rows) + 1]] <- data.table(
      tag = tag, cell = cell, rep = rep, env = env, arm = "A_unrestricted_marker", floor = floor,
      n_eligible_markers = nrow(map), n_eligible_units = nrow(units1),
      n_tests = length(pm_obs), frac_markers_retained = 1, frac_dqtn_covered = if (n_dq) 1 else NA_real_,
      n_significant = length(sig_A), n_regions = nrow(detail_A), TP = TP_A, FP = FP_A,
      n_detectable_qtn = n_dq, n_recovered = length(rec_A), n_qtn_covered_by_eligible = n_dq)

    ## Arm B: floor-filtered marker -- BH RECOMPUTED on the eligible subset
    p_B <- pm_obs[eligible_markers]
    rows[[length(rows) + 1]] <- .arm_row("B_floor_filtered_marker", floor, eligible_markers,
      nrow(eligible_units1), p_B, function(sig_markers) stage2_seed_from_markers(stage1, map, sig_markers))

    ## Arm C: Stage-1 Simes -- subset the floor=1 raw p, re-BH
    p_C <- simes_p1[eligible_units1$core_snp]
    rows[[length(rows) + 1]] <- .arm_row("C_stage1_simes", floor, eligible_markers,
      nrow(eligible_units1), p_C, function(sig_core) stage2_seed_from_units(stage1, map, sig_core))

    ## Arm D: Stage-1 consensus (EMMAX only) -- subset the floor=1 raw p, re-BH
    p_D <- con_p1[eligible_units1$core_snp]
    rows[[length(rows) + 1]] <- .arm_row("D_stage1_consensus", floor, eligible_markers,
      nrow(eligible_units1), p_D, function(sig_core) stage2_seed_from_units(stage1, map, sig_core))
    rows[[length(rows)]][, method := "emmax"]
    rows[[length(rows) - 1]][, method := "emmax"]
    rows[[length(rows) - 2]][, method := "emmax"]
    rows[[length(rows) - 3]][, method := "emmax"]

    if (have_lfmm) {
      qLA <- lf$marker$q; names(qLA) <- lf$marker$marker
      sig_LA <- lf$marker$marker[lf$marker$significant_snp]
      detail_LA <- score_stage2_regions(
        assemble_stage2(stage1, map, GTs, b$LD_decay, stage2_seed_from_markers(stage1, map, sig_LA),
                        REGION_ASSEMBLY$score_threshold, REGION_ASSEMBLY$distance_threshold), qtn_lut_match)
      TP_LA <- if (nrow(detail_LA)) sum(detail_LA$TP) else 0L; FP_LA <- nrow(detail_LA) - TP_LA
      rec_LA <- if (nrow(detail_LA)) unique(unlist(detail_LA$linked_qtn[detail_LA$TP])) else character()
      r <- data.table(
        tag = tag, cell = cell, rep = rep, env = env, arm = "A_unrestricted_marker", floor = floor,
        n_eligible_markers = nrow(map), n_eligible_units = nrow(units1),
        n_tests = length(lp_obs), frac_markers_retained = 1, frac_dqtn_covered = if (n_dq) 1 else NA_real_,
        n_significant = length(sig_LA), n_regions = nrow(detail_LA), TP = TP_LA, FP = FP_LA,
        n_detectable_qtn = n_dq, n_recovered = length(rec_LA), n_qtn_covered_by_eligible = n_dq,
        method = "lfmm")
      rows[[length(rows) + 1]] <- r

      p_LB <- lp_obs[eligible_markers]
      rB <- .arm_row("B_floor_filtered_marker", floor, eligible_markers, nrow(eligible_units1), p_LB,
                     function(sig_markers) stage2_seed_from_markers(stage1, map, sig_markers))
      rB[, method := "lfmm"]; rows[[length(rows) + 1]] <- rB

      p_LC <- lsimes_p1[eligible_units1$core_snp]
      rC <- .arm_row("C_stage1_simes", floor, eligible_markers, nrow(eligible_units1), p_LC,
                     function(sig_core) stage2_seed_from_units(stage1, map, sig_core))
      rC[, method := "lfmm"]; rows[[length(rows) + 1]] <- rC
    }
  }
  rbindlist(rows, fill = TRUE)
}

if (sys.nframe() == 0L) {
  combos <- CJ(tag = TAGS_ALL, cell = CELLS_ALL, rep = REPS_ALL, env = ENVS_ALL)
  say("[1] floor decomposition for %d combos x %d floors x 4 EMMAX + 3 LFMM arms\n",
      nrow(combos), length(FLOOR_GRID))
  t0 <- Sys.time()
  res <- mclapply(seq_len(nrow(combos)), function(i) {
    tryCatch(floor_decomposition_one_combo(combos$tag[i], combos$cell[i], combos$rep[i], combos$env[i]),
             error = function(e) { message(sprintf("[12_floor_decomposition] %s_%s_rep%d_env%d: %s",
                                                    combos$tag[i], combos$cell[i], combos$rep[i], combos$env[i],
                                                    conditionMessage(e))); NULL })
  }, mc.cores = 7)
  errs <- sum(vapply(res, is.null, logical(1)))
  dt <- rbindlist(Filter(Negate(is.null), res), fill = TRUE)
  say("[2] %s rows from %d combos (%d combo errors) in %.1f min\n", format(nrow(dt), big.mark = ","),
      nrow(combos), errs, as.numeric(difftime(Sys.time(), t0, units = "mins")))
  saveRDS(dt, "out_final_v1/12_floor_decomposition_raw.rds", compress = "xz")

  ## ---- validation checks 1, 3 (in-script, before any pooling) --------------
  a1 <- dt[method == "emmax" & arm == "A_unrestricted_marker" & floor == 1L]
  bf1 <- dt[method == "emmax" & arm == "B_floor_filtered_marker" & floor == 1L]
  key <- c("tag", "cell", "rep", "env")
  chk1 <- merge(a1[, .(tag, cell, rep, env, TP, FP, n_significant)],
               bf1[, .(tag, cell, rep, env, TP, FP, n_significant)], by = key, suffixes = c("_A", "_B1"))
  ok1 <- all(chk1$TP_A == chk1$TP_B1 & chk1$FP_A == chk1$FP_B1 & chk1$n_significant_A == chk1$n_significant_B1)
  say("[3] VALIDATION 'floor=1 arm B == arm A' (EMMAX): %s (%d combos compared)\n", ok1, nrow(chk1))
  if (!ok1) stop("floor=1 arm-B/arm-A mismatch -- see chk1")

  say("\n[4] pooling by (tag, cell, method, arm, floor) -- ratio of pooled counts\n")
  pooled <- dt[, .(TP = sum(TP), FP = sum(FP), n_regions = sum(n_regions),
                   n_detectable_qtn = sum(n_detectable_qtn), n_recovered = sum(n_recovered),
                   n_tests = sum(n_tests), n_eligible_markers = sum(n_eligible_markers),
                   n_eligible_units = sum(n_eligible_units),
                   frac_markers_retained = mean(frac_markers_retained),
                   frac_dqtn_covered = mean(frac_dqtn_covered, na.rm = TRUE)),
               by = .(tag, cell, method, arm, floor)]
  pooled[, `:=`(precision = TP / pmax(TP + FP, 1), recall = n_recovered / pmax(n_detectable_qtn, 1),
                realised_fdp = FP / pmax(TP + FP, 1))]

  ## ---- map-cluster bootstrap CIs, methods/arms/floors paired per resample --
  rep_dt <- dt[, .(TP = sum(TP), FP = sum(FP), n_detectable_qtn = sum(n_detectable_qtn),
                   n_recovered = sum(n_recovered)), by = .(tag, cell, method, arm, floor, rep)]
  n_rep <- length(REPS_ALL)
  strata <- unique(pooled[, .(tag, cell, method, arm, floor)])
  ci_rows <- vector("list", nrow(strata))
  for (i in seq_len(nrow(strata))) {
    s <- strata[i]
    rd <- rep_dt[tag == s$tag & cell == s$cell & method == s$method & arm == s$arm & floor == s$floor][match(REPS_ALL, rep)]
    rd[is.na(TP), `:=`(TP = 0L, FP = 0L, n_detectable_qtn = 0L, n_recovered = 0L)]
    mat <- as.matrix(rd[, .(TP, FP, n_detectable_qtn, n_recovered)])
    bs <- bootstrap_rep_matrix(mat, N_BOOTSTRAP, SEEDS[["bootstrap"]])
    prec_b <- bs[, "TP"] / pmax(bs[, "TP"] + bs[, "FP"], 1)
    rec_b <- bs[, "n_recovered"] / pmax(bs[, "n_detectable_qtn"], 1)
    fdp_b <- bs[, "FP"] / pmax(bs[, "TP"] + bs[, "FP"], 1)
    ci_rows[[i]] <- data.table(tag = s$tag, cell = s$cell, method = s$method, arm = s$arm, floor = s$floor,
                               precision_ci_lo = ci_quantile(prec_b)[1], precision_ci_hi = ci_quantile(prec_b)[2],
                               recall_ci_lo = ci_quantile(rec_b)[1], recall_ci_hi = ci_quantile(rec_b)[2],
                               fdp_ci_lo = ci_quantile(fdp_b)[1], fdp_ci_hi = ci_quantile(fdp_b)[2])
  }
  ci_dt <- rbindlist(ci_rows)
  sweep <- merge(pooled, ci_dt, by = c("tag", "cell", "method", "arm", "floor"))
  setorder(sweep, method, arm, tag, cell, floor)
  fwrite(sweep, "results/simulation_floor_sweep.tsv", sep = "\t")
  say("[5] wrote results/simulation_floor_sweep.tsv (%d rows)\n", nrow(sweep))

  ## ---- grand-pooled (across all 14 strata) contrasts, paired bootstrap -----
  say("\n[6] grand-pooled paired contrasts: B-A, C-B, C-A per method x floor\n")
  grand_rep <- dt[, .(TP = sum(TP), FP = sum(FP), n_detectable_qtn = sum(n_detectable_qtn),
                      n_recovered = sum(n_recovered), n_tests = sum(n_tests), n_regions = sum(n_regions)),
                 by = .(method, arm, floor, rep)]
  contrast_pairs <- list(c("B_floor_filtered_marker", "A_unrestricted_marker"),
                         c("C_stage1_simes", "B_floor_filtered_marker"),
                         c("C_stage1_simes", "A_unrestricted_marker"),
                         c("D_stage1_consensus", "B_floor_filtered_marker"),
                         c("D_stage1_consensus", "A_unrestricted_marker"))
  contrast_rows <- list()
  ## [!] `this_floor`, NOT `floor` -- grand_rep has its own `floor` COLUMN, and
  ## `floor == floor` inside a data.table `[` filter resolves the RHS to that
  ## column too (data.table's NSE looks up the column before the loop
  ## variable of the same name), so it silently compared the column to
  ## itself and matched every row regardless of which floor was current.
  ## PK caught this from the contrast table repeating identical rows at
  ## every floor -- confirmed and fixed here, not touched anywhere else
  ## (R/13/R/14 use `..var`-prefixed or `s$floor`-style lookups specifically
  ## to avoid this exact class of bug; this was the one place a bare loop
  ## variable collided with a same-named data.table column).
  for (m in unique(grand_rep$method)) for (this_floor in FLOOR_GRID) for (pr in contrast_pairs) {
    x <- grand_rep[method == m & arm == pr[1] & floor == this_floor][match(REPS_ALL, rep)]
    z <- grand_rep[method == m & arm == pr[2] & floor == this_floor][match(REPS_ALL, rep)]
    if (!nrow(x) || !nrow(z) || all(is.na(x$TP)) || all(is.na(z$TP))) next
    x[is.na(TP), `:=`(TP = 0L, FP = 0L, n_detectable_qtn = 0L, n_recovered = 0L, n_tests = 0L, n_regions = 0L)]
    z[is.na(TP), `:=`(TP = 0L, FP = 0L, n_detectable_qtn = 0L, n_recovered = 0L, n_tests = 0L, n_regions = 0L)]
    bx <- bootstrap_rep_matrix(as.matrix(x[, .(TP, FP, n_detectable_qtn, n_recovered, n_tests, n_regions)]),
                               N_BOOTSTRAP, SEEDS[["bootstrap"]])
    bz <- bootstrap_rep_matrix(as.matrix(z[, .(TP, FP, n_detectable_qtn, n_recovered, n_tests, n_regions)]),
                               N_BOOTSTRAP, SEEDS[["bootstrap"]])
    prec_x <- bx[, "TP"] / pmax(bx[, "TP"] + bx[, "FP"], 1); prec_z <- bz[, "TP"] / pmax(bz[, "TP"] + bz[, "FP"], 1)
    rec_x <- bx[, "n_recovered"] / pmax(bx[, "n_detectable_qtn"], 1); rec_z <- bz[, "n_recovered"] / pmax(bz[, "n_detectable_qtn"], 1)
    fdp_x <- bx[, "FP"] / pmax(bx[, "TP"] + bx[, "FP"], 1); fdp_z <- bz[, "FP"] / pmax(bz[, "TP"] + bz[, "FP"], 1)
    d_prec <- prec_x - prec_z; d_rec <- rec_x - rec_z; d_fdp <- fdp_x - fdp_z
    d_tests <- bx[, "n_tests"] - bz[, "n_tests"]; d_regions <- bx[, "n_regions"] - bz[, "n_regions"]
    px <- sum(x$TP) / pmax(sum(x$TP) + sum(x$FP), 1); pz <- sum(z$TP) / pmax(sum(z$TP) + sum(z$FP), 1)
    rx <- sum(x$n_recovered) / pmax(sum(x$n_detectable_qtn), 1); rz <- sum(z$n_recovered) / pmax(sum(z$n_detectable_qtn), 1)
    fx <- sum(x$FP) / pmax(sum(x$TP) + sum(x$FP), 1); fz <- sum(z$FP) / pmax(sum(z$TP) + sum(z$FP), 1)
    contrast_rows[[length(contrast_rows) + 1]] <- data.table(
      method = m, floor = this_floor, arm_x = pr[1], arm_z = pr[2],
      diff_precision = px - pz, diff_precision_ci_lo = ci_quantile(d_prec)[1], diff_precision_ci_hi = ci_quantile(d_prec)[2],
      diff_recall = rx - rz, diff_recall_ci_lo = ci_quantile(d_rec)[1], diff_recall_ci_hi = ci_quantile(d_rec)[2],
      diff_fdp = fx - fz, diff_fdp_ci_lo = ci_quantile(d_fdp)[1], diff_fdp_ci_hi = ci_quantile(d_fdp)[2],
      diff_n_tests = sum(x$n_tests) - sum(z$n_tests), diff_n_tests_ci_lo = ci_quantile(d_tests)[1], diff_n_tests_ci_hi = ci_quantile(d_tests)[2],
      diff_n_regions = sum(x$n_regions) - sum(z$n_regions), diff_n_regions_ci_lo = ci_quantile(d_regions)[1], diff_n_regions_ci_hi = ci_quantile(d_regions)[2])
  }
  contrasts <- rbindlist(contrast_rows)
  fwrite(contrasts, "results/simulation_floor_contrasts.tsv", sep = "\t")
  say("[7] wrote results/simulation_floor_contrasts.tsv (%d rows)\n", nrow(contrasts))

  ## ---- phenotype-blind floor-choice table (no association testing) --------
  say("\n[8] phenotype-blind floor-choice table (structural only, no test statistics)\n")
  blind_rows <- list()
  for (i in seq_len(nrow(combos))) {
    combo_id <- sprintf("%s_%s_rep%d_env%d", combos$tag[i], combos$cell[i], combos$rep[i], combos$env[i])
    bf <- file.path(stage_dir("02_build_ld_units", combo_id), "ld_units.rds")
    if (!file.exists(bf)) next
    bb <- readRDS(bf)
    u1 <- LDscnR:::.ld_outlier_units(bb$stage1, bb$map, 1L)
    for (floor in FLOOR_GRID) {
      el <- u1[n_markers >= floor]
      mk <- unique(unlist(el$members, use.names = FALSE))
      cov_by_chr <- el[, .(n_units = .N, n_markers = sum(n_markers)), by = Chr]
      blind_rows[[length(blind_rows) + 1]] <- data.table(
        tag = combos$tag[i], cell = combos$cell[i], rep = combos$rep[i], env = combos$env[i], floor = floor,
        n_eligible_units = nrow(el), n_tests = nrow(el), frac_markers_retained = length(mk) / nrow(bb$map),
        n_chr_represented = nrow(cov_by_chr))
    }
  }
  blind <- rbindlist(blind_rows)
  blind_pooled <- blind[, .(mean_n_eligible_units = mean(n_eligible_units), mean_n_tests = mean(n_tests),
                            mean_frac_markers_retained = mean(frac_markers_retained),
                            mean_n_chr_represented = mean(n_chr_represented)), by = floor]
  setorder(blind_pooled, floor)
  fwrite(blind_pooled, "results/simulation_floor_choice_blind.tsv", sep = "\t")
  say("[9] wrote results/simulation_floor_choice_blind.tsv\n")
  print(blind_pooled)

  write_receipt("12_floor_decomposition", inputs = character(),
                params = list(floor_grid = FLOOR_GRID, alpha = ALPHA, size_floor_canonical = SIZE_FLOOR,
                              n_bootstrap = N_BOOTSTRAP, seed_bootstrap = SEEDS[["bootstrap"]]),
                outputs = c("results/simulation_floor_sweep.tsv", "results/simulation_floor_contrasts.tsv",
                           "results/simulation_floor_choice_blind.tsv"))
}
