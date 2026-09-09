## =============================================================================
## final_analysis/R/06_summarise.R
##
## Pool the full 1,400-combination grid's per-combo truth_scores.rds into
## the primary reportable estimates, per CLAUDE_REANALYSIS_INSTRUCTIONS.md's
## "Pooling and uncertainty" section.
##
## Point estimates are RATIOS OF POOLED COUNTS, never means of per-combo
## precision/recall:
##   precision = sum(TP) / (sum(TP) + sum(FP))
##   recall    = recovered unique QTN / detectable unique QTN
## "Unique QTN" is scoped WITHIN each combo -- QTN identity is not
## comparable across reps (each rep is an independent simulated genome with
## its own reference map, rec_map<rep>.rds), so pooling recall across
## combos means summing each combo's own (recovered, detectable) counts,
## never deduplicating QTN identity globally. 05_score_truth.R's saved `FN`
## and `n_detectable_qtn` already give this per combo:
##   n_recovered_this_combo = n_detectable_qtn - FN
##
## UNCERTAINTY: the ten environmental continuations for a map share its
## burn-in (confirmed directly against the .ini files, see module_sim_3sp53/
## R/14_random_removal_control.R's own history) -- NOT ten independent
## replicates. The primary CI is therefore a CLUSTER bootstrap over the ten
## map/burn-in IDs (rep), retaining ALL ten environmental continuations
## belonging to each resampled rep, with all methods paired within a
## resample (so method CONTRASTS, not just each method's own CI, are valid
## paired comparisons). A crossed two-way bootstrap (reps and envs resampled
## INDEPENDENTLY, Cartesian product) is run as a labelled sensitivity check
## on whether treating the ten environmental surfaces as themselves sampled
## changes the picture.
##
## Both bootstraps are vectorised as matrix algebra over PRE-AGGREGATED
## rep-level (or rep x env-level) sums, not by literally re-subsetting rows
## per replicate -- resampling only ever changes how many times a whole
## rep's (or rep,env cell's) already-pooled counts are counted, never the
## counts themselves.
## =============================================================================
suppressMessages({library(data.table)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))

METHODS <- c("emmax_snp", "emmax_snp_nonsingleton", "emmax_simes", "emmax_consensus")
REFERENCE_METHOD <- "emmax_snp"   ## "its unrestricted marker-wise engine" -- the contrast baseline

## LFMM is "retain only as a portability analysis" (instructions) -- a SEPARATE
## arm with its own within-engine baseline (lfmm_simes vs lfmm_snp), not mixed
## into the primary EMMAX methods/contrasts above. Does the same phenotype-
## blind Stage-1 restriction still help when the association engine changes?
LFMM_METHODS <- c("lfmm_snp", "lfmm_simes")
LFMM_REFERENCE_METHOD <- "lfmm_snp"

## Stage-2 assembled-region scoring (05_score_truth.R, added 2026-09-10):
## one WHOLE region = one hypothesis, fixing the unit-level artefact where
## a single reported region can contain both a truth-linked AND a
## non-linked Stage-1 unit, showing mixed TP/FP for what is really one
## call. Reported as a PAIRED contrast against its own unit-level
## counterpart (region vs unit, SAME method, same bootstrap draw) -- a
## different comparison in kind from REFERENCE_METHOD above (different
## methods, same scoring granularity).
REGION_PAIRS <- list(
  c(region = "emmax_simes_region", unit = "emmax_simes"),
  c(region = "emmax_consensus_region", unit = "emmax_consensus"),
  c(region = "lfmm_simes_region", unit = "lfmm_simes")
)

## ---- 1. load every combo's truth scores into one table -----------------------
load_all_scores <- function() {
  rows <- list()
  for (tag in TAGS_ALL) for (cell in CELLS_ALL) for (rep in REPS_ALL) for (envn in ENVS_ALL) {
    combo_id <- sprintf("%s_%s_rep%d_env%d", tag, cell, rep, envn)
    f <- file.path(stage_dir("05_score_truth", combo_id), "truth_scores.rds")
    if (!file.exists(f)) stop("missing truth_scores.rds for ", combo_id, " -- gate 4 must be complete before summarising.")
    rows[[combo_id]] <- readRDS(f)$scores
  }
  dt <- rbindlist(rows)
  dt[, n_recovered := n_detectable_qtn - FN]
  dt[]
}

## ---- 2. point estimates: ratios of POOLED counts, never means ----------------
pooled_point <- function(dt) {
  dt[, .(n_combo = .N, n_tested = sum(n_tested), n_significant = sum(n_significant),
        TP = sum(TP), FP = sum(FP), FN = sum(FN),
        precision = sum(TP) / max(sum(TP) + sum(FP), 1),
        n_detectable_qtn = sum(n_detectable_qtn), n_recovered = sum(n_recovered),
        recall = sum(n_recovered) / max(sum(n_detectable_qtn), 1),
        detectable_qtn_covered = sum(detectable_qtn_covered),
        coverage = sum(detectable_qtn_covered) / max(sum(n_detectable_qtn), 1),
        mean_eligible_marker_fraction = mean(eligible_marker_fraction)),
     by = .(tag, cell, method)]
}

## ---- 3. primary map-cluster bootstrap -----------------------------------------
## Pre-aggregate to (tag, cell, method, rep) -- the unit the bootstrap
## resamples IS the rep (map/burn-in ID); all 10 envs within a rep always
## move together, so summing them first loses nothing the resampling could
## ever separate.
rep_level_stats <- function(dt) {
  dt[, .(TP = sum(TP), FP = sum(FP), n_detectable_qtn = sum(n_detectable_qtn), n_recovered = sum(n_recovered),
        n_tested = sum(n_tested), n_significant = sum(n_significant)),
     by = .(tag, cell, method, rep)]
}

## For ONE (tag, cell) stratum: B bootstrap replicates of the 10 reps (with
## replacement), vectorised as a 10 x B multiplicity matrix, applied via
## matrix multiplication to each method's 10-vector of rep-level sums.
## Returns a list: point (data.table) and boot (long data.table, one row
## per method x replicate, for percentile CIs and paired contrasts).
map_cluster_bootstrap_stratum <- function(rep_dt, B, seed, methods = METHODS) {
  set.seed(seed)
  n_rep <- length(REPS_ALL)
  draws <- matrix(sample.int(n_rep, size = n_rep * B, replace = TRUE), nrow = n_rep, ncol = B)
  mult <- apply(draws, 2, tabulate, nbins = n_rep)   ## n_rep x B

  boot_rows <- list()
  for (m in methods) {
    rm <- rep_dt[method == m][order(rep)]
    stopifnot("rep_level_stats must have exactly REPS_ALL rows per method" = nrow(rm) == n_rep && identical(rm$rep, REPS_ALL))
    TP_b <- as.numeric(crossprod(mult, rm$TP))
    FP_b <- as.numeric(crossprod(mult, rm$FP))
    ndq_b <- as.numeric(crossprod(mult, rm$n_detectable_qtn))
    nrec_b <- as.numeric(crossprod(mult, rm$n_recovered))
    boot_rows[[m]] <- data.table(method = m, b = seq_len(B),
                                 precision = TP_b / pmax(TP_b + FP_b, 1),
                                 recall = nrec_b / pmax(ndq_b, 1))
  }
  rbindlist(boot_rows)
}

ci <- function(x, probs = c(0.025, 0.975)) stats::quantile(x, probs, na.rm = TRUE, names = FALSE)

## ---- 4. crossed two-way bootstrap (sensitivity) -------------------------------
## Reps AND envs resampled INDEPENDENTLY, Cartesian product -- needs
## (rep, env)-level cells, not the rep-collapsed table above.
combo_level_stats <- function(dt) {
  dt[, .(TP = sum(TP), FP = sum(FP), n_detectable_qtn = sum(n_detectable_qtn), n_recovered = sum(n_recovered)),
     by = .(tag, cell, method, rep, env)]
}

crossed_bootstrap_stratum <- function(combo_dt, B, seed) {
  set.seed(seed)
  n_rep <- length(REPS_ALL); n_env <- length(ENVS_ALL)
  rep_draws <- matrix(sample.int(n_rep, size = n_rep * B, replace = TRUE), nrow = n_rep, ncol = B)
  env_draws <- matrix(sample.int(n_env, size = n_env * B, replace = TRUE), nrow = n_env, ncol = B)
  rep_mult <- apply(rep_draws, 2, tabulate, nbins = n_rep)   ## n_rep x B
  env_mult <- apply(env_draws, 2, tabulate, nbins = n_env)   ## n_env x B

  boot_rows <- list()
  for (m in METHODS) {
    cm <- combo_dt[method == m]
    TPmat <- matrix(0, n_rep, n_env); FPmat <- TPmat; ndqmat <- TPmat; nrecmat <- TPmat
    idx <- cbind(match(cm$rep, REPS_ALL), match(cm$env, ENVS_ALL))
    TPmat[idx] <- cm$TP; FPmat[idx] <- cm$FP; ndqmat[idx] <- cm$n_detectable_qtn; nrecmat[idx] <- cm$n_recovered
    prec <- numeric(B); rec <- numeric(B)
    for (b in seq_len(B)) {
      rw <- rep_mult[, b]; ew <- env_mult[, b]
      TP_b <- sum((rw %o% ew) * TPmat); FP_b <- sum((rw %o% ew) * FPmat)
      ndq_b <- sum((rw %o% ew) * ndqmat); nrec_b <- sum((rw %o% ew) * nrecmat)
      prec[b] <- TP_b / max(TP_b + FP_b, 1); rec[b] <- nrec_b / max(ndq_b, 1)
    }
    boot_rows[[m]] <- data.table(method = m, b = seq_len(B), precision = prec, recall = rec)
  }
  rbindlist(boot_rows)
}

## ---- 3b. one arm's map-cluster bootstrap + paired contrasts, all strata -------
## Factored out of the driver so the primary EMMAX arm and the separate LFMM
## portability arm run through identical bootstrap/contrast logic, each with
## its own methods/reference_method (and its own bootstrap seed stream --
## seed_offset keeps the two arms' resampling independent).
bootstrap_arm <- function(dt, point, methods, reference_method, B, seed_offset = 0) {
  rep_dt <- rep_level_stats(dt)
  strata <- unique(dt[, .(tag, cell)])
  perf_rows <- list(); contrast_rows <- list()
  for (i in seq_len(nrow(strata))) {
    tg <- strata$tag[i]; cl <- strata$cell[i]
    say("    %s/%s\n", tg, cl)
    rd <- rep_dt[tag == tg & cell == cl]
    boot <- map_cluster_bootstrap_stratum(rd, B, seed = SEEDS[["bootstrap"]] + seed_offset + i, methods = methods)
    pt <- point[tag == tg & cell == cl]
    for (m in methods) {
      bm <- boot[method == m]
      pr <- pt[method == m]
      perf_rows[[length(perf_rows) + 1]] <- data.table(
        tag = tg, cell = cl, method = m, n_tested = pr$n_tested, n_significant = pr$n_significant,
        TP = pr$TP, FP = pr$FP, FN = pr$FN, n_detectable_qtn = pr$n_detectable_qtn, n_recovered = pr$n_recovered,
        coverage = pr$coverage,
        precision = pr$precision, precision_ci_lo = ci(bm$precision)[1], precision_ci_hi = ci(bm$precision)[2],
        recall = pr$recall, recall_ci_lo = ci(bm$recall)[1], recall_ci_hi = ci(bm$recall)[2])
    }
    ## paired contrasts vs the reference method, from the SAME bootstrap replicates
    ref <- boot[method == reference_method][order(b)]
    for (m in setdiff(methods, reference_method)) {
      bm <- boot[method == m][order(b)]
      d_prec <- bm$precision - ref$precision
      d_rec <- bm$recall - ref$recall
      pt_m <- point[tag == tg & cell == cl & method == m]
      pt_r <- point[tag == tg & cell == cl & method == reference_method]
      contrast_rows[[length(contrast_rows) + 1]] <- data.table(
        tag = tg, cell = cl, method = m, reference_method = reference_method,
        diff_precision = pt_m$precision - pt_r$precision, diff_precision_ci_lo = ci(d_prec)[1], diff_precision_ci_hi = ci(d_prec)[2],
        diff_recall = pt_m$recall - pt_r$recall, diff_recall_ci_lo = ci(d_rec)[1], diff_recall_ci_hi = ci(d_rec)[2])
    }
  }
  list(performance = rbindlist(perf_rows), contrasts = rbindlist(contrast_rows))
}

## ---- 3c. region-vs-unit granularity contrasts (paired, SAME bootstrap draw) ---
## Different in kind from bootstrap_arm()'s contrasts: those compare
## different METHODS at a shared scoring granularity against one common
## reference; this compares the SAME method's two scoring granularities
## (region vs unit) against EACH OTHER, pair by pair -- so there is no
## single shared reference_method, and pairs is a list of (region, unit)
## name pairs instead. Absolute performance rows are produced ONLY for the
## region-level methods (the unit-level ones already have one from
## bootstrap_arm() above); this arm exists for the paired contrast.
granularity_arm <- function(dt, point, pairs, B, seed_offset = 0) {
  methods <- unique(unlist(lapply(pairs, unname)))
  rep_dt <- rep_level_stats(dt)
  strata <- unique(dt[, .(tag, cell)])
  perf_rows <- list(); contrast_rows <- list()
  for (i in seq_len(nrow(strata))) {
    tg <- strata$tag[i]; cl <- strata$cell[i]
    say("    %s/%s\n", tg, cl)
    rd <- rep_dt[tag == tg & cell == cl]
    boot <- map_cluster_bootstrap_stratum(rd, B, seed = SEEDS[["bootstrap"]] + seed_offset + i, methods = methods)
    pt <- point[tag == tg & cell == cl]
    for (m in vapply(pairs, `[[`, character(1), "region")) {
      bm <- boot[method == m]; pr <- pt[method == m]
      perf_rows[[length(perf_rows) + 1]] <- data.table(
        tag = tg, cell = cl, method = m, n_tested = pr$n_tested, n_significant = pr$n_significant,
        TP = pr$TP, FP = pr$FP, FN = pr$FN, n_detectable_qtn = pr$n_detectable_qtn, n_recovered = pr$n_recovered,
        coverage = pr$coverage,
        precision = pr$precision, precision_ci_lo = ci(bm$precision)[1], precision_ci_hi = ci(bm$precision)[2],
        recall = pr$recall, recall_ci_lo = ci(bm$recall)[1], recall_ci_hi = ci(bm$recall)[2])
    }
    for (pair in pairs) {
      m_region <- pair[["region"]]; m_unit <- pair[["unit"]]
      b_region <- boot[method == m_region][order(b)]
      b_unit <- boot[method == m_unit][order(b)]
      d_prec <- b_region$precision - b_unit$precision
      d_rec <- b_region$recall - b_unit$recall
      pt_r <- point[tag == tg & cell == cl & method == m_region]
      pt_u <- point[tag == tg & cell == cl & method == m_unit]
      contrast_rows[[length(contrast_rows) + 1]] <- data.table(
        tag = tg, cell = cl, method = m_region, reference_method = m_unit,
        diff_precision = pt_r$precision - pt_u$precision, diff_precision_ci_lo = ci(d_prec)[1], diff_precision_ci_hi = ci(d_prec)[2],
        diff_recall = pt_r$recall - pt_u$recall, diff_recall_ci_lo = ci(d_rec)[1], diff_recall_ci_hi = ci(d_rec)[2])
    }
  }
  list(performance = rbindlist(perf_rows), contrasts = rbindlist(contrast_rows))
}

## ---- driver --------------------------------------------------------------------
summarise_grid <- function(B = N_BOOTSTRAP, do_crossed_sensitivity = TRUE) {
  say("[1] loading all 1400 combos' truth scores\n")
  dt <- load_all_scores()
  have_lfmm <- all(LFMM_METHODS %in% dt$method)
  region_methods_present <- unique(unlist(lapply(REGION_PAIRS, unname)))
  have_region <- all(region_methods_present %in% dt$method)
  n_methods_seen <- length(unique(dt$method))
  say("    %d rows (%d combos x %d methods: %d EMMAX%s%s)\n", nrow(dt), nrow(dt) / n_methods_seen, n_methods_seen,
      length(METHODS), if (have_lfmm) sprintf(" + %d LFMM", length(LFMM_METHODS)) else "",
      if (have_region) sprintf(" + %d Stage-2-region", length(region_methods_present)) else "")

  say("\n[2] point estimates (pooled counts, one row per cell x tag x method)\n")
  point <- pooled_point(dt)
  print(point[, .(tag, cell, method, n_tested, n_significant, TP, FP, precision, recall, coverage)])

  say("\n[3] primary map-cluster bootstrap (B=%d per stratum)\n", B)
  arm <- bootstrap_arm(dt, point, METHODS, REFERENCE_METHOD, B)
  performance <- arm$performance
  contrasts <- arm$contrasts

  performance_lfmm <- NULL; contrasts_lfmm <- NULL
  if (have_lfmm) {
    say("\n[3b] LFMM portability arm: map-cluster bootstrap (B=%d per stratum)\n", B)
    arm_lfmm <- bootstrap_arm(dt, point, LFMM_METHODS, LFMM_REFERENCE_METHOD, B, seed_offset = 5000)
    performance_lfmm <- arm_lfmm$performance
    contrasts_lfmm <- arm_lfmm$contrasts
  }

  performance_region <- NULL; contrasts_region <- NULL
  if (have_region) {
    say("\n[3c] Stage-2 region-vs-unit granularity contrasts (B=%d per stratum)\n", B)
    arm_region <- granularity_arm(dt, point, REGION_PAIRS, B, seed_offset = 6000)
    performance_region <- arm_region$performance
    contrasts_region <- arm_region$contrasts
  }

  sensitivity <- NULL
  if (do_crossed_sensitivity) {
    say("\n[4] crossed two-way bootstrap sensitivity (grand-pooled across all cells/tags, B=%d)\n", B)
    combo_dt <- combo_level_stats(dt)
    ## grand pooled (not per-stratum, per instructions' "if time permits" framing) --
    ## collapse cell/tag into one pool for this sensitivity pass
    combo_dt_grand <- combo_dt[, .(TP = sum(TP), FP = sum(FP), n_detectable_qtn = sum(n_detectable_qtn), n_recovered = sum(n_recovered)),
                               by = .(method, rep, env)]
    combo_dt_grand[, `:=`(tag = "ALL", cell = "ALL")]
    boot_cross <- crossed_bootstrap_stratum(combo_dt_grand, B, seed = SEEDS[["bootstrap"]] + 9999)
    rep_dt_grand <- dt[, .(TP = sum(TP), FP = sum(FP), n_detectable_qtn = sum(n_detectable_qtn), n_recovered = sum(n_recovered)),
                       by = .(method, rep)]
    boot_primary_grand <- map_cluster_bootstrap_stratum(rep_dt_grand, B, seed = SEEDS[["bootstrap"]] + 8888)
    point_grand <- pooled_point(copy(dt)[, `:=`(tag = "ALL", cell = "ALL")])
    sens_rows <- list()
    for (m in METHODS) {
      pt <- point_grand[method == m]
      bp <- boot_primary_grand[method == m]; bc <- boot_cross[method == m]
      sens_rows[[m]] <- data.table(
        method = m, precision = pt$precision,
        primary_ci_lo = ci(bp$precision)[1], primary_ci_hi = ci(bp$precision)[2],
        crossed_ci_lo = ci(bc$precision)[1], crossed_ci_hi = ci(bc$precision)[2],
        recall = pt$recall,
        primary_recall_ci_lo = ci(bp$recall)[1], primary_recall_ci_hi = ci(bp$recall)[2],
        crossed_recall_ci_lo = ci(bc$recall)[1], crossed_recall_ci_hi = ci(bc$recall)[2])
    }
    sensitivity <- rbindlist(sens_rows)
    say("    grand-pooled: primary (map-cluster, conditional on observed env surfaces) vs crossed (envs also resampled)\n")
    print(sensitivity)
  }

  list(point = point, performance = performance, contrasts = contrasts,
       performance_lfmm = performance_lfmm, contrasts_lfmm = contrasts_lfmm,
       performance_region = performance_region, contrasts_region = contrasts_region,
       sensitivity = sensitivity, raw = dt)
}

if (sys.nframe() == 0L) {
  res <- summarise_grid()
  dir.create("results", showWarnings = FALSE)
  fwrite(res$performance, "results/simulation_performance.tsv", sep = "\t")
  fwrite(res$contrasts, "results/simulation_method_contrasts.tsv", sep = "\t")
  if (!is.null(res$sensitivity)) fwrite(res$sensitivity, "results/simulation_bootstrap_sensitivity.tsv", sep = "\t")
  written <- "results/simulation_performance.tsv, results/simulation_method_contrasts.tsv, results/simulation_bootstrap_sensitivity.tsv"
  if (!is.null(res$performance_lfmm)) {
    fwrite(res$performance_lfmm, "results/simulation_performance_lfmm.tsv", sep = "\t")
    fwrite(res$contrasts_lfmm, "results/simulation_lfmm_portability_contrast.tsv", sep = "\t")
    written <- paste0(written, ", results/simulation_performance_lfmm.tsv, results/simulation_lfmm_portability_contrast.tsv")
  }
  if (!is.null(res$performance_region)) {
    fwrite(res$performance_region, "results/simulation_performance_region.tsv", sep = "\t")
    fwrite(res$contrasts_region, "results/simulation_region_granularity_contrast.tsv", sep = "\t")
    written <- paste0(written, ", results/simulation_performance_region.tsv, results/simulation_region_granularity_contrast.tsv")
  }
  saveRDS(res, "results/simulation_summary_full.rds")
  say("\nwrote %s\n", written)
}
