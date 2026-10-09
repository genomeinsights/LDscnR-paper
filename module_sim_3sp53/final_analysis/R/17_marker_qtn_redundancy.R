## =============================================================================
## final_analysis/R/17_marker_qtn_redundancy.R
##
## PK's question (2026-09-21): what proportion of significant STAGE-1
## (non-clustered, individual-marker) SNPs are false positives vs. true
## positives, where "a QTN can only be counted once" -- i.e. if many
## significant markers are all in LD with the SAME single QTN, that should
## not count as many independent true positives.
##
## [!] Confirmed by code inspection before computing anything: the existing
## marker-level TP in 05_score_truth.R's .score() closure is NOT
## QTN-deduplicated -- `TP <- sum(truth_linked)` counts every significant
## marker linked to >=1 detectable QTN as its own TP, so N markers all
## tagging the same QTN contribute N to TP, not 1. Recall IS already
## QTN-deduplicated (`recovered_qtn <- unique(unlist(...))`), precision is
## not. This script quantifies the gap directly: for every significant
## Stage-1 marker (emmax_snp / lfmm_snp, i.e. the UNCLUSTERED, single-marker
## hypothesis -- not any Stage-2 assembled region), partition into:
##   - false positive: linked to no detectable QTN
##   - QTN-linked: linked to >=1 detectable QTN, further split into
##       - "credited" (dedup TP): the one marker counted per recovered QTN
##       - "redundant": every other significant marker tagging an
##         ALREADY-credited QTN (real signal, but not a distinct discovery)
##
## qtn_linked_fraction = n_QTN_linked / n_significant
## distinct_qtn_per_call = n_unique_QTN_recovered / n_significant
##
## Neither quantity is called marker precision. The first counts repeated
## marker support for the same QTN; the second is a distinct-causal-variant
## yield per marker call and is not a one-to-one precision estimate.
##
## DOWNSTREAM ONLY in spirit -- EMMAX/LFMM significance calls are read from
## the already-cached 03_emmax.R/04_lfmm.R output (no association rerun);
## the QTN-LD lookup (qtn_ld_table()) is the SAME cheap, deterministic
## computation 05_score_truth.R already does for its own scoring, just with
## per-marker linkage detail RETAINED here instead of only the aggregate
## count truth_scores.rds keeps. No Stage-1 clustering, no EMMAX/LFMM rerun.
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(parallel)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "helpers_stage2_truth.R"))   ## neutral_chromosomes()
say("=== 17_marker_qtn_redundancy ===\n\n")

C1_CELLS <- c("V0.5_c1", "V1_c1", "V2_c1")

redundancy_one_combo <- function(tag, cell, rep, env) {
  combo_id <- sprintf("%s_%s_rep%d_env%d", tag, cell, rep, env)
  em <- readRDS(file.path(stage_dir("03_emmax", combo_id), "emmax.rds"))
  lfmm_file <- file.path(stage_dir("04_lfmm", combo_id), "lfmm.rds")
  have_lfmm <- file.exists(lfmm_file)
  lf <- if (have_lfmm) readRDS(lfmm_file) else NULL
  b <- readRDS(file.path(stage_dir("02_build_ld_units", combo_id), "ld_units.rds"))
  GTs <- b$GTs; map <- copy(b$map)

  p <- colSums(GTs) / nrow(GTs) / 2
  map[, p_freq := p[match(marker, colnames(GTs))]]
  map[, Va := 2 * p_freq * (1 - p_freq) * allelic_values^2]
  map <- flag_true_qtns(map, va_col = "Va", maf_col = "MAF", maf_min = MAF_KEEP, p_va_min = VA_SHARE_DETECTABLE)
  detectable_qtn <- map[true_pos_QTN == TRUE, marker]
  thr <- score_thresholds(b$LD_decay$decay_sum, rho_r2 = TRUTH_RHO_R2, rho_d = TRUTH_RHO_D, dmax_cap = TRUTH_DMAX_CAP)
  qtn_lut_match <- if (length(detectable_qtn)) {
    qtn_lut <- qtn_ld_table(GTs, map, candidate_markers = map$marker, cores = 1)
    qtn_lut[r2 > thr$r2min & dist_bp < thr$dmax & qtn_marker %in% detectable_qtn]
  } else data.table(marker = character(), qtn_marker = character())

  .partition <- function(sig_markers) {
    n_sig <- length(sig_markers)
    if (n_sig == 0L) return(data.table(n_significant = 0L, n_qtn_linked_naive = 0L, n_false_positive = 0L,
                                       n_unique_qtn_recovered = 0L, n_redundant = 0L))
    linked <- lapply(sig_markers, function(mm) unique(qtn_lut_match[marker == mm, qtn_marker]))
    n_linked <- lengths(linked)
    qtn_linked_naive <- sum(n_linked > 0)
    fp <- sum(n_linked == 0)
    unique_qtn <- unique(unlist(linked[n_linked > 0]))
    data.table(n_significant = n_sig, n_qtn_linked_naive = qtn_linked_naive, n_false_positive = fp,
              n_unique_qtn_recovered = length(unique_qtn), n_redundant = qtn_linked_naive - length(unique_qtn))
  }

  ## Per-chromosome split (2026-10-09, additive columns): significant markers on
  ## a neutral chromosome (no QTN of non-zero effect) are false by construction.
  ntrl_markers <- map[as.character(Chr) %in% neutral_chromosomes(map), marker]
  .partition_split <- function(sig_markers) {
    base <- .partition(sig_markers)
    on_ntrl <- sig_markers %chin% ntrl_markers
    nt <- .partition(sig_markers[on_ntrl])
    cbind(base, n_significant_ntrl = nt$n_significant, n_qtn_linked_ntrl = nt$n_qtn_linked_naive,
          n_false_positive_ntrl = nt$n_false_positive)
  }
  rows <- list()
  rows[["emmax_snp"]] <- cbind(method = "emmax_snp", .partition_split(em$marker$marker[em$marker$significant_snp]))
  rows[["emmax_snp_nonsingleton"]] <- cbind(method = "emmax_snp_nonsingleton", .partition_split(em$marker$marker[em$marker$significant_snp_nonsingleton]))
  if (have_lfmm) rows[["lfmm_snp"]] <- cbind(method = "lfmm_snp", .partition_split(lf$marker$marker[lf$marker$significant_snp]))

  out <- rbindlist(rows)
  out[, `:=`(tag = tag, cell = cell, rep = rep, env = env, n_detectable_qtn = length(detectable_qtn))]
  out[]
}

if (sys.nframe() == 0L) {
  combos <- CJ(tag = TAGS_ALL, cell = CELLS_ALL, rep = REPS_ALL, env = ENVS_ALL)
  say("[1] marker-QTN redundancy for %d combos\n", nrow(combos))
  t0 <- Sys.time()
  res <- mclapply(seq_len(nrow(combos)), function(i) {
    tryCatch(list(ok = TRUE, data = redundancy_one_combo(combos$tag[i], combos$cell[i], combos$rep[i], combos$env[i])),
             error = function(e) {
               msg <- sprintf("[17_marker_qtn_redundancy] %s_%s_rep%d_env%d: %s",
                              combos$tag[i], combos$cell[i], combos$rep[i], combos$env[i], conditionMessage(e))
               message(msg); list(ok = FALSE, data = NULL, error = msg)
             })
  }, mc.cores = 12)
  ok_flags <- vapply(res, `[[`, logical(1), "ok")
  if (!all(ok_flags))
    stop(sprintf("%d/%d combos failed -- see [17_marker_qtn_redundancy] messages above", sum(!ok_flags), nrow(combos)))
  dt <- rbindlist(lapply(res, `[[`, "data"))
  say("[2] %d/%d combos OK in %.1f min\n", uniqueN(dt[, .(tag, cell, rep, env)]), nrow(combos), as.numeric(difftime(Sys.time(), t0, units = "mins")))

  fwrite(dt, "results/simulation_marker_qtn_redundancy_raw.tsv", sep = "\t")

  ## ---- pooled summary: ratios of pooled counts, full grid and c=1-only ----
  .pool <- function(d, label) {
    p <- d[, .(n_combo = uniqueN(paste(tag, cell, rep, env)),
              n_significant = sum(n_significant), n_qtn_linked_naive = sum(n_qtn_linked_naive),
              n_false_positive = sum(n_false_positive), n_unique_qtn_recovered = sum(n_unique_qtn_recovered),
              n_redundant = sum(n_redundant)), by = method]
    p[, `:=`(scope = label,
            qtn_linked_fraction = n_qtn_linked_naive / pmax(n_significant, 1),
            distinct_qtn_per_call = n_unique_qtn_recovered / pmax(n_significant, 1),
            frac_redundant_among_sig = n_redundant / pmax(n_significant, 1),
            frac_redundant_among_qtn_linked = n_redundant / pmax(n_qtn_linked_naive, 1),
            mean_markers_per_recovered_qtn = n_qtn_linked_naive / pmax(n_unique_qtn_recovered, 1))]
    p[]
  }
  summary_all <- .pool(dt, "full_grid_1400")
  summary_c1 <- .pool(dt[cell %chin% C1_CELLS], "c1_primary_600")
  summary_dt <- rbindlist(list(summary_all, summary_c1))
  setcolorder(summary_dt, c("scope", "method"))
  fwrite(summary_dt, "results/simulation_marker_qtn_redundancy_summary.tsv", sep = "\t")

  ## ---- explicit singleton contrast for primary c=1 simulations -----------
  ## emmax_snp_nonsingleton is a subset of emmax_snp under the SAME original
  ## marker-wise BH decision, so subtraction isolates significant markers in
  ## singleton Stage-1 units without rerunning multiple-testing correction.
  c1 <- dt[cell %chin% C1_CELLS]
  all_m <- c1[method == "emmax_snp",
              .(tag, cell, rep, env, n_all = n_significant,
                linked_all = n_qtn_linked_naive, qtn_all = n_unique_qtn_recovered)]
  multi_m <- c1[method == "emmax_snp_nonsingleton",
                .(tag, cell, rep, env, n_multi = n_significant,
                  linked_multi = n_qtn_linked_naive, qtn_multi = n_unique_qtn_recovered)]
  singleton_raw <- merge(all_m, multi_m, by = c("tag", "cell", "rep", "env"))
  singleton_raw[, `:=`(n_singleton = n_all - n_multi,
                       linked_singleton = linked_all - linked_multi,
                       qtn_recovered_only_by_singletons = qtn_all - qtn_multi)]
  singleton_raw[, `:=`(unmatched_singleton = n_singleton - linked_singleton,
                       unmatched_multi = n_multi - linked_multi)]
  fwrite(singleton_raw, "results/simulation_marker_qtn_singleton_raw.tsv", sep = "\t")

  singleton_summary <- rbindlist(list(
    singleton_raw[, .(group = "all significant EMMAX markers",
                      n_significant = sum(n_all), n_qtn_linked = sum(linked_all),
                      n_unmatched = sum(n_all - linked_all),
                      n_distinct_qtn_recovered = sum(qtn_all))],
    singleton_raw[, .(group = "markers in Stage-1 units >=2",
                      n_significant = sum(n_multi), n_qtn_linked = sum(linked_multi),
                      n_unmatched = sum(unmatched_multi),
                      n_distinct_qtn_recovered = sum(qtn_multi))],
    singleton_raw[, .(group = "markers in singleton Stage-1 units",
                      n_significant = sum(n_singleton), n_qtn_linked = sum(linked_singleton),
                      n_unmatched = sum(unmatched_singleton),
                      n_distinct_qtn_recovered = sum(qtn_recovered_only_by_singletons))]
  ), fill = TRUE)
  singleton_summary[, `:=`(
    qtn_linked_fraction = n_qtn_linked / pmax(n_significant, 1),
    unmatched_fraction = n_unmatched / pmax(n_significant, 1))]
  fwrite(singleton_summary, "results/simulation_marker_qtn_singleton_summary.tsv", sep = "\t")

  set.seed(170921L)
  map_ids <- sort(unique(singleton_raw$rep))
  boot_difference <- replicate(2000L, {
    sampled <- sample(map_ids, length(map_ids), replace = TRUE)
    z <- rbindlist(lapply(sampled, function(id) singleton_raw[rep == id]))
    sum(z$linked_multi) / sum(z$n_multi) -
      sum(z$linked_singleton) / sum(z$n_singleton)
  })
  singleton_contrast <- data.table(
    contrast = "QTN-linked fraction: units >=2 minus singleton units",
    estimate = singleton_summary[group == "markers in Stage-1 units >=2", qtn_linked_fraction] -
      singleton_summary[group == "markers in singleton Stage-1 units", qtn_linked_fraction],
    lo = unname(quantile(boot_difference, 0.025)),
    hi = unname(quantile(boot_difference, 0.975)),
    n_boot = 2000L)
  fwrite(singleton_contrast, "results/simulation_marker_qtn_singleton_contrast.tsv", sep = "\t")
  say("\n[3] wrote results/simulation_marker_qtn_redundancy_raw.tsv, results/simulation_marker_qtn_redundancy_summary.tsv\n")
  print(summary_dt[, .(scope, method, n_significant, n_qtn_linked_naive, n_unique_qtn_recovered, n_redundant,
                       qtn_linked_fraction, distinct_qtn_per_call, mean_markers_per_recovered_qtn)])
  print(singleton_summary)
  print(singleton_contrast)

  write_receipt("17_marker_qtn_redundancy", inputs = character(),
                params = list(rho_r2 = TRUTH_RHO_R2, rho_d = TRUTH_RHO_D, dmax_cap = TRUTH_DMAX_CAP,
                              c1_cells = C1_CELLS),
                outputs = c("results/simulation_marker_qtn_redundancy_raw.tsv",
                            "results/simulation_marker_qtn_redundancy_summary.tsv",
                            "results/simulation_marker_qtn_singleton_raw.tsv",
                            "results/simulation_marker_qtn_singleton_summary.tsv",
                            "results/simulation_marker_qtn_singleton_contrast.tsv"))
  cat("MARKER_QTN_REDUNDANCY_DONE\n")
}
