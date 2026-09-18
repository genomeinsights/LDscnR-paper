## =============================================================================
## final_analysis/R/14_summarise_additional_analyses.R
##
## Analysis 2: among reported Stage-2 regions, are regions supported by fewer
## markers or fewer Stage-1 units more likely to be false positives?
##
## Unit of analysis: one reported Stage-2 region (results/
## simulation_stage2_region_details.rds, R/11) -- NOT a repeat of the
## Stage-1-unit-size analysis module_sim_3sp53/R/12+R/14 already did (see
## ADDITIONAL_ANALYSES_AUDIT.md: that pair used the pre-correction parser
## AND scored at the unit scale, not the reported-region scale).
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(parallel)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "helpers_stage2_truth.R"))
say("=== 14_summarise_additional_analyses ===\n\n")

REGION_METHODS <- c("emmax_snp_region", "emmax_simes_region", "emmax_consensus_region",
                    "lfmm_snp_region", "lfmm_simes_region")

## ---- load region details (R/11) ---------------------------------------------
region_detail_file <- "results/simulation_stage2_region_details.rds"
if (!file.exists(region_detail_file))
  stop("R/11_stage2_region_details.R has not produced: ", region_detail_file, " -- run it first.")
dt <- readRDS(region_detail_file)
dt[, FP := !TP]

## ---- size bins, chosen from the observed distribution (see ---------------
## ADDITIONAL_ANALYSES_AUDIT.md for the counts that justified these edges):
## n_markers is heavily right-skewed (median 6-8, 90th pctile ~30, max 451);
## n_units is extremely concentrated (82-83% of ALL regions are a single
## Stage-1 unit, max 15) -- bins are finer at the small end for n_units.
MARKER_BIN_EDGES <- c(0, 1, 2, 5, 10, 20, 50, Inf)
MARKER_BIN_LABELS <- c("1", "2", "3-5", "6-10", "11-20", "21-50", "51+")
UNIT_BIN_EDGES <- c(0, 1, 2, 4, 9, Inf)
UNIT_BIN_LABELS <- c("1", "2", "3-4", "5-9", "10+")
dt[, marker_bin := cut(n_markers, MARKER_BIN_EDGES, labels = MARKER_BIN_LABELS)]
dt[, unit_bin := cut(n_units, UNIT_BIN_EDGES, labels = UNIT_BIN_LABELS)]

## =============================================================================
## [1] Descriptive: FP proportion by size bin, per method x BGS, with
## map-cluster bootstrap CIs (rep-level counts, same shared primitive as
## R/12/R/13).
## =============================================================================
say("[1] descriptive: FP proportion by size bin\n")

.descriptive_by <- function(size_col_bin) {
  rows <- list()
  for (m in REGION_METHODS) for (tg in TAGS_ALL) {
    sub <- dt[method == m & tag == tg]
    if (!nrow(sub)) next
    for (bn in levels(sub[[size_col_bin]])) {
      b <- sub[get(size_col_bin) == bn]
      n_tot <- nrow(b)
      if (!n_tot) { rows[[length(rows) + 1]] <- data.table(method = m, tag = tg, bin = bn, n_regions = 0L,
                                                           n_TP = 0L, n_FP = 0L, fp_prop = NA_real_,
                                                           fp_prop_ci_lo = NA_real_, fp_prop_ci_hi = NA_real_,
                                                           median_size_TP = NA_real_, median_size_FP = NA_real_); next }
      rep_dt <- b[, .(TP = sum(TP), FP = sum(FP)), by = rep][match(REPS_ALL, rep)]
      rep_dt[is.na(TP), `:=`(TP = 0L, FP = 0L)]
      bs <- bootstrap_rep_matrix(as.matrix(rep_dt[, .(TP, FP)]), N_BOOTSTRAP, SEEDS[["bootstrap"]])
      fp_b <- bs[, "FP"] / pmax(bs[, "TP"] + bs[, "FP"], 1)
      size_col <- if (size_col_bin == "marker_bin") "n_markers" else "n_units"
      rows[[length(rows) + 1]] <- data.table(
        method = m, tag = tg, bin = bn, n_regions = n_tot, n_TP = sum(b$TP), n_FP = sum(b$FP),
        fp_prop = sum(b$FP) / n_tot, fp_prop_ci_lo = ci_quantile(fp_b)[1], fp_prop_ci_hi = ci_quantile(fp_b)[2],
        median_size_TP = stats::median(b[[size_col]][b$TP]), median_size_FP = stats::median(b[[size_col]][b$FP]))
    }
  }
  d <- rbindlist(rows)
  d[, bin := factor(bin, levels = if (size_col_bin == "marker_bin") MARKER_BIN_LABELS else UNIT_BIN_LABELS)]
  setorder(d, method, tag, bin)
  d[]
}
desc_markers <- .descriptive_by("marker_bin")[, size_measure := "n_markers"]
desc_units <- .descriptive_by("unit_bin")[, size_measure := "n_units"]
desc <- rbind(desc_markers, desc_units)
fwrite(desc, "results/simulation_stage2_size_truth.tsv", sep = "\t")
say("[2] wrote results/simulation_stage2_size_truth.tsv (%d rows)\n", nrow(desc))

## =============================================================================
## [2] Model: FP ~ log2(size) + method + cell + BGS, marker-count and
## unit-count fitted separately. Map/burn-in-clustered uncertainty via the
## bootstrap primitive (refit-in-bootstrap-replicate, not a new
## clustered-SE dependency).
## =============================================================================
say("\n[3] logistic model: FP ~ log2(size) + method + cell + BGS\n")

.fit_size_model <- function(size_col) {
  d <- copy(dt)[n_markers > 0]
  d[, log2_size := log2(get(size_col))]
  d[, method := factor(method, levels = REGION_METHODS)]
  d[, cell := factor(cell, levels = CELLS_ALL)]
  d[, tag := factor(tag, levels = TAGS_ALL)]
  fit <- stats::glm(FP ~ log2_size + method + cell + tag, data = d, family = stats::binomial())
  co <- summary(fit)$coefficients
  point <- data.table(term = rownames(co), estimate = co[, "Estimate"])

  ## map/burn-in cluster bootstrap: refit on data resampled by (cell,rep)
  ## cluster. [!] Deliberately NOT (tag,cell,rep) (PK, 2026-09-18 review):
  ## BGS and no-BGS simulations are PAIRED by parameter combination, map and
  ## environment and share a deterministic seed (materials_and_methods.tex) --
  ## they are not independent replicates of each other. Clustering by
  ## (tag,cell,rep) would let a bootstrap draw select one treatment's rows
  ## for a given (cell,rep) without its paired treatment, breaking that
  ## pairing silently despite this file's own header claiming methods/BGS
  ## stay paired. (cell,rep) as the cluster means both tag rows for a
  ## resampled (cell,rep) always move together.
  clusters <- unique(d[, .(cell, rep)])
  n_cl <- nrow(clusters)
  set.seed(SEEDS[["bootstrap"]])
  boot_coef <- matrix(NA_real_, nrow = N_BOOTSTRAP, ncol = nrow(point), dimnames = list(NULL, point$term))
  d[, cl_id := .GRP, by = .(cell, rep)]
  cl_rows <- split(seq_len(nrow(d)), d$cl_id)   ## computed ONCE, not per bootstrap replicate
  for (bb in seq_len(N_BOOTSTRAP)) {
    draw <- sample.int(n_cl, n_cl, replace = TRUE)
    idx <- unlist(cl_rows[as.character(draw)], use.names = FALSE)
    fit_b <- tryCatch(stats::glm(FP ~ log2_size + method + cell + tag, data = d[idx], family = stats::binomial()),
                      error = function(e) NULL, warning = function(w) NULL)
    if (!is.null(fit_b)) {
      cb <- coef(fit_b)
      boot_coef[bb, names(cb)[names(cb) %in% colnames(boot_coef)]] <- cb[names(cb) %in% colnames(boot_coef)]
    }
  }
  point[, `:=`(ci_lo = apply(boot_coef, 2, ci_quantile)[1, ][term],
              ci_hi = apply(boot_coef, 2, ci_quantile)[2, ][term],
              n_boot_converged = colSums(!is.na(boot_coef))[term])]
  list(point = point, n_obs = nrow(d), n_events = sum(d$FP), aic = AIC(fit))
}

model_markers <- .fit_size_model("n_markers")
model_units <- .fit_size_model("n_units")
models <- rbind(model_markers$point[, size_measure := "n_markers"],
                model_units$point[, size_measure := "n_units"])
fwrite(models, "results/simulation_stage2_size_models.tsv", sep = "\t")
say("[4] wrote results/simulation_stage2_size_models.tsv\n")
say("    n_markers model: n=%d, FP events=%d, AIC=%.1f\n", model_markers$n_obs, model_markers$n_events, model_markers$aic)
say("    n_units   model: n=%d, FP events=%d, AIC=%.1f\n", model_units$n_obs, model_units$n_events, model_units$aic)
print(models[term == "log2_size"])

## ---- collinearity check: span and marker density as sensitivity covariates --
say("\n[5] collinearity: log2(n_markers), log2(span+1), log2(density) -- manual VIF (no new dependency)\n")
cov_dt <- dt[n_markers > 0][, `:=`(span = to - from + 1, density = n_markers / pmax(to - from + 1, 1))]
cov_dt <- cov_dt[, .(log2_markers = log2(n_markers), log2_span = log2(span + 1), log2_density = log2(density + 1e-6))]
.vif <- function(d) {
  vapply(names(d), function(v) {
    f <- stats::as.formula(sprintf("%s ~ .", v))
    r2 <- summary(stats::lm(f, data = d))$r.squared
    1 / (1 - r2)
  }, numeric(1))
}
vif_vals <- .vif(cov_dt)
print(vif_vals)
say("    VIF > 5 flags unstable coefficients for that covariate; span/density are NOT included in the\n")
say("    primary log2_size model above for exactly this reason if flagged -- see audit report.\n")

## =============================================================================
## [3] Opportunity-effect sensitivity: is the observed TP-rate-by-size
## gradient stronger than pure genomic coverage would achieve by chance?
## [!] REDESIGNED (PK, 2026-09-18 review, point 2). The previous version drew
## K=30 GENERIC contiguous windows per (combo, chromosome, SIZE BIN) at the
## bin's MEDIAN size -- weighted per stratum, not per region -- and compared
## just two point-estimate slopes with no uncertainty on their difference.
## This version relocates EVERY OBSERVED REGION INDIVIDUALLY: each region's
## own exact constituent-marker INDEX-SPACING PATTERN on its own chromosome
## (not merely its marker count) is preserved via a circular shift of that
## exact pattern to a uniformly random start point on the SAME chromosome
## (standard circular-permutation / GAT-style relocation null; wraps around
## so every shift is equally valid regardless of chromosome-end proximity --
## no separate "near the end" fallback is needed). Each region contributes
## exactly one relocation draw per replicate, repeated for R_RELOC
## independent replicates, giving a FULL null distribution of the
## opportunity-only TP-rate-vs-size slope -- so the slope DIFFERENCE
## (observed - opportunity-null) gets an empirical CI and p-value, not a
## single point comparison against one arbitrarily-drawn null slope.
## =============================================================================
say("\n[6] opportunity-effect sensitivity (exact per-region circular-shift relocation null)\n")
R_RELOC <- 200L
bin_size <- c("1" = 1, "2" = 2, "3-5" = 4, "6-10" = 8, "11-20" = 15, "21-50" = 35, "51+" = 100)

## one combo -> a FIXED-SIZE (length(MARKER_BIN_LABELS) x R_RELOC) long table
## of (marker_bin, replicate, n_hits, n_draws), regardless of how many
## regions that combo has -- keeps the cross-combo rbind small and the
## per-region relocation detail local to this function.
## [!] Safe deterministic seed (PK, 2026-09-18 second review, point 2a):
## strtoi() on a full 8-hex-digit CRC32 string overflows R's signed 32-bit
## integer for any hash >= 0x80000000 -- roughly HALF of all values -- and
## silently returns NA. NA then propagates through the seed arithmetic to
## set.seed(NA), which throws "supplied seed is not a valid integer,"
## caught by the caller's tryCatch and silently dropped as a failed combo.
## Confirmed in the overnight run's log: 619/1,223 combos (50.6%) lost this
## way, while `observed_by_bin` (used for obs_slope) was still built from
## the FULL, unfiltered grid -- comparing a full-grid observed slope against
## a ~half-grid null slope. Fixed by splitting the 8-hex-digit hash into two
## 4-hex-digit halves (each safely < 2^16) and combining them with DOUBLE
## arithmetic (the literal 65536 is a double, so `strtoi(...) * 65536`
## promotes before it could overflow), then taking %% within double
## precision before ever coercing to R's 32-bit integer for set.seed().
.safe_combo_seed <- function(combo_id, base_seed) {
  hex <- digest::digest(combo_id, algo = "crc32")
  crc_val <- strtoi(substr(hex, 1, 4), base = 16L) * 65536 + strtoi(substr(hex, 5, 8), base = 16L)
  (base_seed + crc_val) %% (.Machine$integer.max - 1L)
}

opportunity_one_combo <- function(tag, cell, rep, env) {
  combo_id <- sprintf("%s_%s_rep%d_env%d", tag, cell, rep, env)
  set.seed(.safe_combo_seed(combo_id, SEEDS[["bootstrap"]]))
  b <- readRDS(file.path(stage_dir("02_build_ld_units", combo_id), "ld_units.rds"))
  GTs <- b$GTs; map <- copy(b$map)
  p <- colSums(GTs) / nrow(GTs) / 2
  map[, p_freq := p[match(marker, colnames(GTs))]]
  map[, Va := 2 * p_freq * (1 - p_freq) * allelic_values^2]
  map <- flag_true_qtns(map, va_col = "Va", maf_col = "MAF", maf_min = MAF_KEEP, p_va_min = VA_SHARE_DETECTABLE)
  detectable_qtn <- map[true_pos_QTN == TRUE, marker]
  ## [!] Zero-opportunity combos must remain, not be dropped (PK, 2026-09-18
  ## second review, point 2b): a combo with no detectable QTN at all, or no
  ## marker passing the LD/distance criteria, has observed regions that are
  ## necessarily false positives. `is_qtn_linked` staying EMPTY (rather than
  ## returning NULL here) makes every relocation draw below correctly score
  ## as a miss (n_hits=0) -- a real, informative zero contribution to the
  ## null denominator -- instead of silently removing those regions from
  ## the comparison while they remain in the observed-side totals.
  is_qtn_linked <- character(0)
  if (length(detectable_qtn)) {
    thr <- score_thresholds(b$LD_decay$decay_sum, rho_r2 = TRUTH_RHO_R2, rho_d = TRUTH_RHO_D, dmax_cap = TRUTH_DMAX_CAP)
    qtn_lut <- qtn_ld_table(GTs, map, candidate_markers = map$marker, cores = 1)
    qtn_lut_match <- qtn_lut[r2 > thr$r2min & dist_bp < thr$dmax & qtn_marker %in% detectable_qtn]
    is_qtn_linked <- unique(qtn_lut_match$marker)
  }

  ## copy args to distinctly-named locals before filtering `dt` -- `dt` has
  ## its own tag/cell/rep/env COLUMNS with the same names as this function's
  ## arguments; ..env alone silently failed data.table's ..var lookup here
  ## previously (root cause not chased down), so this sidesteps it entirely,
  ## same fix as R/12's `this_floor`.
  q_tag <- tag; q_cell <- cell; q_rep <- rep; q_env <- env
  sub <- dt[tag == q_tag & cell == q_cell & rep == q_rep & env == q_env & method %in% REGION_METHODS]
  if (!nrow(sub)) return(NULL)   ## genuinely no regions for this combo -- nothing to relocate

  setkey(map, Chr)
  chr_cache <- new.env(parent = emptyenv())
  get_chr_markers <- function(ch) {
    key <- as.character(ch)
    v <- chr_cache[[key]]
    if (is.null(v)) { v <- map[Chr == ch][order(Pos), marker]; chr_cache[[key]] <- v }
    v
  }

  hit_mat <- matrix(0L, nrow = nrow(sub), ncol = R_RELOC)
  for (i in seq_len(nrow(sub))) {
    chr_markers <- get_chr_markers(sub$Chr[i])
    n_chr <- length(chr_markers)
    idx <- match(sub$member_markers[[i]], chr_markers)
    idx <- idx[!is.na(idx)]
    if (!length(idx) || n_chr < 2) next
    delta <- idx - min(idx)   ## this region's own exact index-spacing pattern, anchored at 0
    shifts <- sample.int(n_chr, R_RELOC, replace = TRUE) - 1L
    hit_mat[i, ] <- vapply(shifts, function(sh) {
      reloc_idx <- ((delta + sh) %% n_chr) + 1L   ## circular wraparound
      as.integer(any(chr_markers[reloc_idx] %chin% is_qtn_linked))
    }, integer(1))
  }

  bin_of <- factor(as.character(sub$marker_bin), levels = MARKER_BIN_LABELS)
  agg <- rowsum(hit_mat, group = bin_of, reorder = FALSE)   ## n_present_bins x R_RELOC
  n_draws_by_bin <- table(bin_of)
  out_hits <- matrix(0L, nrow = length(MARKER_BIN_LABELS), ncol = R_RELOC,
                     dimnames = list(MARKER_BIN_LABELS, NULL))
  out_hits[rownames(agg), ] <- agg
  out_draws <- stats::setNames(rep(0L, length(MARKER_BIN_LABELS)), MARKER_BIN_LABELS)
  out_draws[names(n_draws_by_bin)] <- as.integer(n_draws_by_bin)

  data.table(marker_bin = rep(MARKER_BIN_LABELS, R_RELOC),
            replicate = rep(seq_len(R_RELOC), each = length(MARKER_BIN_LABELS)),
            n_hits = as.integer(out_hits), n_draws = rep(out_draws, R_RELOC))
}

combos_present <- unique(dt[method %in% REGION_METHODS, .(tag, cell, rep, env)])
t0 <- Sys.time()
opp_res <- mclapply(seq_len(nrow(combos_present)), function(i) {
  tryCatch(opportunity_one_combo(combos_present$tag[i], combos_present$cell[i],
                                 combos_present$rep[i], combos_present$env[i]),
           error = function(e) { message(sprintf("[14 opportunity] %s_%s_rep%d_env%d: %s",
                                                 combos_present$tag[i], combos_present$cell[i],
                                                 combos_present$rep[i], combos_present$env[i], conditionMessage(e))); NULL })
}, mc.cores = 12)   ## bumped from 7 (PK, 2026-09-18): mini has 14 physical cores, 2 left for the system
opp_dt <- rbindlist(Filter(Negate(is.null), opp_res))
## [!] Hard completeness stop (PK, 2026-09-18 second review): each
## successful combo contributes exactly length(MARKER_BIN_LABELS) * R_RELOC
## rows (a fixed shape, by construction of opportunity_one_combo()), so the
## number of combos actually represented in opp_dt is recoverable and
## checkable against combos_present. The previous run silently pooled a
## ~half-complete opp_dt as if it were the full grid (see the seed-overflow
## fix above); this now halts instead. A combo legitimately contributes
## zero HITS (via the empty-is_qtn_linked path above) but never legitimately
## contributes zero ROWS while being in combos_present -- every combo in
## combos_present has >=1 region row in `dt` by construction, so `sub` is
## never empty for it and the function cannot fall through to its `NULL`
## "no regions" branch.
n_combo_rows <- length(MARKER_BIN_LABELS) * R_RELOC
n_combos_succeeded <- nrow(opp_dt) / n_combo_rows
if (n_combos_succeeded != nrow(combos_present))
  stop(sprintf("opportunity relocation incomplete: %s/%d combos succeeded -- see [14 opportunity] messages above for which ones failed; refusing to pool a partial null against the full observed grid",
               format(n_combos_succeeded), nrow(combos_present)))
say("[7] opportunity-null relocations drawn for %d/%d combos (R=%d replicates/region) in %.1f min\n",
    n_combos_succeeded, nrow(combos_present), R_RELOC, as.numeric(difftime(Sys.time(), t0, units = "mins")))

if (nrow(opp_dt)) {
  ## pool hit/draw counts across ALL combos, keeping (marker_bin, replicate)
  ## separate -- this is what lets us fit one slope PER REPLICATE below.
  opp_pooled_r <- opp_dt[, .(n_hits = sum(n_hits), n_draws = sum(n_draws)), by = .(marker_bin, replicate)]
  opp_pooled_r[, log2_size := log2(bin_size[as.character(marker_bin)])]

  ## null distribution of the opportunity-only slope: one slope per replicate.
  null_slopes <- vapply(seq_len(R_RELOC), function(r) {
    dd <- opp_pooled_r[replicate == r & n_draws > 0]
    if (nrow(dd) < 3) return(NA_real_)
    tryCatch(coef(glm(cbind(n_hits, n_draws - n_hits) ~ log2_size, data = dd, family = binomial()))[2],
             error = function(e) NA_real_)
  }, numeric(1))
  null_slopes <- null_slopes[is.finite(null_slopes)]

  observed_by_bin <- dt[method %in% REGION_METHODS, .(n_regions = .N, n_TP = sum(TP)), by = marker_bin]
  observed_by_bin[, observed_tp_rate := n_TP / n_regions]
  observed_by_bin[, marker_bin := factor(marker_bin, levels = MARKER_BIN_LABELS)]
  observed_by_bin[, log2_size := log2(bin_size[as.character(marker_bin)])]
  setorder(observed_by_bin, marker_bin)
  obs_slope <- tryCatch(coef(glm(cbind(n_TP, n_regions - n_TP) ~ log2_size, data = observed_by_bin, family = binomial()))[2],
                        error = function(e) NA_real_)

  ## descriptive table: mean opportunity hit-rate per bin, pooled over ALL
  ## replicates (same shape as the old opp_compare TSV, now built from exact
  ## per-region relocations rather than generic bin-median windows).
  opp_mean <- opp_pooled_r[, .(n_hits = sum(n_hits), n_draws = sum(n_draws)), by = marker_bin]
  opp_mean[, opportunity_hit_rate := n_hits / n_draws]
  opp_mean[, marker_bin := factor(marker_bin, levels = MARKER_BIN_LABELS)]
  setorder(opp_mean, marker_bin)
  opp_compare <- merge(observed_by_bin[, .(marker_bin, n_regions, n_TP, observed_tp_rate)],
                       opp_mean[, .(marker_bin, opportunity_hit_rate)], by = "marker_bin", all = TRUE)
  opp_compare[, excess_ratio := observed_tp_rate / opportunity_hit_rate]
  fwrite(opp_compare, "results/simulation_stage2_opportunity_null.tsv", sep = "\t")
  say("[8] wrote results/simulation_stage2_opportunity_null.tsv\n")
  print(opp_compare)
  say("    excess ratio (observed/opportunity) by bin: %s\n",
      paste(sprintf("%s=%.1fx", opp_compare$marker_bin, opp_compare$excess_ratio), collapse = ", "))

  ## slope-DIFFERENCE uncertainty from the full null distribution (point 2):
  ## obs_slope is fixed (real data, not itself resampled here); null_slopes
  ## varies across the R_RELOC relocation replicates, so diff_dist is the
  ## empirical distribution of (observed - opportunity-null) under repeated
  ## relocation.
  ## [!] Terminology (PK, 2026-09-18 second review): the *_ci_lo/_ci_hi
  ## columns below (kept as-is for continuity with ci_quantile()'s naming
  ## used elsewhere in this pipeline for map-cluster BOOTSTRAP intervals)
  ## are NOT a sampling-uncertainty confidence interval across simulated
  ## maps/replicates -- obs_slope is a single fixed value computed once from
  ## the real, observed region set. They describe variation among R_RELOC
  ## independent relocations OF THAT SAME fixed observed set, i.e. a 95%
  ## relocation interval conditional on the observed regions, not a CI on
  ## the observed slope itself. Call it a "relocation interval" in any
  ## write-up, never a generic "confidence interval."
  diff_dist <- obs_slope - null_slopes
  slope_diff_summary <- data.table(
    obs_slope = obs_slope, n_reloc_replicates = length(null_slopes),
    null_slope_mean = mean(null_slopes), null_slope_median = stats::median(null_slopes),
    null_slope_ci_lo = ci_quantile(null_slopes)[1], null_slope_ci_hi = ci_quantile(null_slopes)[2],
    slope_diff = obs_slope - mean(null_slopes),
    slope_diff_ci_lo = ci_quantile(diff_dist)[1], slope_diff_ci_hi = ci_quantile(diff_dist)[2],
    p_obs_le_null = mean(null_slopes >= obs_slope))
  fwrite(slope_diff_summary, "results/simulation_stage2_opportunity_slope.tsv", sep = "\t")
  say("[9] wrote results/simulation_stage2_opportunity_slope.tsv\n")
  print(slope_diff_summary)
  say("\n    observed log-odds slope vs log2(size): %.4f ; opportunity-null slope: %.4f [95%% relocation interval %.4f, %.4f] (%d replicates)\n",
      slope_diff_summary$obs_slope, slope_diff_summary$null_slope_mean,
      slope_diff_summary$null_slope_ci_lo, slope_diff_summary$null_slope_ci_hi, slope_diff_summary$n_reloc_replicates)
  say("    slope difference (observed - opportunity-null): %.4f [95%% relocation interval %.4f, %.4f], P(null slope >= observed) = %.3f\n",
      slope_diff_summary$slope_diff, slope_diff_summary$slope_diff_ci_lo, slope_diff_summary$slope_diff_ci_hi,
      slope_diff_summary$p_obs_le_null)
  ## [!] Three-way classification (PK, 2026-09-18 second review, spotted while
  ## inspecting the validation run): the previous version only branched on
  ## "CI excludes zero, positive" vs. an else covering both "includes zero"
  ## AND "excludes zero, negative" -- so a significantly SHALLOWER-than-null
  ## observed slope (exactly this run's result: CI entirely negative) was
  ## misreported as "does NOT exclude zero." Numeric outputs (the TSV) were
  ## never affected -- only this printed narrative.
  if (is.na(obs_slope) || !length(null_slopes)) {
    say("    slope comparison unavailable -- see the ratio-by-bin line above instead.\n")
  } else if (slope_diff_summary$slope_diff_ci_lo > 0) {
    say("    the slope-difference CI excludes zero (entirely positive): evidence the size-truth relationship is\n")
    say("    STRONGER than pure genomic-coverage opportunity predicts.\n")
  } else if (slope_diff_summary$slope_diff_ci_hi < 0) {
    say("    the slope-difference CI excludes zero (entirely negative): the size gradient is SHALLOWER than pure\n")
    say("    genomic-coverage opportunity predicts, even though observed TP rate is substantially higher than\n")
    say("    opportunity at every bin (a large but size-DEcreasing excess, per the ratio line above) -- do not\n")
    say("    claim the size gradient itself is selection-driven; if anything this argues the opposite direction.\n")
    say("    See ADDITIONAL_ANALYSES_AUDIT.md.\n")
  } else {
    say("    the slope-difference CI includes zero: no evidence the size gradient's steepness differs from pure\n")
    say("    genomic-coverage opportunity, in either direction.\n")
  }
} else {
  say("[8] opportunity-null produced no usable draws -- falling back to descriptive-only reporting\n")
  say("    per the spec's own escape hatch; do not give a mechanistic interpretation of the size effect.\n")
}

write_receipt("14_summarise_additional_analyses", inputs = region_detail_file,
              params = list(marker_bin_edges = MARKER_BIN_EDGES, unit_bin_edges = UNIT_BIN_EDGES,
                            r_reloc = R_RELOC, n_bootstrap = N_BOOTSTRAP, seed_bootstrap = SEEDS[["bootstrap"]]),
              outputs = c("results/simulation_stage2_size_truth.tsv", "results/simulation_stage2_size_models.tsv",
                         "results/simulation_stage2_opportunity_null.tsv", "results/simulation_stage2_opportunity_slope.tsv"))
say("\nDONE\n")
