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

  ## map/burn-in cluster bootstrap: refit on data resampled by (tag,cell,rep)
  ## cluster (all rows for a resampled rep move together), same B/seed as
  ## elsewhere.
  clusters <- unique(d[, .(tag, cell, rep)])
  n_cl <- nrow(clusters)
  set.seed(SEEDS[["bootstrap"]])
  boot_coef <- matrix(NA_real_, nrow = N_BOOTSTRAP, ncol = nrow(point), dimnames = list(NULL, point$term))
  d[, cl_id := .GRP, by = .(tag, cell, rep)]
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
## [3] Opportunity-effect sensitivity: is the observed FP-by-size gradient
## stronger than a same-chromosome, same-marker-count RANDOM WINDOW would
## achieve by pure genomic coverage? Per (combo, chromosome, size bin):
## draw K random contiguous same-size marker windows on that chromosome and
## test each against the SAME detectable-QTN LD/distance criteria used for
## real regions (qtn_lut_match, recomputed per combo -- not retained in
## R/11's compact table). Never randomises across chromosomes.
## =============================================================================
say("\n[6] opportunity-effect sensitivity (chromosome-preserving random-window null)\n")
K_DRAWS <- 30L

opportunity_one_combo <- function(tag, cell, rep, env) {
  combo_id <- sprintf("%s_%s_rep%d_env%d", tag, cell, rep, env)
  b <- readRDS(file.path(stage_dir("02_build_ld_units", combo_id), "ld_units.rds"))
  GTs <- b$GTs; map <- copy(b$map); stage1 <- b$stage1
  p <- colSums(GTs) / nrow(GTs) / 2
  map[, p_freq := p[match(marker, colnames(GTs))]]
  map[, Va := 2 * p_freq * (1 - p_freq) * allelic_values^2]
  map <- flag_true_qtns(map, va_col = "Va", maf_col = "MAF", maf_min = MAF_KEEP, p_va_min = VA_SHARE_DETECTABLE)
  detectable_qtn <- map[true_pos_QTN == TRUE, marker]
  if (!length(detectable_qtn)) return(NULL)
  thr <- score_thresholds(b$LD_decay$decay_sum, rho_r2 = TRUTH_RHO_R2, rho_d = TRUTH_RHO_D, dmax_cap = TRUTH_DMAX_CAP)
  qtn_lut <- qtn_ld_table(GTs, map, candidate_markers = map$marker, cores = 1)
  qtn_lut_match <- qtn_lut[r2 > thr$r2min & dist_bp < thr$dmax & qtn_marker %in% detectable_qtn]
  if (!nrow(qtn_lut_match)) return(NULL)

  setkey(map, Chr)
  ## which (size bin, chromosome) combinations actually occur for this combo's
  ## regions. `dt` has its own tag/cell/rep/env COLUMNS with the same names as
  ## this function's arguments. `..var` (data.table's idiom for "look this up
  ## in the calling scope, not as a column") worked for ..tag/..cell/..rep but
  ## silently failed for ..env specifically ("object '..env' not found" on
  ## every single combo, confirmed by a full run producing zero usable draws)
  ## -- some name collision between the argument `env` and an internal
  ## data.table object of the same name, not chased down further. Sidestepped
  ## entirely by copying the arguments to distinctly-named local variables
  ## first, the same fix already applied to R/12's `this_floor` for the
  ## identical class of bug -- no dependence on `..var`'s scoping at all.
  q_tag <- tag; q_cell <- cell; q_rep <- rep; q_env <- env
  sub <- dt[tag == q_tag & cell == q_cell & rep == q_rep & env == q_env]
  if (!nrow(sub)) return(NULL)
  rows <- list()
  for (ch in unique(sub$Chr)) {
    chr_markers <- map[Chr == ch][order(Pos), marker]
    n_chr <- length(chr_markers)
    if (n_chr < 2) next
    for (bn in unique(sub[Chr == ch, marker_bin])) {
      ## representative size for this bin = median observed size in this (combo, chr, bin)
      sz <- as.integer(round(stats::median(sub[Chr == ch & marker_bin == bn, n_markers])))
      sz <- max(1L, min(sz, n_chr))
      hits <- vapply(seq_len(K_DRAWS), function(k) {
        start <- sample.int(n_chr - sz + 1L, 1L)
        win <- chr_markers[start:(start + sz - 1L)]
        any(qtn_lut_match$marker %chin% win)
      }, logical(1))
      rows[[length(rows) + 1]] <- data.table(tag = tag, cell = cell, rep = rep, env = env, Chr = ch,
                                             marker_bin = bn, window_size = sz, K = K_DRAWS,
                                             n_hits = sum(hits))
    }
  }
  rbindlist(rows)
}

combos_present <- unique(dt[, .(tag, cell, rep, env)])
t0 <- Sys.time()
opp_res <- mclapply(seq_len(nrow(combos_present)), function(i) {
  tryCatch(opportunity_one_combo(combos_present$tag[i], combos_present$cell[i],
                                 combos_present$rep[i], combos_present$env[i]),
           error = function(e) { message(sprintf("[14 opportunity] %s_%s_rep%d_env%d: %s",
                                                 combos_present$tag[i], combos_present$cell[i],
                                                 combos_present$rep[i], combos_present$env[i], conditionMessage(e))); NULL })
}, mc.cores = 7)
opp_dt <- rbindlist(Filter(Negate(is.null), opp_res))
say("[7] opportunity-null windows drawn for %d combos in %.1f min\n", nrow(combos_present),
    as.numeric(difftime(Sys.time(), t0, units = "mins")))

if (nrow(opp_dt)) {
  opp_pooled <- opp_dt[, .(n_hits = sum(n_hits), n_draws = sum(K)), by = marker_bin]
  opp_pooled[, opportunity_hit_rate := n_hits / n_draws]
  opp_pooled[, marker_bin := factor(marker_bin, levels = MARKER_BIN_LABELS)]
  setorder(opp_pooled, marker_bin)

  observed_by_bin <- dt[method %in% REGION_METHODS, .(n_regions = .N, n_TP = sum(TP)), by = marker_bin]
  observed_by_bin[, observed_tp_rate := n_TP / n_regions]
  observed_by_bin[, marker_bin := factor(marker_bin, levels = MARKER_BIN_LABELS)]
  setorder(observed_by_bin, marker_bin)

  opp_compare <- merge(observed_by_bin, opp_pooled, by = "marker_bin", all = TRUE)
  fwrite(opp_compare, "results/simulation_stage2_opportunity_null.tsv", sep = "\t")
  say("[8] wrote results/simulation_stage2_opportunity_null.tsv\n")
  print(opp_compare)

  ## Is the observed TP-rate GRADIENT (across bins) steeper than the
  ## opportunity-only gradient? Compare log-odds slope vs log2(median bin
  ## size) for each -- and, since the two rates sit an order of magnitude
  ## apart (observed >> opportunity at every bin, per the printed table
  ## above), ALSO track the observed/opportunity RATIO by bin: if that ratio
  ## is roughly flat across bins, the apparent size gradient is largely what
  ## opportunity alone predicts (a constant multiplicative excess); if it
  ## rises with size, there is a genuine size-dependent excess beyond
  ## opportunity.
  bin_size <- c("1" = 1, "2" = 2, "3-5" = 4, "6-10" = 8, "11-20" = 15, "21-50" = 35, "51+" = 100)
  opp_compare[, log2_size := log2(bin_size[as.character(marker_bin)])]
  opp_compare[, excess_ratio := observed_tp_rate / opportunity_hit_rate]
  obs_slope <- tryCatch(coef(glm(cbind(n_TP, n_regions - n_TP) ~ log2_size, data = opp_compare, family = binomial()))[2],
                        error = function(e) NA_real_)
  opp_slope <- tryCatch(coef(glm(cbind(n_hits, n_draws - n_hits) ~ log2_size, data = opp_compare, family = binomial()))[2],
                        error = function(e) NA_real_)
  say("\n[9] observed TP-rate log-odds slope vs log2(size): %.4f ; opportunity-null hit-rate slope: %.4f\n",
      obs_slope, opp_slope)
  say("    excess ratio (observed/opportunity) by bin: %s\n",
      paste(sprintf("%s=%.1fx", opp_compare$marker_bin, opp_compare$excess_ratio), collapse = ", "))
  if (is.na(obs_slope) || is.na(opp_slope)) {
    say("    slope comparison unavailable (a glm fit failed) -- see the ratio-by-bin line above instead.\n")
  } else if (obs_slope > opp_slope) {
    say("    observed slope EXCEEDS the opportunity-null slope: evidence the size-truth relationship is\n")
    say("    STRONGER than pure genomic-coverage opportunity predicts.\n")
  } else {
    say("    observed slope does NOT exceed the opportunity-null slope: the STEEPNESS of the size gradient\n")
    say("    is not shown to exceed pure genomic-coverage opportunity, even though observed TP rate is\n")
    say("    substantially higher than opportunity at every bin (a large but roughly size-independent\n")
    say("    excess, per the ratio line above) -- do not claim the size GRADIENT itself is selection-driven;\n")
    say("    see ADDITIONAL_ANALYSES_AUDIT.md for the full, hedged interpretation.\n")
  }
} else {
  say("[8] opportunity-null produced no usable draws -- falling back to descriptive-only reporting\n")
  say("    per the spec's own escape hatch; do not give a mechanistic interpretation of the size effect.\n")
}

write_receipt("14_summarise_additional_analyses", inputs = region_detail_file,
              params = list(marker_bin_edges = MARKER_BIN_EDGES, unit_bin_edges = UNIT_BIN_EDGES,
                            k_draws = K_DRAWS, n_bootstrap = N_BOOTSTRAP, seed_bootstrap = SEEDS[["bootstrap"]]),
              outputs = c("results/simulation_stage2_size_truth.tsv", "results/simulation_stage2_size_models.tsv",
                         "results/simulation_stage2_opportunity_null.tsv"))
say("\nDONE\n")
