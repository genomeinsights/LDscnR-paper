## =============================================================================
## final_analysis/R/15_env_structure_alignment.R
##
## Tests PK's hypothesis (2026-09-18, see env-structure-alignment-hypothesis.md
## memory): false positives should become common when environmental variation
## follows the same spatial pattern as genetic relatedness -- ENVIRONMENT-
## GENETIC-STRUCTURE ALIGNMENT, not simple collinearity, since relatedness is
## a matrix, not a single variable.
##
## Per combo (tag, cell, rep, env) -- i.e. PER ENVIRONMENTAL CONTINUATION, not
## pooled across the 10 per map/burn-in history, since pooling would average
## away exactly the variation this analysis needs:
##  1. Population-level relationship matrix: block-mean the individual-level
##     GRM within/between populations (env is already population-constant in
##     this module's bundles -- confirmed: 80 pops, env constant within pop).
##  2. PRIMARY alignment measure: R^2 of population-level env ~ top 5
##     eigenvectors of the population-level relationship matrix (5 axes,
##     matching the 5 spatial sampling groups already used elsewhere in this
##     module, e.g. R/13's group-scheme null).
##  3. SENSITIVITY measure: a Mantel-style correlation between the flattened
##     population-level relationship matrix and a flattened population-level
##     environmental-similarity matrix (Gaussian kernel on |env_i - env_j|,
##     bandwidth = median pairwise env difference -- same bandwidth
##     convention as R/13's gen_spatial() null, applied to env differences
##     instead of spatial coordinates). Avoids dependence on the choice of 5
##     axes.
##
## No association models or null permutations are rerun here -- GRM and env
## already sit in each combo's 02_build_ld_units bundle; observed FP
## proportion comes from R/11's region_details.rds; the spatial-null ratio
## comes from R/13's per-env output (results/simulation_null_truth_
## calibration_by_env.tsv, added alongside this analysis since it didn't
## exist before -- see that script's own header).
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(parallel)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "helpers_stage2_truth.R"))   ## ci_quantile()
say("=== 15_env_structure_alignment ===\n\n")

N_AXES <- 5L

alignment_one_combo <- function(tag, cell, rep, env) {
  combo_id <- sprintf("%s_%s_rep%d_env%d", tag, cell, rep, env)
  b <- readRDS(file.path(stage_dir("02_build_ld_units", combo_id), "ld_units.rds"))
  env_dt <- b$env
  GRM <- b$GRM

  ## population-level block-mean relationship matrix: Ind[k,p] = 1/n_p if
  ## individual k is in population p, else 0 -- t(Ind) %*% GRM %*% Ind then
  ## gives, entry-wise, the mean of the corresponding individual-level block,
  ## i.e. exactly the population-level relatedness matrix. Positional
  ## correspondence between GRM rows/cols and env_dt rows is relied on
  ## throughout this bundle already (no dimnames anywhere in it; R/13's own
  ## null generators make the same assumption against the same bundle).
  pops <- unique(env_dt$pop)
  n_pop <- length(pops)
  pop_idx <- match(env_dt$pop, pops)
  n_ind <- nrow(env_dt)
  Ind <- matrix(0, n_ind, n_pop)
  Ind[cbind(seq_len(n_ind), pop_idx)] <- 1
  pop_size <- colSums(Ind)
  Ind <- Ind %*% diag(1 / pop_size, n_pop, n_pop)
  R_pop <- t(Ind) %*% GRM %*% Ind
  R_pop <- (R_pop + t(R_pop)) / 2   ## symmetrise away float asymmetry before eigen()

  env_pop <- env_dt$env[match(pops, env_dt$pop)]   ## one env value per population (already pop-constant)

  ## primary measure: R^2 of env ~ top N_AXES eigenvectors of R_pop
  ax <- min(N_AXES, n_pop - 1L)
  eg <- eigen(R_pop, symmetric = TRUE)
  V <- eg$vectors[, seq_len(ax), drop = FALSE]
  fit <- stats::lm(env_pop ~ V)
  r2_axes <- summary(fit)$r.squared

  ## sensitivity measure: Mantel-style correlation between the flattened
  ## relationship matrix and a flattened env-similarity matrix (Gaussian
  ## kernel, same bandwidth convention as R/13's gen_spatial()).
  env_diff <- as.matrix(stats::dist(env_pop))
  bw <- stats::median(env_diff[lower.tri(env_diff)])
  env_sim <- if (bw > 0) exp(-0.5 * (env_diff / bw)^2) else matrix(1, n_pop, n_pop)
  lt <- lower.tri(R_pop)
  mantel_r <- suppressWarnings(stats::cor(R_pop[lt], env_sim[lt]))

  data.table(tag = tag, cell = cell, rep = rep, env = env, n_pop = n_pop,
            r2_axes = r2_axes, mantel_r = mantel_r)
}

if (sys.nframe() == 0L) {
  by_env <- fread("results/simulation_null_truth_calibration_by_env.tsv")
  combos <- unique(by_env[, .(tag, cell, rep, env)])
  say("[1] computing environment-structure alignment for %d combos\n", nrow(combos))
  t0 <- Sys.time()
  res <- mclapply(seq_len(nrow(combos)), function(i) {
    tryCatch(list(ok = TRUE, data = alignment_one_combo(combos$tag[i], combos$cell[i], combos$rep[i], combos$env[i])),
             error = function(e) {
               msg <- sprintf("[15_env_structure_alignment] %s_%s_rep%d_env%d: %s",
                              combos$tag[i], combos$cell[i], combos$rep[i], combos$env[i], conditionMessage(e))
               message(msg)
               list(ok = FALSE, data = NULL, error = msg)
             })
  }, mc.cores = 12)
  ok_flags <- vapply(res, `[[`, logical(1), "ok")
  if (!all(ok_flags))
    stop(sprintf("%d/%d combos failed -- see [15_env_structure_alignment] messages above; refusing to pool a partial result set", sum(!ok_flags), nrow(combos)))
  align <- rbindlist(lapply(res, `[[`, "data"))
  say("[2] %d/%d combos OK in %.1f min\n", nrow(align), nrow(combos), as.numeric(difftime(Sys.time(), t0, units = "mins")))

  ## ---- join against observed FP proportion (per method, from R/11) and the spatial-null ratio (from R/13) ----
  region_detail_file <- "results/simulation_stage2_region_details.rds"
  if (!file.exists(region_detail_file)) stop("R/11 has not produced ", region_detail_file, " -- run it first.")
  rd <- readRDS(region_detail_file)
  method_map <- c(emmax_consensus = "emmax_consensus_region", emmax_simes = "emmax_simes_region")
  obs <- rd[method %in% method_map, .(n_regions = .N, n_FP = sum(!TP)), by = .(tag, cell, rep, env, method)]
  obs[, method := names(method_map)[match(method, method_map)]]
  obs[, fp_prop := n_FP / n_regions]

  null_spatial <- by_env[scheme == "spatial", .(tag, cell, rep, env, method, ratio_null_obs)]

  gp <- merge(align, obs, by = c("tag", "cell", "rep", "env"), all.x = TRUE, allow.cartesian = TRUE)
  gp <- merge(gp, null_spatial, by = c("tag", "cell", "rep", "env", "method"), all.x = TRUE)
  gp <- gp[!is.na(method)]   ## drop combos with no observed regions for a method -- fp_prop undefined, not a false zero

  fwrite(gp, "results/simulation_env_structure_alignment.tsv", sep = "\t")
  say("[3] wrote results/simulation_env_structure_alignment.tsv (%d rows)\n", nrow(gp))
  print(gp[, .(n = .N, median_r2 = median(r2_axes), median_mantel = median(mantel_r)), by = cell])

  ## =============================================================================
  ## [3b] Complementary full-grid table (PK, 2026-09-18 second review, point
  ## 3): `gp` above is conditioned on >=1 reported region existing for that
  ## (combo, method) -- fp_prop is undefined otherwise, so that conditioning
  ## is correct for fp_prop itself, but it means `gp` silently drops the
  ## 516/1,400 (combo,method) pairs with ZERO reported regions, and cannot
  ## answer "does alignment make a false positive more likely to occur AT
  ## ALL." `full_grid` below keeps every (combo, method) pair (700 combos x
  ## 2 methods = 1,400 rows), zero-filling n_regions/n_FP where R/11 has no
  ## row for that method (a REAL zero, not a missing value -- same
  ## discipline as elsewhere in this pipeline, never max(x,1)).
  ## =============================================================================
  REGION_METHODS_SHORT <- c("emmax_consensus", "emmax_simes")
  full_grid <- align[, .(tag, cell, rep, env, r2_axes, mantel_r, n_pop)][rep(seq_len(.N), each = length(REGION_METHODS_SHORT))]
  full_grid[, method := rep(REGION_METHODS_SHORT, times = nrow(align))]
  full_grid <- merge(full_grid, obs, by = c("tag", "cell", "rep", "env", "method"), all.x = TRUE)
  full_grid[is.na(n_regions), `:=`(n_regions = 0L, n_FP = 0L, fp_prop = NA_real_)]
  full_grid[, any_fp := as.integer(n_FP > 0)]
  ## mean_null (expected spatial-null region count) is defined regardless of
  ## whether any region was actually observed -- the null draws ran either
  ## way -- so this merge should never introduce NAs; asserted, not assumed.
  null_mean <- by_env[scheme == "spatial", .(tag, cell, rep, env, method, mean_null)]
  full_grid <- merge(full_grid, null_mean, by = c("tag", "cell", "rep", "env", "method"), all.x = TRUE)
  if (anyNA(full_grid$mean_null)) stop("full_grid: mean_null missing for some (combo,method) -- R/13's by-env export should cover every one")
  fwrite(full_grid, "results/simulation_env_structure_alignment_full_grid.tsv", sep = "\t")
  say("[3b] wrote results/simulation_env_structure_alignment_full_grid.tsv (%d rows, all combos incl. zero-region ones)\n", nrow(full_grid))

  ## =============================================================================
  ## [4] Slope models. Map-cluster bootstrap (refit per replicate), clustered
  ## by (cell, rep) -- 10 env continuations of a (cell,rep) share a map/
  ## burn-in history and must be resampled together, paired across tag (same
  ## convention as R/14's .fit_size_model()).
  ##
  ## [!] BUG FIXED (PK, 2026-09-18 third review): the previous version built
  ## `cl_rows`/`n_cl` ONCE from the full `gp` table and reused them for the
  ## ratio model too, which is fit on `ratio_d` (`gp` filtered to
  ## ratio_null_obs>0, a SMALLER table with its OWN 1..nrow(ratio_d) row
  ## numbering). A bootstrap index built from gp's numbering routinely
  ## exceeded ratio_d's row count; data.table silently returns an all-NA row
  ## for an out-of-range index rather than erroring, corrupting that
  ## replicate's resample instead of failing loudly. `.boot_slope()` below
  ## now ALWAYS builds its own cluster info fresh from whatever `data` it is
  ## actually given, so this class of mismatch cannot recur regardless of
  ## which subset is passed in (unconditional/adjusted/within-cell/
  ## full-grid all reuse the identical helper now, replacing the previous
  ## separate, duplicated .boot_slope_c()).
  ## =============================================================================
  say("\n[4] slope models: alignment vs. FP occurrence/burden and null calibration\n")
  .prep <- function(d) {
    d <- copy(d)
    d[, cell_f := factor(cell, levels = CELLS_ALL)]
    d[, tag_f := factor(tag, levels = TAGS_ALL)]
    d[, method_f := factor(method, levels = c("emmax_consensus", "emmax_simes"))]
    d
  }
  gp <- .prep(gp)
  full_grid <- .prep(full_grid)

  .boot_slope <- function(fit_fn, extract_fn, data) {
    data <- copy(data)
    clusters <- unique(data[, .(cell, rep)])
    n_cl <- nrow(clusters)
    data[, cl_id := .GRP, by = .(cell, rep)]
    cl_rows <- split(seq_len(nrow(data)), data$cl_id)
    fit0 <- tryCatch(fit_fn(data), error = function(e) NULL)
    if (is.null(fit0)) return(data.table(estimate = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, n = nrow(data), n_boot_converged = 0L))
    point <- extract_fn(fit0)
    set.seed(SEEDS[["bootstrap"]])
    boot <- rep(NA_real_, N_BOOTSTRAP)
    for (bb in seq_len(N_BOOTSTRAP)) {
      draw <- sample.int(n_cl, n_cl, replace = TRUE)
      ## list-index by the draw vector directly (NOT intersect(), which
      ## would deduplicate repeated draws and break resampling-with-
      ## replacement) -- indexing a list by a vector of names naturally
      ## repeats an element as many times as its name appears in the index.
      idx <- unlist(cl_rows[as.character(draw)], use.names = FALSE)
      fit_b <- tryCatch(fit_fn(data[idx]), error = function(e) NULL, warning = function(w) NULL)
      if (!is.null(fit_b)) boot[bb] <- tryCatch(extract_fn(fit_b), error = function(e) NA_real_)
    }
    boot <- boot[is.finite(boot)]
    data.table(estimate = point, ci_lo = if (length(boot)) ci_quantile(boot)[1] else NA_real_,
              ci_hi = if (length(boot)) ci_quantile(boot)[2] else NA_real_,
              n = nrow(data), n_boot_converged = length(boot))
  }

  ratio_d <- gp[ratio_null_obs > 0]
  ratio_d[, log2_ratio := log2(ratio_null_obs)]
  full_grid[, log2_null := log2(mean_null + 1)]   ## log1p-style transform -- mean_null can be exactly 0 (a real, common outcome, not an error)

  ## outcome -> (fit formula template, family/method, data, extractor), one
  ## row per (outcome x alignment measure) so the same machinery covers both
  ## the primary axis-based measure and the Mantel-style sensitivity check
  ## (point 4 of the review: previously computed but never modelled).
  make_specs <- function(align_var) {
    list(
      fp_prop = list(outcome = "fp_prop",
                     fit = function(d) stats::glm(stats::as.formula(sprintf("cbind(n_FP, n_regions - n_FP) ~ %s + method_f", align_var)), data = d, family = stats::binomial()),
                     fit_adj = function(d) stats::glm(stats::as.formula(sprintf("cbind(n_FP, n_regions - n_FP) ~ %s + method_f + cell_f + tag_f", align_var)), data = d, family = stats::binomial()),
                     data = gp),
      null_ratio_log2 = list(outcome = "null_ratio_log2",
                             fit = function(d) stats::lm(stats::as.formula(sprintf("log2_ratio ~ %s + method_f", align_var)), data = d),
                             fit_adj = function(d) stats::lm(stats::as.formula(sprintf("log2_ratio ~ %s + method_f + cell_f + tag_f", align_var)), data = d),
                             data = ratio_d),
      ## complementary, full-grid outcomes (point 3): does alignment predict
      ## whether a false positive occurs at all, how many occur, and the
      ## expected null burden -- all WITHOUT conditioning on >=1 discovery.
      any_fp = list(outcome = "any_fp",
                   fit = function(d) stats::glm(stats::as.formula(sprintf("any_fp ~ %s + method_f", align_var)), data = d, family = stats::binomial()),
                   fit_adj = function(d) stats::glm(stats::as.formula(sprintf("any_fp ~ %s + method_f + cell_f + tag_f", align_var)), data = d, family = stats::binomial()),
                   data = full_grid),
      n_fp_count = list(outcome = "n_fp_count",
                       fit = function(d) stats::glm(stats::as.formula(sprintf("n_FP ~ %s + method_f", align_var)), data = d, family = stats::poisson()),
                       fit_adj = function(d) stats::glm(stats::as.formula(sprintf("n_FP ~ %s + method_f + cell_f + tag_f", align_var)), data = d, family = stats::poisson()),
                       data = full_grid),
      null_mean_log2 = list(outcome = "null_mean_log2",
                            fit = function(d) stats::lm(stats::as.formula(sprintf("log2_null ~ %s + method_f", align_var)), data = d),
                            fit_adj = function(d) stats::lm(stats::as.formula(sprintf("log2_null ~ %s + method_f + cell_f + tag_f", align_var)), data = d),
                            data = full_grid)
    )
  }
  extract_for <- function(align_var) function(fit) unname(coef(fit)[align_var])

  model_rows <- list()
  for (align_var in c("r2_axes", "mantel_r")) {
    specs <- make_specs(align_var)
    ext <- extract_for(align_var)
    for (spec in specs) {
      model_rows[[length(model_rows) + 1]] <- cbind(outcome = spec$outcome, align_measure = align_var, variant = "unconditional", cell = NA_character_,
                                                     .boot_slope(spec$fit, ext, spec$data))
      model_rows[[length(model_rows) + 1]] <- cbind(outcome = spec$outcome, align_measure = align_var, variant = "adjusted_cell_tag_method", cell = NA_character_,
                                                     .boot_slope(spec$fit_adj, ext, spec$data))
    }
  }

  ## within-cell: r2_axes only, fp_prop and null_ratio only (the two
  ## outcomes PK's original spec asked to break out by cell; the
  ## complementary full-grid outcomes and the Mantel measure are reported
  ## unconditional/adjusted only, to keep this already-large model table
  ## bounded -- noted explicitly in the audit, not silently scoped down).
  ext_r2 <- extract_for("r2_axes")
  for (cc in CELLS_ALL) {
    d_fp <- gp[cell == cc]
    if (nrow(d_fp) >= 10) {
      fp_fit_fn <- function(d) stats::glm(cbind(n_FP, n_regions - n_FP) ~ r2_axes + method_f, data = d, family = stats::binomial())
      model_rows[[length(model_rows) + 1]] <- cbind(outcome = "fp_prop", align_measure = "r2_axes", variant = "within_cell", cell = cc,
                                                     .boot_slope(fp_fit_fn, ext_r2, d_fp))
      d_ratio <- ratio_d[cell == cc]
      if (nrow(d_ratio) >= 10) {
        ratio_fit_fn <- function(d) stats::lm(log2_ratio ~ r2_axes + method_f, data = d)
        model_rows[[length(model_rows) + 1]] <- cbind(outcome = "null_ratio_log2", align_measure = "r2_axes", variant = "within_cell", cell = cc,
                                                       .boot_slope(ratio_fit_fn, ext_r2, d_ratio))
      }
    }
  }
  model_dt <- rbindlist(model_rows)
  fwrite(model_dt, "results/simulation_env_structure_alignment_models.tsv", sep = "\t")
  say("[5] wrote results/simulation_env_structure_alignment_models.tsv\n")
  print(model_dt)
  say("\n[6] within-cell CIs excluding zero (flagging, not hiding, any that do):\n")
  excl <- model_dt[variant == "within_cell" & !is.na(ci_lo) & !is.na(ci_hi) & (ci_lo > 0 | ci_hi < 0)]
  if (nrow(excl)) print(excl[, .(outcome, cell, estimate, ci_lo, ci_hi)]) else say("    none\n")

  write_receipt("15_env_structure_alignment", inputs = c("results/simulation_null_truth_calibration_by_env.tsv", region_detail_file),
                params = list(n_axes = N_AXES, n_bootstrap = N_BOOTSTRAP, seed_bootstrap = SEEDS[["bootstrap"]]),
                outputs = c("results/simulation_env_structure_alignment.tsv", "results/simulation_env_structure_alignment_full_grid.tsv",
                           "results/simulation_env_structure_alignment_models.tsv"))
  cat("ENV_ALIGNMENT_DONE\n")
}
