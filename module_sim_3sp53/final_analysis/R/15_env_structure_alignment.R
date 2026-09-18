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
  ## [4] Slope models: does known FP proportion / the spatial-null-to-observed
  ## ratio increase with environment-structure alignment (r2_axes)? Three
  ## variants each, per PK's spec: unconditional, adjusted for (cell, tag,
  ## method) -- cell alone encodes both V and c jointly, so this is
  ## "adjusted for c, V, BGS, method" without a collinear separate V/c term
  ## -- and within-cell (one fit per cell). Map-cluster bootstrap (refit per
  ## replicate), clustered by (cell, rep) -- the same convention as R/14's
  ## .fit_size_model(), for the same reason: the 10 env continuations of a
  ## (cell, rep) share a map/burn-in history and must be resampled together,
  ## paired across tag.
  ## =============================================================================
  say("\n[4] slope models: FP proportion and null-ratio vs. environment-structure alignment\n")
  gp[, cell_f := factor(cell, levels = CELLS_ALL)]
  gp[, tag_f := factor(tag, levels = TAGS_ALL)]
  gp[, method_f := factor(method, levels = c("emmax_consensus", "emmax_simes"))]

  clusters <- unique(gp[, .(cell, rep)])
  n_cl <- nrow(clusters)
  gp[, cl_id := .GRP, by = .(cell, rep)]
  cl_rows <- split(seq_len(nrow(gp)), gp$cl_id)

  .boot_slope <- function(fit_fn, extract_fn, data) {
    fit0 <- tryCatch(fit_fn(data), error = function(e) NULL)
    if (is.null(fit0)) return(data.table(estimate = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, n = nrow(data), n_boot_converged = 0L))
    point <- extract_fn(fit0)
    set.seed(SEEDS[["bootstrap"]])
    boot <- rep(NA_real_, N_BOOTSTRAP)
    for (bb in seq_len(N_BOOTSTRAP)) {
      draw <- sample.int(n_cl, n_cl, replace = TRUE)
      ## list-index by the draw vector directly (NOT intersect(), which would
      ## deduplicate repeated draws and silently break resampling-with-
      ## replacement) -- indexing a list by a character vector naturally
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

  fp_fit_fn <- function(d) stats::glm(cbind(n_FP, n_regions - n_FP) ~ r2_axes + method_f, data = d, family = stats::binomial())
  fp_fit_adj_fn <- function(d) stats::glm(cbind(n_FP, n_regions - n_FP) ~ r2_axes + method_f + cell_f + tag_f, data = d, family = stats::binomial())
  fp_extract <- function(fit) unname(coef(fit)["r2_axes"])

  ratio_d <- gp[ratio_null_obs > 0]
  ratio_d[, log2_ratio := log2(ratio_null_obs)]
  ratio_fit_fn <- function(d) stats::lm(log2_ratio ~ r2_axes + method_f, data = d)
  ratio_fit_adj_fn <- function(d) stats::lm(log2_ratio ~ r2_axes + method_f + cell_f + tag_f, data = d)
  ratio_extract <- function(fit) unname(coef(fit)["r2_axes"])

  model_rows <- list()
  model_rows[["fp_unconditional"]] <- cbind(outcome = "fp_prop", variant = "unconditional", cell = NA_character_, .boot_slope(fp_fit_fn, fp_extract, gp))
  model_rows[["fp_adjusted"]] <- cbind(outcome = "fp_prop", variant = "adjusted_cell_tag_method", cell = NA_character_, .boot_slope(fp_fit_adj_fn, fp_extract, gp))
  model_rows[["ratio_unconditional"]] <- cbind(outcome = "null_ratio_log2", variant = "unconditional", cell = NA_character_, .boot_slope(ratio_fit_fn, ratio_extract, ratio_d))
  model_rows[["ratio_adjusted"]] <- cbind(outcome = "null_ratio_log2", variant = "adjusted_cell_tag_method", cell = NA_character_, .boot_slope(ratio_fit_adj_fn, ratio_extract, ratio_d))

  ## within-cell: same (cell,rep) cluster-bootstrap machinery, restricted to
  ## one cell's rows at a time (n_cl there is that cell's own cluster count).
  for (cc in CELLS_ALL) {
    d_fp <- gp[cell == cc]
    if (nrow(d_fp) >= 10) {
      clusters_c <- unique(d_fp[, .(cell, rep)]); n_cl_c <- nrow(clusters_c)
      d_fp[, cl_id := .GRP, by = rep]
      cl_rows_c <- split(seq_len(nrow(d_fp)), d_fp$cl_id)
      .boot_slope_c <- function(fit_fn, extract_fn, data, n_cl_l, cl_rows_l) {
        fit0 <- tryCatch(fit_fn(data), error = function(e) NULL)
        if (is.null(fit0)) return(data.table(estimate = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, n = nrow(data), n_boot_converged = 0L))
        point <- extract_fn(fit0)
        set.seed(SEEDS[["bootstrap"]])
        boot <- rep(NA_real_, N_BOOTSTRAP)
        for (bb in seq_len(N_BOOTSTRAP)) {
          draw <- sample.int(n_cl_l, n_cl_l, replace = TRUE)
          idx <- unlist(cl_rows_l[as.character(draw)], use.names = FALSE)   ## see .boot_slope()'s comment -- never intersect()
          fit_b <- tryCatch(fit_fn(data[idx]), error = function(e) NULL, warning = function(w) NULL)
          if (!is.null(fit_b)) boot[bb] <- tryCatch(extract_fn(fit_b), error = function(e) NA_real_)
        }
        boot <- boot[is.finite(boot)]
        data.table(estimate = point, ci_lo = if (length(boot)) ci_quantile(boot)[1] else NA_real_,
                  ci_hi = if (length(boot)) ci_quantile(boot)[2] else NA_real_, n = nrow(data), n_boot_converged = length(boot))
      }
      model_rows[[paste0("fp_", cc)]] <- cbind(outcome = "fp_prop", variant = "within_cell", cell = cc,
                                               .boot_slope_c(fp_fit_fn, fp_extract, d_fp, n_cl_c, cl_rows_c))
      d_ratio <- ratio_d[cell == cc]
      if (nrow(d_ratio) >= 10) {
        d_ratio[, cl_id := .GRP, by = rep]
        cl_rows_r <- split(seq_len(nrow(d_ratio)), d_ratio$cl_id)
        model_rows[[paste0("ratio_", cc)]] <- cbind(outcome = "null_ratio_log2", variant = "within_cell", cell = cc,
                                                    .boot_slope_c(ratio_fit_fn, ratio_extract, d_ratio, n_cl_c, cl_rows_r))
      }
    }
  }
  model_dt <- rbindlist(model_rows)
  fwrite(model_dt, "results/simulation_env_structure_alignment_models.tsv", sep = "\t")
  say("[5] wrote results/simulation_env_structure_alignment_models.tsv\n")
  print(model_dt)

  write_receipt("15_env_structure_alignment", inputs = c("results/simulation_null_truth_calibration_by_env.tsv", region_detail_file),
                params = list(n_axes = N_AXES, n_bootstrap = N_BOOTSTRAP, seed_bootstrap = SEEDS[["bootstrap"]]),
                outputs = c("results/simulation_env_structure_alignment.tsv", "results/simulation_env_structure_alignment_models.tsv"))
  cat("ENV_ALIGNMENT_DONE\n")
}
