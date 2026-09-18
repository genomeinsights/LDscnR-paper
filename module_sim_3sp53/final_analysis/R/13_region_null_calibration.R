## =============================================================================
## final_analysis/R/13_region_null_calibration.R
##
## Analysis 1: does a structure-aware null's discovery burden predict the
## known false-positive burden -- AT THE STAGE-2 REGION SCALE (the exploratory
## module_sim_3sp53/R/12_structured_null.R counted significant STAGE-1 UNITS;
## this redo counts assembled Stage-2 regions, on corrected-parser inputs, and
## is the reason this file exists rather than reusing that one -- see
## ADDITIONAL_ANALYSES_AUDIT.md).
##
## Three null constructions (recipes ported from module_sim_3sp53/R/12, NOT
## its numeric outputs -- those are pre-correction-parser and not reusable;
## see that file's own header for the validated recipe this reproduces):
##   group   : permute the population-level environmental value within each
##             of 5 spatial sampling groups (kmeans(k=5) on population (x,y)
##             -- 00_config.R's LFMM_K=5 justification already confirms this
##             recovers the true by-design 4-corners+centre grouping).
##   mvn     : kinship-matched continuous null, s ~ MVN(0, GRM), orthogonalised
##             against the observed phenotype.
##   spatial : spatially autocorrelated continuous null, Gaussian-kernel MVN
##             over individual (x,y), same orthogonalisation.
##
## Per-null-replicate procedure (never approximated by a Stage-1-unit count):
## EMMAX p-values (emmax_fast() on the ALREADY-FITTED emmax_setup() prep --
## no refit) -> BH at ALPHA -> assemble_stage2() -> count Stage-2 regions.
## EMMAX consensus and EMMAX Simes only (LFMM excluded per spec: a valid LFMM
## null needs the full latent-factor model refit per replicate).
##
## Observed Stage-2 region counts are READ from results/
## simulation_stage2_region_details.rds (R/11), never recomputed here --
## guarantees "observed Stage-2 counts in the null analysis equal the
## corresponding observed truth-scoring counts" by construction, not by a
## second implementation that could drift.
##
## max(n_observed, 1) is NEVER used to manufacture a ratio -- zero-observed
## strata get NA and are counted, not silently dropped or floored (the exact
## failure mode module_sim_3sp53/R/12 itself had to be corrected for; see
## that module's README, "2026-09-09 CORRECTED").
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(parallel)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "helpers_stage2_truth.R"))
say("=== 13_region_null_calibration ===\n\n")

NULL_METHODS <- c("emmax_consensus", "emmax_simes")
NULL_SCHEMES <- c("group", "mvn", "spatial")

## ---- null-phenotype generators, one build per combo, reused across B draws --
## Recipe verbatim from module_sim_3sp53/R/12_structured_null.R (validated
## there to recover the true 4-corners+centre design grouping); only the
## INPUTS (this combo's corrected-parser GTs/map/GRM/env) are new.
.build_null_generators <- function(GTs, GRM, env_dt, y) {
  n_ind <- nrow(GTs)

  pos <- unique(env_dt[, .(pop, x, y)])
  set.seed(1L)
  km <- stats::kmeans(pos[, .(x, y)], centers = 5, nstart = 10)
  pos[, pop_group := km$cluster]
  POPT <- merge(unique(env_dt[, .(pop, env)]), pos[, .(pop, pop_group)], by = "pop")
  ## individual -> POPT row index, computed ONCE: POPT$pop's identity/order is
  ## fixed once built (copy()+`:=` below never reorders rows), only $ep is
  ## reshuffled per replicate -- re-matching env_dt$pop against pt$pop inside
  ## gen_group() would re-hash pt$pop from scratch on every one of the B
  ## calls for no reason (same class of issue as .genomic_sort_clusters()'s
  ## fix in helpers_stage2_truth.R, flagged the same day).
  pop_idx <- match(env_dt$pop, POPT$pop)
  gen_group <- function(seed) {
    set.seed(seed)
    pt <- copy(POPT)[, ep := sample(env), by = pop_group]
    pt$ep[pop_idx]
  }

  eK <- eigen(GRM, symmetric = TRUE)
  Lv <- pmax(eK$values, 0); Vk <- eK$vectors
  gen_mvn <- function(seed) {
    set.seed(seed)
    s <- as.numeric(Vk %*% (sqrt(Lv) * stats::rnorm(n_ind)))
    as.numeric(stats::resid(stats::lm(s ~ y)))
  }

  coords <- as.matrix(env_dt[, .(x, y)])
  Dm <- as.matrix(stats::dist(coords)); l_bw <- stats::median(Dm[lower.tri(Dm)])
  eK_sp <- eigen(exp(-0.5 * (Dm / l_bw)^2), symmetric = TRUE)
  Lv_sp <- pmax(eK_sp$values, 0); Vk_sp <- eK_sp$vectors
  gen_spatial <- function(seed) {
    set.seed(seed)
    s <- as.numeric(Vk_sp %*% (sqrt(Lv_sp) * stats::rnorm(n_ind)))
    as.numeric(stats::resid(stats::lm(s ~ y)))
  }

  list(group = gen_group, mvn = gen_mvn, spatial = gen_spatial)
}

## ---- one (combo, method, scheme) -> B Stage-2 region counts -----------------
## Deterministic seed per (b, scheme-offset), matching module_sim_3sp53/R/12's
## own convention (group: seed=b; mvn: seed=b+1e6; spatial: seed=b+2e6) so a
## given b draws genuinely different randomness across schemes.
SCHEME_SEED_OFFSET <- c(group = 0, mvn = 1e6, spatial = 2e6)

null_region_counts <- function(tag, cell, rep, env, method, scheme, B) {
  combo_id <- sprintf("%s_%s_rep%d_env%d", tag, cell, rep, env)
  b <- readRDS(file.path(stage_dir("02_build_ld_units", combo_id), "ld_units.rds"))
  GTs <- b$GTs; map <- copy(b$map); stage1 <- b$stage1; GRM <- b$GRM; env_dt <- b$env
  y <- env_dt$env
  units_base <- LDscnR:::.ld_outlier_units(stage1, map, SIZE_FLOOR)
  gens <- .build_null_generators(GTs, GRM, env_dt, y)
  gen <- gens[[scheme]]
  off <- SCHEME_SEED_OFFSET[[scheme]]

  if (method == "emmax_consensus") {
    um <- ld_unit_matrix(GTs, stage1, map, size_floor = SIZE_FLOOR, repr = "consensus")
    stopifnot(identical(colnames(um), as.character(units_base$unit_id)))
    P <- emmax_setup(um, GRM)
    one <- function(bb) {
      p_null <- emmax_fast(P, gen(bb + off))
      q <- stats::p.adjust(p_null, "BH")
      sig_core <- units_base$core_snp[!is.na(q) & q <= ALPHA]
      cl_sub <- stage2_seed_from_units(stage1, map, sig_core)
      nrow(assemble_stage2(stage1, map, GTs, b$LD_decay, cl_sub,
                           REGION_ASSEMBLY$score_threshold, REGION_ASSEMBLY$distance_threshold))
    }
  } else {  ## emmax_simes
    Pm <- emmax_setup(GTs, GRM)
    marker_names <- map$marker
    one <- function(bb) {
      p_null <- emmax_fast(Pm, gen(bb + off))
      names(p_null) <- marker_names
      simes_p <- vapply(units_base$members, function(mm) LDscnR:::.simes(p_null[mm]), numeric(1))
      q <- stats::p.adjust(simes_p, "BH")
      sig_core <- units_base$core_snp[!is.na(q) & q <= ALPHA]
      cl_sub <- stage2_seed_from_units(stage1, map, sig_core)
      nrow(assemble_stage2(stage1, map, GTs, b$LD_decay, cl_sub,
                           REGION_ASSEMBLY$score_threshold, REGION_ASSEMBLY$distance_threshold))
    }
  }
  vapply(seq_len(B), one, integer(1))
}

## ---- resumable per-(combo,method,scheme) cache ------------------------------
## [!] BUG FIXED (PK, 2026-09-18 second review): the cache-hit branch used to
## return the whole saved list (tag/cell/.../counts), not just $counts, while
## the fresh-compute branch returned counts directly -- a shape mismatch the
## caller's mean()/max() calls didn't catch cleanly: mean() on a list just
## warns and returns NA, but max() throws "invalid 'type' (list) of
## argument", caught by the caller's tryCatch and silently dropped as a
## failed job. Confirmed in the overnight run's log: 1,134/4,200 jobs lost
## this way. NULL_CALIB_VERSION bumped to force every job to recompute under
## the fixed reader (the on-disk $counts values themselves were never
## corrupted -- only this accessor was broken -- but recomputing from a
## clean slate removes any doubt, per PK's request).
NULL_CALIB_VERSION <- 2L
run_one <- function(tag, cell, rep, env, method, scheme, B, force = FALSE) {
  STAGE <- "13_region_null_calibration"
  combo_id <- sprintf("%s_%s_rep%d_env%d", tag, cell, rep, env)
  target <- sprintf("%s__%s__%s", combo_id, method, scheme)
  ld_units_file <- file.path(stage_dir("02_build_ld_units", combo_id), "ld_units.rds")
  PARAMS <- list(method = method, scheme = scheme, B = B, alpha = ALPHA, size_floor = SIZE_FLOOR,
                 region_assembly = REGION_ASSEMBLY, seed_offset = unname(SCHEME_SEED_OFFSET[scheme]),
                 null_calib_version = NULL_CALIB_VERSION)
  if (!force && !stage_stale(STAGE, ld_units_file, PARAMS, target = target)) {
    return(readRDS(file.path(stage_dir(STAGE, target), "null_counts.rds"))$counts)
  }
  counts <- null_region_counts(tag, cell, rep, env, method, scheme, B)
  OUT_DIR <- stage_dir(STAGE, target); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  OUT <- file.path(OUT_DIR, "null_counts.rds")
  saveRDS(list(tag = tag, cell = cell, rep = rep, env = env, method = method, scheme = scheme,
              B = B, counts = counts), OUT)
  write_receipt(STAGE, inputs = ld_units_file, params = PARAMS, outputs = OUT, target = target)
  counts
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  MODE <- if (length(args) >= 1) args[1] else "smoke"   ## "smoke" | "full" | "balanced"
  B <- if (length(args) >= 2) as.integer(args[2]) else if (MODE == "smoke") 30L else 200L

  if (MODE == "smoke") {
    combos <- data.table(tag = c("nobgs", "bgs"), cell = c("V0.5_c1", "V1_c1"), rep = c(1L, 2L), env = c(1L, 3L))
  } else if (MODE == "balanced") {
    ## Documented balanced design (used only if the full grid proves
    ## impractical -- see ADDITIONAL_ANALYSES_AUDIT.md for the timing that
    ## justified this choice): N_REPS_PER_CELL map/burn-in reps per cell, ALL
    ## 10 environmental continuations for each selected rep (required for the
    ## map-level pooling the spec asks for), both BGS treatments, all 7
    ## cells.
    ## [!] PAIRED MAP IDENTIFIERS (PK, 2026-09-18 review, point 3): reps are
    ## sampled ONCE PER CELL and shared across BOTH tag values, not drawn
    ## independently per (tag, cell) as the previous design did. BGS and
    ## no-BGS simulations for the same (cell, rep) are the same map/burn-in
    ## history and env draw by construction (same seed lineage) -- sampling
    ## reps independently per tag meant the two treatments were almost never
    ## actually comparing the same underlying map, undermining any
    ## cell/map-level pairing in the downstream bootstrap or interpretation.
    ## N_REPS_PER_CELL raised from 2 to 5 ("at least five paired maps" per
    ## the review) now that the match()-vectorisation fix (2026-09-xx) has
    ## cut per-replicate cost -- retimed on the mini before the full launch,
    ## see ADDITIONAL_ANALYSES_AUDIT.md. 7 cells x 5 reps x 2 tags x 10 envs
    ## = 700 combos.
    N_REPS_PER_CELL <- 5L
    set.seed(SEEDS[["bootstrap"]])
    reps_by_cell <- data.table(cell = CELLS_ALL)[, .(rep = sample(REPS_ALL, N_REPS_PER_CELL)), by = cell]
    combos <- reps_by_cell[, CJ(tag = TAGS_ALL, env = ENVS_ALL), by = .(cell, rep)]

    ## [!] Hard design-completeness assertion (PK, 2026-09-18 second review):
    ## verify the NOMINAL design table itself -- N_REPS_PER_CELL distinct
    ## reps per cell, BOTH tags present for every selected (cell, rep), and
    ## all 10 envs per (tag, cell, rep) -- before spending any compute on it.
    ## True by construction of the CJ() above, but asserted explicitly rather
    ## than assumed, per the review.
    reps_per_cell_n <- combos[, uniqueN(rep), by = cell]
    tags_per_cellrep_n <- combos[, uniqueN(tag), by = .(cell, rep)]
    envs_per_tagcellrep_n <- combos[, uniqueN(env), by = .(tag, cell, rep)]
    if (!all(reps_per_cell_n$V1 == N_REPS_PER_CELL))
      stop("balanced design malformed: not every cell has exactly ", N_REPS_PER_CELL, " reps")
    if (!all(tags_per_cellrep_n$V1 == length(TAGS_ALL)))
      stop("balanced design malformed: not every (cell, rep) has both tags")
    if (!all(envs_per_tagcellrep_n$V1 == length(ENVS_ALL)))
      stop("balanced design malformed: not every (tag, cell, rep) has all ", length(ENVS_ALL), " envs")
    say("[0] design check OK: %d cells x %d reps/cell x %d tags x %d envs = %d combos\n",
        uniqueN(combos$cell), N_REPS_PER_CELL, length(TAGS_ALL), length(ENVS_ALL), nrow(combos))
  } else {
    combos <- CJ(tag = TAGS_ALL, cell = CELLS_ALL, rep = REPS_ALL, env = ENVS_ALL)
  }
  jobs <- CJ(i = seq_len(nrow(combos)), method = NULL_METHODS, scheme = NULL_SCHEMES)
  say("[1] MODE=%s : %d combos x %d methods x %d schemes = %d jobs, B=%d\n",
      MODE, nrow(combos), length(NULL_METHODS), length(NULL_SCHEMES), nrow(jobs), B)

  t0 <- Sys.time()
  res <- mclapply(seq_len(nrow(jobs)), function(k) {
    i <- jobs$i[k]
    tryCatch({
      counts <- run_one(combos$tag[i], combos$cell[i], combos$rep[i], combos$env[i], jobs$method[k], jobs$scheme[k], B)
      ## [!] Payload shape/content validation (PK, 2026-09-18 second
      ## review): the global nrow(dt)==nrow(jobs) completeness check catches
      ## a MISSING row, but not a malformed one -- a corrupt cache file
      ## could still return e.g. a wrong-length or NA-containing vector and
      ## produce one silently-bad row that passes the row-count check.
      ## Throwing here routes it through the same tryCatch/completeness-stop
      ## path as any other job failure, rather than letting a bad value
      ## reach mean()/max() and pool silently.
      if (!is.numeric(counts) || length(counts) != B ||
          !all(is.finite(counts)) || any(counts < 0) || any(counts != round(counts)))
        stop(sprintf("malformed counts payload: class=%s length=%d (expected %d)", class(counts)[1], length(counts), B))
      data.table(tag = combos$tag[i], cell = combos$cell[i], rep = combos$rep[i], env = combos$env[i],
                method = jobs$method[k], scheme = jobs$scheme[k], B = B, mean_null = mean(counts), max_null = max(counts))
    }, error = function(e) {
      message(sprintf("[13_region_null_calibration] %s/%s/rep%d/env%d %s/%s: %s",
                      combos$tag[i], combos$cell[i], combos$rep[i], combos$env[i],
                      jobs$method[k], jobs$scheme[k], conditionMessage(e)))
      NULL
    })
  }, mc.cores = 12)   ## bumped from 7 (PK, 2026-09-18): mini has 14 physical cores, 2 left for the system
  dt <- rbindlist(Filter(Negate(is.null), res))
  ## [!] Hard completeness stop (PK, 2026-09-18 second review): a partial
  ## result set was previously silently pooled as if complete (the
  ## cache-shape bug fixed above lost 1,134/4,200 jobs this way without the
  ## pipeline ever stopping). Any job failure now halts the run rather than
  ## being averaged over.
  if (nrow(dt) != nrow(jobs))
    stop(sprintf("only %d/%d jobs succeeded -- see [13_region_null_calibration] messages above for which ones failed; refusing to pool a partial result set", nrow(dt), nrow(jobs)))
  say("[2] %d/%d jobs succeeded in %.1f min (%.3fs/job)\n", nrow(dt), nrow(jobs),
      as.numeric(difftime(Sys.time(), t0, units = "mins")),
      as.numeric(difftime(Sys.time(), t0, units = "secs")) / nrow(jobs))
  print(dt[, .(mean_null = mean(mean_null), max_null = max(max_null), n = .N), by = .(method, scheme)])
  saveRDS(dt, sprintf("out_final_v1/13_null_calib_%s_perjob.rds", MODE))

  ## =============================================================================
  ## [3] Pooling: join per-(combo,method,scheme) E[N_null] against the
  ## OBSERVED Stage-2 region counts and known TP/FP, READ from
  ## results/simulation_stage2_region_details.rds (R/11) -- never
  ## recomputed here (satisfies "observed Stage-2 counts in the null
  ## analysis equal the corresponding observed truth-scoring counts" by
  ## construction, not by a second implementation).
  ## =============================================================================
  say("\n[3] pooling against observed Stage-2 truth (R/11's region details)\n")
  region_detail_file <- "results/simulation_stage2_region_details.rds"
  if (!file.exists(region_detail_file))
    stop("R/11_stage2_region_details.R has not produced: ", region_detail_file, " -- run it first.")
  rd <- readRDS(region_detail_file)
  method_map <- c(emmax_consensus = "emmax_consensus_region", emmax_simes = "emmax_simes_region")
  obs <- rd[method %in% method_map, .(n_obs_regions = .N, TP = sum(TP), FP = sum(!TP)),
           by = .(tag, cell, rep, env, method)]
  obs[, method := names(method_map)[match(method, method_map)]]

  ## every combo x method x scheme this run actually computed
  full <- merge(dt, obs, by = c("tag", "cell", "rep", "env", "method"), all.x = TRUE)
  full[is.na(n_obs_regions), `:=`(n_obs_regions = 0L, TP = 0L, FP = 0L)]   ## combo scored 0 regions for this method -- a real zero, not missing

  ## per-env (NOT pooled across the 10 envs) copy, kept for analyses that need
  ## per-env variation rather than the map/burn-in-level pooled summary below
  ## (e.g. environment-structure-alignment, PK 2026-09-18 -- see
  ## env-structure-alignment-hypothesis.md; pooling across envs would average
  ## away exactly the per-env variation that analysis needs).
  full[, ratio_null_obs := ifelse(n_obs_regions > 0, mean_null / n_obs_regions, NA_real_)]
  full[, FDP_truth := ifelse((TP + FP) > 0, FP / (TP + FP), NA_real_)]
  fwrite(full, "results/simulation_null_truth_calibration_by_env.tsv", sep = "\t")
  say("[3b] wrote results/simulation_null_truth_calibration_by_env.tsv (%d rows, per-env not pooled)\n", nrow(full))

  ## ---- map/burn-in level: pool the (up to) 10 environmental continuations ----
  say("[4] map/burn-in-level pooling (%d combos retained %d envs each on average)\n",
      uniqueN(full[, .(tag, cell, rep)]), round(nrow(full) / uniqueN(full[, .(tag, cell, rep, method, scheme)])))
  map_level <- full[, .(sum_E_null = sum(mean_null), sum_n_obs = sum(n_obs_regions),
                        sum_TP = sum(TP), sum_FP = sum(FP), n_envs = .N),
                   by = .(tag, cell, rep, method, scheme)]
  ## R_null/obs, r -- NEVER max(n_obs,1): a zero-observed map/burn-in gets NA,
  ## retained and counted, not floored or dropped (mandatory checks 7-8).
  map_level[, R_null_obs := ifelse(sum_n_obs > 0, sum_E_null / sum_n_obs, NA_real_)]
  map_level[, FDP_truth := ifelse((sum_TP + sum_FP) > 0, sum_FP / (sum_TP + sum_FP), NA_real_)]
  frac_zero_obs <- map_level[, .(frac_zero_obs_regions = mean(sum_n_obs == 0)), by = .(method, scheme)]
  say("    zero-observed-region map/burn-in groups: %s\n",
      paste(sprintf("%s/%s=%.1f%%", frac_zero_obs$method, frac_zero_obs$scheme,
                    100 * frac_zero_obs$frac_zero_obs_regions), collapse = "; "))

  fwrite(map_level, "results/simulation_null_truth_calibration.tsv", sep = "\t")
  say("[5] wrote results/simulation_null_truth_calibration.tsv (%d rows)\n", nrow(map_level))

  ## ---- summary per (method, scheme): calibration diagnostics + map-cluster bootstrap --
  say("\n[6] calibration summary per (method, scheme)\n")
  summary_rows <- list()
  for (m in NULL_METHODS) for (s in NULL_SCHEMES) {
    ml <- map_level[method == m & scheme == s]
    complete <- ml[!is.na(R_null_obs) & !is.na(FDP_truth)]
    rho <- if (nrow(complete) >= 3) stats::cor(complete$R_null_obs, complete$FDP_truth, method = "spearman") else NA_real_
    fit <- if (nrow(complete) >= 3) stats::lm(FDP_truth ~ R_null_obs, data = complete) else NULL
    slope <- if (!is.null(fit)) unname(coef(fit)[2]) else NA_real_
    intercept <- if (!is.null(fit)) unname(coef(fit)[1]) else NA_real_
    mean_signed_diff <- if (nrow(complete)) mean(complete$R_null_obs - complete$FDP_truth) else NA_real_
    median_abs_diff <- if (nrow(complete)) stats::median(abs(complete$R_null_obs - complete$FDP_truth)) else NA_real_

    ## cell(+tag)-adjusted correlation (PK, 2026-09-18 review, point 3):
    ## residualise RANKS of R_null_obs and FDP_truth against demographic
    ## cell and BGS treatment before correlating, so an apparent association
    ## driven purely by baseline differences between cells (or between BGS
    ## treatments) doesn't get reported as evidence of null/truth
    ## calibration. NA when too few cells vary to fit the adjustment.
    adj_rho <- NA_real_
    if (nrow(complete) >= 3 && uniqueN(complete$cell) >= 2) {
      rx <- rank(complete$R_null_obs); ry <- rank(complete$FDP_truth)
      rx_resid <- stats::resid(stats::lm(rx ~ factor(cell) + factor(tag), data = complete))
      ry_resid <- stats::resid(stats::lm(ry ~ factor(cell) + factor(tag), data = complete))
      adj_rho <- suppressWarnings(stats::cor(rx_resid, ry_resid))
    }

    ## map-cluster bootstrap CI on the grand-pooled R_null/obs and FDP_truth.
    ## [!] Clustered by (cell, rep) ONLY, NOT (tag, cell, rep) (PK,
    ## 2026-09-18 review, point 3/5 -- same principle as R/14_summarise_
    ## additional_analyses.R's .fit_size_model() fix): BGS and no-BGS rows
    ## for the same (cell, rep) are a PAIRED map/burn-in draw and must be
    ## resampled together, never independently. Each cluster's bgs/no-bgs
    ## rows are pre-summed into one row before resampling -- equivalent to
    ## always keeping both tag rows in or out of a bootstrap draw together.
    cluster_level <- ml[, .(sum_E_null = sum(sum_E_null), sum_n_obs = sum(sum_n_obs),
                            sum_TP = sum(sum_TP), sum_FP = sum(sum_FP)), by = .(cell, rep)]
    n_r <- nrow(cluster_level)
    mat <- as.matrix(cluster_level[, .(sum_E_null, sum_n_obs, sum_TP, sum_FP)])
    bs <- if (n_r >= 2) bootstrap_rep_matrix(mat, N_BOOTSTRAP, SEEDS[["bootstrap"]]) else NULL
    ratio_b <- if (!is.null(bs)) bs[, "sum_E_null"] / bs[, "sum_n_obs"] else NA_real_
    fdp_b <- if (!is.null(bs)) bs[, "sum_FP"] / (bs[, "sum_TP"] + bs[, "sum_FP"]) else NA_real_
    ratio_b <- ratio_b[is.finite(ratio_b)]; fdp_b <- fdp_b[is.finite(fdp_b)]

    ## null_role: MVN is a NEGATIVE CONTROL (matched kinship structure only,
    ## no spatial/environmental signal -- it should show weak/no calibration
    ## if the pipeline is sound) not a candidate FDP estimator; group/spatial
    ## are alternative BIOLOGICAL null recipes, not interchangeable
    ## estimates of the true FDP (PK, 2026-09-18 review, point 3).
    null_role <- c(group = "candidate_null", mvn = "negative_control", spatial = "candidate_null")[[s]]

    summary_rows[[length(summary_rows) + 1]] <- data.table(
      method = m, scheme = s, null_role = null_role,
      n_map_groups = nrow(ml), n_complete = nrow(complete),
      frac_zero_obs_regions = mean(ml$sum_n_obs == 0),
      spearman_rho = rho, cell_tag_adjusted_rho = adj_rho, slope = slope, intercept = intercept,
      mean_signed_diff = mean_signed_diff, median_abs_diff = median_abs_diff,
      pooled_E_null = sum(ml$sum_E_null), pooled_n_obs = sum(ml$sum_n_obs),
      pooled_TP = sum(ml$sum_TP), pooled_FP = sum(ml$sum_FP),
      pooled_R_null_obs = if (sum(ml$sum_n_obs) > 0) sum(ml$sum_E_null) / sum(ml$sum_n_obs) else NA_real_,
      pooled_R_null_obs_ci_lo = if (length(ratio_b)) ci_quantile(ratio_b)[1] else NA_real_,
      pooled_R_null_obs_ci_hi = if (length(ratio_b)) ci_quantile(ratio_b)[2] else NA_real_,
      pooled_FDP_truth = if (sum(ml$sum_TP) + sum(ml$sum_FP) > 0) sum(ml$sum_FP) / (sum(ml$sum_TP) + sum(ml$sum_FP)) else NA_real_,
      pooled_FDP_truth_ci_lo = if (length(fdp_b)) ci_quantile(fdp_b)[1] else NA_real_,
      pooled_FDP_truth_ci_hi = if (length(fdp_b)) ci_quantile(fdp_b)[2] else NA_real_)
  }
  summary_dt <- rbindlist(summary_rows)
  fwrite(summary_dt, "results/simulation_null_truth_summary.tsv", sep = "\t")
  say("[7] wrote results/simulation_null_truth_summary.tsv\n")
  print(summary_dt[, .(method, scheme, null_role, n_map_groups, frac_zero_obs_regions,
                       spearman_rho, cell_tag_adjusted_rho, slope, pooled_R_null_obs, pooled_FDP_truth)])

  write_receipt("13_region_null_calibration", inputs = region_detail_file,
                params = list(mode = MODE, B = B, null_methods = NULL_METHODS, null_schemes = NULL_SCHEMES,
                              alpha = ALPHA, size_floor = SIZE_FLOOR, n_bootstrap = N_BOOTSTRAP,
                              seed_bootstrap = SEEDS[["bootstrap"]], seed_offsets = SCHEME_SEED_OFFSET,
                              null_calib_version = NULL_CALIB_VERSION),
                outputs = c("results/simulation_null_truth_calibration.tsv", "results/simulation_null_truth_summary.tsv",
                           "results/simulation_null_truth_calibration_by_env.tsv"))
  cat("NULLCALIB_DONE\n")
}
