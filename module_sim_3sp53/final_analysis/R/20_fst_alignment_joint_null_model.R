## =============================================================================
## final_analysis/R/20_fst_alignment_joint_null_model.R
##
## EXPLORATORY -- not wired into the receipt/stage_stale pipeline (see header
## caveat at the bottom). Follow-up to R/19's pooled Fst/LD correlations
## (PK, 2026-09-26 review):
##
##   "One provenance concern needs resolving before promotion: R/19 reads LD
##   estimates from the older module_sim bundles, whereas the null outcomes
##   and alignment estimates use the corrected final-analysis bundles...
##   I would not interpret those LD correlations as final until a_pred and b
##   are recomputed from the same corrected bundles. My suggested next step
##   is one joint, map-cluster-bootstrapped analysis of expected spatial-
##   null count using FST and environment-structure alignment together,
##   adjusted for cell, BGS treatment and method. Report c=1 separately and
##   treat c=1.5/2 as stress tests. That would tell us whether genomic
##   differentiation adds information beyond alignment, rather than adding
##   another strong-looking pooled correlation to the manuscript."
##
## This script does exactly that, and deliberately drops LD-decay-rate/
## background-LD entirely -- PK's own suggested next step does not ask for
## them, and R/19's LD numbers are flagged as not-provenance-clean (see that
## script's updated header). Fst is retained: unlike LD-decay, Fst depends
## only on genotype dosage x population label, not marker sort-order, so
## the pre-correction-parser popgen_summary.rds is safe to join onto
## corrected-parser combo IDs (established in R/19's header; independently,
## module_sim_bgs5/README.md reports Fst r=1.00 between pipeline
## generations on shared cells).
##
## Inputs, BOTH already computed from the CORRECTED final-analysis pipeline
## (no provenance mismatch, unlike R/19's LD predictors):
##   - results/simulation_env_structure_alignment_full_grid.tsv (R/15) --
##     r2_axes/mantel_r (alignment) and mean_null (expected SPATIAL-null
##     region count; already scheme = "spatial" specific, already zero-
##     filled for combos with no reported regions -- a real zero, not a
##     missing value).
##   - module_sim_3sp53/results/popgen_summary.rds$per_replicate -- Fst,
##     joined by (tag, cell, rep, env). Full 7-cell grid.
##
## Model: log2(mean_null + 1) ~ r2_axes + Fst [+ method_f + cell_f + tag_f],
## fit three ways per scope (Fst-only, alignment-only, joint) so the joint
## model's own coefficients can be read against their marginal counterparts
## -- if Fst's coefficient collapses toward 0 (and its CI widens to include
## 0) once r2_axes is in the model, Fst was mostly a proxy for alignment;
## if it survives, Fst adds independent information.
##
## Scope, per PK's request: "c1" (the 3 high-dispersal cells, V0.5_c1/
## V1_c1/V2_c1 -- CELLS_ALL's own "_c1" suffix) is the PRIMARY result;
## "c1.5_c2" (the remaining 4 cells) is reported as a STRESS TEST, not
## pooled into the headline; "pooled_all_cells" is reported too, for direct
## comparison against R/19's pooled numbers, but is not the primary result.
##
## [!] FIXED (PK, 2026-09-26 second review) -- TWO issues in the first
## version of this script:
##
## 1. The nested F-test (previously here, `anova(align_only, joint,
##    test="F")`) treats every row as an independent observation. It is not:
##    rows are repeated measures across 10 environmental continuations, 2
##    methods, and (see point 2) a shared map/burn-in across cells. Those
##    p-values were invalid and have been REMOVED, not just re-labelled --
##    replaced below with a map-cluster-bootstrapped delta-R^2 (joint model
##    R^2 minus alignment-only R^2), which answers the same "does Fst add
##    explanatory power" question but with a percentile CI built from the
##    same resampling scheme as the coefficient CIs, rather than a
##    parametric test that assumes independence this design doesn't have.
##
## 2. `.boot_slope()` previously clustered by (cell, rep). That is too fine:
##    per `~/gitlab/LDscnR-NEMO/make_prod.sh`'s STEP 1, the genetic maps /
##    QTN positions / environmental values / dispersal template (10 sets,
##    one per `rep`) are built ONCE and reused identically across ALL 7
##    cells and both tags -- confirmed directly in that script, not
##    inferred -- so rows sharing a `rep` value are not independent even
##    when they come from different cells (selection intensity V is applied
##    downstream of this shared template; only the per-cell/per-tag adapt-
##    phase RNG seed differs, per that script's own seed formula). The
##    analysis-plan document's own "Pooling and uncertainty" section
##    confirms this framing ("cluster-bootstrap the ten map/burn-in IDs...
##    retain all ten environmental continuations... keep all methods
##    paired"): TEN clusters, not up to seventy. `.boot_slope()` below now
##    clusters by `rep` alone. Checked empirically before trusting this:
##    within the c1_primary scope, a given `rep` value is NOT always present
##    in all 3 cells (e.g. rep=8/9/10 appear only in V1_c1 in this balanced
##    5-reps-per-cell dataset) -- clusters are uneven in size, which a
##    cluster bootstrap handles correctly, but it means "10 clusters" is an
##    upper bound, not always achieved (9 distinct rep values survive in
##    c1_primary here). PK independently spot-checked this fix by
##    resampling by rep across the c=1 cells and got a still-positive but
##    wider Fst interval (~5.23-20.37) than the (cell,rep)-clustered
##    version's [3.97, 19.83] -- consistent with fewer, coarser independent
##    units, and the primary finding (Fst survives, alignment doesn't, in
##    c1_primary) does not reverse.
##
## `.boot_slope()`/`.prep()` are sourced from R/15_env_structure_alignment.R
## (with the clustering fix above applied), where both now (2026-09-26) live
## at top level for exactly this kind of reuse; R/15 itself and its saved
## models.tsv have since been refreshed with the same fix (R/15b_refresh_
## alignment_models_after_cluster_fix.R).
##
## 3. [!] FIXED (PK, 2026-09-26 third review) -- the delta-R^2 in fix note 1
##    above is NOT valid evidence that Fst adds information, and is
##    downgraded here to a descriptive-only statistic. Both nested models
##    were fit to the SAME (bootstrap-resampled) observations in every
##    replicate: adding a predictor to an OLS fit cannot lower its IN-SAMPLE
##    R^2, so a delta-R^2 interval excluding zero is close to guaranteed
##    even for a completely uninformative added predictor -- it was never a
##    test of anything, however clustering-aware its resampling was. Kept
##    below (renamed `delta_r2_insample`) purely as a descriptive fit-
##    improvement measure, printed alongside an explicit "not evidence"
##    label, NOT as a result to interpret. The actual evidence for "does
##    Fst add information beyond alignment" is the ADJUSTED Fst coefficient
##    CI from fix note 2's model_dt (scope c1_primary, variant "joint",
##    predictor "Fst") -- that CI can, in principle, include zero, which is
##    what makes it a real test. A stronger check than either -- genuinely
##    OUT-OF-SAMPLE, holding out whole map IDs so a spurious predictor CAN
##    show up worse, not just no-better -- is added below as leave-rep-out
##    cross-validation (`.loro_delta_r2()`): refit align_only/joint on all
##    reps but one, predict the held-out rep's rows, pool held-out residuals
##    across all folds, and compare OOS R^2 between the two models. A CI on
##    the OOS delta-R^2 comes from resampling WHICH folds' already-computed
##    held-out residuals get pooled (fast -- no per-bootstrap-replicate
##    refitting needed, since the n_reps fold-level fits are each computed
##    exactly once).
## =============================================================================
suppressMessages({library(data.table)})
MODULE_ROOT <- path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis")
source(file.path(MODULE_ROOT, "R", "00_config.R"))
source(file.path(MODULE_ROOT, "R", "helpers_stage2_truth.R"))   ## ci_quantile()
source(file.path(MODULE_ROOT, "R", "15_env_structure_alignment.R"))   ## .prep()/.boot_slope() -- now top-level there (2026-09-26), does not rerun 15's own guarded analysis
say("=== 20_fst_alignment_joint_null_model (EXPLORATORY) ===\n\n")

## ---- descriptive-only in-sample delta-R^2 (NOT evidence -- see header fix --
## note 3). Kept for completeness/comparison to what an unwary reader might
## have computed, printed with an explicit disclaimer, never as a test.
.boot_delta_r2_insample <- function(fit_fn_a, fit_fn_b, data) {
  data <- copy(data)
  clusters <- sort(unique(data$rep))
  n_cl <- length(clusters)
  cl_rows <- split(seq_len(nrow(data)), data$rep)
  r2 <- function(fit) summary(fit)$r.squared
  fit_a0 <- tryCatch(fit_fn_a(data), error = function(e) NULL)
  fit_b0 <- tryCatch(fit_fn_b(data), error = function(e) NULL)
  if (is.null(fit_a0) || is.null(fit_b0))
    return(data.table(r2_a = NA_real_, r2_b = NA_real_, delta_r2 = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, n = nrow(data), n_boot_converged = 0L))
  r2_a0 <- r2(fit_a0); r2_b0 <- r2(fit_b0); point <- r2_b0 - r2_a0
  set.seed(SEEDS[["bootstrap"]])
  boot <- rep(NA_real_, N_BOOTSTRAP)
  for (bb in seq_len(N_BOOTSTRAP)) {
    draw <- sample(clusters, n_cl, replace = TRUE)
    idx <- unlist(cl_rows[as.character(draw)], use.names = FALSE)
    d_b <- data[idx]
    fa <- tryCatch(fit_fn_a(d_b), error = function(e) NULL, warning = function(w) NULL)
    fb <- tryCatch(fit_fn_b(d_b), error = function(e) NULL, warning = function(w) NULL)
    if (!is.null(fa) && !is.null(fb)) boot[bb] <- tryCatch(r2(fb) - r2(fa), error = function(e) NA_real_)
  }
  boot <- boot[is.finite(boot)]
  data.table(r2_a = r2_a0, r2_b = r2_b0, delta_r2 = point,
            ci_lo = if (length(boot)) ci_quantile(boot)[1] else NA_real_,
            ci_hi = if (length(boot)) ci_quantile(boot)[2] else NA_real_,
            n = nrow(data), n_boot_converged = length(boot))
}

## ---- leave-rep-out (LORO) out-of-sample delta-R^2 -- the ACTUAL evidence --
## for "does Fst add predictive information beyond alignment" (header fix
## note 3). Each rep's rows are held out exactly once; align_only/joint are
## refit on every OTHER rep and used to predict the held-out rep -- a
## spurious predictor CAN and typically does make out-of-sample R^2 worse,
## unlike in-sample R^2, which is why this is a real test and
## .boot_delta_r2_insample() above is not. Each fold's model is fit exactly
## once (n_reps refits per model, not thousands); the CI on the pooled OOS
## delta-R^2 comes from resampling WHICH folds' already-computed held-out
## residuals get pooled, not from refitting per bootstrap draw.
.loro_delta_r2 <- function(fit_fn_a, fit_fn_b, data) {
  data <- copy(data)
  reps <- sort(unique(data$rep))
  y_bar <- mean(data$log2_null)   ## fixed baseline for SS_tot, shared by both models and all folds
  fold_ss <- vector("list", length(reps))
  for (i in seq_along(reps)) {
    held_out <- data$rep == reps[i]
    train <- data[!held_out]; test <- data[held_out]
    fa <- tryCatch(fit_fn_a(train), error = function(e) NULL)
    fb <- tryCatch(fit_fn_b(train), error = function(e) NULL)
    if (is.null(fa) || is.null(fb) || !nrow(test)) { fold_ss[[i]] <- NULL; next }
    pa <- tryCatch(predict(fa, newdata = test), error = function(e) rep(NA_real_, nrow(test)))
    pb <- tryCatch(predict(fb, newdata = test), error = function(e) rep(NA_real_, nrow(test)))
    ## a fold where prediction fails outright (e.g. an unseen factor level --
    ## checked, does not occur in practice here, but guarded regardless) must
    ## be DROPPED, not scored with na.rm=TRUE -- that would silently count a
    ## fully-failed prediction as zero squared error, inflating OOS R^2.
    if (anyNA(pa) || anyNA(pb)) { fold_ss[[i]] <- NULL; next }
    fold_ss[[i]] <- data.table(rep = reps[i],
                               sse_a = sum((test$log2_null - pa)^2),
                               sse_b = sum((test$log2_null - pb)^2),
                               sst = sum((test$log2_null - y_bar)^2), n = nrow(test))
  }
  fold_dt <- rbindlist(fold_ss)
  if (!nrow(fold_dt) || sum(fold_dt$sst) == 0)
    return(data.table(oos_r2_a = NA_real_, oos_r2_b = NA_real_, oos_delta_r2 = NA_real_,
                      ci_lo = NA_real_, ci_hi = NA_real_, n_folds = nrow(fold_dt)))
  pooled_r2 <- function(fd) data.table(r2_a = 1 - sum(fd$sse_a) / sum(fd$sst), r2_b = 1 - sum(fd$sse_b) / sum(fd$sst))
  point <- pooled_r2(fold_dt)
  n_fo <- nrow(fold_dt)
  set.seed(SEEDS[["bootstrap"]])
  boot <- rep(NA_real_, N_BOOTSTRAP)
  for (bb in seq_len(N_BOOTSTRAP)) {
    draw <- sample.int(n_fo, n_fo, replace = TRUE)
    pr <- pooled_r2(fold_dt[draw])
    boot[bb] <- pr$r2_b - pr$r2_a
  }
  boot <- boot[is.finite(boot)]
  data.table(oos_r2_a = point$r2_a, oos_r2_b = point$r2_b, oos_delta_r2 = point$r2_b - point$r2_a,
            ci_lo = if (length(boot)) ci_quantile(boot)[1] else NA_real_,
            ci_hi = if (length(boot)) ci_quantile(boot)[2] else NA_real_, n_folds = n_fo)
}

C1_CELLS      <- c("V0.5_c1", "V1_c1", "V2_c1")
C1_5_C2_CELLS <- setdiff(CELLS_ALL, C1_CELLS)

## ---- 1. load already-corrected-pipeline alignment + null-burden table ------
full_grid <- fread(file.path(MODULE_ROOT, "results", "simulation_env_structure_alignment_full_grid.tsv"))
say("[1] loaded %d rows from R/15's full_grid (%d combos x %d methods)\n",
    nrow(full_grid), nrow(unique(full_grid[, .(tag, cell, rep, env)])), uniqueN(full_grid$method))

## ---- 2. join Fst (full 7-cell grid, pre-correction-parser but genotype/pop- --
## label-based -- safe; see header) -------------------------------------------
popgen <- readRDS(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53"), "results", "popgen_summary.rds"))
fst_dt <- popgen$per_replicate[, .(tag, cell, rep, env, Fst)]
dt <- merge(full_grid, fst_dt, by = c("tag", "cell", "rep", "env"), all.x = TRUE)
if (anyNA(dt$Fst)) stop("Fst missing for some full_grid rows -- popgen_summary.rds should cover the full 7-cell grid")
dt[, log2_null := log2(mean_null + 1)]
dt <- .prep(dt)   ## adds cell_f/tag_f/method_f, from R/15
say("[2] joined Fst (no missing values) -> %d rows\n\n", nrow(dt))

## ---- 3. joint vs. marginal models, per scope --------------------------------
FIT_SPECS <- list(
  fst_only   = function(d) stats::lm(log2_null ~ Fst + method_f + cell_f + tag_f, data = d),
  align_only = function(d) stats::lm(log2_null ~ r2_axes + method_f + cell_f + tag_f, data = d),
  joint      = function(d) stats::lm(log2_null ~ r2_axes + Fst + method_f + cell_f + tag_f, data = d)
)
## within a single-cell scope, cell_f has one level and drops out on its own
## (model.matrix silently omits a constant column) -- no special-casing needed.

model_rows <- list()
insample_rows <- list()
loro_rows <- list()
for (scope_name in c("c1_primary", "c1.5_c2_stress_test", "pooled_all_cells")) {
  cells_here <- switch(scope_name, c1_primary = C1_CELLS, c1.5_c2_stress_test = C1_5_C2_CELLS, pooled_all_cells = CELLS_ALL)
  d_scope <- dt[cell %in% cells_here]
  d_scope[, cell_f := droplevels(cell_f)]

  fits <- lapply(FIT_SPECS, function(f) f(d_scope))
  names(fits) <- names(FIT_SPECS)

  ## bootstrapped coefficient CIs -- Fst's coefficient from BOTH the
  ## Fst-only fit (marginal) and the joint fit (adjusted for r2_axes); same
  ## for r2_axes. Read the pair side by side: shrinkage from marginal to
  ## joint is the "does it survive adjustment" answer, per predictor. THIS
  ## is the real evidence for the coefficient-level question (can include
  ## zero -- a genuine test), alongside the OOS delta-R^2 below.
  for (variant in c("fst_only", "align_only", "joint")) {
    fit_fn <- FIT_SPECS[[variant]]
    if ("Fst" %in% names(coef(fits[[variant]])))
      model_rows[[length(model_rows) + 1]] <- cbind(
        scope = scope_name, variant = variant, predictor = "Fst",
        .boot_slope(fit_fn, function(fit) unname(coef(fit)["Fst"]), d_scope))
    if ("r2_axes" %in% names(coef(fits[[variant]])))
      model_rows[[length(model_rows) + 1]] <- cbind(
        scope = scope_name, variant = variant, predictor = "r2_axes",
        .boot_slope(fit_fn, function(fit) unname(coef(fit)["r2_axes"]), d_scope))
  }

  ## descriptive-only in-sample delta-R^2 -- NOT evidence (header fix note 3).
  insample_rows[[length(insample_rows) + 1]] <- cbind(
    scope = scope_name, comparison = "joint_vs_align_only",
    .boot_delta_r2_insample(FIT_SPECS$align_only, FIT_SPECS$joint, d_scope))

  ## the actual evidence for "does Fst add information beyond alignment" --
  ## leave-rep-out out-of-sample delta-R^2.
  loro_rows[[length(loro_rows) + 1]] <- cbind(
    scope = scope_name, comparison = "joint_vs_align_only",
    .loro_delta_r2(FIT_SPECS$align_only, FIT_SPECS$joint, d_scope))
}
model_dt <- rbindlist(model_rows)
setcolorder(model_dt, c("scope", "variant", "predictor", "estimate", "ci_lo", "ci_hi", "n", "n_boot_converged"))
insample_dt <- rbindlist(insample_rows)
setnames(insample_dt, c("r2_a", "r2_b"), c("r2_align_only", "r2_joint"))
loro_dt <- rbindlist(loro_rows)

fwrite(model_dt, file.path(MODULE_ROOT, "results", "exploratory_fst_alignment_joint_models.tsv"), sep = "\t")
fwrite(insample_dt, file.path(MODULE_ROOT, "results", "exploratory_fst_alignment_joint_delta_r2_insample.tsv"), sep = "\t")
fwrite(loro_dt, file.path(MODULE_ROOT, "results", "exploratory_fst_alignment_joint_delta_r2_oos.tsv"), sep = "\t")
say("[3] wrote results/exploratory_fst_alignment_joint_models.tsv, _delta_r2_insample.tsv, _delta_r2_oos.tsv\n\n")

say("--- coefficient estimates (log2(mean_null+1) scale), marginal vs. joint ---\n")
say("--- (map-cluster bootstrap, clustered by rep -- the shared map/burn-in ID) ---\n")
print(model_dt[order(scope, predictor, variant)])

say("\n--- IN-SAMPLE delta-R^2 -- DESCRIPTIVE ONLY, NOT EVIDENCE (adding a predictor cannot ---\n")
say("--- lower in-sample R^2 by construction; a positive interval here is expected regardless ---\n")
say("--- of whether Fst is truly informative -- see this script's header, fix note 3) ---\n")
print(insample_dt[order(scope), .(scope, comparison, r2_align_only, r2_joint, delta_r2, ci_lo, ci_hi, n, n_boot_converged)])

say("\n--- OUT-OF-SAMPLE delta-R^2 (leave-rep-out) -- the real evidence for whether Fst adds ---\n")
say("--- information beyond alignment: CAN be negative for an uninformative predictor ---\n")
print(loro_dt[order(scope), .(scope, comparison, oos_r2_a, oos_r2_b, oos_delta_r2, ci_lo, ci_hi, n_folds)])

say("\n--- read: scope should be judged from c1_primary first; c1.5_c2 is a stress test, ---\n")
say("--- not pooled into the headline; pooled_all_cells is for comparison to R/19 only ---\n")

cat("JOINT_MODEL_DONE\n")
