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
## if it survives, Fst adds independent information. A likelihood-ratio F
## test (joint vs. alignment-only) on the POINT-ESTIMATE fit answers the
## same question directly, reported once per scope, not bootstrapped.
##
## Scope, per PK's request: "c1" (the 3 high-dispersal cells, V0.5_c1/
## V1_c1/V2_c1 -- CELLS_ALL's own "_c1" suffix) is the PRIMARY result;
## "c1.5_c2" (the remaining 4 cells) is reported as a STRESS TEST, not
## pooled into the headline; "pooled_all_cells" is reported too, for direct
## comparison against R/19's pooled numbers, but is not the primary result.
##
## Map-cluster bootstrap (resample by (cell,rep), bgs/nobgs paired within a
## draw) copies R/15's own `.boot_slope()`/`.prep()` verbatim rather than
## sourcing them -- both are defined INSIDE R/15's `if (sys.nframe() == 0L)`
## guard (its unreachable-bundle-dependent main block), not at that file's
## top level, so `source()`-ing R/15 does not actually expose them (checked
## directly rather than assumed: it does not). `ci_quantile()` alone is
## genuinely reusable -- it lives at helpers_stage2_truth.R's top level --
## so that one is sourced, not copied.
## =============================================================================
suppressMessages({library(data.table)})
MODULE_ROOT <- path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis")
source(file.path(MODULE_ROOT, "R", "00_config.R"))
source(file.path(MODULE_ROOT, "R", "helpers_stage2_truth.R"))   ## ci_quantile()
say("=== 20_fst_alignment_joint_null_model (EXPLORATORY) ===\n\n")

## ---- copied verbatim from R/15_env_structure_alignment.R (see note above) --
.prep <- function(d) {
  d <- copy(d)
  d[, cell_f := factor(cell, levels = CELLS_ALL)]
  d[, tag_f := factor(tag, levels = TAGS_ALL)]
  d[, method_f := factor(method, levels = c("emmax_consensus", "emmax_simes"))]
  d
}
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
    idx <- unlist(cl_rows[as.character(draw)], use.names = FALSE)
    fit_b <- tryCatch(fit_fn(data[idx]), error = function(e) NULL, warning = function(w) NULL)
    if (!is.null(fit_b)) boot[bb] <- tryCatch(extract_fn(fit_b), error = function(e) NA_real_)
  }
  boot <- boot[is.finite(boot)]
  data.table(estimate = point, ci_lo = if (length(boot)) ci_quantile(boot)[1] else NA_real_,
            ci_hi = if (length(boot)) ci_quantile(boot)[2] else NA_real_,
            n = nrow(data), n_boot_converged = length(boot))
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
lrt_rows   <- list()
for (scope_name in c("c1_primary", "c1.5_c2_stress_test", "pooled_all_cells")) {
  cells_here <- switch(scope_name, c1_primary = C1_CELLS, c1.5_c2_stress_test = C1_5_C2_CELLS, pooled_all_cells = CELLS_ALL)
  d_scope <- dt[cell %in% cells_here]
  d_scope[, cell_f := droplevels(cell_f)]

  fits <- lapply(FIT_SPECS, function(f) f(d_scope))
  names(fits) <- names(FIT_SPECS)

  ## bootstrapped coefficient CIs -- Fst's coefficient from BOTH the
  ## Fst-only fit (marginal) and the joint fit (adjusted for r2_axes); same
  ## for r2_axes. Read the pair side by side: shrinkage from marginal to
  ## joint is the "does it survive adjustment" answer, per predictor.
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

  ## direct answer to "does Fst add information beyond alignment" -- nested F
  ## test on the POINT-ESTIMATE fits (not bootstrapped; a single well-defined
  ## test on the observed data, reported alongside the bootstrap CIs above).
  lrt <- stats::anova(fits$align_only, fits$joint, test = "F")
  lrt_rows[[length(lrt_rows) + 1]] <- data.table(
    scope = scope_name, comparison = "joint_vs_align_only",
    df = lrt$Df[2], F_stat = lrt$F[2], p_value = lrt[["Pr(>F)"]][2],
    r2_align_only = summary(fits$align_only)$r.squared, r2_joint = summary(fits$joint)$r.squared,
    n = nrow(d_scope))
}
model_dt <- rbindlist(model_rows)
setcolorder(model_dt, c("scope", "variant", "predictor", "estimate", "ci_lo", "ci_hi", "n", "n_boot_converged"))
lrt_dt <- rbindlist(lrt_rows)

fwrite(model_dt, file.path(MODULE_ROOT, "results", "exploratory_fst_alignment_joint_models.tsv"), sep = "\t")
fwrite(lrt_dt, file.path(MODULE_ROOT, "results", "exploratory_fst_alignment_joint_lrt.tsv"), sep = "\t")
say("[3] wrote results/exploratory_fst_alignment_joint_models.tsv, results/exploratory_fst_alignment_joint_lrt.tsv\n\n")

say("--- coefficient estimates (log2(mean_null+1) scale), marginal vs. joint ---\n")
print(model_dt[order(scope, predictor, variant)])

say("\n--- does Fst add information beyond alignment? nested F test, joint vs. align_only ---\n")
print(lrt_dt[order(scope)])

say("\n--- read: scope should be judged from c1_primary first; c1.5_c2 is a stress test, ---\n")
say("--- not pooled into the headline; pooled_all_cells is for comparison to R/19 only ---\n")

cat("JOINT_MODEL_DONE\n")
