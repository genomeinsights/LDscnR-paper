## =============================================================================
## final_analysis/R/15b_refresh_alignment_models_after_cluster_fix.R
##
## One-off refresh, not a new pipeline stage: regenerates results/
## simulation_env_structure_alignment_models.tsv using R/15's FIXED
## `.boot_slope()` (clustered by `rep` alone, not (cell, rep) -- see that
## file's header/top-level note, PK's 2026-09-26 review), without rerunning
## R/15's own combo-level computation (`alignment_one_combo()`), which needs
## the 02_build_ld_units bundles -- NOT reachable on this host right now
## (external volume unmounted; see R/00_config.R's PATHS$parsed).
##
## Reads `gp` and `full_grid` back from the two TSVs R/15 already wrote
## (results/simulation_env_structure_alignment.tsv,
## _full_grid.tsv) -- both bundle-independent, already on disk, unaffected
## by this fix (the fix is only in how the bootstrap resamples rows it
## already has, not in how those rows were computed) -- and re-runs R/15's
## own modelling section [4]-[6] verbatim (copied, not reimplemented) against
## the now-fixed `.boot_slope()`/`.prep()`.
##
## When the bundles ARE reachable again (e.g. run on the mini), prefer
## rerunning R/15 itself end to end over this script -- this one exists only
## to verify and apply the clustering fix without a full recompute.
## =============================================================================
suppressMessages({library(data.table)})
MODULE_ROOT <- path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis")
source(file.path(MODULE_ROOT, "R", "00_config.R"))
source(file.path(MODULE_ROOT, "R", "helpers_stage2_truth.R"))   ## ci_quantile()
source(file.path(MODULE_ROOT, "R", "15_env_structure_alignment.R"))   ## fixed .prep()/.boot_slope() (top-level; guarded main block does not run)
say("=== 15b_refresh_alignment_models_after_cluster_fix ===\n\n")

gp        <- fread(file.path(MODULE_ROOT, "results", "simulation_env_structure_alignment.tsv"))
full_grid <- fread(file.path(MODULE_ROOT, "results", "simulation_env_structure_alignment_full_grid.tsv"))
say("[1] loaded cached gp (%d rows) and full_grid (%d rows), unchanged by this fix\n", nrow(gp), nrow(full_grid))

## ---- [4] slope models (copied verbatim from R/15, only .boot_slope's own --
## internals differ now) -------------------------------------------------------
gp <- .prep(gp)
full_grid <- .prep(full_grid)

ratio_d <- gp[ratio_null_obs > 0]
ratio_d[, log2_ratio := log2(ratio_null_obs)]
full_grid[, log2_null := log2(mean_null + 1)]

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

## within-cell rows are UNAFFECTED by this fix -- restricting to one cell
## already makes (cell,rep) and rep-alone identical clustering (cell is
## constant within the subset) -- recomputed anyway so the file is complete
## and internally consistent, not because they were expected to change.
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
fwrite(model_dt, file.path(MODULE_ROOT, "results", "simulation_env_structure_alignment_models.tsv"), sep = "\t")
say("[2] overwrote results/simulation_env_structure_alignment_models.tsv (rep-clustered bootstrap)\n")
print(model_dt)
say("\n[3] within-cell CIs excluding zero (flagging, not hiding, any that do):\n")
excl <- model_dt[variant == "within_cell" & !is.na(ci_lo) & !is.na(ci_hi) & (ci_lo > 0 | ci_hi < 0)]
if (nrow(excl)) print(excl[, .(outcome, cell, estimate, ci_lo, ci_hi)]) else say("    none\n")

cat("REFRESH_DONE\n")
