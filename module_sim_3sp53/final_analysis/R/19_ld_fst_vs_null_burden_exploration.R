## =============================================================================
## final_analysis/R/19_ld_fst_vs_null_burden_exploration.R
##
## EXPLORATORY -- not wired into the receipt/stage_stale pipeline, and not
## cited anywhere yet. PK asked (2026-09-26): does background LD / the
## LD-decay rate, and separately mean pairwise Fst, predict how easily a
## structured null produces outlier regions? Pooled, and per cell. This
## script is a first look; keep or promote to a numbered, receipted stage
## only once it's judged useful enough to matter for the manuscript.
##
## [!] PROVENANCE FLAG (PK, 2026-09-26 review) -- DO NOT CITE THE LD-DECAY-
## RATE / BACKGROUND-LD RESULTS BELOW AS FINAL. This script's `ld_decay_rate`
## and `background_ld` columns come from module_sim/'s OLDER, PRE-CORRECTION-
## PARSER bundles (see section 3 below), while every other input here --
## `mean_null`, `ratio_null_obs`, `FDP_truth`, and Fst's join target -- comes
## from the CORRECTED final-analysis pipeline. Per final_analysis/README.md's
## offset_propagation_check section (qc/offset_propagation_check.R, line
## ~74-88), the parser fix moves QTN/marker sort-order, and
## `ld_complexity_reduction()`'s LD-decay windowing is explicitly ORDER-
## dependent -- unlike Fst, which depends only on genotype dosage x
## population label and is safe to mix across parser generations (see
## section 2's own note). So the `ld_decay_rate`/`background_ld` rows in
## this script's correlation table are a genuine, unresolved provenance
## mismatch, not just an untidy label -- do not promote, cite, or build on
## them until `a_pred`/`b` are recomputed from final_analysis's OWN
## corrected per-combo bundles (needs either the external volume PATHS
## $parsed points at, or a fresh 01_parse_nemo.R + LD-decay-only rerun
## against /Volumes/Nemo's raw archives). The Fst results in this same
## script are NOT affected by this flag. See R/20_fst_alignment_joint_
## null_model.R for the follow-up analysis PK asked for instead, which
## drops LD-decay/background-LD entirely rather than trying to patch this.
##
## "How easily the null produces outlier regions" is read literally as the
## null's own EXPECTED region count per surrogate (`mean_null` in
## R/13_region_null_calibration.R's output) -- NOT `ratio_null_obs`
## (null/observed), because that ratio's denominator is the OBSERVED count,
## which is a property of the real signal, not of the null. A predictor that
## makes the null noisier should show up in `mean_null` directly. Both
## `ratio_null_obs` and `FDP_truth` (empirical false-discovery proportion
## against real QTN truth) are reported too, as secondary framings of the
## same question, with the same "never impute a zero denominator to 1"
## convention documented in R/13 and module_sim_3sp53/README.md's 2026-09-09
## correction -- rows with a zero denominator are NA, not floored, and how
## often that happens is reported per stratum.
##
## Data sources and important COVERAGE ASYMMETRY between the two predictors:
##   - Null-discovery outcomes: results/simulation_null_truth_calibration_by_env.tsv
##     (R/13) -- all 7 cells x 2 tags, 5 reps/cell (the "balanced" run), 10
##     envs, emmax_consensus/emmax_simes x group/mvn/spatial. `mvn` is a
##     NEGATIVE CONTROL (matched kinship, no real spatial/environmental
##     signal) -- reported separately, never pooled with group/spatial.
##   - Fst: module_sim_3sp53/results/popgen_summary.rds$per_replicate --
##     FULL 1400-combo grid (all 7 cells), computed by the pre-correction-
##     parser pipeline, but genotype dosages and population labels (all Fst
##     depends on) are untouched by that parser's marker-position bug --
##     see final_analysis/README.md's offset-propagation-check section, and
##     module_sim_bgs5/README.md's independent Fst r=1.00 cross-check. Safe
##     to join onto final_analysis's corrected-parser combo IDs.
##   - LD-decay rate (`a_pred`) and background LD (`b`): final_analysis's OWN
##     per-combo LD_decay objects are NOT reachable on this host right now
##     (PATHS$parsed lives on an unmounted external volume; out_final_v1/
##     only has 18_bgs_marker_noise/ materialised). The only per-combo
##     LD_decay actually on disk here is module_sim/out/02_bundle/'s bundles
##     (symlinked to /Volumes/Nemo, which IS mounted) -- but module_sim only
##     ever built 4 of CELLS_ALL's 7 cells (V0.5_c1, V0.5_c2, V1_c1.5,
##     V2_c1; see module_sim/R/00_config.R:74). So the LD-decay-rate/
##     background-LD analysis below covers only those 4 cells, while the
##     Fst analysis covers all 7 -- reported explicitly, not silently
##     truncated to match.
##
## Known caveats NOT modelled here (exploratory, so kept simple -- revisit
## if this gets promoted):
##   - bgs/nobgs share the same map/burn-in draw per (cell,rep) (R/13's own
##     bootstrap resamples them together for this reason) -- `tag` is kept
##     as a column throughout so this is inspectable, but the correlations
##     below pool bgs+nobgs rather than block-resampling the pair.
##   - Method (emmax_consensus/emmax_simes) and scheme (group/mvn/spatial)
##     are kept as separate strata rather than averaged together, since
##     they are genuinely different null recipes/statistics, not repeat
##     measures of the same thing.
##   - Plain Spearman correlation, no map-cluster bootstrap CI (contrast
##     with R/13's own summary table) -- a first look, not an inferential
##     claim.
## =============================================================================
suppressMessages({library(data.table); library(ggplot2)})
MODULE_ROOT <- path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis")
source(file.path(MODULE_ROOT, "R", "00_config.R"))
say("=== 19_ld_fst_vs_null_burden_exploration (EXPLORATORY) ===\n\n")

BUNDLE_DIR   <- path.expand("~/gitlab/LDscnR-paper/module_sim/out/02_bundle")
LD_DECAY_CELLS <- c("V0.5_c1", "V0.5_c2", "V1_c1.5", "V2_c1")

## ---- 1. null-discovery outcomes ---------------------------------------------
null_dt <- fread(file.path(MODULE_ROOT, "results", "simulation_null_truth_calibration_by_env.tsv"))
say("[1] loaded %d null-calibration rows (%d distinct tag/cell/rep/env combos)\n",
    nrow(null_dt), nrow(unique(null_dt[, .(tag, cell, rep, env)])))

## ---- 2. Fst (full 7-cell grid) ----------------------------------------------
popgen <- readRDS(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53"), "results", "popgen_summary.rds"))
fst_dt <- popgen$per_replicate[, .(tag, cell, rep, env, Fst)]
say("[2] loaded Fst for %d combos (all %d cells)\n", nrow(fst_dt), uniqueN(fst_dt$cell))

## ---- 3. LD-decay rate + background LD, read once per distinct combo --------
## Only for combos this script actually needs (present in null_dt), restricted
## to the 4 cells module_sim actually built -- see header.
say("[!] PROVENANCE WARNING: ld_decay_rate/background_ld below come from the\n")
say("    OLDER pre-correction-parser module_sim bundles, not final_analysis's\n")
say("    own corrected bundles -- do not cite as final. See this script's\n")
say("    header and R/20_fst_alignment_joint_null_model.R.\n")
combos <- unique(null_dt[cell %in% LD_DECAY_CELLS, .(tag, cell, rep, env)])
say("[3] reading LD_decay from %d module_sim bundles (%d cells covered: %s)...\n",
    nrow(combos), length(LD_DECAY_CELLS), paste(LD_DECAY_CELLS, collapse = ", "))

.read_decay <- function(tag, cell, rep, env) {
  f <- file.path(BUNDLE_DIR, sprintf("bundle_%s_rep%d_%s_env%d.rds", tag, rep, cell, env))
  if (!file.exists(f)) return(data.table(a_pred = NA_real_, b = NA_real_))
  b <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(b) || is.null(b$LD_decay$decay_sum)) return(data.table(a_pred = NA_real_, b = NA_real_))
  ds <- b$LD_decay$decay_sum
  ## a_pred/b are already identical across a combo's chromosome rows (a_pred
  ## is compute_LD_decay()'s own cross-chromosome-robust rate; b is a single
  ## genome-wide background-LD constant) -- confirmed directly, not assumed.
  data.table(a_pred = ds$a_pred[1], b = ds$b[1])
}
decay_rows <- combos[, .read_decay(tag, cell, rep, env), by = .(tag, cell, rep, env)]
n_missing <- sum(is.na(decay_rows$a_pred))
say("    %d/%d combos missing a bundle or LD_decay (excluded, not imputed)\n", n_missing, nrow(decay_rows))

## ---- 4. join everything at (tag, cell, rep, env) grain ----------------------
dt <- merge(null_dt, fst_dt, by = c("tag", "cell", "rep", "env"), all.x = TRUE)
dt <- merge(dt, decay_rows, by = c("tag", "cell", "rep", "env"), all.x = TRUE)
setnames(dt, c("a_pred", "b"), c("ld_decay_rate", "background_ld"))
fwrite(dt, file.path(MODULE_ROOT, "results", "exploratory_ld_fst_null_burden_combos.tsv"), sep = "\t")
say("[4] joined table: %d rows -> results/exploratory_ld_fst_null_burden_combos.tsv\n\n", nrow(dt))

## ---- 5. Spearman correlations, pooled and per cell ---------------------------
## outcome, then predictor, kept separate by method x scheme (never pooling
## mvn -- a negative control -- with group/spatial, and never averaging
## consensus/simes together; see header).
OUTCOMES   <- c("mean_null", "ratio_null_obs", "FDP_truth")
PREDICTORS <- list(ld_decay_rate = "ld_decay_rate", background_ld = "background_ld", Fst = "Fst")

.cor_row <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  n_valid <- sum(ok)
  if (n_valid < 5) return(data.table(rho = NA_real_, p_value = NA_real_, n_valid = n_valid))
  ct <- tryCatch(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE),
                error = function(e) NULL)
  if (is.null(ct)) return(data.table(rho = NA_real_, p_value = NA_real_, n_valid = n_valid))
  data.table(rho = unname(ct$estimate), p_value = ct$p.value, n_valid = n_valid)
}

results_list <- list()
for (pred_name in names(PREDICTORS)) {
  pred_col <- PREDICTORS[[pred_name]]
  cells_here <- if (pred_col %in% c("ld_decay_rate", "background_ld")) LD_DECAY_CELLS else sort(unique(dt$cell))
  scopes <- c("ALL", cells_here)
  for (outcome in OUTCOMES) {
    for (m in unique(dt$method)) {
      for (s in unique(dt$scheme)) {
        for (scope in scopes) {
          sub <- dt[method == m & scheme == s]
          if (scope != "ALL") sub <- sub[cell == scope]
          if (!nrow(sub)) next
          row <- .cor_row(sub[[pred_col]], sub[[outcome]])
          row[, `:=`(predictor = pred_name, outcome = outcome, method = m, scheme = s, scope = scope,
                     n_rows = nrow(sub), frac_outcome_na = mean(!is.finite(sub[[outcome]])))]
          results_list[[length(results_list) + 1]] <- row
        }
      }
    }
  }
}
res <- rbindlist(results_list)
setcolorder(res, c("predictor", "outcome", "method", "scheme", "scope", "n_rows", "n_valid",
                   "frac_outcome_na", "rho", "p_value"))
fwrite(res, file.path(MODULE_ROOT, "results", "exploratory_ld_fst_null_burden_correlations.tsv"), sep = "\t")
say("[5] wrote %d correlation rows -> results/exploratory_ld_fst_null_burden_correlations.tsv\n\n", nrow(res))

## ---- 6. console highlight: pooled ("ALL"), primary outcome, candidate nulls --
say("--- pooled (scope = ALL), outcome = mean_null, candidate nulls only ---\n")
print(res[scope == "ALL" & outcome == "mean_null" & scheme != "mvn",
          .(predictor, outcome, method, scheme, n_rows, n_valid, rho, p_value)][order(predictor, method, scheme)])

say("\n--- pooled (scope = ALL), outcome = mean_null, mvn NEGATIVE CONTROL ---\n")
print(res[scope == "ALL" & outcome == "mean_null" & scheme == "mvn",
          .(predictor, outcome, method, scheme, n_rows, n_valid, rho, p_value)][order(predictor, method)])

say("\n--- per cell, outcome = mean_null, scheme = group, method = emmax_consensus ---\n")
print(res[scope != "ALL" & outcome == "mean_null" & scheme == "group" & method == "emmax_consensus",
          .(predictor, scope, n_rows, n_valid, rho, p_value)][order(predictor, scope)])

## ---- 7. two quick diagnostic plots (candidate nulls, mean_null, per cell) ---
plot_dt <- melt(dt[scheme != "mvn"],
                measure.vars = c("ld_decay_rate", "background_ld", "Fst"),
                variable.name = "predictor", value.name = "predictor_value")
plot_dt <- plot_dt[is.finite(predictor_value) & is.finite(mean_null)]

p_ld <- ggplot(plot_dt[predictor %in% c("ld_decay_rate", "background_ld")],
               aes(predictor_value, mean_null, colour = scheme)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
  facet_grid(method ~ cell + predictor, scales = "free_x") +
  theme_bw(11) + theme(strip.background = element_blank()) +
  labs(x = NULL, y = "null-produced regions per surrogate (mean_null)",
       title = "Exploratory: LD-decay rate / background LD vs null-region burden, by cell")
ggsave(file.path(MODULE_ROOT, "figures", "exploratory_ld_decay_vs_null_burden.pdf"),
      p_ld, width = 14, height = 6)

p_fst <- ggplot(plot_dt[predictor == "Fst"],
                aes(predictor_value, mean_null, colour = scheme)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
  facet_grid(method ~ cell, scales = "free_x") +
  theme_bw(11) + theme(strip.background = element_blank()) +
  labs(x = "mean pairwise Fst", y = "null-produced regions per surrogate (mean_null)",
       title = "Exploratory: mean pairwise Fst vs null-region burden, by cell")
ggsave(file.path(MODULE_ROOT, "figures", "exploratory_fst_vs_null_burden.pdf"),
      p_fst, width = 14, height = 5)
say("[7] wrote figures/exploratory_ld_decay_vs_null_burden.pdf, figures/exploratory_fst_vs_null_burden.pdf\n")

cat("EXPLORATION_DONE\n")
