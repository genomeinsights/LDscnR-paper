## =============================================================================
## final_analysis/R/21_floor_stability_select_floors.R
##
## EXPLORATORY, additive analysis (PK, 2026-09-27) -- does NOT touch
## R/12_floor_decomposition.R, its outputs, or any canonical/manuscript
## number. New results/figures only, under distinct `floor_stability_*`
## filenames.
##
## Question this whole analysis (R/21-R/23) answers: in the 600 primary c=1
## simulations, does requiring a reported region to persist across floors
## increase precision, and how much recall is lost? Hypothesis, not an
## expected result to confirm.
##
## Step 1 (this script): PHENOTYPE-BLIND floor selection, per combo -- pick
## the integer Stage-1 size floor closest to a 99.5%/99.7%/99.9% reduction in
## test count relative to ALL ASSAYED MARKERS, using cluster sizes alone,
## before any association result or truth label is examined. See
## R/helpers_floor_stability.R's select_stability_floors() for the exact
## selection rule (tie -> smaller floor) and its own worked-example smoke
## test (not part of this script's own run).
##
## Scope: the 600 c=1 combos only (V0.5_c1, V1_c1, V2_c1; both tags; all 10
## reps x 10 envs) -- CELLS_ALL's own 3 "_c1" cells, the full grid for them
## (not a balanced/subsampled subset).
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(parallel)})
MODULE_ROOT <- path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis")
source(file.path(MODULE_ROOT, "R", "00_config.R"))
source(file.path(MODULE_ROOT, "R", "helpers_floor_stability.R"))
say("=== 21_floor_stability_select_floors ===\n\n")

C1_CELLS <- c("V0.5_c1", "V1_c1", "V2_c1")
TARGETS <- c(0.995, 0.997, 0.999)
CANONICAL_TARGET <- 0.997

select_floors_one_combo <- function(tag, cell, rep, env) {
  combo_id <- sprintf("%s_%s_rep%d_env%d", tag, cell, rep, env)
  bf <- file.path(stage_dir("02_build_ld_units", combo_id), "ld_units.rds")
  if (!file.exists(bf)) return(list(ok = FALSE, data = NULL, error = sprintf("missing ld_units.rds: %s", combo_id)))
  b <- readRDS(bf)
  units1 <- LDscnR:::.ld_outlier_units(b$stage1, b$map, 1L)
  n_total <- nrow(b$map)
  sel <- if (FLOOR_SCHEME == "pct") select_stability_floors(units1, n_total, targets = TARGETS)
         else select_density_floors(units1, b$map)
  sel[, `:=`(tag = tag, cell = cell, rep = rep, env = env, combo_id = combo_id, n_total_markers = n_total)]
  list(ok = TRUE, data = sel, error = NA_character_)
}

if (sys.nframe() == 0L) {
  combos <- CJ(tag = TAGS_ALL, cell = C1_CELLS, rep = REPS_ALL, env = ENVS_ALL)
  say("[1] phenotype-blind floor selection for %d c=1 combos x %d targets\n", nrow(combos), length(TARGETS))
  t0 <- Sys.time()
  res <- mclapply(seq_len(nrow(combos)), function(i) {
    tryCatch(select_floors_one_combo(combos$tag[i], combos$cell[i], combos$rep[i], combos$env[i]),
             error = function(e) list(ok = FALSE, data = NULL,
                                      error = sprintf("%s_%s_rep%d_env%d: %s", combos$tag[i], combos$cell[i],
                                                      combos$rep[i], combos$env[i], conditionMessage(e))))
  }, mc.cores = 12)
  ok_flags <- vapply(res, `[[`, logical(1), "ok")
  if (!all(ok_flags)) {
    errs <- vapply(res[!ok_flags], `[[`, character(1), "error")
    message(paste(errs, collapse = "\n"))
    stop(sprintf("%d/%d combos failed floor selection -- see messages above; refusing to pool a partial result set",
                 sum(!ok_flags), nrow(combos)))
  }
  sel <- rbindlist(lapply(res, `[[`, "data"))
  say("[2] %d rows (%d combos x %d targets) -- 0 failures\n", nrow(sel), nrow(combos), length(TARGETS))

  fwrite(sel, file.path(MODULE_ROOT, "results", fs_name("floor_selection.tsv")), sep = "\t")
  say("[3] wrote results/floor_stability_floor_selection.tsv\n")

  ## ---- collapse diagnostics: per combo, do 99.5/99.7/99.9 pick the SAME ----
  ## integer floor? A collapsed pair means that "recovery" at the collapsed
  ## side is trivial (same region set by construction), not an independent
  ## confirmation -- R/23 reads this table to gate its tier logic, never
  ## re-derives collapse status itself.
  wide <- dcast(sel, tag + cell + rep + env + combo_id + n_total_markers ~ target, value.var = "floor")
  setnames(wide, c("0.995", "0.997", "0.999"), c("floor_995", "floor_997", "floor_999"))
  wide[, `:=`(collapsed_995_997 = floor_995 == floor_997,
             collapsed_997_999 = floor_997 == floor_999,
             collapsed_995_999 = floor_995 == floor_999)]
  wide[, n_distinct_floors := vapply(seq_len(.N), function(i) uniqueN(c(floor_995[i], floor_997[i], floor_999[i])), integer(1))]
  fwrite(wide, file.path(MODULE_ROOT, "results", fs_name("floor_selection_wide.tsv")), sep = "\t")
  say("[4] wrote results/floor_stability_floor_selection_wide.tsv\n")

  say("\n[5] collapse frequency (of %d combos):\n", nrow(wide))
  say("    99.5%% == 99.7%%: %d (%.1f%%)\n", sum(wide$collapsed_995_997), 100 * mean(wide$collapsed_995_997))
  say("    99.7%% == 99.9%%: %d (%.1f%%)\n", sum(wide$collapsed_997_999), 100 * mean(wide$collapsed_997_999))
  say("    99.5%% == 99.9%%: %d (%.1f%%)\n", sum(wide$collapsed_995_999), 100 * mean(wide$collapsed_995_999))
  say("    n_distinct_floors distribution: %s\n", paste(capture.output(print(table(wide$n_distinct_floors))), collapse = " "))

  say("\n[6] achieved-reduction summary by target (should sit close to the nominal target if floor granularity allows):\n")
  print(sel[, .(mean_achieved = mean(achieved_reduction), min_achieved = min(achieved_reduction),
               max_achieved = max(achieved_reduction), mean_n_eligible_units = mean(n_eligible_units)), by = target])

  write_receipt(fs_stage("21_floor_stability_select_floors"), inputs = character(),
                params = list(targets = TARGETS, canonical_target = CANONICAL_TARGET, cells = C1_CELLS),
                outputs = c(paste0("results/", fs_name("floor_selection.tsv")), paste0("results/", fs_name("floor_selection_wide.tsv"))))
  cat("FLOOR_SELECTION_DONE\n")
}
