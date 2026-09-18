## =============================================================================
## final_analysis/R/rerun_score_truth_grid.R
##
## Full-1,400-combo re-run of score_truth() (05_score_truth.R) after the
## Stage-2 ordering-unification fix (PK, 2026-09-18 review, point 1;
## region_scoring_version 3 -> 4 in 05_score_truth.R's PARAMS). Stages 01-04
## are untouched by this fix, so run_combo() is deliberately NOT used here --
## calling score_truth() directly lets its own stage_stale() decide (it will
## correctly detect the PARAMS change via region_scoring_version and
## recompute), without re-touching parse/build/EMMAX/LFMM at all.
##
## [!] run_full_grid.sh is NOT used for this rerun: its shell-level
## "ALREADY_DONE" check tests only truth_scores.rds's EXISTENCE, bypassing
## score_truth()'s own stage_stale()/PARAMS-based invalidation entirely --
## it would skip all 1,400 combos since they already exist from the prior
## (version=3) run.
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(SNPRelate); library(LEA); library(parallel)})
MOD <- path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis")
source(file.path(MOD, "R", "00_config.R"))
source(file.path(MOD, "R", "01_parse_nemo.R"))
source(file.path(MOD, "R", "02_build_ld_units.R"))
source(file.path(MOD, "R", "03_emmax.R"))
source(file.path(MOD, "R", "04_lfmm.R"))
source(file.path(MOD, "R", "05_score_truth.R"))

combos <- CJ(tag = TAGS_ALL, cell = CELLS_ALL, rep = REPS_ALL, env = ENVS_ALL)
say("[1] rescoring %d combos (region_scoring_version bump only -- stages 01-04 short-circuit via their own stage_stale())\n", nrow(combos))
t0 <- Sys.time()
res <- mclapply(seq_len(nrow(combos)), function(i) {
  tryCatch({
    score_truth(combos$tag[i], combos$cell[i], combos$rep[i], combos$env[i])
    "OK"
  }, error = function(e) {
    message(sprintf("[rerun_score_truth_grid] %s_%s_rep%d_env%d: %s",
                    combos$tag[i], combos$cell[i], combos$rep[i], combos$env[i], conditionMessage(e)))
    "FAILED"
  })
}, mc.cores = 7)
status <- unlist(res)
say("[2] done in %.1f min: %s\n", as.numeric(difftime(Sys.time(), t0, units = "mins")),
    paste(sprintf("%s=%d", names(table(status)), table(status)), collapse = ", "))
if (any(status == "FAILED")) {
  say("    %d combo(s) FAILED -- see messages above.\n", sum(status == "FAILED"))
}
cat("RESCORE_GRID_DONE\n")
