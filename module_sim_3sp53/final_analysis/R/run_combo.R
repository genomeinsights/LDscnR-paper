## final_analysis/R/run_combo.R
##
## Chains parse -> Stage 1/GRM -> EMMAX -> truth-scoring for one
## (tag, cell, rep, env) combination. Used by both the gate-3 14-
## combination grid and (later) the full 1,400-combination driver.
##
## Usage: Rscript R/run_combo.R <tag> <cell> <rep> <env>
suppressMessages({library(data.table); library(LDscnR); library(SNPRelate)})
MOD <- path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis")
source(file.path(MOD, "R", "01_parse_nemo.R"))
source(file.path(MOD, "R", "02_build_ld_units.R"))
source(file.path(MOD, "R", "03_emmax.R"))
source(file.path(MOD, "R", "05_score_truth.R"))

run_combo <- function(tag, cell, rep, env, force = FALSE) {
  combo_id <- sprintf("%s_%s_rep%d_env%d", tag, cell, rep, env)
  say("\n########## %s ##########\n", combo_id)

  parsed_file <- file.path(PATHS$parsed, sprintf("nemo_%s_rep%d_%s_env%d.rds", tag, rep, cell, env))
  if (force || !file.exists(parsed_file)) {
    res <- parse_nemo_run(tag, cell, rep, env, offset = 0L)
    saveRDS(res, parsed_file)
  } else {
    say("[parse] %s exists -- skipping (pass force=TRUE to rebuild)\n", basename(parsed_file))
  }

  invisible(build_ld_units(tag, cell, rep, env, force = force))
  invisible(run_emmax(tag, cell, rep, env, force = force))
  scores <- score_truth(tag, cell, rep, env, force = force)
  scores
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 4) stop("Usage: Rscript R/run_combo.R <tag> <cell> <rep> <env>")
  invisible(run_combo(args[1], args[2], as.integer(args[3]), as.integer(args[4])))
}
