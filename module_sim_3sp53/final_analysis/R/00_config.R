## =============================================================================
## final_analysis/R/00_config.R
##
## Every result-affecting value for the clean pipeline, centralised here per
## CLAUDE_REANALYSIS_INSTRUCTIONS.md's Phase 2 requirement -- Phase 1 (this
## file, as it stands) only needs the parser's own paths and design constants;
## later phases add to this file rather than duplicating a value in a stage
## script.
##
## Host-dependent paths are NOT assumed to exist -- see .check_paths() at the
## bottom, called once, which fails with a clear message naming the missing
## path rather than a downstream file-not-found deep in a stage script.
## =============================================================================
suppressMessages({library(data.table); library(digest)})

MODULE_ROOT <- path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis")

## ---- 1. WHERE THINGS ARE (raw = read-only; everything else = NEW roots, ----
## none overlapping module_sim_3sp53's own out/ or parsed/ so no old output
## can ever satisfy a receipt or get silently reused) --------------------------
PATHS <- list(
  ## raw NEMO output -- READ-ONLY. Same source module_sim_3sp53 uses (both
  ## tags live in the one directory); see ../CLAUDE_REANALYSIS_INSTRUCTIONS.md
  ## and ~/gitlab/LDscnR-NEMO/make_prod.sh for provenance (1400 adapt runs +
  ## 140 shared burn-ins, verified complete, 0 failures).
  raw_nemo       = "/Volumes/Large_storage/prod_out",
  raw_recmap_dir = path.expand("~/LDscnR-NEMO/params_3spC_53cM/rds"),
  raw_env_dir    = path.expand("~/LDscnR-NEMO/params_3spC_53cM"),

  ## Parsed bundles -- NEW root, "_final_v1" suffix, same volume as the old
  ## (large) root but never the same directory.
  parsed = "/Volumes/Large_storage/module_sim_3sp53_parsed_final_v1",

  ## Everything else this subtree produces (QA reports now; stage outputs,
  ## once later phases exist) -- kept inside the git-tracked final_analysis/
  ## tree rather than the external volume, since these are small and meant
  ## to be reviewed directly.
  qc  = file.path(MODULE_ROOT, "qc"),
  out = file.path(MODULE_ROOT, "out_final_v1"),

  module = MODULE_ROOT
)

## ---- 2. THE FULL DESIGN GRID -------------------------------------------------
## 7 cells x 2 tags x 10 reps (map IDs -- NEMO's raw run-directory naming
## calls this axis "chr<N>", which is NOT the 2 genomic chromosomes; see
## make_prod.sh's "10 chromosomes" == 10 map IDs) x 10 envs (environmental
## continuations sharing their rep's single burn-in -- NOT independent
## replicates; see CLAUDE_REANALYSIS_INSTRUCTIONS.md's "Pooling and
## uncertainty" section) = 1,400 combinations.
CELLS_ALL <- c("V0.5_c1", "V0.5_c1.5", "V0.5_c2", "V1_c1", "V1_c1.5", "V2_c1", "V2_c1.5")
TAGS_ALL  <- c("nobgs", "bgs")
REPS_ALL  <- 1:10
ENVS_ALL  <- 1:10

## ---- 3. PARSER-LEVEL CONSTANTS -----------------------------------------------
## NEMO's own design (make_prod.sh: BGS_MAP_RES=3.3e-6, quanti_loci=101,
## delet_loci=2000). Sample size: 320 raw individuals per (tag,cell,rep,env)
## before subsampling, 160 after step=2 (see Stop points -- assert, don't
## assume).
SUBSAMPLE_STEP <- 2L
MAF_KEEP       <- 0.10   ## manuscript default (materials_and_methods.tex); Stage-1/2 use, not the parser itself, but centralised here per the "no second hard-coded copy" rule
GRID_SIDE      <- 48L    ## environmental surface is a 48x48 spatial grid

## ---- 4. PROVENANCE (captured once per run, saved with every QA report) ------
ldscnr_sha <- function() {
  tryCatch(system2("git", c("-C", path.expand("~/gitlab/LDscnR"), "rev-parse", "--short", "HEAD"),
                   stdout = TRUE, stderr = FALSE), error = function(e) NA_character_)
}
r_version <- function() paste(R.version$major, R.version$minor, sep = ".")

## ---- 5. FAIL LOUDLY ON A MISSING HOST PATH -----------------------------------
.check_paths <- function() {
  need <- c(PATHS$raw_nemo, PATHS$raw_recmap_dir, PATHS$raw_env_dir)
  missing <- need[!dir.exists(need)]
  if (length(missing)) {
    stop("final_analysis/R/00_config.R: required path(s) not found on this host:\n  ",
         paste(missing, collapse = "\n  "),
         "\nThis pipeline expects to run on the machine holding the raw NEMO output ",
         "(currently: Petri's mini). Check PATHS in 00_config.R if that has changed.",
         call. = FALSE)
  }
  for (p in c(PATHS$parsed, PATHS$qc, PATHS$out)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
}
.check_paths()

say <- function(fmt, ...) cat(sprintf(fmt, ...))
