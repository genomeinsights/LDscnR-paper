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

## ---- 3b. STAGE 1 (LD complexity reduction) + KINSHIP --------------------------
## Manuscript defaults (CLAUDE_REANALYSIS_INSTRUCTIONS.md, "Primary EMMAX
## comparison"): rho=0.50, MAF>0.10 (already MAF_KEEP above), minimum unit
## size 2. Same LD-decay/GRM settings module_sim_3sp53's own R/00_config.R
## used (n_win_decay=20 canonical on both halves of this project; GRM basis =
## stage-1 pruned representatives, GCTA estimator -- confirmed there to move
## discoveries materially vs a centred-product estimator on identical
## markers, 0.843x at fixed tau, p=8e-5).
DECAY_ARGS <- list(
  min_maf_decay = 0.1, q = 0.95, n_sub_bg = 5000, n_win_decay = 20, overlap = 0.5,
  max_SNPs_decay = Inf, prob_robust = 0.95, max_pairs = 5000, ld_method = "corr",
  n_strata = 20, keep_el = FALSE, slide = 500, rho_targets = 0.99, cores = 1
)
RHO_GRID   <- c(seq(0.05, 0.95, by = 0.05), 0.99)
CR_RHO     <- 0.50    ## Stage-1 rho -- manuscript default, per instructions
SIZE_FLOOR <- 2L       ## minimum Stage-1 unit size -- manuscript default, per instructions
GRM_BASIS  <- "stage1_pruned"
GRM_METHOD <- "GCTA"

## ---- 3c. ASSOCIATION TESTING ---------------------------------------------------
ALPHA <- 0.05
UNIT_REPR <- "consensus_dosage"   ## emmax_consensus's per-unit variable
## LFMM_K justified 2026-09-09 (PK): "The individuals sampled from the 48x48
## grid are clustered in the corners and one population in the middle,
## that's where K=5 comes from. This was deliberate to avoid the discussion
## of what K to use if the samples were... randomly sampled across the
## landscape." I.e. K=5 is a property of the SAMPLING DESIGN (5 discrete
## spatial groups the 80 populations fall into, by construction), verifiable
## from population (x,y) coordinates ALONE -- independent of genotypes or
## phenotype, exactly the "structure diagnostic independent of association
## truth" the instructions ask for, not a value chosen by maximising QTN
## recovery. Already confirmed independently during this project's
## structured-null work (kmeans(k=5) on population coordinates recovers the
## same 4-corners-plus-centre grouping -- see module_sim_3sp53/R/12_
## structured_null.R's own "5 spatial population groups" section). Not
## re-verified again here with fresh code since it would be the identical
## check on the identical (fixed, design-level) input.
LFMM_K <- 5L
## LFMM itself remains out of scope for this checkpoint (Phase 1 + primary
## EMMAX gate only) -- this resolves the K justification, not the decision
## to build/run the LFMM stage.

## ---- 3d. REGION ASSEMBLY (post hoc; cannot affect a p-value or BH decision) ---
## Same fixed package behaviour as module_sim_3sp53/R/00_config.R documents:
## ld_outlier_test()'s stage2_discovered branch hardcodes ld_w_threshold=0 and
## min_r2_rho is resolved from stage1$params$rho at call time, not passed.
DISTANCE_THRESHOLD <- 1e5   ## module_3sp's Table 2 decision, harmonised across the paper's modules
SCORE_THRESHOLD    <- 0.80
REGION_ASSEMBLY <- list(ld_w_threshold = 0, min_r2_rho = NA_real_,
                        score_threshold = SCORE_THRESHOLD, distance_threshold = DISTANCE_THRESHOLD,
                        merge_over = "discovered")

## ---- 3e. TRUTH DEFINITION AND SCORING ------------------------------------------
## Va_j = 2 p_j (1-p_j) a_j^2, computed on the analysed (post-subsample,
## post-MAF) individuals. Primary detectability: MAF>0.10 (MAF_KEEP) AND at
## least this share of the total QTN additive variance ON THAT CHROMOSOME.
VA_SHARE_DETECTABLE <- 0.05
VA_SHARE_SENSITIVITY_GRID <- c(0.01, 0.02, 0.05, 0.10, 0.20)   ## includes the primary value

## Decay-relative truth-matching thresholds -- LDscnR::score_thresholds()'s
## OWN defaults already match the instructions' intended values exactly
## (checked directly against R/ld_benchmark.R: rho_r2=0.75, rho_d=0.95); set
## explicitly here anyway, per "no second hard-coded copy" / central config.
TRUTH_RHO_R2 <- 0.75
TRUTH_RHO_D  <- 0.95
TRUTH_DMAX_CAP <- DISTANCE_THRESHOLD   ## "the manuscript distance cap from central config"

## ---- 3f. UNCERTAINTY -----------------------------------------------------------
N_BOOTSTRAP <- 2000L   ## cluster-bootstrap replicates (map/burn-in IDs primary; crossed two-way sensitivity)

## ---- 3g. SEEDS -------------------------------------------------------------------
SEEDS <- c(bundle = 1L, clusters = 11L, bootstrap = 41L)

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

say <- function(fmt, ...) { cat(sprintf(fmt, ...)); flush(stdout()) }

## ---- 6. RECEIPT / CACHING MACHINERY -------------------------------------------
## Ported from module_sim_3sp53/R/00_config.R (proven infrastructure, not a
## modelling choice -- same reasoning that file gives for copying it from
## module_3sp rather than re-inventing it). Satisfies "Cache expensive
## phenotype-blind objects and association statistics separately... rescorable
## without recomputing LD, the GRM, EMMAX or LFMM": each stage gets its own
## receipt keyed by combo_id, fingerprinting inputs (by content hash) and
## params (by identity) so a later stage can be rerun (e.g. after a truth-
## threshold change) without invalidating an earlier one's cache.
##
## NOT ported: the old check_ldscnr() hard pin against a specific LDscnR
## commit/source-hash -- this session's own work already moved that source
## (ld_region_stability.R fix, outlier-scan commit c29023f), so a hard pin
## copied from the old config would immediately and spuriously fail. Records
## the SHA for provenance instead of gating on it; re-introduce a pin once the
## final_analysis pipeline's own LDscnR dependency stabilises.
stage_dir <- function(stage, target = "") file.path(PATHS$out, stage, target)
receipt_path <- function(stage, target = "") file.path(stage_dir(stage, target), "_receipt.rds")
sha <- function(f) if (file.exists(f)) digest::digest(f, algo = "sha256", file = TRUE) else NA_character_

write_receipt <- function(stage, inputs = character(), params = list(), outputs = character(), target = "") {
  dir.create(stage_dir(stage, target), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(stage = stage, when = Sys.time(), git = git_sha(), ldscnr_sha = ldscnr_sha(),
               r_version = r_version(),
               inputs = data.table(path = inputs, sha256 = vapply(inputs, sha, "")),
               params = params, outputs = outputs), receipt_path(stage, target))
  invisible(TRUE)
}

stage_stale <- function(stage, inputs = character(), params = list(), target = "") {
  rp <- receipt_path(stage, target)
  if (!file.exists(rp)) { message("  [", stage, "] no receipt -- will run"); return(TRUE) }
  r <- readRDS(rp)
  if (!identical(params, r$params)) { message("  [", stage, "] parameters changed -- will run"); return(TRUE) }
  now <- vapply(inputs, sha, "")
  old <- stats::setNames(r$inputs$sha256, r$inputs$path)
  ch <- names(now)[is.na(old[names(now)]) | old[names(now)] != now]
  if (length(ch)) { message("  [", stage, "] inputs changed: ", paste(basename(ch), collapse = ", "),
                            " -- will run"); return(TRUE) }
  message("  [", stage, "] up to date (", format(r$when, "%Y-%m-%d %H:%M"), ")"); FALSE
}

git_sha <- function() tryCatch(system2("git", c("-C", PATHS$module, "rev-parse", "--short", "HEAD"),
                                       stdout = TRUE, stderr = FALSE), error = function(e) NA_character_)
