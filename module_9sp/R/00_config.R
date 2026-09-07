## =============================================================================
## module_9sp/R/00_config.R
##
## EVERY PARAMETER OF THE 9sp PIPELINE, SET HERE AND NOWHERE ELSE.
##
## Same contract as module_3sp/R/00_config.R, same reason: no stage script may
## contain a numeric constant that affects a result. This module mirrors
## module_3sp's structure and conventions deliberately (PK: "exactly the same
## analyses as for 3sp") -- built from scratch against the current LDscnR
## package, not a port of the legacy LD-scaling-genome-scans/empirical_data/
## R/9sp_sticlebacks.R C-score pipeline (which this supersedes for 9sp the
## same way module_3sp superseded 3sp_sticklebacks.R).
##
## 9sp = nine-spined stickleback (Pungitius pungitius) -- a DIFFERENT SPECIES
## from 3sp's three-spined stickleback (Gasterosteus aculeatus), not a bigger
## sample of the same one. Real differences from 3sp, not just scale:
##   - 149 individuals, ~1.36M markers pre-MAF-filter (3sp: 117, ~882k).
##   - 30 populations, 6 pop_locality regions, 4 lineages (Admixed/EL/WA/WL) --
##     3sp's phenotype table has none of these; population structure here is
##     both finer-grained AND stronger.
##   - The legacy 9sp script's EMMAX call includes `Covar = pheno$lineage`;
##     3sp's canonical model has no covariate beyond the GRM. NOT YET DECIDED
##     for this rebuild -- see EMMAX_COVAR below (PK, 2026-09-07: "flag it,
##     decide once we see the GRM/structure").
##   - NO EcoPeaks-equivalent external validation reference was found for
##     nine-spined stickleback anywhere in this codebase (searched
##     2026-09-07). PK: skip that piece for now rather than block on it --
##     so this config has no ECOPEAK_BEDS/recmap entries at all, not empty
##     placeholders for them.
##
## Sourced by every stage. Defines parameters, paths, and the receipt machinery
## that makes "what needs rerunning" a computed answer instead of a judgement
## call -- identical machinery to module_3sp, copied rather than re-invented,
## because it is infrastructure, not a modelling choice.
## =============================================================================
suppressMessages({library(data.table); library(digest)})

## ---- 1. WHERE THINGS ARE ------------------------------------------------------
PATHS <- list(
  module   = path.expand("~/gitlab/LDscnR-paper/module_9sp"),
  out      = path.expand("~/gitlab/LDscnR-paper/module_9sp/out"),
  figures  = path.expand("~/gitlab/LDscnR-paper/module_9sp/figures"),
  tables   = path.expand("~/gitlab/LDscnR-paper/module_9sp/tables"),
  ## RAW INPUT -- copied in from LD-scaling-genome-scans (byte-verified, see
  ## data/PROVENANCE.csv), same reasoning as module_3sp's raw_3sp: 290 MB, over
  ## GitHub's 100 MB limit either way, so self-contained on disk, not in git.
  ## Contains GT (149 x ~1.36M dosage matrix), map (Chr/Pos/marker/maf),
  ## pheno (population/lineage/locality/GPS), ecotype_bin, and a GRM whose
  ## basis is NOT YET CHARACTERISED (the legacy script overwrites it with its
  ## own LD-pruned-basis GRM before use -- see 00_config.R header note above;
  ## do not assume the bundled GRM is usable as-is).
  raw_9sp  = path.expand("~/gitlab/LDscnR-paper/module_9sp/data/9sp_data.rds"),
  ## LFMM F-values, precomputed. Same INHERITED-not-computed status as 3sp's
  ## lfmm_F.rds -- see LFMM_SOURCE below.
  lfmm_F   = path.expand("~/gitlab/LDscnR-paper/module_9sp/data/lfmm_F.rds"),
  provenance = path.expand("~/gitlab/LDscnR-paper/module_9sp/data/PROVENANCE.csv"),
  cache    = path.expand("~/gitlab/LDscnR-paper/module_9sp/cache")
)
dir.create(PATHS$figures, recursive = TRUE, showWarnings = FALSE)
dir.create(PATHS$tables, recursive = TRUE, showWarnings = FALSE)

## ---- 2. SEEDS -----------------------------------------------------------------
SEEDS <- c(bundle = 1L, clusters = 11L, nulls = 41L, sensitivity = 41L)

## ---- 3. STAGE 02: THE BUNDLE ---------------------------------------------------
## Starting values MIRROR module_3sp's canonical settings, NOT yet re-derived for
## this species/panel. Nothing here is validated against 9sp's own LD-decay
## behaviour -- that validation is what running 02_bundle.R once actually is.
## Flagged explicitly rather than presented as settled, matching 00_config.R's
## own contract: a number that "just happens" to match 3sp's is still a claim
## about 9sp until checked.
DECAY_ARGS <- list(
  min_maf_decay  = 0.1,
  q              = 0.95,
  n_sub_bg       = 5000,
  n_win_decay    = 20,        # CANONICAL in 3sp; unverified starting point here
  overlap        = 0.5,
  max_SNPs_decay = Inf,
  prob_robust    = 0.95,
  max_pairs      = 5000,
  ld_method      = "corr",
  n_strata       = 20,
  keep_el        = TRUE,
  ## [!] 3sp uses 500 ("LD decays faster so 500 is enough") -- an empirical
  ## claim about 3sp, not a package default. The legacy 9sp script itself used
  ## slide_win_ld = 1000. Left at 3sp's 500 as the starting value pending 9sp's
  ## own decay fit; revisit if the fit looks under-resolved (e.g. background LD
  ## not reached within the window).
  slide          = 500,
  rho_targets    = 0.99,
  cores          = 1
)
RHO_GRID <- c(seq(0.05, 0.95, by = 0.05), 0.99)

MAF_KEEP <- 0.1   # matches the legacy 9sp script's own filter (map_9sp$maf > 0.1)

## LFMM. Same status as module_3sp: INHERITED, not computed by this pipeline.
## lfmm_F.rds here is over 9sp's own FULL pre-MAF map (verify the exact
## pre/post-MAF alignment in 01_inputs.R/02_bundle.R the same way module_3sp's
## 04_lfmm.R does -- do not assume the same row-order logic transfers without
## checking against THIS map).
LFMM_SOURCE <- "inherit"

## ---- GRM: BASIS AND ESTIMATOR --------------------------------------------------
## STARTING POINT ONLY, inherited from 3sp's reasoning (the same operation that
## defines the test units selects the kinship markers), NOT yet independently
## verified for 9sp. 3sp's own basis/estimator sweep (stage1_pruned vs ld_w<=0.1
## vs greedy vs none; GCTA vs centred) is exactly the kind of check this panel
## has not had yet -- population structure here is stronger (30 pops, 4
## lineages) so the basis/estimator choice may matter more, not less, than it
## did for 3sp.
GRM_BASIS      <- "stage1_pruned"
GRM_BASIS_SUPP <- c("greedy", "none")
GRM_GREEDY     <- list(ld.threshold = 0.2, slide.max.bp = 5e5)
GRM_METHOD     <- "GCTA"
GRM_LDW_THRESHOLD <- 0.1
GRM_LDW_OP        <- "<"

## ---- COVARIATE: DECIDED (PK, 2026-09-07) ----------------------------------------
## SET to "lineage", confirming the legacy 9sp script's own choice, on the strength
## of a direct check: LD-pruned PCA + a lineage x ecotype crosstab (149 individuals)
## shows WL (39 Freshwater / 0 Marine) and WA (4 / 0) are ENTIRELY monomorphic for
## ecotype, EL is 90% Freshwater (54/6), Admixed is 91% Marine (4/42). 43 of 149
## individuals (WL+WA) carry NO within-lineage ecotype contrast at all -- ecotype is
## almost a deterministic function of lineage over most of the sample, not
## independently replicated across genetic backgrounds the way 3sp's regions are.
## Without a covariate the GRM (built from genome-wide relatedness, itself
## correlated with lineage) and the ecotype fixed effect are close to collinear over
## most of the sample -- residualising ecotype on lineage before testing (the same
## mechanism emmax()'s own Covar argument uses: `Y <- resid(lm(Y ~ Covar))`,
## replicated manually here since emmax_fast() -- the 25x-faster scan this module
## uses throughout, see emmax_fast.R -- has no Covar parameter) is the standard fix.
## PC1/PC2 do NOT cleanly separate lineages as dominant axes (5%/4.2% variance
## explained, wide overlapping spread) -- the confound is specifically with ecotype
## SAMPLING within lineage, not raw genetic distance, so this is not "the GRM
## already handles it via population structure" territory.
EMMAX_COVAR <- "lineage"

## ---- 4. STAGE 03: CLUSTERING ----------------------------------------------------
## Same starting values as 3sp, same "not yet re-derived" caveat as DECAY_ARGS.
CR_RHO             <- 0.5
DISTANCE_THRESHOLD <- 1e5
MIN_R2_RHO         <- 0.5
SCORE_THRESHOLD    <- 0.80
LDW_FLAG           <- 0.05

## ---- 5. STAGES 04+: TESTING ------------------------------------------------------
ALPHA <- 0.05
## SIZE_FLOOR = 10, SET (PK, 2026-09-07) -- NOT the mechanical "2x median
## stage-1 cluster size" rule 3sp's 8 came from. 02_bundle.R's first run
## reported 753,625 Stage-1 clusters, median size 1.00 (72.3% singletons), so
## that rule would give floor=2 here -- excluding only true singletons, no real
## multiplicity reduction (flagged, not applied; same situation as
## module_sim's SIZE_FLOOR=2). PK's call instead: floor=10 is the value
## consistently used across the Fang et al. papers and is the more obviously
## canonical choice on its own terms -- 3sp's 8 came from a derivation ("2x
## 4.11") that is no longer how this project sets the floor, not a value to
## propagate to a new panel just because it was 3sp's.
##
## [!] NOT SPECIFIC TO LD-COMPLEXITY REDUCTION (PK): this floor is a general
## claim that fewer than 10 supporting markers is too little evidence to trust
## a locus, and applies AT LEAST as much to single-SNP analyses as to Stage-1
## clusters -- arguably more, since a lone SNP has no within-cluster
## corroboration at all. Any single-marker arm this module reports (marker-
## wise EMMAX/LFMM BH counts, matching module_3sp's) should be read with this
## in mind, and a marker-density filter (>=10 markers in the same local
## neighbourhood) is a live candidate for the single-SNP arms too, not only for
## which Stage-1 clusters enter testing. Not yet implemented as a filter
## anywhere -- recorded here so it is not lost before stage 03/04 need it.
SIZE_FLOOR <- 10L

REGION_ASSEMBLY <- list(
  ld_w_threshold     = 0,
  min_r2_rho         = MIN_R2_RHO,
  score_threshold    = SCORE_THRESHOLD,
  distance_threshold = DISTANCE_THRESHOLD,
  merge_over         = "discovered"
)
STATISTICS <- c("consensus", "Simes")
ENGINES    <- c("EMMAX", "LFMM")

UNIT_REPR      <- "consensus_dosage"
UNIT_REPR_ALL  <- c("consensus_dosage", "eMLG", "representative", "best_snp")
EMLG_ARGS      <- list(input = "auto", cor_th = 0.8, l_min = 10)
BEST_SNP_ARGS  <- list(fill = TRUE, round_fill = TRUE)

## ---- 6. STAGE 06: NULLS -----------------------------------------------------------
## DECIDED (PK, 2026-09-07), matched to EMMAX_COVAR rather than chosen
## independently: permute ecotype WITHIN LINEAGE, at the POPULATION level (same
## granularity as 3sp's perm_regional() -- individuals within one population
## always share one permuted label, never permuted individually), analogous to
## 3sp's "regional" scheme but conditioned on lineage instead of locality.
##
## WL and WA are entirely monomorphic for ecotype (see EMMAX_COVAR note), so
## sample(ecotype) within those strata deterministically reproduces the observed
## labels on every draw -- ZERO permutation variability from 43 of 149
## individuals. This is not a bug to work around: those individuals cannot
## validate an ecotype-specific signal via permutation regardless of scheme,
## because there is no within-lineage ecotype contrast to permute. The null is
## still valid; it is just correctly uninformative about that portion of the
## sample. A sensitivity check restricting to EL + Admixed only (the two
## lineages with real within-lineage contrast: 54/6 and 4/42) is worth running
## alongside the full-sample result once stage 03 exists, to see whether the
## reported signal holds up on the "cleaner" subset.
NPERM_CONSENSUS <- 1000L
NPERM_SIMES     <- 200L
PERM_SCHEMES    <- "lineage"

## ---- 7. STAGE 08: SENSITIVITY -------------------------------------------------
SWEEP <- list(
  rho           = c(0.35, 0.50, 0.65),
  size_floor    = NULL,      # depends on SIZE_FLOOR above being set first
  grm_basis     = c("stage1_pruned", "greedy", "none"),
  grm_estimator = c("GCTA", "centred"),
  unit_repr     = UNIT_REPR_ALL
)

## ---- 8. EXTERNAL VALIDATION: SKIPPED (PK, 2026-09-07) -------------------------
## No EcoPeaks-equivalent reference found for nine-spined stickleback (searched
## this codebase 2026-09-07). No ECOPEAK_BEDS, no recmap, no rotation-null
## stage for this module until/unless a reference is identified -- deliberately
## absent rather than stubbed out, so a stage referencing them fails loudly
## instead of silently no-op'ing.

## =============================================================================
## THE LDscnR VERSION IS PART OF THE PIPELINE. Same pin as module_3sp -- both
## modules run against the same package checkout, so a mismatch here would mean
## the two panels are not comparable "the same analyses" in the first place.
## =============================================================================
LDSCNR_PIN <- list(
  repo    = path.expand("~/gitlab/LDscnR"),
  sha     = "011a165c8c99",
  branch  = "outlier-scan",
  version = "0.0.0.9000",
  src_sha = "362892c1ee4c6628",
  installed_after = NULL
)

check_ldscnr <- function(stop_on_fail = !nzchar(Sys.getenv("LDSCNR_LAX"))) {
  v <- as.character(utils::packageVersion("LDscnR"))
  d <- utils::packageDescription("LDscnR")
  built_raw <- if (is.null(d$Built)) NA_character_ else sub(";.*$", "", sub("^[^;]*;[^;]*; *", "", d$Built))
  bt <- suppressWarnings(as.POSIXct(built_raw, tz = "UTC"))
  if (length(bt) != 1L) bt <- as.POSIXct(NA_character_, tz = "UTC")
  g <- function(...) tryCatch(system2("git", c("-C", LDSCNR_PIN$repo, ...),
                                      stdout = TRUE, stderr = FALSE), error = function(e) character())
  head_sha <- substr(paste(g("rev-parse","HEAD"), collapse = ""), 1, 12)
  dirty <- length(g("status","--porcelain","--untracked-files=no")) > 0
  src <- sort(g("ls-files","R/"))
  cur <- if (!length(src)) NA_character_ else {
    fp <- file.path(LDSCNR_PIN$repo, src)
    substr(digest::digest(paste(vapply(fp[file.exists(fp)],
             function(f) paste(readLines(f, warn = FALSE), collapse = "\n"), ""), collapse = "\n"),
             algo = "sha256", serialize = FALSE), 1, 16) }
  built_ok <- TRUE
  ok <- identical(v, LDSCNR_PIN$version) && !dirty &&
        identical(cur, LDSCNR_PIN$src_sha) && built_ok
  cat(sprintf("  LDscnR %s | commit %s (%s)%s%s\n", v, head_sha, LDSCNR_PIN$branch,
              if (dirty) " [TRACKED CHANGES]" else "",
              if (identical(head_sha, LDSCNR_PIN$sha)) ""
              else sprintf(" [moved from pinned %s -- no R/ change, so the install stands]",
                           LDSCNR_PIN$sha)))
  cat(sprintf("  source hash %s %s pin %s | built %s\n", cur,
              if (identical(cur, LDSCNR_PIN$src_sha)) "==" else "!=", LDSCNR_PIN$src_sha,
              if (is.na(bt)) "unknown (loaded via devtools::load_all(), not an install)"
              else format(bt, "%Y-%m-%d %H:%M")))
  if (!ok) { m <- paste("LDscnR does not match the pin.",
      sprintf("Reinstall and update the pin:\n    R CMD INSTALL %s\n", LDSCNR_PIN$repo),
      "  Set LDSCNR_LAX=1 to proceed anyway.")
    if (stop_on_fail) stop(m) else warning(m) }
  invisible(list(version = v, sha = head_sha, branch = LDSCNR_PIN$branch,
                 built = as.character(bt), dirty = dirty, src_sha = cur, ok = ok))
}

## =============================================================================
## RECEIPT MACHINERY -- identical to module_3sp's, copied not re-invented.
## =============================================================================
stage_dir <- function(stage) file.path(PATHS$out, stage)
receipt_path <- function(stage) file.path(stage_dir(stage), "_receipt.rds")

sha <- function(f) if (file.exists(f)) digest(f, algo = "sha256", file = TRUE) else NA_character_

git_sha <- function() tryCatch(system2("git", c("-C", PATHS$module, "rev-parse", "--short", "HEAD"),
                                       stdout = TRUE, stderr = FALSE), error = function(e) NA_character_)

write_receipt <- function(stage, inputs = character(), params = list(), outputs = character()) {
  dir.create(stage_dir(stage), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(stage = stage, when = Sys.time(), git = git_sha(),
               ldscnr = tryCatch(check_ldscnr(stop_on_fail = FALSE), error = function(e) NA),
               inputs = data.table(path = inputs, sha256 = vapply(inputs, sha, "")),
               params = params, outputs = outputs), receipt_path(stage))
  invisible(TRUE)
}

stage_stale <- function(stage, inputs = character(), params = list()) {
  rp <- receipt_path(stage)
  if (!file.exists(rp)) { message("  [", stage, "] no receipt -- will run"); return(TRUE) }
  r <- readRDS(rp)
  if (!identical(params, r$params)) { message("  [", stage, "] parameters changed -- will run"); return(TRUE) }
  now <- vapply(inputs, sha, "")
  old <- setNames(r$inputs$sha256, r$inputs$path)
  ch <- names(now)[is.na(old[names(now)]) | old[names(now)] != now]
  if (length(ch)) { message("  [", stage, "] inputs changed: ", paste(basename(ch), collapse=", "),
                            " -- will run"); return(TRUE) }
  message("  [", stage, "] up to date (", format(r$when, "%Y-%m-%d %H:%M"), ")"); FALSE
}

say <- function(...) { cat(sprintf(...)); flush(stdout()) }
