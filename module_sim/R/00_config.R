## =============================================================================
## module_sim/R/00_config.R
##
## EVERY PARAMETER OF THE SIMULATION PIPELINE, SET HERE AND NOWHERE ELSE.
##
## Same contract as module_3sp/R/00_config.R, same reason: no stage script may
## contain a numeric constant that affects a result. This file was cleared and
## rebuilt from scratch on 2026-09-04 (PK) -- module_sim previously held the
## superseded C-score / single-tau RMSC chain, archived in git history rather
## than carried forward. Nothing here inherits that pipeline's parameters;
## everything is re-derived or re-justified.
##
## Sourced by every stage. Defines parameters, paths, and the receipt machinery
## that makes "what needs rerunning" a computed answer instead of a judgement
## call -- identical machinery to module_3sp, copied rather than re-invented,
## because it is infrastructure, not a modelling choice.
## =============================================================================
suppressMessages({library(data.table); library(digest)})

## ---- 1. WHERE THINGS ARE ------------------------------------------------------
PATHS <- list(
  module = path.expand("~/gitlab/LDscnR-paper/module_sim"),
  out    = path.expand("~/gitlab/LDscnR-paper/module_sim/out"),
  ## RAW SIMULATION OUTPUT -- never written by this pipeline, and NOT copied in.
  ##
  ## [!] UNLIKE module_3sp, this is a DELIBERATE departure from "copy raw inputs
  ## into the module," and the reason is size, not principle: the full grid
  ## across (V, c) cells, environments and both BGS arms is far past anything
  ## that can be duplicated onto this disk or into git. The external-volume
  ## dependency this creates is exactly what module_3sp's own header names as a
  ## defect in the pipeline it replaced ("pointed this at /Volumes/Nemo... one
  ## reason its bundle could not be rebuilt on demand"). It is accepted here
  ## because there is no alternative, not because it stopped mattering -- 01
  ## verifies presence and hashes what it reads, same discipline as module_3sp,
  ## so the dependency is at least VISIBLE and CHECKED rather than assumed.
  raw_root = "/Volumes/Nemo/Nemo_sim",
  ## Per-(V,c,env,chr) bundles: raw genotypes, map (architecture + ground truth),
  ## phenotype. TAGS below selects which arm.
  raw_nobgs = "/Volumes/Nemo/Nemo_sim/regen_sim_data_nobgs",
  raw_bgs5  = "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5",
  cache = file.path(path.expand("~/gitlab/LDscnR-paper/module_sim"), "cache")
)
PATHS$el_dir <- file.path(PATHS$cache, "edge_lists")

## ---- 2. WHICH SLICE OF THE GRID THIS BUILD TARGETS -----------------------------
## [!] OPEN DECISION, FLAGGED RATHER THAN GUESSED. The full grid is 4 selection
## cells (V0.5_c1, V0.5_c2, V1_c1.5, V2_c1) x 2 BGS arms x 10 environments x 10
## chromosomes -- module_3sp has exactly one dataset; this module does not.
## Building the whole grid through every stage before anything is verified would
## repeat today's module_3sp lesson (rebuild small, verify, then scale) in the
## worst possible way -- hours of compute before a single number can be checked.
##
## Set to the SMALLEST cell that lets every stage be verified end to end. Widen
## CELLS/TAGS/ENVS once 02-04 are confirmed correct on this one, not before.
CELLS <- "V2_c1"          # widen to c("V0.5_c1","V0.5_c2","V1_c1.5","V2_c1") once verified
TAGS  <- "nobgs"          # widen to c("nobgs","bgs5") once verified -- see RAW_SCORING below
ENVS  <- 1L                # widen once verified
CHRS  <- 1:10               # all 10 chromosome files pool into one per-(cell,tag,env) unit

## ---- 3. SEEDS -------------------------------------------------------------------
SEEDS <- c(bundle = 1L, clusters = 11L, nulls = 41L, sensitivity = 41L)

## ---- 4. RAW SCORING: REBUILD, NOT INHERIT --------------------------------------
## [!] THE SINGLE BIGGEST ASSUMPTION IN THIS FILE, stated plainly so PK can
## override it in one line.
##
## Every bundle under regen_sim_data_* already carries emx_p/emx_F, lfmm_p/lfmm_F,
## a GRM, grm_markers, and a complexity_reduction object -- method OUTPUTS baked
## in by the superseded regen_sim_data.R pipeline at parse time. Checked directly
## (2026-09-04): that GRM's grm_method is "complexity_chain" over 13,963 of
## 30,922 markers, NOT the stage-1-representatives basis module_3sp settled on,
## and complexity_reduction$params$rho is EMPTY -- unseeded, unrecorded, exactly
## the "cannot be reproduced even in principle" problem module_3sp's own header
## describes for its superseded bundle.
##
## RAW_SCORING = "rebuild": treat the bundle's GTs, map's ARCHITECTURE/TRUTH
## columns (Chr, Pos, marker, type, chr_type, true_QTN, focal_QTN,
## max_LD_with_QTN, MAF, Va) and env as the raw substrate -- the sim-side
## equivalent of module_3sp's raw_3sp.RData -- and recompute decay, ld_w,
## stage-1 clustering, the kinship basis, and both engine scans fresh, in this
## module's own stages, under this file's own seeds. Discard emx_p, emx_F,
## lfmm_p (see below), GRM, grm_markers, complexity_reduction from the bundle.
## The alternative, "inherit", would carry the same unreproducibility into every
## downstream number that module_3sp was rebuilt specifically to remove.
RAW_SCORING <- "rebuild"

## LFMM. [!] NOT YET CONFIRMED reproducible from raw genotypes within this
## codebase -- unlike EMMAX, no function here computes it, matching module_3sp's
## LFMM_SOURCE situation exactly. UNLIKE module_3sp, this data is entirely ours
## (not an external panel from another lab), so "was the bundle's lfmm_p/lfmm_F
## computed reproducibly by regen_sim_data.R itself" is an answerable question,
## not an inherited unknown -- it has not been answered yet. Defaulting to
## "inherit" (module_3sp's honest default) until that is checked; do not read
## this as a claim that it IS reproducible.
LFMM_SOURCE <- "inherit"

## ---- 5. STAGE 02: THE BUNDLE (decay, ld_w, stage-1 clustering, kinship) --------
## n_win_decay = 20 is canonical independently on BOTH halves of this project --
## established here from window-size testing early in the sims work, and
## independently the value module_3sp pins for the same reason. Not a
## coincidence worth re-litigating; carried forward as agreed infrastructure.
DECAY_ARGS <- list(
  min_maf_decay  = 0.1,
  q              = 0.95,
  n_sub_bg       = 5000,
  n_win_decay    = 20,
  overlap        = 0.5,
  max_SNPs_decay = Inf,
  prob_robust    = 0.95,
  max_pairs      = 5000,
  ld_method      = "corr",
  n_strata       = 20,
  keep_el        = TRUE,
  slide          = 500,
  rho_targets    = 0.99,
  cores          = 1
)
RHO_GRID <- c(seq(0.05, 0.95, by = 0.05), 0.99)

MAF_KEEP <- 0.1     # same convention as module_3sp; applied to the raw sim map before anything else

## GRM basis and estimator. SAME reasoning as module_3sp, because the reasoning
## was never dataset-specific: the same operation that defines the test units
## (stage-1 clustering) selects the kinship markers, so there is one decision
## instead of two. And separately: GCTA vs a centred product on IDENTICAL
## markers is now KNOWN to move discoveries materially (module_3sp CR; confirmed
## on this simulated data too via grm_comparison.R's chain vs chain_centred arm,
## 0.843x at fixed tau, p=8e-5) -- so the estimator is reported explicitly, not
## left implicit.
GRM_BASIS  <- "stage1_pruned"
GRM_METHOD <- "GCTA"

## ---- 6. STAGE 03: CLUSTERING ----------------------------------------------------
## rho = 0.5: the standing decision for this project (not re-opened here).
## distance_threshold = 1e5: module_3sp's Table 2 decision (PK, 2026-09-04),
## carried here because module_sim's OWN git history already harmonised to this
## value once (commit 8dbb09a, "Harmonise the distance cap to 1e5 across the
## paper modules") before the module was cleared -- so this is not a new choice,
## it is restoring one already made. [!] Config header there also records why
## the simulations cannot themselves justify this number: LD blocks here are
## short relative to either 1e5 or 5e5, so the cap is nearly inert on this data
## and decisive only on the stickleback panel (Eda). Recorded here so a future
## reader does not mistake this for a value the sims independently support.
CR_RHO             <- 0.5
DISTANCE_THRESHOLD <- 1e5
SCORE_THRESHOLD    <- 0.80

## ---- 7. STAGES 04+: TESTING ------------------------------------------------------
ALPHA <- 0.05
## [!] SIZE_FLOOR IS NOT SET HERE, DELIBERATELY. module_3sp's floor of 8 is
## "2x the median stage-1 cluster size of 4.11" for ONE dataset -- a derived
## quantity, not a constant, and this module's median cluster size will differ
## by (V, c) cell (dispersal kernel c drives background LD, which drives
## cluster size -- established earlier this project). Stage 02 must COMPUTE and
## RECORD the median per cell it builds; stage 04 reads that recorded value
## rather than a number typed in here. A single SIZE_FLOOR hardcoded across
## cells with different LD structure would silently re-introduce the coupling
## bug this file's whole contract exists to prevent.

STATISTICS <- c("consensus", "Simes")
ENGINES    <- c("EMMAX", "LFMM")
UNIT_REPR      <- "consensus_dosage"
UNIT_REPR_ALL  <- c("consensus_dosage", "eMLG", "representative", "best_snp")
EMLG_ARGS      <- list(input = "auto", cor_th = 0.8, l_min = 10)
BEST_SNP_ARGS  <- list(fill = TRUE, round_fill = TRUE)

## ---- REGION ASSEMBLY (post hoc; cannot affect a p-value) -------------------------
## Same construction as module_3sp's REGION_ASSEMBLY, because ld_outlier_test()
## is the same package function and its stage2_discovered branch hardcodes
## ld_w_threshold=0 and merge-over-significant-clusters-only internally --
## verified by reading R/ld_outlier_test.R directly (2026-09-04) rather than
## assumed from module_3sp's comment. min_r2_rho is likewise pulled from
## stage1$params$rho at call time, not passed. Kept as a list here purely as
## documentation of fixed package behaviour, matching module_3sp's own caveat:
## do not expect changing these fields to change anything.
REGION_ASSEMBLY <- list(
  ld_w_threshold     = 0,
  min_r2_rho         = NA_real_,             # resolved from stage1$params$rho at call time
  score_threshold    = SCORE_THRESHOLD,
  distance_threshold = DISTANCE_THRESHOLD,
  merge_over         = "discovered"
)

## ---- 8. GROUND TRUTH -------------------------------------------------------------
## [!] REPLACES module_3sp's EcoPeak rotation null -- there is no external
## validation set here, there is KNOWN truth. Recall/precision against
## true_QTN/max_LD_with_QTN is the sim-side analogue of the rotation p-value,
## and the neutral-region convention (chr_type == "ntrl": the marker-level flag
## for "not within the QTN-affected span of its chromosome," NOT a whole neutral
## chromosome -- checked directly against real bundles 2026-09-04, correcting an
## earlier looser description of this convention) is the convention-free
## false-positive control this project already validated once
## (neutral_chr_control.R, sec_sims claim: cluster-level BH at nominal 5% running
## at roughly 70% FDP). Both should be reported for every scan, not one or the
## other.
TRUTH_DIST_MAX <- NULL   # derived from score_thresholds(decay_sum, ...) at call time, per module_3sp precedent

## ---- 9. NULLS ----------------------------------------------------------------------
## [!] NOT YET DESIGNED. module_3sp's rotation null has no sim-side counterpart
## (see section 8); a permutation null here needs a structure-aware surrogate,
## and this project's own earlier finding is that the obvious one -- MVN(0,
## sigma_g^2 K + sigma_e^2 I) -- is anti-conservative by 5.9x in the tail where
## BH operates, and that orthogonalising against the covariate only rescues some
## bases (genetic: 90x -> 1x chance; spatial: locked at 224x either way). Do not
## default to the naive MVN null for this pipeline's reported number without
## re-stating that finding next to it.
NPERM <- NULL   # TBD once the null basis is chosen

## ---- 10. SENSITIVITY -----------------------------------------------------------------
## Same axes module_3sp sweeps, because they are dataset-agnostic method
## questions the panel side is already testing -- kept parallel deliberately, so
## a finding on one half has a natural counterpart to check on the other
## (module_3sp CR/CS/CY-CZ this week all took exactly this shape).
SWEEP <- list(
  rho           = c(0.35, 0.50, 0.65),
  grm_basis     = c("stage1_pruned", "greedy", "none"),
  grm_estimator = c("GCTA", "centred"),
  unit_repr     = UNIT_REPR_ALL
)

## =============================================================================
## LDSCNR VERSION PIN -- identical machinery to module_3sp, same package, same
## branch, reused rather than duplicated with different logic that could drift.
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
## RECEIPT MACHINERY -- identical to module_3sp, reused verbatim.
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
