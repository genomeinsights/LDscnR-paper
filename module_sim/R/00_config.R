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
## [!] NEMO_ROOT ADDED 2026-09-04, RUNNING ON A SECOND MACHINE (PK: "4 on mini2
## and 4 here"). Everything under raw_root was hardcoded to /Volumes/Nemo/
## Nemo_sim -- a LOCAL disk on this machine (confirmed: /dev/disk5s1, apfs,
## local), not reachable from mini2 at all. Rather than mount it there (would
## need either physically moving the drive or enabling network file sharing --
## a system-settings change, not done without asking), the exact subset of
## data this grid needs (7 cells x 10 reps x env1: ~279 MB of archives + 104 MB
## of reference maps/env -- 383 MB total, verified by size before deciding this
## was practical) is copied to mini2 under ~/Nemo_data instead of /Volumes/Nemo
## (mini2's /Volumes is not writable without sudo; checked directly, not
## assumed). NEMO_ROOT makes that a config value, not a code change, on either
## machine: unset, everything resolves exactly as before.
NEMO_ROOT <- Sys.getenv("SIM_NEMO_ROOT", "/Volumes/Nemo/Nemo_sim")
PATHS <- list(
  module = path.expand("~/gitlab/LDscnR-paper/module_sim"),
  out    = path.expand("~/gitlab/LDscnR-paper/module_sim/out"),
  ## RAW SIMULATION OUTPUT -- never written by this pipeline, and NOT copied
  ## into the repository (only the subset above may be copied to a SECOND
  ## machine's disk, not into git; git remains untouched by NEMO_ROOT).
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
  raw_root = NEMO_ROOT,
  ## [!] CORRECTED 2026-09-04. Originally pointed at regen_sim_data_{nobgs,bgs5}
  ## -- traced those on PK's instruction and they are NOT raw: every bundle
  ## already carries emx_p/emx_F/lfmm_p/lfmm_F, a GRM and grm_markers built via
  ## GRM_METHOD="complexity_chain" (regen_sim_data.R:43) -- which is
  ## ld_prune_and_eMLG()'s stage-2 group representatives, the SAME set
  ## module_3sp's 02_bundle.R explicitly calls out as "a DIFFERENT set...
  ## deliberately not used" for the kinship basis. Their own upstream,
  ## parsed_sim_data2, is ALSO already scored (same columns). The true raw
  ## source is NEMO's own native per-file output, one tgz per (tag,chr,V,c,env):
  raw_nemo_nobgs = file.path(NEMO_ROOT, "Nemo_out_nobgs"),
  raw_nemo_bgs5  = file.path(NEMO_ROOT, "Nemo_out_bgs"),
  ## Reference chromosome maps (position, type, allelic_values -- the genome
  ## MODEL, shared across every simulation replicate) and the spatial
  ## environment surfaces (one file per env index, shared across V/c/chr).
  ## Neither is a simulation OUTPUT; both are simulation INPUT that every
  ## replicate reads. Verified present for chr1/env1 (2026-09-04).
  raw_recmap_dir = file.path(NEMO_ROOT, "maps_500kb_with_allelic_values", "chromosome_maps_500kb_rds"),
  raw_env_dir    = file.path(NEMO_ROOT, "env"),
  ## R_parsing/'s output: the clean, unscored bundle (GTs + map architecture/
  ## truth + env), one file per (tag,chr,V,c,env) -- the sim-side equivalent of
  ## module_3sp's raw_3sp.RData, and what 02_bundle.R in R/ will read as ITS raw
  ## input. Same external-volume reasoning as raw_root: too large for git or for
  ## a full local copy once this widens past one file.
  parsed = file.path(NEMO_ROOT, "module_sim_parsed"),
  cache  = file.path(path.expand("~/gitlab/LDscnR-paper/module_sim"), "cache")
)
PATHS$el_dir <- file.path(PATHS$cache, "edge_lists")
PATHS$untar  <- file.path(PATHS$cache, "untar")   # scratch for unpacked .tgz; not an output

## ---- 2. WHICH SLICE OF THE GRID THIS BUILD TARGETS -----------------------------
## [!] OPEN DECISION, FLAGGED RATHER THAN GUESSED. The full grid is 4 selection
## cells (V0.5_c1, V0.5_c2, V1_c1.5, V2_c1) x 2 BGS arms x 10 environments x 10
## REPLICATES -- module_3sp has exactly one dataset; this module does not.
## Building the whole grid through every stage before anything is verified would
## repeat today's module_3sp lesson (rebuild small, verify, then scale) in the
## worst possible way -- hours of compute before a single number can be checked.
##
## Set to the SMALLEST cell that lets every stage be verified end to end. Widen
## CELLS/TAGS/ENVS once 02-04 are confirmed correct on this one, not before.
CELLS <- "V2_c1"          # widen to c("V0.5_c1","V0.5_c2","V1_c1.5","V2_c1") once verified
TAGS  <- "nobgs"          # widen to c("nobgs","bgs5") once verified -- see RAW_SCORING below
ENVS  <- 1L                # widen once verified
## [!] RENAMED FROM "CHRS" 2026-09-04 (PK). "chr1".."chr10" in NEMO's own
## filenames are not ten genomic chromosomes to pool -- they are ten
## INDEPENDENT SIMULATION REPLICATES, each already containing its own two
## chromosomes (sharing one recombination map; confirmed directly, rec_map1.rds
## has rows Chr %in% c(1,2), same bp range for both, and rec_map2.rds is a
## genuinely different map -- one per replicate, not one per chromosome). A
## replicate's parsed file therefore already IS a complete two-chromosome
## bundle; there is no cross-replicate pooling to do before 02_bundle.R.
## Analyses run PER REPLICATE (PK); only the final PR/recall scoring pools
## across all ten (20 chromosomes total) -- that pooling is a later stage, not
## R_parsing/ or 02_bundle.R.
## WIDENED 2026-09-04 (PK) -- verified end to end on 21 combinations (7 cells x
## 3 reps, zero failures, ~23s/replicate) before widening further, matching
## this file's own "verify small, then scale" rule.
REPS <- 1:10

## ---- 3. SEEDS -------------------------------------------------------------------
SEEDS <- c(bundle = 1L, clusters = 11L, nulls = 41L, sensitivity = 41L)

## ---- 4. RAW SCORING: REBUILD, NOT INHERIT --------------------------------------
## [!] THE SINGLE BIGGEST ASSUMPTION IN THIS FILE, stated plainly so PK can
## override it in one line.
##
## Traced one level further than the first version of this file (2026-09-04,
## PK): regen_sim_data_* is not raw (previous paragraph in git history), and
## neither is ITS OWN input parsed_sim_data2 -- checked directly, that also
## already carries emx_p/emx_F/lfmm_p/lfmm_F/ld_w_095 baked into the map. The
## true raw substrate is NEMO's own native output (raw_nemo_nobgs/raw_nemo_bgs5
## above), one .tgz per (tag,chr,V,c,env) holding a .map + .snp_geno pair.
## R_parsing/ reads THAT directly. Confirmed by extracting one archive: exactly
## a .map, a .snp_geno, a run log and a stats file, nothing else.
##
## RAW_SCORING = "rebuild": R_parsing/ produces a bundle carrying only what NEMO
## actually simulated -- GTs, Chr/Pos/marker, per-marker type (ntrl/QTN/delet),
## allelic_values, true_QTN, MAF, and env (spatial position + environmental
## value) -- and R/ (02_bundle.R onward) recomputes decay, ld_w, stage-1
## clustering, the kinship basis and both engine scans fresh from that, under
## this file's own seeds. Nothing scored is inherited from any prior parse.
##
## [!] chr_type / max_LD_with_QTN / focal_QTN / bp_to_focal_QTN are DELIBERATELY
## NOT part of the raw parse, though the old Parse_sim_data.R computed something
## with those names at parse time. Read directly (2026-09-04): that script
## assigned chr_type by `ifelse(Chr==1,"QTN","ntrl")` -- a hardcoded per-FILE
## label, wrong on its face once a chromosome other than 1 carries a QTN, which
## checked bundles show every chromosome does (chr1..chr10 each have 1-3 QTN
## markers in V2_c1/env1). Whatever produced the CURRENT bundles' genuinely
## per-marker chr_type (QTN and ntrl values within every chromosome, matching
## the "53.4% of markers are neutral-region" finding this project already
## validated) is not this script, or not this version of it. Rather than
## inherit an unverified derivation, these are left as method-scored quantities
## for R/ to compute properly (LD-distance-to-nearest-QTN, a real choice with a
## threshold, same status as dmax/r2min) -- not raw truth.
RAW_SCORING <- "rebuild"

## Individual subsampling. keep_inds <- seq(1, 320, by = 2) in the old parser --
## every other individual, 160 of 320. PK has confirmed this is deliberate
## project policy, not a parsing artefact: "The reason I'm subsampling is
## because in a situation where the signal is clear, any method will do." Kept
## as the SAME rule, moved here so it is a config value rather than a number
## inline in a parsing script.
SUBSAMPLE_STEP <- 2L    # keep_inds <- seq(1, 320, by = SUBSAMPLE_STEP)

## MAF filter. [!] CONFLICT, NOT RESOLVED. This file inherited MAF_KEEP = 0.1
## from module_3sp on the assumption "same convention" -- never independently
## confirmed for sims. The old Parse_sim_data.R used min_maf = 0.05. Left at 0.1
## below because that is what was already committed, but flagged here rather
## than silently kept: this needs a sim-specific decision, not an inherited one.

## LFMM. [!] UPDATED 2026-09-04 (PK): unlike module_3sp, this data is entirely
## ours, so LFMM does not need to be inherited -- it can be computed fresh, in
## the same run as EMMAX, from the same GTs and env. The bundle this pipeline
## parses from (NEMO's own output) has no lfmm_p, and there is nothing to
## inherit from any prior stage; the call pattern below is reused, not
## reinvented, from the superseded module_sim_LDscnR/regen_sim_data.R, which
## already computed LFMM the same way for the same kind of data:
##   LEA::write.lfmm(GTs, ...); LEA::write.env(y, ...)
##   proj <- LEA::lfmm2(geno, env, K = LFMM_K)
##   pv   <- LEA::lfmm2.test(proj, geno, env, genomic.control = TRUE, full = TRUE)
##   lfmm_p <- pv$pvalues; lfmm_F <- pv$fscores / pv$gif
## K = 5 unchanged from that script; not re-derived here, and worth revisiting
## once more than one replicate is in scope.
LFMM_SOURCE <- "compute"
LFMM_K <- 5L

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
## SIZE_FLOOR, mechanically applying module_3sp's own rule (2x the median
## stage-1 cluster size) now that 02_bundle.R has actually reported that number
## for this replicate: median cluster size 1.00 (rep 1, V2_c1, env1, nobgs;
## 63.8% singleton clusters). 2 x 1 = 2.
##
## [!] FLAGGED, NOT A CONFIRMED PK DECISION -- this is the rule applied
## mechanically to unblock the scan stage, not a fresh judgement call, and it
## is a much thinner floor than 3sp's 8: at floor 2, almost every non-singleton
## cluster clears it, so this excludes only singletons rather than doing the
## real multiplicity-reduction work a floor is meant for. Median cluster size
## will differ by (V, c) cell (dispersal kernel c drives background LD, which
## drives cluster size) and plausibly by replicate within one cell -- this
## value is specific to rep 1 and should be re-checked, not silently reused,
## once more replicates or cells are in scope.
SIZE_FLOOR <- 2L

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
## RECEIPT MACHINERY -- from module_3sp, with one addition: an optional
## `target` subdirectory.
##
## [!] ADDED 2026-09-04, RUNNING GENUINE PARALLELISM (PK: "4 on mini2 and 4
## here"). module_3sp has exactly one target per stage, so a single shared
## receipt per stage is correct there and this is left able to reproduce that
## exact behaviour: target = "" (the default) resolves to the same path as
## before, so nothing about module_3sp or module_sim's own earlier single-
## target runs changes. module_sim now has up to 70 (cell, rep) combinations
## that can run CONCURRENTLY, and a single shared receipt.rds per stage would
## have multiple processes racing to overwrite the same file -- at best each
## write clobbers the last, silently discarding provenance for every
## combination but whichever wrote last; at worst two processes' writes
## interleave and the file becomes unreadable, crashing a THIRD, unrelated
## process's readRDS() on its own turn. Passing target = combo_id below gives
## every combination its own receipt subdirectory instead.
stage_dir <- function(stage, target = "") file.path(PATHS$out, stage, target)
receipt_path <- function(stage, target = "") file.path(stage_dir(stage, target), "_receipt.rds")

sha <- function(f) if (file.exists(f)) digest(f, algo = "sha256", file = TRUE) else NA_character_

git_sha <- function() tryCatch(system2("git", c("-C", PATHS$module, "rev-parse", "--short", "HEAD"),
                                       stdout = TRUE, stderr = FALSE), error = function(e) NA_character_)

write_receipt <- function(stage, inputs = character(), params = list(), outputs = character(), target = "") {
  dir.create(stage_dir(stage, target), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(stage = stage, when = Sys.time(), git = git_sha(),
               ldscnr = tryCatch(check_ldscnr(stop_on_fail = FALSE), error = function(e) NA),
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
  old <- setNames(r$inputs$sha256, r$inputs$path)
  ch <- names(now)[is.na(old[names(now)]) | old[names(now)] != now]
  if (length(ch)) { message("  [", stage, "] inputs changed: ", paste(basename(ch), collapse=", "),
                            " -- will run"); return(TRUE) }
  message("  [", stage, "] up to date (", format(r$when, "%Y-%m-%d %H:%M"), ")"); FALSE
}

say <- function(...) { cat(sprintf(...)); flush(stdout()) }
