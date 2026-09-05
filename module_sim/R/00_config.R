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
  ## source is NEMO's own native per-file output, one tgz per (tag,chr,V,c,env).
  ##
  ## [!] CORRECTED AGAIN, MORE SERIOUSLY, 2026-09-05 (PK: "We are at bgs5").
  ## This pointed at Nemo_out_nobgs/Nemo_out_bgs from 2026-09-04 until now --
  ## checked directly, Nemo_out_bgs's bgs arm has 3 deleterious loci in the
  ## WHOLE GENOME (adapt_bgs_chr1_V0.5_c1_env1), far below even the
  ## already-documented-unmeasurable bgs2 parameterization (200 loci, memory
  ## nemo-bgs-unmeasurable-settings). Every module_sim result built against it
  ## is archived, frozen, at ../module_sim_bgs2/ -- a numerically correct
  ## measurement of an arm with no real BGS signal to find, not a code bug.
  ## The correctly-parameterized dataset ("scenario B", 1720 deleterious loci
  ## confirmed in the same spot-check) is bgs5 -- ONE directory holding BOTH
  ## tags together (unlike the old separate Nemo_out_nobgs/Nemo_out_bgs), so
  ## both PATHS entries below now point at the same place. bgs5 currently
  ## covers 4 of CELLS_ALL's 7 cells (V0.5_c1, V0.5_c2, V1_c1.5, V2_c1) --
  ## the other 3 (V0.5_c1.5, V1_c1, V2_c1.5) are still being simulated (PK);
  ## grid drivers check archive existence and SKIP missing cells rather than
  ## fail, so dropping the remaining archives into this same bgs5/ directory
  ## needs no code change here or anywhere else.
  raw_nemo_nobgs = file.path(NEMO_ROOT, "bgs5"),
  raw_nemo_bgs5  = file.path(NEMO_ROOT, "bgs5"),
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
  ##
  ## [!] RENAMED 2026-09-05, moving to bgs5 as the raw source (see raw_nemo_bgs5
  ## above) -- a fresh directory name rather than reusing "module_sim_parsed"
  ## (now deleted) so a stray old file from the wrong raw source can never be
  ## silently read; every parsed file here is bgs5-sourced by construction.
  parsed = file.path(NEMO_ROOT, "module_sim_parsed_bgs5"),
  cache  = file.path(path.expand("~/gitlab/LDscnR-paper/module_sim"), "cache")
)
## [!] PATHS$el_dir REMOVED 2026-09-05. Was a per-combination path for
## compute_LD_decay's edge-list output -- nothing downstream ever read it
## (02_bundle.R passes gds directly to ld_complexity_reduction()), and across
## the first 1400-combination grid the edge lists it enabled accumulated to
## 382 GB. Fixed at the source: 02_bundle.R no longer passes el_data_folder
## or overrides keep_el (see DECAY_ARGS below), so nothing is written at all.
PATHS$untar  <- file.path(PATHS$cache, "untar")   # scratch for unpacked .tgz; not an output

## ---- 2. WHICH SLICE OF THE GRID THIS BUILD TARGETS -----------------------------
## Target grid: 7 cells (3 selection intensities x 3 dispersal levels, minus
## V1_c2/V2_c2 -- never simulated) x 2 BGS arms x 10 environments x 10
## REPLICATES. bgs5 currently HAS 4 of those 7 (V0.5_c1, V0.5_c2, V1_c1.5,
## V2_c1); the other 3 (V0.5_c1.5, V1_c1, V2_c1.5) are still being simulated
## (PK, 2026-09-05) -- grid drivers check archive existence and skip missing
## cells rather than fail, so this list does not need editing when they land,
## only the archives dropped into bgs5/.
##
## CELLS below is a SINGLE-CELL DEFAULT for the per-stage scripts' own
## Sys.getenv() fallback (smallest slice that lets one stage be verified end
## to end without a grid driver) -- the actual grid drivers (run_grid.sh etc.)
## carry their own CELLS_ALL array with all 7.
CELLS <- "V2_c1"          # widen to c("V0.5_c1","V0.5_c2","V1_c1.5","V2_c1") once verified
TAGS  <- "nobgs"          # other value is "bgs" (NOT "bgs5" -- that's just PATHS$raw_nemo_bgs5's
                            # own name; archive filenames and TARGET_TAG use "bgs", see
                            # R_parsing/01_parse_nemo.R's adapt_%s_chr%d_... pattern)
## [!] REP IS NOT THE STATISTICAL REPLICATE AXIS -- ENV IS. Corrected
## 2026-09-05 after PK asked directly: "is not each replicate a different
## map, and the environments are the true replicates?" Verified by untarring
## chr1_V1_c1_env1 and chr1_V1_c1_env2 and diffing them -- SAME rec_map1.rds
## (rep only varies the recombination map, shared across every env for that
## rep) but COMPLETELY DIFFERENT .snp_geno/.map contents (different md5,
## different sizes) between env1 and env2 of the same rep. So:
##   REP (chr1..chr10)  = a different recombination map/genomic architecture
##                         each time -- NOT an exchangeable replicate; pooling
##                         across reps mixes different maps.
##   ENV (env1..env10)  = same map, independent population/genotype
##                         realization (different QTN placement, different
##                         drift) -- THIS is the true replicate axis for
##                         replicate-averaging.
## R/05_pool.R replicate-averages over ENV within each REP (not across REPS)
## for exactly this reason. REP is still widened to the full grid below --
## PK: "fully cross reps x envs (10x10)" -- so results can also show whether
## precision/recall replicate ACROSS different maps, which is a genuinely
## different, useful question from within-map replicate variance.
ENVS <- 1:10
## [!] RENAMED FROM "CHRS" 2026-09-04 (PK). "chr1".."chr10" in NEMO's own
## filenames are not ten genomic chromosomes to pool -- they are ten
## independent recombination maps (see the ENVS comment above for why this
## makes REP the map axis, not the replicate axis). Each replicate already
## contains its own two chromosomes (sharing one recombination map; confirmed
## directly, rec_map1.rds has rows Chr %in% c(1,2), same bp range for both,
## and rec_map2.rds is a genuinely different map -- one per replicate, not one
## per chromosome). A replicate's parsed file therefore already IS a complete
## two-chromosome bundle; there is no cross-replicate pooling to do before
## 02_bundle.R. Analyses run PER (rep, env) COMBINATION; only the final
## PR/recall scoring pools across them -- that pooling is a later stage, not
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
  keep_el        = FALSE,  # [!] was TRUE; see 02_bundle.R's comment -- caused a 382 GB leak, fixed 2026-09-05
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
## SIZE_FLOOR = 2L, CONFIRMED (PK, 2026-09-03/04): checked across 5 sampled
## V0.5_c1 bundles (different reps/envs) -- median stage-1 cluster size pins
## at 1.00 everywhere sampled, so floor 2 consistently excludes only true
## singleton clusters, not the general multiplicity-reduction work a floor is
## meant for (unlike 3sp's floor of 8). PK's framing: SIZE_FLOOR is a
## convenience post-hoc knob for how small a cluster to still trust, not a
## value to mechanically re-derive per cell -- a smaller floor trades known
## lower precision for an unknown recall gain, and floor 2 (singletons-only)
## is the accepted operating point for this simulation grid. Not swept per
## (V,c) cell/replicate.
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
  ## radix sort, NOT the default locale-collated one -- R's sort() otherwise
  ## orders filenames by LC_COLLATE, which differs across machines (this
  ## Mac's "C" vs mini2's "C.UTF-8") and silently changes the concatenation
  ## order the hash below is computed over, aliasing a byte-identical R/
  ## directory into two different pins. Found 2026-09-04 provisioning mini2.
  src <- sort(g("ls-files","R/"), method = "radix")
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
