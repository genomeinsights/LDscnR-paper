## =============================================================================
## module_3sp/R/00_config.R
##
## EVERY PARAMETER OF THE 3sp PIPELINE, SET HERE AND NOWHERE ELSE.
##
## The contract this file exists to enforce: no stage script may contain a numeric
## constant that affects a result. If you find one, it is a bug in that stage, not a
## convenience. The reason is concrete rather than stylistic -- in the pipeline this
## replaces, `ld_w <= 0.1` in one script and `ld_w < 0.1` in another silently agreed
## on this dataset and would not have on the next, and nothing could have detected it
## because neither script knew the other existed.
##
## Sourced by every stage. Defines parameters, paths, and the receipt machinery that
## makes "what needs rerunning" a computed answer instead of a judgement call.
## =============================================================================
suppressMessages({library(data.table); library(digest)})

## ---- 1. WHERE THINGS ARE ----------------------------------------------------
## Five locations, two of which hold data that cannot be committed. Named here so a
## stage never hardcodes one, and so a machine with a different layout needs one edit.
PATHS <- list(
  module   = path.expand("~/gitlab/LDscnR-paper/module_3sp"),
  out      = path.expand("~/gitlab/LDscnR-paper/module_3sp/out"),
  ## RAW INPUTS -- never written by this pipeline.
  ##
  ## COPIED IN so the module is self-contained on this filesystem (PK). Both are
  ## byte-verified copies of the originals under LD-scaling-genome-scans, checked by
  ## SHA-256 at copy time and re-checked by 01_inputs.R on every run. They are NOT in
  ## git -- 3sp_data.RData is 159 MB, over GitHub's 100 MB hard limit, and bzip2 only
  ## reaches 103 MB -- so this is self-containment on disk, not in the repository. Those
  ## are different properties and only the second survives a clone.
  raw_3sp  = path.expand("~/gitlab/LDscnR-paper/module_3sp/data/3sp_data.RData"),
  ## LFMM F-values over the FULL pre-MAF map. See LFMM_SOURCE below -- inherited, not
  ## computed by anything in any repository.
  lfmm_F   = path.expand("~/gitlab/LDscnR-paper/module_3sp/data/lfmm_F.rds"),
  ## NOT copied, deliberately. These are TRACKED in this repository already and OWNED by
  ## other sessions -- the BEDs by ldscnr-fe (kingman2021), the map by ldscnr-d4. A second
  ## copy inside module_3sp would be a copy that can silently drift from the maintained
  ## one, and the failure would look like a result changing for no reason. Self-contained
  ## does not mean duplicated.
  ecopeaks = path.expand("~/gitlab/LDscnR-paper/kingman2021/data/liftover"),
  recmap   = path.expand("~/gitlab/LDscnR-paper/3sp_data/rec_maps/stickleback_recomb_3crosses_gasAcu1.tsv"),
  ## Provenance of the two copied files: where they came from and their SHA-256 at copy
  ## time, as a RECORD rather than a live path. An earlier draft kept origin_* pointing at
  ## LD-scaling-genome-scans so the copies could be re-verified against them, which
  ## reintroduced exactly the external dependency the copy removed (PK). Checking against
  ## a recorded hash catches a corrupted or replaced local copy and needs no other
  ## repository to exist.
  provenance = path.expand("~/gitlab/LDscnR-paper/module_3sp/data/PROVENANCE.csv"),
  ## CACHE -- large, regenerable, gitignored. Inside the module (PK), so the whole
  ## pipeline lives under one root and nothing resolves to another repository. It will
  ## reach GB scale: the GDS, the edge lists and the bundle all land here.
  cache    = path.expand("~/gitlab/LDscnR-paper/module_3sp/cache")
)
## Edge lists from compute_LD_decay. The pipeline being replaced pointed this at
## /Volumes/Nemo, an external volume, which is one reason its bundle could not be
## rebuilt on demand. Kept local and regenerable.
PATHS$el_dir <- file.path(PATHS$cache, "edge_lists")

## The EcoPeak BEDs, named explicitly rather than globbed: a glob would silently
## widen if the kingman session adds a cohort, and the overlap denominators would
## change with no diff to show it.
ECOPEAK_BEDS <- c("pv_c155.specific.bed", "pv_c150.specific.bed")

## ---- 2. SEEDS ---------------------------------------------------------------
## One seed per stage, not one global seed. A global seed makes every stage's stream
## depend on how much randomness the stages before it consumed, so adding a
## diagnostic to stage 05 would silently change stage 06's nulls.
SEEDS <- c(bundle = 1L, clusters = 11L, nulls = 41L, sensitivity = 41L, external = 41L)

## ---- 3. STAGE 02: THE BUNDLE ------------------------------------------------
## compute_LD_decay subsamples in three places (background, per-chromosome thinning,
## stratified pair sampling). The bundle being replaced was fitted before seeding
## existed, so its ld_w -- which sets BOTH the GRM basis and the stage-1 clustering --
## cannot be reproduced even in principle. That is the single strongest reason to
## rebuild rather than inherit, and it means the new numbers may differ slightly from
## the published ones. Expected, not a failure.
DECAY_ARGS <- list(
  min_maf_decay  = 0.1,
  q              = 0.95,
  n_sub_bg       = 5000,
  n_win_decay    = 20,        # CANONICAL. The sweep at 10/50 is stage 09, not here.
  overlap        = 0.5,
  max_SNPs_decay = Inf,      # must be uncapped except if you are only using this for LD-decay estimation
  prob_robust    = 0.95,
  max_pairs      = 5000,
  ld_method      = "corr",
  n_strata       = 20,
  keep_el        = TRUE,
  slide          = 500,     # LD-decays faster so 500 is enough
  rho_targets    = 0.99,
  cores          = 1
)
RHO_GRID <- c(seq(0.05, 0.95, by = 0.05), 0.99)   # ld_w columns computed in place, only 0.95 is used.

## MAF filter, applied to the raw panel before anything else. It was NOT in my first
## draft of this file and it should have been: it sets how many markers exist, so it
## sets the multiplicity of every test downstream. Exactly the kind of constant that
## was living inside a script rather than in a config.
MAF_KEEP <- 0.1                   # markers with maf > MAF_KEEP are retained

## LFMM. THE PIPELINE DOES NOT COMPUTE THIS AND CANNOT. regen_3sp_data.R reads
## lfmm_F.rds -- F-values over the full pre-MAF map, dated Sep 2025 -- and no script in
## any repository produces it. So the LFMM engine is INHERITED, and a "from scratch"
## rebuild is from scratch for everything except this.
##   "inherit"   use lfmm_F.rds, hash it, and say so in the receipt (the honest default)
##   "omit"      drop the LFMM engine entirely; EMMAX-only pipeline, fully reproducible
## Set to "omit" only if PK decides the second engine is not worth an unreproducible root.
LFMM_SOURCE <- "inherit"

## GRM. The estimator is load-bearing and was discovered to be so only by accident:
## GCTA and a plain centred product over the IDENTICAL markers give 79 and 59
## discoveries from kinships correlating 0.9898. Stage 08 tests both; this is the one
## the pipeline reports, and it is a choice rather than a default.
## ---- GRM BASIS: which markers define the kinship ---------------------------
## STAGE-1 PRUNING IS THE REPORTED BASIS, and the reason is that it is not a separate
## decision at all. THE SAME OPERATION THAT DEFINES THE TEST UNITS SELECTS THE KINSHIP
## MARKERS (PK) -- one choice serving two purposes. An `ld_w` threshold, which is what
## the superseded bundle used, is a second arbitrary number existing only to serve the
## kinship, and it has to be justified on its own.
##
## THE MEASUREMENTS THAT MADE THIS THE EASY CALL (DA-era sweeps, all centred rebuilds so
## the estimator is held fixed):
##
##   basis                              markers   disc   regions   on peak
##   stage-1 representatives            616,497     57        21         9
##   ld_w <= 0.1                        742,858     59        22        10
##   greedy, SNPRelate r2<0.2 / 500 kb   31,070     58        --        --
##   NO PRUNING                         790,446      2         2         0   <- 52% FDP, p = 0.14
##
## So WHICH pruning is nearly irrelevant and WHETHER you prune is decisive. The bases
## also agree on the kinship itself, not merely on the answer: a fifth of the markers
## differ between stage-1 representatives and ld_w <= 0.1, and their off-diagonal
## correlation is r = 0.9994. With 117 individuals, ~600k markers is far past the point
## where swapping markers moves a relatedness estimate.
##
## The manuscript therefore compares THREE arms and no more: stage-1 pruning as the
## principled basis, with greedy pruning and no pruning in the supplement (PK).
GRM_BASIS      <- "stage1_pruned"          # reported: representatives of the stage-1 clusters
GRM_BASIS_SUPP <- c("greedy", "none")      # supplementary comparators only
GRM_GREEDY     <- list(ld.threshold = 0.2, slide.max.bp = 5e5)   # SNPRelate snpgdsLDpruning

## THE ESTIMATOR IS THE PART THAT ACTUALLY MOVES THINGS, and it was discovered by
## accident: GCTA and a plain centred product over the IDENTICAL markers give 79 and 59
## discoveries from kinships correlating 0.9898. Reported as GCTA, tested in the sweep.
GRM_METHOD <- "GCTA"

## `ld_w` no longer selects the GRM basis. It survives only as the flagging input to
## stage-1 pruning itself; the threshold below is kept because ld_prune_and_eMLG takes it,
## NOT because it defines the kinship.
GRM_LDW_THRESHOLD <- 0.1
GRM_LDW_OP        <- "<"    # STRICT. The two conventions agreed on this dataset by luck.

## [!] ORDERING CONSEQUENCE. With the basis being stage-1 representatives, the STAGE-1
## CLUSTERING MUST RUN BEFORE THE KINSHIP. regen_3sp_data.R already computes
## ld_complexity_reduction before the GRM and simply does not use it for the basis, so the
## order is available -- but stage 02 and stage 03 are no longer cleanly separable and
## 02_bundle.R has to produce the clustering it prunes from. The sim side already works
## this way (regen_sim_data.R takes ld_prune_and_eMLG(...)$pruned when GRM_METHOD is not
## "ld_w_threshold"), so this aligns the two halves rather than diverging them.

## ---- 4. STAGE 03: CLUSTERING ------------------------------------------------
CR_RHO             <- 0.5      # ld_complexity_reduction, stage-1
DISTANCE_THRESHOLD <- 1e5      # stage-2 grouping, by distance between clusters
MIN_R2_RHO         <- 0.5      # derived per chromosome from the decay fit
SCORE_THRESHOLD    <- 0.80
LDW_FLAG           <- 0.05     # ours; fe's assembly uses 0.025, run alongside in stage 07

## ---- 5. STAGES 04-07: TESTING -----------------------------------------------
ALPHA      <- 0.05
SIZE_FLOOR <- 8L      # 2 x the median stage-1 cluster size of 4.11 (PK)

## ---- REGION ASSEMBLY (post hoc; cannot affect a p-value) --------------------
## THE STAGE-2 MACHINERY MERGES THE SIGNIFICANT STAGE-1 CLUSTERS. Not a fixed physical
## distance, which is what this line used to be.
##
## PK's point, and it is what makes this free: used only to merge already-significant
## clusters, `ld_w_threshold` IS NOT A GATE. It exists to skip pairwise comparisons that
## could not have clustered anyway, so it can be 0 and the assembly introduces NO new
## parameters -- it inherits stage-1's own rho and the canonical score and distance
## thresholds. My objection to LD-aware assembly was that it dragged four parameters into
## the region definition; used this way it drags none.
##
## MERGE OVER THE DISCOVERED CLUSTERS ONLY. This is the whole difference between this and
## the other LD-aware assembly in the literature of this project: running stage-2 over the
## FULL partition and keeping groups that contain a discovery lets a region inherit the
## extent of a group built mostly from UNDISCOVERED clusters, which is how a 19.11 Mb
## region arises from clusters whose own spans are kilobases. Merging over the discovered
## clusters means a region can only span discovered signal.
##
## Measured on this panel (DA): 34 regions, 14 on EcoPeaks, fold 2.61x, span-preserving
## rotation p = 0.0002, median occupancy 0.649, widest 1.14 Mb.
REGION_ASSEMBLY <- list(
  ld_w_threshold     = 0,                    # NOT a gate here -- 0 filters nothing
  min_r2_rho         = MIN_R2_RHO,           # the rho that defined stage-1
  score_threshold    = SCORE_THRESHOLD,      # 0.80
  distance_threshold = DISTANCE_THRESHOLD,   # 1e5
  merge_over         = "discovered"          # NOT the full partition
)

## The old physical merge, RETAINED AS A CHECK rather than as the report. The two agreed
## closely on this panel -- 33 regions against 34, fold 2.60x against 2.61x -- and an
## agreement between an arbitrary rule and a motivated one is worth continuing to measure,
## because it is the evidence that the reported regions are not an artefact of either.
REGION_GAP_CHECK <- 3e5
STATISTICS <- c("consensus", "Simes")
ENGINES    <- c("EMMAX", "LFMM")

## ---- HOW A CLUSTER IS REPRESENTED AS ONE VARIABLE ---------------------------
## A SEPARATE AXIS FROM `STATISTICS` above, and they are easy to confuse. `STATISTICS`
## is how a cluster's evidence is combined (one test on a summary variable, or Simes over
## its members' p-values). THIS is what that summary variable IS. Four candidates, all to
## be tested (PK):
##
##   eMLG           make_eMLGs(GTs, map_cl, cor_th, l_min)$eMLG -- the package's own
##                  block consensus, individuals x clusters. Ideal when you want one value
##                  per LD block.
##   representative the cluster's core SNP. Chosen for CENTRALITY (highest median r2 to
##                  the rest of the cluster), not for signal.
##   best_snp       eMLG_best_snp(result, GTs, fill = TRUE) -- the member SNP most
##                  correlated with the block consensus, with its missing calls filled
##                  FROM the consensus. The middle option: a real SNP, but picked by
##                  agreement with the block rather than by centrality.
##   consensus_dosage  LDscnR:::consensus_dosage -- polarise members to a common allele,
##                  then row-mean. THIS IS WHAT EVERY REPORTED NUMBER SO FAR USED, and it
##                  is NOT obviously the same object as make_eMLGs()'s eMLG. Included so
##                  that question gets answered rather than assumed.
##
## WHY IT MATTERS RATHER THAN BEING A DETAIL: the eMLG averages a block's markers, which
## dilutes signal that is genuinely SNP-specific -- differing even between strongly linked
## markers. Where that happens a single SNP is better; where it does not, averaging is
## better estimated. Which regime this panel is in is an empirical question we have never
## put. All four go in the sensitivity sweep and we decide afterwards what to keep.
UNIT_REPR      <- "consensus_dosage"      # reported, for continuity with published numbers
UNIT_REPR_ALL  <- c("consensus_dosage", "eMLG", "representative", "best_snp")
EMLG_ARGS      <- list(input = "auto", cor_th = 0.8, l_min = 10)   # make_eMLGs
BEST_SNP_ARGS  <- list(fill = TRUE, round_fill = TRUE)             # eMLG_best_snp
##
## [!] l_min = 10 is make_eMLGs' own default and is LARGER than our SIZE_FLOOR of 8, so
## at the default it would silently drop clusters we test. Stage 04 must either set
## l_min <= SIZE_FLOOR or report how many units it loses; do not leave that implicit.

## ---- 6. STAGE 06: NULLS -----------------------------------------------------
## Two budgets because the costs differ by ~100x: a consensus surrogate is one scan
## over ~1,300 columns, a Simes surrogate rescans all ~790,000 markers.
NPERM_CONSENSUS <- 1000L
NPERM_SIMES     <- 200L
N_ROTATIONS     <- 10000L
ROTATION_SCHEME <- "within"   # nominated in advance; EcoPeaks are non-uniform among chromosomes
PERM_SCHEMES    <- c("regional", "population")   # 'individual' is invalid here, reported only as the reference that overstates

## ---- 7. STAGE 08: SENSITIVITY -----------------------------------------------
## grm_basis replaces the old grm_ldw axis: the ld_w series was four points on a
## parameter that no longer selects the basis, and it was also the series that mixed
## estimators and produced the "57-79 across 24-fold" error. Three named arms instead.
SWEEP <- list(
  rho           = c(0.35, 0.50, 0.65),
  size_floor    = c(4L, 6L, 8L, 12L, 20L),
  grm_basis     = c("stage1_pruned", "greedy", "none"),
  grm_estimator = c("GCTA", "centred"),
  unit_repr     = UNIT_REPR_ALL
)

## ---- 8. STAGE 09: EXTERNAL VALIDATION ---------------------------------------
## 20 is canonical (PK, section 3). 10 and 50 are the sensitivity arms.
## 5 is DROPPED: it was the value the superseded bundle used, and keeping it would invite
## comparing the new pipeline against a fit that also differs in seeding -- two changes at
## once. Its result is on record in 3SP_RESULTS.md and does not need recomputing.
## At 50 the windows are ~362 kb, FINER than the map's own 489 kb median bin, so a fall
## there is the reference running out rather than the proxy failing.
DECAY_NWIN_SWEEP <- c(10, 20, 50)
## WHAT "CONCORDANT ONLY" MEANS. A subset of the map's 880 bins where the three F2
## crosses AGREE about the local rate. From the map's own README (ldscnr-d4):
##
##   disagree    rank_sd ranked WITHIN rate decile; uniform on [0,1] and orthogonal to
##               rate. 0 = the crosses agree about this bin relative to others of
##               similar rate.
##   concordant  >= 2 crosses AND disagree <= 0.75      (655 of 880 bins)
##
## Restricting to it tests the decay proxy against BETTER-MEASURED reference values, so a
## real correlation should hold up at least as well there. It does: marginally stronger at
## every window count, and 19/19 chromosomes positive against 17/19 at n_win 5.
##
## [!] THE HISTORY MATTERS, because the flag changed under us. `concordant` was
## previously built on `rel_sd`, which correlates +0.635 (Pearson) / +0.716 (Spearman)
## with rate -- so it retained 100% of the lowest-rate decile and 25% of the highest. That
## is a RATE FILTER WEARING A QUALITY LABEL, and under it the correlation attenuated,
## which I misread as evidence against the claim. d4 rebuilt the flag on `disagree`,
## orthogonal to rate by construction; 206 of 880 bins flipped. `rel_sd` is now marked in
## the map's README as "descriptive only -- do not filter on it". Never condition on
## rel_sd, and check cor(flag, predictor) before reading any subset comparison.
MAP_SUBSETS      <- c("all bins", "concordant only")


## =============================================================================
## THE LDscnR VERSION IS PART OF THE PIPELINE, so it is pinned and checked.
##
## `Version` alone cannot do this: the package is 0.0.0.9000 on every commit, so two
## builds three hours and one substantive commit apart are indistinguishable by version.
## The pin is therefore the SOURCE COMMIT, and the check is whether the INSTALLED build
## could possibly contain it.
##
## THIS WAS NOT HYPOTHETICAL WHEN IT WAS WRITTEN. The installed build was dated
## 2026-09-02 15:59:40 while R/ld_group_map.R had been modified at 18:58 and HEAD
## committed at 19:00 -- so every result produced that day, including the stage-2 region
## assembly comparison, ran against source three hours older than the repository. Nothing
## reported it because nothing looked.
##
## The branch situation is now resolved: pvalue-api was 27 commits ahead of main and has
## been merged as a FAST-FORWARD, so main and pvalue-api are the same commit and the
## pipeline no longer depends on an unmerged branch. R CMD check on the built tarball is
## Status: OK (0 errors, 0 warnings, 0 notes) and the test suite passes.
##
## The pin below is the merge result. It differs from the commit the package was BUILT
## from, and that is fine and deliberately tolerated: the only later commit touches
## README.md, so no file under R/ is newer than the build. That is exactly why the check
## tests R/ mtimes against the build time rather than comparing SHAs alone -- a SHA
## comparison would demand a pointless reinstall for a documentation commit.
LDSCNR_PIN <- list(
  repo    = path.expand("~/gitlab/LDscnR"),
  sha     = "011a165c8c99",
  branch  = "outlier-scan", # the promoted API is on this branch, not yet merged to main
  version = "0.0.0.9000",
  ## SHA-256 over the CONCATENATED CONTENT of every tracked file under R/, taken at the
  ## moment the package was installed from this commit. This is the check that works.
  src_sha = "362892c1ee4c6628",
  installed_after = NULL   # informational only; set at each pin refresh, never checked
)

## WHY CONTENT AND NOT MTIMES. The first version of this check compared each tracked R/
## file's mtime against the installed build time. It was defeated immediately by an
## ordinary `git checkout`: switching to main rewrote the mtime of all 16 changed files to
## the present, so a package installed from exactly that source was reported as stale
## against 16 files. mtime records when git touched a file, not what is in it.
##
## Content hashing is immune to that, and it is the same discipline as inputs_3sp: hash
## the bytes, never trust a timestamp. The pin is asserted by HAVING INSTALLED from this
## commit -- `installed_after` records when -- and the hash then detects any later drift
## of the source away from what was installed.
check_ldscnr <- function(stop_on_fail = !nzchar(Sys.getenv("LDSCNR_LAX"))) {
  v <- as.character(utils::packageVersion("LDscnR"))
  d <- utils::packageDescription("LDscnR")
  ## d$Built can be ABSENT, not merely malformed: devtools::load_all() -- which
  ## 02_bundle.R and every module_3sp scan script call, deliberately, per PK, since
  ## development is still active -- attaches the source tree directly and
  ## packageDescription() then returns a DESCRIPTION with no Built field at all (that
  ## field is stamped only by R CMD INSTALL/build). sub() on NULL silently returns
  ## character(0), as.POSIXct(character(0)) is a length-0 POSIXct, and a length-0
  ## condition in if() is a hard error in current R -- found by running 03_EMMAX.R for
  ## real, not by the interactive checks that verified this function earlier, none of
  ## which had also called load_all() first.
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
  ## `installed_after` is INFORMATIONAL, not a gate -- it once was, and failed on a
  ## six-second gap between capturing wall-clock time and R's Built field (which
  ## truncates to the minute). The src_sha comparison already proves the install
  ## reflects this exact source; a timestamp cannot prove anything the hash does not.
  ## Already established as INFORMATIONAL, not a gate (a real six-second-gap
  ## failure earlier today) -- always TRUE. Absent under load_all() is a second,
  ## separate reason it must never gate: load_all() runs whatever source is on
  ## disk right now, so the question this timestamp exists to answer (does the
  ## install reflect the current source?) is moot by construction. src_sha above
  ## is the check that still matters, and it is unaffected by any of this.
  built_ok <- TRUE
  ## WHAT COUNTS AS A FAILURE, and this took three attempts to get right.
  ## The invariant that matters is "the installed code IS the current code", and that is
  ## the SOURCE HASH -- not the mtimes (rewritten by any checkout) and not the commit id
  ## (changed by a README edit that cannot affect the install). The commit is recorded as
  ## provenance and its drift is REPORTED, not fatal.
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
## RECEIPT MACHINERY
##
## Each stage writes a receipt recording its parameters, the SHA-256 of every input,
## the git SHA and the time. A stage is STALE if any of those changed. This is what
## makes "what needs rerunning" computable rather than remembered -- which matters
## because the honest answer today, across a chain that crosses two repositories, is
## that nobody knows.
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

## TRUE when the stage must run: no receipt, a changed parameter, or a changed input.
## Deliberately conservative -- it reports WHY, so a surprising rerun is explainable
## rather than mysterious.
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
