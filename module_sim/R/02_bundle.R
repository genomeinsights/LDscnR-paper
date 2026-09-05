## =============================================================================
## module_sim/R/02_bundle.R
##
## BUILD THE DATA BUNDLE FOR ONE REPLICATE: GDS, LD decay, ld_w, stage-1
## clustering, kinship. Mirrors module_3sp/R/02_bundle.R's structure and its
## checkpointing machinery (copied, not shared, since that machinery is not yet
## factored into 00_config.R on either side).
##
## [!] CORRECTED 2026-09-04, after the first version's design turned out to be
## wrong. It POOLED across CHRS (a loop reading and cbind/rbind-ing several
## "chromosome files" into one panel), modelled on module_3sp's single-panel
## genome. That model does not apply here: "chr1".."chr10" are ten INDEPENDENT
## REPLICATES (PK), each already a complete two-chromosome bundle once
## R_parsing/01_parse_nemo.R parses it correctly (that script had its own,
## related bug -- collapsing both chromosomes into one label -- fixed the same
## day). There is nothing to pool at this stage: one replicate's one parsed
## file already contains both its chromosomes. Analyses run PER REPLICATE;
## pooling across all ten (20 chromosomes) happens only at the final PR/recall
## scoring stage, not here.
##
## Runs no association test -- same separation of concerns as module_3sp:
## the scan is a later stage, because it depends on the kinship this stage
## produces, and conflating the two is exactly what made "rebuild the kinship"
## and "rerun the scan" impossible to decide independently in every pipeline
## this project has since replaced.
##
## THE KINSHIP BASIS IS THE REASON THIS STAGE EXISTS. PK asked directly whether
## the sim bundles used module_3sp's canonical basis; they did not (checked
## 2026-09-04, see 00_config.R's RAW_SCORING note and the R_parsing/01
## commit message) -- the sim GRM_METHOD="complexity_chain" routes through
## ld_prune_and_eMLG()$pruned, the stage-2 group representatives module_3sp's
## own 02_bundle.R explicitly calls "a DIFFERENT set... deliberately not used".
## This stage builds the CORRECT basis: unique(na.omit(stage1$pruned)) --
## clusters$core_snp straight off ld_complexity_reduction(), same as 3sp,
## verified against the package source rather than assumed.
##
## SCOPE: one replicate (REPS currently 1L in 00_config.R). Widen REPS and rerun
## this stage per replicate once this one is verified end to end -- there is no
## per-replicate loop to write here, each replicate's bundle is independent all
## the way through the scan stage.
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(SNPRelate); library(digest)})
devtools::load_all("~/gitlab/LDscnR")
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim"), "R", "00_config.R"))
STAGE <- "02_bundle"

## ---- CHECKPOINTS for the two expensive steps, same rationale as module_3sp:
## the decay fit and the clustering are what this stage can die partway
## through, and a fingerprinted cache means a crash or a downstream-bug rerun
## does not repeat either from scratch. stage1's fingerprint folds decay_fp in,
## so any decay-parameter change invalidates the clustering too, transitively.
.cache_step <- function(name, fp, compute) {
  f <- file.path(PATHS$cache, paste0(name, ".rds"))
  if (file.exists(f) && !nzchar(Sys.getenv("FORCE"))) {
    x <- readRDS(f)
    if (identical(x$fp, fp)) {
      say("    [cache hit] %s (computed %s) -- skipping\n", name,
          format(x$when, "%Y-%m-%d %H:%M"))
      return(x$value)
    }
    say("    [cache stale] %s -- fingerprint changed, recomputing\n", name)
  }
  t0 <- Sys.time()
  value <- compute()
  dir.create(PATHS$cache, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(fp = fp, value = value, when = Sys.time()), f)
  say("    [cached] %s -> %s (%.1f min)\n", name, f,
      as.numeric(difftime(Sys.time(), t0, units = "mins")))
  value
}

say("=== %s ===\n\n", STAGE)
invisible(check_ldscnr())

TARGET_TAG <- Sys.getenv("SIM_TAG", TAGS[1]); TARGET_CELL <- Sys.getenv("SIM_CELL", CELLS[1]); TARGET_ENV <- as.integer(Sys.getenv("SIM_ENV", ENVS[1])); TARGET_REP <- as.integer(Sys.getenv("SIM_REP", REPS[1]))
## combo_id keys every per-combination path below (GDS, edge lists, decay/
## clustering cache, receipt) -- defined here, before the staleness check,
## since that check itself needs it as the receipt target.
combo_id <- sprintf("%s_%s_rep%d_env%d", TARGET_TAG, TARGET_CELL, TARGET_REP, TARGET_ENV)
parsed_file <- file.path(PATHS$parsed,
  sprintf("nemo_%s_rep%d_%s_env%d.rds", TARGET_TAG, TARGET_REP, TARGET_CELL, TARGET_ENV))
if (!file.exists(parsed_file)) stop("R_parsing/01_parse_nemo.R has not produced: ",
  basename(parsed_file), " -- run it for this replicate first.")

## Staleness check BEFORE the expensive work, matching module_3sp's placement
## (computed immediately after resolving the raw input, not at the end next to
## write_receipt()) -- a rerun with nothing changed should not pay for the
## whole decay fit before finding out it didn't need to.
INPUTS <- unname(parsed_file)
PARAMS <- list(rep = TARGET_REP, decay_args = DECAY_ARGS, rho_grid = RHO_GRID,
               cr_rho = CR_RHO, grm_basis = GRM_BASIS, grm_method = GRM_METHOD,
               seed_bundle = SEEDS[["bundle"]], seed_clusters = SEEDS[["clusters"]])
if (!stage_stale(STAGE, INPUTS, PARAMS, target = combo_id) && !nzchar(Sys.getenv("FORCE"))) {
  say("\nNothing to do. Set FORCE=1 to rerun anyway.\n"); quit(save = "no")
}

## ---- 1. read the replicate's parsed bundle ------------------------------------
## Already has both chromosomes -- no pooling. See header for why the first
## version's cross-file pool was the wrong model.
say("[1] reading %s\n", basename(parsed_file))
b <- readRDS(parsed_file)
GTs <- b$GTs; map <- b$map; env <- b$env
stopifnot("GTs columns and map rows must correspond 1:1" = identical(colnames(GTs), map$marker),
          "GTs rows and env rows must correspond 1:1" = nrow(GTs) == nrow(env))
say("    %d individuals x %s markers ; %d QTN ; chromosomes: %s\n", nrow(GTs),
    format(ncol(GTs), big.mark = ","), sum(map$true_QTN),
    paste(sprintf("%s (%s markers, %d QTN)", names(table(map$Chr)),
          format(table(map$Chr), big.mark = ","),
          sapply(split(map$type, map$Chr), function(t) sum(t == "QTN"))), collapse = ", "))

## ---- 2. GDS ------------------------------------------------------------------
## [!] PER-COMBINATION PATH, NOT PATHS$cache/sim.gds. That single shared path
## was safe under the sequential driver (deleted and rebuilt fresh every run,
## nothing to alias) but is NOT safe once mini2 and this machine (or several
## processes on one machine) build different combinations' GDS files at the
## same time -- concurrent delete+rebuild on one path would corrupt whichever
## process loses the race. Fixed 2026-09-04, widening to real parallelism (PK:
## "4 on mini2 and 4 here").
gds_path <- file.path(PATHS$cache, paste0("sim_", combo_id, ".gds"))
dir.create(dirname(gds_path), recursive = TRUE, showWarnings = FALSE)
if (file.exists(gds_path)) { say("\n[2] removing stale %s\n", basename(gds_path)); unlink(gds_path) }
say("\n[2] GDS -> %s\n", gds_path)
gds <- create_gds_from_geno(geno = GTs, map = map, gds_path)
on.exit({ try(snpgdsClose(gds), silent = TRUE); unlink(gds_path) }, add = TRUE)
stopifnot(file.exists(gds_path), file.size(gds_path) > 0)
say("    %.1f MB\n", file.size(gds_path) / 1e6)

## ---- 3. LD decay, SEEDED, ld_w in place ---------------------------------------
say("\n[3] LD decay: n_win_decay = %d, seed %d\n", DECAY_ARGS$n_win_decay, SEEDS[["bundle"]])
## [!] tag/cell ADDED 2026-09-04, widening REPS to a real grid rather than one
## cell. The fingerprint previously omitted them, relying only on rep + marker
## count -- with a single cell ever tested that was invisible, but two
## different cells could in principle land on the same marker count after
## their own MAF filtering and silently alias a cached decay/clustering. The
## cache FILENAME is now also keyed per (tag, cell, rep) below, which is the
## primary fix; this closes the same gap in the fingerprint itself.
decay_fp <- digest(list(decay_args = DECAY_ARGS, rho_grid = RHO_GRID,
                        tag = TARGET_TAG, cell = TARGET_CELL, rep = TARGET_REP,
                        seed = SEEDS[["bundle"]], markers = ncol(GTs)), algo = "sha256")
## el_data_folder ALSO PER-COMBINATION, same reasoning as gds_path above --
## compute_LD_decay writes "Chr1.el"/"Chr2.el" into this folder on every call;
## a shared el_dir would let concurrent combinations overwrite each other's
## edge lists mid-write. Nothing downstream reads these back (confirmed
## 2026-09-04: ld_complexity_reduction below is passed gds directly, not
## el_data_folder), so this is pure side-effect output, but a write race is
## still worth avoiding rather than shrugging off because it happens not to be
## read.
el_dir <- file.path(PATHS$el_dir, combo_id)
LD_decay <- .cache_step(paste0("ld_decay_", combo_id), decay_fp, function() {
  set.seed(SEEDS[["bundle"]])
  do.call(compute_LD_decay,
         c(list(gds = gds, el_data_folder = el_dir, ld_w_rho = RHO_GRID,
                seed = SEEDS[["bundle"]]),
           DECAY_ARGS)) })
## [!] CLEANED UP IMMEDIATELY. Found 2026-09-05 rebuilding for bgs5: nothing
## ever removed these, and across the first 1400-combination grid they
## (never read back downstream, per this section's own comment above)
## accumulated to 382 GB before anyone noticed. el_data_folder still writes
## them (needed internally during this call), but there is no reason to keep
## the result once compute_LD_decay has returned.
##
## [!] el_dir IS A FILE PREFIX, NOT A DIRECTORY. compute_LD_decay writes
## "<el_dir>Chr1.el"/"<el_dir>Chr2.el" literally concatenated (no path
## separator) -- unlink(el_dir, recursive=TRUE) is a silent no-op against a
## path that was never actually created, which is why the first attempt at
## this cleanup (same day) did not work: measured 560 MB surviving one
## combination's smoke test that should have logged this line. Glob instead.
unlink(Sys.glob(paste0(el_dir, "*.el")))
ld_ws <- LD_decay$ld_ws[map$marker, , drop = FALSE]
ld95  <- if ("rho_0.95" %in% colnames(ld_ws)) "rho_0.95" else "0.95"
map[, ld_w_095 := ld_ws[, ld95]]
say("    %d chromosome(s) ; ld_w matrix %s x %d\n", nrow(LD_decay$decay_sum),
    format(nrow(ld_ws), big.mark = ","), ncol(ld_ws))

## ---- 4. stage-1 clustering -- BEFORE the kinship ------------------------------
say("\n[4] stage-1 clustering (ld_complexity_reduction, rho = %.2f)\n", CR_RHO)
stage1_fp <- digest(list(decay_fp = decay_fp, cr_rho = CR_RHO), algo = "sha256")
stage1 <- .cache_step(paste0("stage1_", combo_id), stage1_fp, function() {
  set.seed(SEEDS[["clusters"]])
  ld_complexity_reduction(map = map, LD_decay = LD_decay, rho = CR_RHO, gds = gds) })
cl <- as.data.table(stage1$clusters)
nl <- if ("n_loci" %in% names(cl)) cl$n_loci else lengths(cl$members)
say("    %s clusters ; median size %.2f -- SIZE_FLOOR is not yet set (00_config.R),\n",
    format(nrow(cl), big.mark = ","), median(nl))
say("    this is the number that should decide it, per module_3sp's own 2x-median convention.\n")

## ---- 5. the kinship basis IS the stage-1 representatives ----------------------
## unique(na.omit(stage1$pruned)) -- clusters$core_snp, straight off the
## ld_complexity_reduction() object. NOT ld_prune_and_eMLG()$pruned (the
## superseded sim bundles' basis) -- verified this is a materially different,
## coarser set (00_config.R's RAW_SCORING note; R_parsing/01's commit message).
say("\n[5] kinship basis: %s (module_3sp's canonical basis, not the superseded\n", GRM_BASIS)
say("    sim bundles' complexity_chain/ld_prune_and_eMLG basis)\n")
grm_markers <- switch(GRM_BASIS,
  stage1_pruned = unique(na.omit(stage1$pruned)),
  greedy = { set.seed(SEEDS[["bundle"]])
    unlist(snpgdsLDpruning(gds, ld.threshold = 0.2, slide.max.bp = 5e5,
                           autosome.only = FALSE, verbose = FALSE), use.names = FALSE) },
  none = map$marker,
  stop("unknown GRM_BASIS: ", GRM_BASIS))
say("    %s of %s markers (%.1f%%)\n", format(length(grm_markers), big.mark = ","),
    format(nrow(map), big.mark = ","), 100 * length(grm_markers) / nrow(map))

say("    GRM: snpgdsGRM(method = \"%s\")\n", GRM_METHOD)
t0 <- Sys.time()
GRM <- snpgdsGRM(gds, snp.id = grm_markers, method = GRM_METHOD,
                 verbose = FALSE, autosome.only = FALSE)$grm
ut <- upper.tri(GRM)
say("    %d x %d ; mean diagonal %.4f ; off-diagonal mean %+.4f sd %.4f ; %.2f min\n",
    nrow(GRM), ncol(GRM), mean(diag(GRM)), mean(GRM[ut]), sd(GRM[ut]),
    as.numeric(difftime(Sys.time(), t0, units = "mins")))

## ---- 6. save -----------------------------------------------------------------
## No emx, no lfmm: those are a later stage. Ground truth (type, true_QTN,
## allelic_values, MAF) carried through unchanged from R_parsing/ -- this stage
## does not touch truth, only method choices.
OUT <- file.path(stage_dir(STAGE), sprintf("bundle_%s_rep%d_%s_env%d.rds", TARGET_TAG, TARGET_REP, TARGET_CELL, TARGET_ENV))
dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)   # OUTPUT dir; safe to share -- dir.create is idempotent, and combo-specific FILENAMES already disambiguate the .rds itself
saveRDS(list(
  GTs = GTs, map = map, env = env,
  ld_ws = ld_ws, LD_decay = LD_decay,
  stage1 = stage1,
  GRM = GRM, grm_markers = grm_markers,
  settings = list(tag = TARGET_TAG, cell = TARGET_CELL, env = TARGET_ENV, rep = TARGET_REP,
                  decay_args = DECAY_ARGS, rho_grid = RHO_GRID, cr_rho = CR_RHO,
                  grm_basis = GRM_BASIS, grm_method = GRM_METHOD,
                  seed_bundle = SEEDS[["bundle"]], seed_clusters = SEEDS[["clusters"]])
), OUT)

write_receipt(STAGE, inputs = INPUTS, params = PARAMS, outputs = OUT, target = combo_id)
say("\n[6] wrote %s (%.1f MB)\n    receipt: %s\n", OUT, file.size(OUT) / 1e6, receipt_path(STAGE, combo_id))
say("\n    Next: a scan stage (EMMAX, matching module_3sp/R/03_EMMAX.R) -- not written yet.\n")
