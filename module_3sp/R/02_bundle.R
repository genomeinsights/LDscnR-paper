## =============================================================================
## module_3sp/R/02_bundle.R
##
## BUILD THE DATA BUNDLE: genotypes, map, LD-decay, ld_w, stage-1 clustering, kinship.
## Runs no association test -- the scans are stage 05, because they depend on the kinship
## this stage produces and conflating the two is why "rebuild the kinship" and "rerun the
## scan" could not previously be decided independently.
##
## THE ORDERING IS THE POINT OF THIS STAGE. The kinship basis is the STAGE-1
## REPRESENTATIVES, so the clustering must run BEFORE the kinship:
##
##     decay  ->  ld_w  ->  stage-1 clustering  ->  representatives  ->  GRM
##
## The superseded pipeline computed the clustering here and then did not use it, taking
## the basis from an `ld_w` threshold instead -- a second arbitrary number serving only
## the kinship. Stage-1 pruning does BOTH JOBS AT ONCE (PK): the clusters are the test
## units and their representatives are the kinship markers, so there is one decision
## where there were two. Measured cost of the switch: 57 discoveries against 59, with the
## kinships correlating r = 0.9994. Nothing is being traded away.
##
## WHY THE RESULT MAY NOT MATCH THE PUBLISHED NUMBERS, and this is expected rather than a
## defect. Two things change deliberately:
##   1. n_win_decay is 20, not the 5 the old bundle used (PK).
##   2. compute_LD_decay is SEEDED here. The old fit predates seeding and subsamples in
##      three places, so its ld_w -- which set both the GRM basis and the clustering --
##      cannot be reproduced even in principle.
## Both propagate to every downstream number. A small shift is the honest consequence of
## making the root reproducible.
##
## COST: hours. The decay fit over the panel's 20 chromosomes dominates (Chr19 is
## absent from this panel entirely -- not dropped, never present; verified against the raw
## map rather than assumed). Writes nothing outside
## out/02_bundle/ and module_3sp/cache/.
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(SNPRelate); library(digest)})
devtools::load_all("~/gitlab/LDscnR")
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_3sp"), "R", "00_config.R"))
STAGE <- "02_bundle"

## ---- CHECKPOINTS for the two expensive steps --------------------------------
## The decay fit and the clustering are the steps this stage can die partway through --
## hours for the first, minutes for the second -- and until now a crash meant starting
## the whole bundle over. Each gets its own cache file in module_3sp/cache/, keyed by a
## fingerprint of everything it depends on, so a rerun after a crash (or a `source()`
## after fixing an unrelated bug downstream) reuses what already finished rather than
## recomputing it.
##
## THE FINGERPRINT IS THE CORRECTNESS GUARANTEE, not the file's mtime or presence.
## stage1's fingerprint folds in decay's fingerprint, so a change to ANY decay parameter
## invalidates the clustering too, transitively, without stage1's own code needing to
## know what decay depends on.
##
## FORCE=1 busts both caches, same as it busts the stage-level receipt above.
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

INPUTS <- c(PATHS$raw_3sp, PATHS$provenance)
PARAMS <- list(maf_keep = MAF_KEEP, decay_args = DECAY_ARGS, rho_grid = RHO_GRID,
               cr_rho = CR_RHO, grm_basis = GRM_BASIS, grm_method = GRM_METHOD,
               seed = SEEDS[["bundle"]])

say("=== %s ===\n\n", STAGE)
if (!stage_stale(STAGE, INPUTS, PARAMS) && !nzchar(Sys.getenv("FORCE"))) {
  say("\nNothing to do. Set FORCE=1 to rebuild anyway.\n"); quit(save = "no")
}
dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)
dir.create(PATHS$el_dir,     recursive = TRUE, showWarnings = FALSE)
t_all <- Sys.time()

## ---- 1. raw genotypes, MAF filter -------------------------------------------
say("[1] raw genotypes\n")
e <- new.env(); load(PATHS$raw_3sp, envir = e)
GTs <- e$GTs_3sp
colnames(GTs) <- e$map_3sp$marker          # consensus_dosage and make_eMLGs index by name
map <- as.data.table(e$map_3sp)
keep <- map$maf > MAF_KEEP
GTs  <- GTs[, keep]; map <- map[keep]
eco  <- as.integer(as.factor(e$pheno_3sp$ecotype)) - 1L   # Freshwater 0, Marine 1
pheno <- as.data.table(e$pheno_3sp)
n <- nrow(GTs)
say("    %d individuals x %s markers (maf > %.2f, %s dropped) ; Marine = %d\n",
    n, format(ncol(GTs), big.mark=","), MAF_KEEP, format(sum(!keep), big.mark=","), sum(eco))
say("    %d populations in %d localities\n", uniqueN(pheno$pop_ID), uniqueN(pheno$pop_locality))

## ---- 2. GDS ------------------------------------------------------------------
## Built here, in module_3sp/cache/, and nowhere else. The old 3sp_data/3sp.gds is a
## different file belonging to the superseded pipeline; this stage neither reads it nor
## touches it. Any stale copy at OUR path is removed first, so a rerun cannot silently
## fit the decay against a GDS built from a different MAF filter.
gds_path <- file.path(PATHS$cache, "3sp.gds")
dir.create(dirname(gds_path), recursive = TRUE, showWarnings = FALSE)
if (file.exists(gds_path)) { say("\n[2] removing stale %s\n", basename(gds_path))
                             unlink(gds_path) }
say("\n[2] GDS -> %s\n", gds_path)
## create_gds_from_geno RETURNS AN OPEN HANDLE (it ends in snpgdsOpen). Closing it is not
## housekeeping: while it is open, a later rerun's unlink() above can fail and the stage
## would then fit against the old file. Registered with on.exit rather than closed at the
## end, because the decay fit is the step most likely to die hours in -- and that is
## exactly when the handle would otherwise be left dangling.
gds <- create_gds_from_geno(geno = GTs, map = map, gds_path)
on.exit(try(snpgdsClose(gds), silent = TRUE), add = TRUE)
stopifnot(file.exists(gds_path), file.size(gds_path) > 1e6)
say("    %.0f MB\n", file.size(gds_path)/1e6)

## ---- 3. LD decay, SEEDED, and ld_w in place ---------------------------------
## SEEDED VIA THE FUNCTION'S OWN `seed` ARGUMENT. R_3sp/119 had to seed in the caller
## because the package installed then had no `seed` formal; the installed version does, and
## an explicit argument beats a set.seed() whose scope a reader has to infer. The
## caller-side set.seed stays as a belt-and-braces for any draw taken outside the seeded
## region -- it costs nothing and removes the need to verify that claim.
##
## This is the step whose absence made the old bundle unreproducible: compute_LD_decay
## subsamples the background, thins per chromosome, and samples pairs within strata, and
## the old fit fixed none of it.
stopifnot("seed" %in% names(formals(compute_LD_decay)))
say("\n[3] LD decay: n_win_decay = %d, seed %d\n", DECAY_ARGS$n_win_decay, SEEDS[["bundle"]])
decay_fp <- digest(list(decay_args = DECAY_ARGS, rho_grid = RHO_GRID, maf_keep = MAF_KEEP,
                        seed = SEEDS[["bundle"]], markers = ncol(GTs)), algo = "sha256")
LD_decay <- .cache_step("ld_decay", decay_fp, function() {
  set.seed(SEEDS[["bundle"]])
  do.call(compute_LD_decay,
         c(list(gds = gds, el_data_folder = PATHS$el_dir, ld_w_rho = RHO_GRID,
                seed = SEEDS[["bundle"]]),
           DECAY_ARGS)) })
ld_ws <- LD_decay$ld_ws[map$marker, , drop = FALSE]
ld95  <- if ("rho_0.95" %in% colnames(ld_ws)) "rho_0.95" else "0.95"
map[, ld_w_095 := ld_ws[, ld95]]
say("    %d chromosomes ; ld_w matrix %s x %d\n", nrow(LD_decay$decay_sum),
    format(nrow(ld_ws), big.mark=","), ncol(ld_ws))

## ---- 4. stage-1 clustering -- BEFORE the kinship ----------------------------
say("\n[4] stage-1 clustering (ld_complexity_reduction, rho = %.2f)\n", CR_RHO)
## Folds decay_fp IN rather than recomputing an independent one: if decay's inputs change,
## stage1's fingerprint changes too, even though stage1's own code never touches DECAY_ARGS.
stage1_fp <- digest(list(decay_fp = decay_fp, cr_rho = CR_RHO), algo = "sha256")
stage1 <- .cache_step("stage1", stage1_fp, function() {
  set.seed(SEEDS[["clusters"]])
  ## gds = gds IS REQUIRED HERE, not optional. compute_LD_decay() was called with
  ## el_data_folder = PATHS$el_dir, which writes edge lists to disk rather than
  ## keeping them in the returned object -- by_chr[[ch]]$el is NULL on the object
  ## itself. ld_complexity_reduction() needs edges for the merge step, so without
  ## gds it has no way to get them and errors: "LD_decay$by_chr[['Chr1']]$el is
  ## NULL -- either run compute_LD_decay() with keep_el = TRUE ... or pass gds =".
  ld_complexity_reduction(map = map, LD_decay = LD_decay, rho = CR_RHO, gds = gds) })
cl <- as.data.table(stage1$clusters)
nl <- if ("n_loci" %in% names(cl)) cl$n_loci else lengths(cl$members)
say("    %s clusters ; %s with >= %d markers ; median size %.2f\n",
    format(nrow(cl), big.mark=","), format(sum(nl >= SIZE_FLOOR), big.mark=","),
    SIZE_FLOOR, median(nl))

## ---- 5. the kinship basis IS the stage-1 representatives --------------------
## stage1$pruned is one core SNP per cluster and is identical to clusters$core_snp
## (verified). ld_prune_and_eMLG()$pruned is a DIFFERENT set -- stage-2 group
## representatives -- and is deliberately not used here.
say("\n[5] kinship basis: %s\n", GRM_BASIS)
grm_markers <- switch(GRM_BASIS,
  stage1_pruned = unique(na.omit(stage1$pruned)),
  greedy = { set.seed(SEEDS[["bundle"]])
    unlist(snpgdsLDpruning(gds, ld.threshold = GRM_GREEDY$ld.threshold,
                           slide.max.bp = GRM_GREEDY$slide.max.bp,
                           autosome.only = FALSE, verbose = FALSE), use.names = FALSE) },
  none = map$marker,
  stop("unknown GRM_BASIS: ", GRM_BASIS))
say("    %s of %s markers (%.1f%%)\n", format(length(grm_markers), big.mark=","),
    format(nrow(map), big.mark=","), 100*length(grm_markers)/nrow(map))

say("    GRM: snpgdsGRM(method = \"%s\")\n", GRM_METHOD)
t0 <- Sys.time()
GRM <- snpgdsGRM(gds, snp.id = grm_markers, method = GRM_METHOD,
                 verbose = FALSE, autosome.only = FALSE)$grm
ut <- upper.tri(GRM)
say("    %d x %d ; mean diagonal %.4f ; off-diagonal mean %+.4f sd %.4f ; %.1f min\n",
    nrow(GRM), ncol(GRM), mean(diag(GRM)), mean(GRM[ut]), sd(GRM[ut]),
    as.numeric(difftime(Sys.time(), t0, units="mins")))

## ---- 6. save -----------------------------------------------------------------
## No emx, no lfmm: those are stage 05. A bundle that carries its own association scan
## cannot be rebuilt without re-deciding the scan, which is the coupling this pipeline
## exists to break.
OUT <- file.path(stage_dir(STAGE), "bundle.rds")
saveRDS(list(
  GTs = GTs, map = map, eco = eco, pheno = pheno,
  ld_ws = ld_ws, LD_decay = LD_decay,
  stage1 = stage1,
  GRM = GRM, grm_markers = grm_markers,
  settings = list(maf_keep = MAF_KEEP, decay_args = DECAY_ARGS, rho_grid = RHO_GRID,
                  cr_rho = CR_RHO, grm_basis = GRM_BASIS, grm_method = GRM_METHOD,
                  seed_bundle = SEEDS[["bundle"]], seed_clusters = SEEDS[["clusters"]])
), OUT)
write_receipt(STAGE, inputs = INPUTS, params = PARAMS, outputs = OUT)
say("\n[6] wrote %s (%.0f MB) in %.1f min total\n", OUT, file.size(OUT)/1e6,
    as.numeric(difftime(Sys.time(), t_all, units="mins")))
say("    receipt: %s\n", receipt_path(STAGE))
say("\n    Next: 03_clusters.R (stage-2 grouping). The scans are 05, not here.\n")
