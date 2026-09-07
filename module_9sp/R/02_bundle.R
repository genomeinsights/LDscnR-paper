## =============================================================================
## module_9sp/R/02_bundle.R
##
## BUILD THE DATA BUNDLE: genotypes, map, LD-decay, ld_w, stage-1 clustering, kinship.
## Same structure and ordering as module_3sp/R/02_bundle.R -- copied, not
## re-invented: decay -> ld_w -> stage-1 clustering -> representatives -> GRM,
## same reason (the kinship basis IS the stage-1 representatives, so clustering
## must run before the kinship).
##
## RAW FORMAT DIFFERS FROM 3sp: 9sp_data.rds is a plain list (readRDS), not an
## .RData environment -- $GT/$map/$pheno/$ecotype_bin/$GRM. The bundled $GRM is
## NOT used here (its basis is uncharacterised, see 00_config.R); this stage
## builds its own from the stage-1 representatives, same as 3sp.
##
## SIZE_FLOOR is NA in 00_config.R -- deliberately, since 3sp's value (8) was
## derived FROM 3sp's own Stage-1 median cluster size and has no basis for this
## panel yet. This stage reports the median cluster size but does not gate
## anything on SIZE_FLOOR, unlike 3sp's version.
##
## COST: hours. The decay fit over ~1.36M markers, 20 chromosomes, dominates.
## Writes nothing outside out/02_bundle/ and module_9sp/cache/.
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(SNPRelate); library(digest)})
devtools::load_all("~/gitlab/LDscnR")
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_9sp"), "R", "00_config.R"))
STAGE <- "02_bundle"
invisible(check_ldscnr())

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

PATHS$el_dir <- file.path(PATHS$cache, "edge_lists")
INPUTS <- c(PATHS$raw_9sp, PATHS$provenance)
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
d <- readRDS(PATHS$raw_9sp)
GTs <- d$GT
map <- as.data.table(d$map)
keep <- map$maf > MAF_KEEP
GTs  <- GTs[, keep]; map <- map[keep]
colnames(GTs) <- map$marker     # consensus_dosage and make_eMLGs index by name
eco   <- as.integer(d$ecotype_bin)   # already 0/1: Freshwater 0, Marine 1 (verified against $pheno$ecotype)
pheno <- as.data.table(d$pheno)
## GTs, eco and pheno are aligned by ROW POSITION ONLY (raw file carries no sample IDs,
## same as 3sp) -- same identity caveat, same explicit count check rather than trusting
## silently.
stopifnot(nrow(GTs) == length(eco), nrow(GTs) == nrow(pheno))
n <- nrow(GTs)
say("    %d individuals x %s markers (maf > %.2f, %s dropped) ; Marine = %d\n",
    n, format(ncol(GTs), big.mark=","), MAF_KEEP, format(sum(!keep), big.mark=","), sum(eco))
say("    %d populations, %d localities, %d lineages\n",
    uniqueN(pheno$pop_ID), uniqueN(pheno$pop_locality), uniqueN(pheno$lineage))

## ---- 2. GDS ------------------------------------------------------------------
gds_path <- file.path(PATHS$cache, "9sp.gds")
dir.create(dirname(gds_path), recursive = TRUE, showWarnings = FALSE)
if (file.exists(gds_path)) { say("\n[2] removing stale %s\n", basename(gds_path))
                             unlink(gds_path) }
say("\n[2] GDS -> %s\n", gds_path)
gds <- create_gds_from_geno(geno = GTs, map = map, gds_path)
on.exit(try(snpgdsClose(gds), silent = TRUE), add = TRUE)
stopifnot(file.exists(gds_path), file.size(gds_path) > 1e6)
say("    %.0f MB\n", file.size(gds_path)/1e6)

## ---- 3. LD decay, SEEDED, and ld_w in place ---------------------------------
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
stage1_fp <- digest(list(decay_fp = decay_fp, cr_rho = CR_RHO), algo = "sha256")
stage1 <- .cache_step("stage1", stage1_fp, function() {
  set.seed(SEEDS[["clusters"]])
  ld_complexity_reduction(map = map, LD_decay = LD_decay, rho = CR_RHO, gds = gds) })
cl <- as.data.table(stage1$clusters)
nl <- if ("n_loci" %in% names(cl)) cl$n_loci else lengths(cl$members)
## SIZE_FLOOR is NA -- this is the measurement that will set it (2x median, matching 3sp's
## own derivation), not a gate applied here. Reported plainly so it can be read off and set
## in 00_config.R once, not re-derived ad hoc by a later stage.
say("    %s clusters ; median size %.2f (2x median = %.0f, for SIZE_FLOOR)\n",
    format(nrow(cl), big.mark=","), median(nl), 2*median(nl))
say("    size distribution: singletons %s (%.1f%%) ; >=2 markers %s\n",
    format(sum(nl==1), big.mark=","), 100*mean(nl==1), format(sum(nl>=2), big.mark=","))

## ---- 5. the kinship basis IS the stage-1 representatives --------------------
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
stopifnot(nrow(GRM) == nrow(GTs), ncol(GRM) == nrow(GTs))
ut <- upper.tri(GRM)
say("    %d x %d ; mean diagonal %.4f ; off-diagonal mean %+.4f sd %.4f ; %.1f min\n",
    nrow(GRM), ncol(GRM), mean(diag(GRM)), mean(GRM[ut]), sd(GRM[ut]),
    as.numeric(difftime(Sys.time(), t0, units="mins")))

## ---- 6. save -----------------------------------------------------------------
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
say("\n    Next: derive SIZE_FLOOR from the median cluster size above, set it in\n")
say("    00_config.R, then 03_EMMAX.R.\n")
