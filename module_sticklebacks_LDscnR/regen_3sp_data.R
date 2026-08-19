## =====================================================================
## module_sticklebacks_LDscnR / regen_3sp_data.R
##
## Regenerate the 3sp input data FROM the raw parsed genotypes, using the current
## LDscnR machinery (compute_LD_decay / compute_ld_w / emmax), replacing the
## outdated ld_decay() machinery in the earlier exploratory 3sp code.
##
## GRM basis: LD-INDEPENDENT markers (ld_w_095 < GRM_LDW_THRESHOLD). We first tried
## the sim complexity-reduction chain (ld_complexity_reduction -> ld_prune_and_eMLG,
## as in module_sim) for method-identity with the sims, but on 3sp it keeps ~96% of
## markers (97.8% singleton clusters) and OVER-corrects -> the C-score collapses. A
## direct local-LD threshold restores a moderate independent set; see the
## GRM_LDW_THRESHOLD note below. sim<->3sp consistency is in the principle, not the
## selector.
##
## EMMAX engine note: the observed C-score and the structure-aware null both run
## through emmax_fast() on ONE shared GRM (the observed phenotype is scored
## alongside the surrogate "permuted" phenotypes by the same fast engine). So the
## kinship is what must be saved and kept identical -- which this script does.
## The single-SNP emmax() result saved here (F-test p + genomic control) is the
## conventional-association track for the -log10(q) Manhattan; LFMM is unchanged
## and reused from the previous run.
##
## Outputs one self-contained bundle to module_sticklebacks_LDscnR/data/ that the
## analysis driver (run_3sp_LDscnR.R) consumes.
##
## Run from the LDscnR-paper root (heavy; review before running):
##   Rscript module_sticklebacks_LDscnR/regen_3sp_data.R
## =====================================================================

suppressMessages({
  library(data.table); library(SNPRelate)
  devtools::load_all("/Users/petrikem/gitlab/LDscnR", quiet = TRUE)   # or library(LDscnR)
})

## ---- paths -----------------------------------------------------------
RAW_DATA  <- "~/gitlab/LD-scaling-genome-scans/empirical_data/3sp/3sp_data.RData"  # GTs_3sp, map_3sp, pheno_3sp
LFMM_F    <- "~/gitlab/LD-scaling-genome-scans/empirical_data/3sp/lfmm_F.rds"       # old LFMM F (full map)
OUT_DIR   <- "module_sticklebacks_LDscnR/data"

OUT_FILE  <- file.path(OUT_DIR, "3sp_LDscnR_data.rds")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)


## ---- canonical settings (identical to module_sim/R/regen_stats.R) ----
RHO_GRID   <- c(seq(0.05, 0.95, by = 0.05), 0.99)          # ld_w columns
## UNCAPPED decay (max_SNPs_decay = Inf, n_win_decay = 5, as in the sims): the
## earlier caps (10000 / 20) cut C-score sensitivity ~4x and lost Chr4-Eda + Chr7.
## ld_w is computed IN PLACE (ld_w_rho) from the edge lists built here, which are
## then dropped (keep_el = FALSE, no el_data_folder) -- no large edge lists saved.
DECAY_ARGS <- list(min_maf_decay = 0.1, q = 0.95, n_sub_bg = 5000, n_win_decay = 5,
                   overlap = 0.5, max_SNPs_decay = Inf, prob_robust = 0.95,
                   max_pairs = 5000, ld_method = "corr", n_strata = 20, keep_el = FALSE,
                   slide = 1000, rho_targets = c(0.99), cores = 1, ld_w_rho = RHO_GRID)
CR_RHO     <- 0.5                                           # ld_complexity_reduction rho (stage1, reference only)
## GRM basis: LD-INDEPENDENT markers via a direct local-LD threshold, NOT the sim
## complexity-reduction chain. Checked on 3sp: the chain keeps ~96% of markers
## (97.8% singleton clusters at rho=0.5) and OVER-corrects -> the C-score collapses
## (72 SNPs C>0, only Chr17 survives). ld_w_095 < 0.1 keeps a moderate independent
## set (~24%), gives a well-calibrated GRM (gif ~1.10, no GC) and recovers the known
## regions (Chr1, Chr20, Chr17, Chr5). The sim<->3sp consistency is in the PRINCIPLE
## (independent-marker GRM), not the selector -- the same chain behaves oppositely
## on the two datasets' LD structures.
GRM_LDW_THRESHOLD <- 0.1

## ---- 1. load raw genotypes + MAF>0.1 filter --------------------------
e <- new.env(); load(path.expand(RAW_DATA), envir = e)
GTs <- e$GTs_3sp
colnames(GTs) <- e$map_3sp$marker                          # consensus_dosage indexes GTs by marker name
map <- data.table::as.data.table(e$map_3sp)                # Chr, Pos, marker, maf
keep <- map$maf > 0.1                                       # canonical MAF filter
GTs <- GTs[, keep]; map <- map[keep]
eco <- as.integer(as.factor(e$pheno_3sp$ecotype)) - 1L      # Freshwater=0, Marine=1
n <- nrow(GTs)
cat(sprintf("[1] %d individuals x %d markers (MAF>0.1) ; Marine = %d\n", n, ncol(GTs), sum(eco)))

gds_path <- file.path(OUT_DIR, "3sp.gds")
if (file.exists(gds_path)) unlink(gds_path)
gds <- create_gds_from_geno(geno = GTs, map = map, gds_path)

## ---- 2. LD-decay (current machinery) + per-SNP ld_w ------------------
cat("[2] estimating LD-decay (compute_LD_decay) ...\n")
LD_decay <- do.call(compute_LD_decay, c(list(gds = gds), DECAY_ARGS))   # ld_w computed in place
ld_ws <- LD_decay$ld_ws[map$marker, , drop = FALSE]                     # rows named by marker
ld95_col <- if ("rho_0.95" %in% colnames(ld_ws)) "rho_0.95" else "0.95"
map[, ld_w_095 := ld_ws[, ld95_col]]                       # the local-LD support used for the GRM basis
cat(sprintf("[2] LD-decay over %d chromosomes ; ld_w matrix %d x %d\n",
            nrow(LD_decay$decay_sum), nrow(ld_ws), ncol(ld_ws)))

## ---- 3. GRM from LD-independent markers (the shared kinship) ----------
## stage1 (complexity reduction) is still computed + saved as the LD-complexity
## reduction object, but is NOT the GRM basis (see GRM_LDW_THRESHOLD note above).
cat("[3] complexity reduction (reference) + LD-independent GRM ...\n")
stage1 <- ld_complexity_reduction(map = map, LD_decay = LD_decay, rho = CR_RHO)
grm_markers <- map$marker[which(map$ld_w_095 < GRM_LDW_THRESHOLD)]
cat(sprintf("[3] GRM markers (ld_w_095 < %.2f): %d / %d (%.1f%%)\n",
            GRM_LDW_THRESHOLD, length(grm_markers), nrow(map), 100 * length(grm_markers) / nrow(map)))
GRM <- snpgdsGRM(gds, snp.id = grm_markers, method = "GCTA",
                 verbose = FALSE, autosome.only = FALSE)$grm

## ---- 4. single-SNP EMMAX (F-test p + genomic control) ----------------
## conventional-association track for the Manhattan; the C-score + structured
## null re-run the SAME GRM through emmax_fast downstream (run_3sp_LDscnR.R).
cat("[4] EMMAX scan ...\n")
emx <- emmax(eco, GTs, K = GRM)
map[, emx_F := emx$F][, emx_p := emx$pval]
gif <- stats::median(map$emx_F) / stats::qf(0.5, 1, n - 2, lower.tail = FALSE)
if (gif > 1.1) {                                           # genomic control
  map[, emx_F := emx_F / gif]
  map[, emx_p := stats::pf(emx_F, 1, n - 2, lower.tail = FALSE)]
}
map[, emx_q := stats::p.adjust(emx_p, "fdr")]
cat(sprintf("[4] EMMAX gif = %.3f ; min emx_q = %.3f ; emx_q<=0.05 : %d\n",
            gif, min(map$emx_q), sum(map$emx_q <= 0.05)))

## ---- 5. LFMM (unchanged; reuse previous F-values) --------------------
lfmm_F_full <- readRDS(path.expand(LFMM_F))                # length = full (pre-MAF) map
stopifnot(length(lfmm_F_full) == length(keep))
map[, lfmm_F := lfmm_F_full[keep]]
map[, lfmm_p := stats::pf(lfmm_F, 1, n - 2, lower.tail = FALSE)]
map[, lfmm_q := stats::p.adjust(lfmm_p, "fdr")]
cat(sprintf("[5] LFMM reused ; min lfmm_q = %.3f ; lfmm_q<=0.05 : %d\n",
            min(map$lfmm_q), sum(map$lfmm_q <= 0.05)))

## ---- 6. save one self-contained bundle -------------------------------
snpgdsClose(gds)
saveRDS(list(
  GTs   = GTs,                     # individuals x SNPs (MAF>0.1), cols = markers
  map   = map,                     # Chr, Pos, marker, maf, ld_w_095, emx_*, lfmm_*
  eco   = eco,                     # phenotype (Marine = 1)
  ld_ws = ld_ws,                   # SNP x rho local-LD matrix (rows = markers)
  LD_decay = LD_decay,             # compute_LD_decay object (decay_sum + ld_ws; edge lists dropped)
  complexity_reduction = list(     # the LD-complexity-reduction object (reference; not the GRM basis)
    stage1 = stage1),
  GRM   = GRM,                     # shared kinship (from the LD-independent markers)
  grm_markers = grm_markers,       # the ld_w_095 < threshold marker set used for the GRM
  grm_ld_w_threshold = GRM_LDW_THRESHOLD,
  emx   = emx,                     # raw EMMAX (F, pval, Rsq)
  emx_gif = gif,
  settings = list(rho_grid = RHO_GRID, decay_args = DECAY_ARGS, cr_rho = CR_RHO,
                  grm_ld_w_threshold = GRM_LDW_THRESHOLD, maf_keep = keep)
), OUT_FILE)
cat(sprintf("[6] wrote %s\n", OUT_FILE))
