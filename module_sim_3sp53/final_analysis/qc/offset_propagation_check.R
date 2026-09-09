## final_analysis/qc/offset_propagation_check.R
##
## Does the locus-index fix require rerunning Stage-1/GRM/EMMAX, or only
## remapping truth (Chr/Pos/type/allelic_values) onto already-computed
## p-values? PK's question, tested directly on bgs/V0.5_c2/rep1/env1.
##
## ANSWER (2026-09-09): NO -- the full rerun is required, not just a
## truth-remap. Result (see qc/offset_propagation_check.log for the full
## run): GRM differs (max abs diff 1.15e-02; 6534 vs 6569 markers selected
## for kinship); emmax_snp's p-value differs for ALL 17,713 markers common
## to both parses (max abs diff 3.15e-02), not just markers near a QTN.
##
## Mechanism: genotypes ARE unaffected (confirmed -- this is not a
## re-litigation of that), but the bug still reaches association results
## through Stage-1 clustering. A QTN's position shifts by a whole QTN-
## spacing (~100kb+) under the bug, which moves its SORT-ORDER position
## relative to potentially hundreds of neighbouring markers (map is sorted
## by Chr,Pos before anything downstream runs). ld_complexity_reduction()'s
## LD-decay windowing uses that sort order, so the shift reshuffles cluster
## membership for a broad neighbourhood around each QTN, not just the QTN
## itself. GRM_BASIS="stage1_pruned" means different clustering -> different
## representative markers selected for kinship -> different GRM. Because
## EMMAX's kinship correction is a GLOBAL adjustment applied to every test,
## a changed GRM perturbs every marker's p-value, including markers on the
## chromosome with no QTN at all.
##
## Practical consequence: any combination already parsed/analysed under the
## OLD offset (there are none in final_analysis/ -- Phase 1 only ever wrote
## corrected bundles -- but this matters for module_sim_3sp53's OWN archived
## results, which used the old parser throughout) cannot be salvaged by
## re-labelling truth alone. Stage 1 through the association stage all need
## rerunning against the corrected map.
##
## Builds a manual OLD-OFFSET bundle (bypassing parse_nemo_run's strict
## assertions -- deliberately reproducing known-buggy behaviour for
## comparison, not claiming it's correct) with the SAME genotype columns as
## the corrected bundle, then runs Stage 1 + GRM + emmax_snp on BOTH and
## compares: (1) GRM, (2) per-marker p-values matched by underlying NEMO
## marker identity (not by Chr:Pos label, which differs), (3) Stage-1 unit
## membership for markers near a QTN.
suppressMessages({library(data.table); library(LDscnR); library(SNPRelate)})
source("R/01_parse_nemo.R")   ## also sources 00_config.R

TAG <- "bgs"; CELL <- "V0.5_c2"; REP <- 1L; ENVN <- 1L

## ---- corrected bundle (already built for gate 2) -----------------------------
corrected <- readRDS(file.path(PATHS$parsed, sprintf("nemo_%s_rep%d_%s_env%d.rds", TAG, REP, CELL, ENVN)))
## corrected$map does NOT carry nemo_marker (parse_nemo_run's final select
## drops it) -- recover it by re-deriving the corrected join directly (same
## .join_nemo_to_refmap(offset=0) the real parser uses) and matching on the
## Chr:Pos string, which IS what corrected$map$marker was built from.
.recover_nemo_marker <- function(map_dt, nemo_map_local, refmap0_local, offset) {
  j <- .join_nemo_to_refmap(nemo_map_local, refmap0_local, offset)$refmap
  j[, mk := paste0("Chr", Chr, ":", bp)]
  j$nemo_marker[match(map_dt$marker, j$mk)]
}

## ---- manual OLD-OFFSET bundle, same shape, bypassing parse_nemo_run's ---------
## strict per-type range assertion (delet hits it here; irrelevant to this
## comparison, which only cares about ntrl/quant-driven GRM/EMMAX behaviour).
bf <- file.path(PATHS$raw_nemo, sprintf("adapt_%s_chr%d_%s_env%d", TAG, REP, CELL, ENVN), "GENO")
files <- list.files(bf, full.names = TRUE)
map_nemo <- fread(files[grepl("\\.map$", files)])
GTs_raw <- fread(files[grepl("snp_geno", files)])
nemo_map <- data.table(marker = map_nemo$trait.locus, do.call(rbind, strsplit(map_nemo$trait.locus, ".", fixed = TRUE)))
setnames(nemo_map, c("V1", "V2"), c("type", "idx_raw")); nemo_map[, idx_raw := as.numeric(idx_raw)]
sample_info <- GTs_raw[, .(pop, ID)]
GTs_raw_mat <- as.matrix(GTs_raw[, 6:ncol(GTs_raw), with = FALSE])
refmap0 <- readRDS(file.path(PATHS$raw_recmap_dir, sprintf("rec_map%d.rds", REP)))

j1 <- .join_nemo_to_refmap(nemo_map, refmap0, 1L)$refmap   ## old offset
map_old <- j1[nemo_marker %in% colnames(GTs_raw_mat)]
GTs_old <- GTs_raw_mat[, map_old$nemo_marker]
map_old[, Pos := bp]; map_old[, Chr := paste0("Chr", Chr)]
map_old[, marker := paste0(Chr, ":", Pos)]
colnames(GTs_old) <- map_old$marker
ord <- order(map_old$Chr, map_old$bp); map_old <- map_old[ord]; GTs_old <- GTs_old[, ord]
dup <- duplicated(map_old$marker)
if (any(dup)) { map_old <- map_old[!dup]; GTs_old <- GTs_old[, !dup] }
keep_inds <- seq(1, nrow(GTs_old), by = SUBSAMPLE_STEP)
GTs_old <- GTs_old[keep_inds, ]
maf <- colSums(GTs_old) / nrow(GTs_old) / 2
map_old[, MAF := pmin(maf, 1 - maf)]
keep_snps <- map_old$MAF > MAF_KEEP
map_old <- map_old[keep_snps]; GTs_old <- GTs_old[, keep_snps]
map_old <- map_old[, .(Chr, Pos, marker, type, allelic_values, MAF, nemo_marker)]

say("corrected: %d markers (%d QTN) ; old-offset: %d markers (%d QTN)\n",
    nrow(corrected$map), sum(corrected$map$true_QTN), nrow(map_old), sum(map_old$type == "QTN"))

## ---- run Stage 1 + GRM + emmax_snp on BOTH ------------------------------------
run_stage12 <- function(GTs, map, env_dt, label) {
  say("\n--- %s ---\n", label)
  gds_path <- file.path(PATHS$out, "gds_tmp", paste0("cmp_", label, ".gds"))
  dir.create(dirname(gds_path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(gds_path)) unlink(gds_path)
  gds <- create_gds_from_geno(geno = GTs, map = map, gds_path)
  on.exit({ try(SNPRelate::snpgdsClose(gds), silent = TRUE); unlink(gds_path) }, add = TRUE)
  set.seed(SEEDS[["bundle"]])
  LD_decay <- do.call(compute_LD_decay, c(list(gds = gds, ld_w_rho = RHO_GRID, seed = SEEDS[["bundle"]]), DECAY_ARGS))
  set.seed(SEEDS[["clusters"]])
  stage1 <- ld_complexity_reduction(map = map, LD_decay = LD_decay, rho = CR_RHO, gds = gds)
  grm_markers <- unique(na.omit(stage1$pruned))
  GRM <- SNPRelate::snpgdsGRM(gds, snp.id = grm_markers, method = GRM_METHOD, verbose = FALSE, autosome.only = FALSE)$grm
  Pm <- emmax_setup(GTs, GRM)
  pm_obs <- emmax_fast(Pm, env_dt$env)
  list(stage1 = stage1, grm_markers = grm_markers, GRM = GRM, p = stats::setNames(pm_obs, map$marker), map = map)
}

res_new <- run_stage12(corrected$GTs, corrected$map, corrected$env, "corrected")
res_old <- run_stage12(GTs_old, map_old, corrected$env, "old_offset")

say("\n=== COMPARISON ===\n")
say("GRM identical: %s (max abs diff %.2e)\n", isTRUE(all.equal(res_new$GRM, res_old$GRM)),
    max(abs(res_new$GRM - res_old$GRM)))
say("n grm_markers: corrected=%d, old=%d, identical SET (by nemo_marker): %s\n",
    length(res_new$grm_markers), length(res_old$grm_markers),
    setequal(corrected$map[match(res_new$grm_markers, corrected$map$marker)]$nemo_marker,
             map_old[match(res_old$grm_markers, map_old$marker)]$nemo_marker))

## match p-values by underlying nemo_marker identity, not by (offset-dependent) Chr:Pos label
corrected_nemo <- .recover_nemo_marker(corrected$map, nemo_map, refmap0, 0L)
p_new_by_nemo <- stats::setNames(res_new$p, corrected_nemo[match(names(res_new$p), corrected$map$marker)])
p_old_by_nemo <- stats::setNames(res_old$p, map_old$nemo_marker[match(names(res_old$p), map_old$marker)])
common_nemo <- intersect(names(p_new_by_nemo), names(p_old_by_nemo))
say("markers in both (by nemo identity): %d of %d (new) / %d (old)\n", length(common_nemo), length(p_new_by_nemo), length(p_old_by_nemo))
diffs <- abs(p_new_by_nemo[common_nemo] - p_old_by_nemo[common_nemo])
say("emmax_snp p-value identical for common markers: %s (max abs diff %.2e, %d differ by >1e-10)\n",
    all(diffs < 1e-10, na.rm = TRUE), max(diffs, na.rm = TRUE), sum(diffs > 1e-10, na.rm = TRUE))

## Stage-1 unit membership near QTN: for each QTN (by nemo identity), which
## unit does it belong to in each parse, and do the OTHER members match (by
## nemo identity)?
corrected_qtn_nemo <- corrected_nemo[corrected$map$true_QTN]
qtn_nemo <- intersect(corrected_qtn_nemo, map_old[type == "QTN", nemo_marker])
say("\nQTN present as QTN-type in BOTH parses (by nemo identity): %d\n", length(qtn_nemo))
cl_new <- as.data.table(res_new$stage1$clusters); cl_old <- as.data.table(res_old$stage1$clusters)
for (qn in qtn_nemo) {
  new_marker <- corrected$map$marker[corrected_nemo == qn][1]
  old_marker <- map_old[nemo_marker == qn, marker]
  unit_new <- cl_new[sapply(members, function(m) new_marker %in% m)]
  unit_old <- cl_old[sapply(members, function(m) old_marker %in% m)]
  members_new_nemo <- if (nrow(unit_new)) sort(corrected_nemo[match(unlist(unit_new$members[1]), corrected$map$marker)]) else character(0)
  members_old_nemo <- if (nrow(unit_old)) sort(map_old$nemo_marker[match(unlist(unit_old$members[1]), map_old$marker)]) else character(0)
  say("  %s: unit size new=%d old=%d, SAME MEMBER SET (by nemo identity): %s\n",
      qn, length(members_new_nemo), length(members_old_nemo), setequal(members_new_nemo, members_old_nemo))
}
