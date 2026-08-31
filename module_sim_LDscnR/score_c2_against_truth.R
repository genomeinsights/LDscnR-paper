## =====================================================================
## module_sim_LDscnR / score_c2_against_truth.R
##
## Score the second-tier C-score against truth. C2 is a RANKING, not a test --
## the framework is explicit that it controls no error rate -- so the question is
## whether it orders genuinely causal regions above spurious ones, and whether it
## does so better than the simpler orderings already available.
##
## This matters for empirical use. On real data q_R routinely pins every region
## at the p-floor 1/(1+B) and cannot separate them; C2 is what the framework
## assigns the ranking job to. If C2 ranks well, the operating-point question
## (which tau, which l_min) becomes secondary, because C2 integrates over the
## grid instead of choosing a cell. If it ranks no better than region mass, then
## the extra machinery is not earning its place.
##
## Regions are labelled with classify_ors(), which returns is_TP per region --
## the same focal-QTN assignment and dedup evaluate_ors() uses in aggregate.
##
## Competing rankings, all available without truth:
##   c2     fraction of usable (tau, l_min) cells in which the region is significant
##   s_R    summed C-mass of the region
##   size   number of markers
##   maxC   the region's highest single-marker C
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/score_c2_against_truth.R <scan_dir> [outdir]
## =====================================================================
suppressMessages({ library(data.table); library(LDscnR) })
`%||%` <- function(a, b) if (is.null(a)) b else a

a <- commandArgs(trailingOnly = TRUE)
if (!length(a)) stop("usage: score_c2_against_truth.R <scan_dir> [outdir]")
SCAN_DIR <- a[1]; OUT <- if (length(a) >= 2) a[2] else SCAN_DIR
PANEL_DIR <- Sys.getenv("PANEL_DIR", "/Volumes/Nemo/Nemo_sim/analysis_inputs")
ENGINE <- Sys.getenv("ENGINE", "emmax"); BASIS <- Sys.getenv("BASIS", "env_orth")
RHO_LD <- 0.75; RHO_D <- 0.95; DCAP <- 1e5   # matches the bundles' stage-2 cap
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

## rank-based AUC with proper tie handling (C2 is discrete, so ties are common
## and Wilcoxon/rank AUC is the right form rather than a step-counting ROC)
auc_of <- function(score, label) {
  if (length(unique(label)) < 2) return(NA_real_)
  r <- rank(score, ties.method = "average")
  n1 <- sum(label); n0 <- sum(!label)
  (sum(r[label]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

scans <- list.files(SCAN_DIR, pattern = sprintf("^scan_.*_%s_%s_.*[.]rds$", ENGINE, BASIS), full.names = TRUE)
cat(sprintf("[1] %d genome(s), engine %s, basis %s\n", length(scans), ENGINE, BASIS)); flush.console()

per_region <- list(); per_genome <- list(); k <- 0L; g <- 0L
for (f in scans) {
  cell <- sub(sprintf("^scan_(V[0-9.]+_c[0-9.]+_env[0-9]+)_%s_.*$", ENGINE), "\\1", basename(f))
  pf <- file.path(PANEL_DIR, sprintf("panel_%s.rds", cell)); if (!file.exists(pf)) next
  panel <- readRDS(pf); map <- as.data.table(panel$map)
  if (!"true_pos_QTN" %in% names(map)) map <- flag_true_qtns(map)
  if (!sum(map$true_pos_QTN %in% TRUE)) next
  x <- readRDS(f)
  if (is.null(x$c2) || !nrow(x$c2$regions)) { cat(sprintf("  [skip] %s: no C2\n", cell)); next }
  C <- x$null$C_obs; if (!any(C > 0)) next
  th <- score_thresholds(as.data.table(panel$decay_sum), rho_r2 = RHO_LD, rho_d = RHO_D, dmax_cap = DCAP)

  ## rebuild the anchor regions so membership is available for labelling
  tau <- x$params$tau; L <- x$params$l_min
  r <- ld_regions(names(C)[which(C >= tau)], x$edges); r <- r[lengths(r) >= L]
  if (!length(r)) next
  pos <- stats::setNames(map$Pos, map$marker); chrv <- stats::setNames(as.character(map$Chr), map$marker)
  sp <- data.table(chr = unname(chrv[vapply(r, `[`, character(1), 1L)]),
                   lo = vapply(r, function(z) min(pos[z]), numeric(1)),
                   hi = vapply(r, function(z) max(pos[z]), numeric(1)),
                   size = lengths(r),
                   s_R = vapply(r, function(z) sum(C[z]), numeric(1)),
                   maxC = vapply(r, function(z) max(C[z]), numeric(1)))
  qtab <- qtn_ld_table(panel$GTs, map, names(C)[C > 0], 2e6, cores = 1)
  cls <- classify_ors(r, map, qtab, th$r2min, th$dmax)
  sp[, is_TP := cls$is_TP[match(seq_len(.N), cls$CL_id)]]

  ## attach C2 by span (the anchor cell is the same, so spans match exactly)
  c2t <- as.data.table(x$c2$regions)[, .(chr, lo, hi, c2)]
  sp <- merge(sp, c2t, by = c("chr", "lo", "hi"), all.x = TRUE)
  sp[is.na(c2), c2 := 0]
  sp[, cell := cell]
  k <- k + 1L; per_region[[k]] <- sp

  g <- g + 1L
  per_genome[[g]] <- data.table(cell, n_regions = nrow(sp), n_TP = sum(sp$is_TP),
    auc_c2 = auc_of(sp$c2, sp$is_TP), auc_sR = auc_of(sp$s_R, sp$is_TP),
    auc_size = auc_of(sp$size, sp$is_TP), auc_maxC = auc_of(sp$maxC, sp$is_TP),
    n_usable_cells = x$c2$n_usable, n_cells = x$c2$n_cells,
    distinct_c2 = uniqueN(sp$c2), distinct_qR = uniqueN(x$regions$regions$q_R))
  cat(sprintf("  %-14s %3d regions, %2d TP | AUC c2 %.3f | s_R %.3f | size %.3f | maxC %.3f | distinct c2 %d vs q_R %d\n",
              cell, nrow(sp), sum(sp$is_TP), per_genome[[g]]$auc_c2, per_genome[[g]]$auc_sR,
              per_genome[[g]]$auc_size, per_genome[[g]]$auc_maxC,
              per_genome[[g]]$distinct_c2, per_genome[[g]]$distinct_qR)); flush.console()
}
pg <- rbindlist(per_genome); prg <- rbindlist(per_region)
fwrite(pg, file.path(OUT, "c2_vs_truth_by_genome.csv")); fwrite(prg, file.path(OUT, "c2_vs_truth_regions.csv"))

cat(sprintf("\n=== ranking AUC, replicate-averaged over %d genomes (0.5 = no better than chance) ===\n", nrow(pg)))
print(pg[, .(ranking = c("C2","s_R","size","maxC"),
             mean_AUC = round(c(mean(auc_c2, na.rm=TRUE), mean(auc_sR, na.rm=TRUE),
                                mean(auc_size, na.rm=TRUE), mean(auc_maxC, na.rm=TRUE)), 3),
             SE = round(c(sd(auc_c2,na.rm=TRUE), sd(auc_sR,na.rm=TRUE),
                          sd(auc_size,na.rm=TRUE), sd(auc_maxC,na.rm=TRUE))/sqrt(.N), 3))])
cat(sprintf("\n  C2 beats s_R in %d of %d genomes | paired mean diff %+.3f\n",
            sum(pg$auc_c2 > pg$auc_sR, na.rm=TRUE), nrow(pg), mean(pg$auc_c2 - pg$auc_sR, na.rm=TRUE)))
cat(sprintf("  resolution: C2 gives %.1f distinct values per genome, q_R gives %.1f\n",
            mean(pg$distinct_c2), mean(pg$distinct_qR)))
