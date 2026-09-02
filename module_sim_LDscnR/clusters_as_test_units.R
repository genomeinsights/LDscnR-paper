## =====================================================================
## module_sim_LDscnR / clusters_as_test_units.R
##
## An alternative design: use the LD CLUSTERS as the test units all along,
## instead of clustering the outliers after the fact.
##
##   current  association -> select outliers -> cluster the outliers -> score
##   here     cluster ALL markers (genotype only) -> a cluster is CALLED if it
##            contains a significant SNP -> score the clusters
##
## Three things this fixes:
##   - units are phenotype-INDEPENDENT, so both methods are scored on exactly
##     the same fixed set. The current design lets each method's selection
##     define its own units, which is why C's larger candidate set produced more
##     regions and why over/under-splitting was a confound at all.
##   - satellite clusters cannot exist: a QTN's neighbourhood is ONE unit.
##   - the multiple-testing burden is a known number of units, fixed in advance
##     rather than created by the selection.
##
## CAVEAT, and it is the mirror image of over-splitting: a cluster built at a
## fixed threshold can MERGE two distinct QTN into one unit, capping recall at
## the design level. Reported as `merged_units`.
##
## STAGE 2, RECOMPUTED. The bundles persist complexity_reduction$stage1 and the
## pruned marker vector, but not the stage-2 group membership -- the parse ran
## ld_prune_and_eMLG() to build grm_markers and kept only $pruned. Recomputing
## costs ~0.6 s per bundle file with compute_unflagged_eMLG = FALSE, and it is
## verifiable: the recomputed $pruned must equal the stored grm_markers, which
## confirms the parse arguments (min_r2_rho 0.5, distance_threshold 1e5,
## ld_w_threshold 0.025, score_threshold 0.80). Any file where it does not match
## is reported rather than silently used.
##
## Stage 2 differs little from stage 1 here -- 13,437 groups against 13,456
## clusters on a test file, same 60% singletons -- because the merge step only
## touches clusters flagged by ld_w > 0.025. So this is a caveat removed, not a
## result changed. Set UNITS=stage1 to compare directly.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/clusters_as_test_units.R
## Env: SIM_DATA, OUT, CELLS, MIN_LOCI (unit size floor), TAU
## =====================================================================
suppressMessages({library(data.table); library(LDscnR)})
## scoring distance cap. 1e5 matches the bundles' clustering distance_threshold
## (commit 8dbb09a harmonised these); earlier runs used a stale 5e5, which scored
## regions built at 100 kb against a 500 kb truth window.
DCAP <- as.numeric(Sys.getenv("DCAP", "1e5"))
SIM   <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
OUT   <- Sys.getenv("OUT", "module_sim_LDscnR/results/bgs5_headtohead")
CELLS <- strsplit(Sys.getenv("CELLS", "V0.5_c1,V1_c1.5,V2_c1"), ",")[[1]]
MINL  <- as.integer(strsplit(Sys.getenv("MIN_LOCI", "1,2,3,5"), ",")[[1]])
TAU   <- as.numeric(Sys.getenv("TAU", "0.05"))
UNITS <- Sys.getenv("UNITS", "stage2")
QSTAR <- seq(0, 0.95, by = 0.05)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
OUTF <- file.path(OUT, sprintf("clusters_as_test_units_%s.csv", Sys.getenv("UNITS", "stage2")))

## one bundle file = 2 chromosomes, which is the unit stage-1 clustering was run on
one <- function(cell, tag, env, i) {
  f <- sprintf("%s/adapt_%s_chr%d_%s_env%d.rds", SIM, tag, i, cell, env)
  if (!file.exists(f)) return(NULL)
  x <- readRDS(f); m <- flag_true_qtns(as.data.table(x$map))
  if (UNITS == "stage2") {
    pr <- ld_prune_and_eMLG(GTs = x$GTs, stage1 = x$complexity_reduction$stage1,
            LD_decay = x$LD_decay, ld_w_col = "ld_w_095", ld_w_threshold = 0.025,
            score_threshold = 0.80, min_r2_rho = 0.5, distance_threshold = 1e5,
            compute_unflagged_eMLG = FALSE, cores = 1)
    ## the recomputation must reproduce the stored GRM marker set, or the parse
    ## arguments are wrong and every unit below would be built on the wrong
    ## partition -- fail loudly rather than proceed
    if (!identical(sort(pr$pruned), sort(x$grm_markers)))
      return(data.table(cell, tag, env, file = i, min_loci = NA_integer_,
                        mismatch = TRUE))
    g <- as.data.table(pr$groups)
    ms <- data.table(marker = unlist(g$members, use.names = FALSE), CL_id = rep.int(g$group_id, lengths(g$members)), n_loci = rep.int(g$n_loci, lengths(g$members)))
  } else {
    ms <- as.data.table(x$complexity_reduction$stage1$map_snp)[, .(marker, CL_id, n_loci)]
  }
  m <- merge(m, ms, by = "marker", all.x = TRUE)
  m <- m[!is.na(CL_id)]
  qtn <- m[true_pos_QTN %in% TRUE]
  if (!nrow(qtn) || !nrow(m)) return(NULL)
  ds <- as.data.table(x$LD_decay$decay_sum)
  th <- score_thresholds(ds, rho_r2 = 0.75, rho_d = 0.95, dmax_cap = DCAP)
  p <- m$emx_p; if (is.null(p) || all(is.na(p))) return(NULL)
  C <- ld_cscore(p, x$ld_ws[m$marker, , drop = FALSE], alpha = 0.05,
                 rho = colnames(x$ld_ws), qstar = QSTAR)
  cv <- rep(0, nrow(m)); names(cv) <- m$marker; cv[names(C)] <- C
  m[, `:=`(C = cv, q = p.adjust(p, "BH"))]

  ## TRUTH ON UNITS: a cluster tags QTN j if any member has r2 >= r2min with it
  ## AND lies within dmax. Dedup exactly as the OR rule does -- one unit per QTN,
  ## the one with the highest r2 -- so a QTN split across units is one discovery.
  link <- rbindlist(lapply(seq_len(nrow(qtn)), function(j) {
    ch <- as.character(qtn$Chr[j])
    near <- m[as.character(Chr) == ch & abs(Pos - qtn$Pos[j]) < th$dmax]
    if (!nrow(near)) return(NULL)
    r2 <- suppressWarnings(cor(x$GTs[, qtn$marker[j]], x$GTs[, near$marker],
                               use = "pairwise.complete.obs")^2)
    ok <- which(is.finite(r2) & r2 >= th$r2min)
    if (!length(ok)) return(NULL)
    data.table(qtn = j, CL_id = near$CL_id[ok], r2 = as.numeric(r2)[ok])
  }))
  rbindlist(lapply(MINL, function(ML) {
    keep <- m[n_loci >= ML]
    if (!nrow(keep)) return(NULL)
    units <- unique(keep$CL_id)
    lk <- if (is.null(link) || !nrow(link)) NULL else link[CL_id %in% units]
    ## best unit per QTN, then one QTN per unit -- a unit holding two QTN is one TP
    true_units <- integer(0); merged <- 0L
    if (!is.null(lk) && nrow(lk)) {
      best <- lk[, .(CL_id = CL_id[which.max(r2)]), by = qtn]
      merged <- sum(duplicated(best$CL_id))
      true_units <- unique(best$CL_id)
    }
    called <- function(sel) unique(keep[sel, CL_id])
    sc <- function(cu) {
      TP <- length(intersect(cu, true_units)); FP <- length(setdiff(cu, true_units))
      FN <- length(setdiff(true_units, cu))
      pr <- if (TP + FP > 0) TP/(TP+FP) else NA_real_
      rc <- if (TP + FN > 0) TP/(TP+FN) else NA_real_
      c(length(cu), TP, FP, FN, pr, rc, if (is.finite(pr) && is.finite(rc)) pr*rc else NA_real_) }
    a <- sc(called(keep$q < 0.05)); cc <- sc(called(keep$C >= TAU))
    data.table(cell, tag, env, file = i, min_loci = ML,
               n_units = length(units), n_true_units = length(true_units),
               n_qtn = nrow(qtn), merged_units = merged,
               a_called = a[1], a_TP = a[2], a_prec = a[5], a_rec = a[6], a_PR = a[7],
               C_called = cc[1], C_TP = cc[2], C_prec = cc[5], C_rec = cc[6], C_PR = cc[7])
  }))
}
grid <- CJ(cell = CELLS, tag = c("bgs","nobgs"), env = 1:10, file = 1:10)
res <- rbindlist(lapply(seq_len(nrow(grid)), function(r)
  one(grid$cell[r], grid$tag[r], grid$env[r], grid$file[r])), fill = TRUE)
fwrite(res, OUTF)
if ("mismatch" %in% names(res) && any(res$mismatch %in% TRUE)) {
  cat(sprintf("*** %d files failed the grm_markers reproduction check ***\n",
              sum(res$mismatch %in% TRUE)))
  print(res[mismatch %in% TRUE, .(cell, tag, env, file)])
  res <- res[!(mismatch %in% TRUE)] }
res <- res[!is.na(min_loci)]
cat(sprintf("  units: %s | %d bundle-files, %d cells\n\n", UNITS,
            uniqueN(res[, .(cell,tag,env,file)]), uniqueN(res$cell)))
cat("=== clusters as test units, by unit-size floor ===\n")
print(res[, .(n = .N, units = round(mean(n_units)), true_units = round(mean(n_true_units), 2),
              QTN = round(mean(n_qtn), 2), merged = round(mean(merged_units), 2),
              a_called = round(mean(a_called), 1), a_PR = round(mean(a_PR, na.rm = TRUE), 3),
              C_called = round(mean(C_called), 1), C_PR = round(mean(C_PR, na.rm = TRUE), 3),
              dPR = round(mean(C_PR - a_PR, na.rm = TRUE), 3)), by = min_loci][order(min_loci)])
cat("\n=== per cell at min_loci = 2 ===\n")
print(res[min_loci == 2, .(n = .N, a_prec = round(mean(a_prec, na.rm=TRUE), 2),
          a_rec = round(mean(a_rec, na.rm=TRUE), 2), a_PR = round(mean(a_PR, na.rm=TRUE), 3),
          C_prec = round(mean(C_prec, na.rm=TRUE), 2), C_rec = round(mean(C_rec, na.rm=TRUE), 2),
          C_PR = round(mean(C_PR, na.rm=TRUE), 3)), by = cell][order(cell)])
