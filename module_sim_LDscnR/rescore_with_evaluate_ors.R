## =====================================================================
## module_sim_LDscnR / rescore_with_evaluate_ors.R
##
## Rescore saved ld_scan fits with the module's canonical truth matcher,
## LDscnR::evaluate_ors(), instead of the positional-containment test that
## analyse_one_dataset.R's evaluation section uses.
##
## WHY. Containment asks "does the region's span cover a QTN's base pair". That
## is not what an association scan detects: it detects markers IN LD WITH a
## causal variant, which need not lie inside the region. evaluate_ors() instead
## assigns each region a focal QTN by LD and distance (r2 > r2_match,
## dist_bp < d_match, thresholds from score_thresholds() on the panel's own decay
## fit), then DEDUPLICATES -- if several regions claim one QTN, the
## highest-evidence region keeps it and the rest are dropped. Dropped regions
## that point at a real QTN count as NEITHER TP nor FP.
##
## It also scores against `true_pos_QTN`, not `true_QTN`: flag_true_qtns() keeps
## only QTNs with MAF > 0.1 that carry >= 5% of their chromosome's additive
## variance. Across these ten cells that is 135 of 181 QTN loci -- the other 46
## are undetectable in principle, and counting them in the denominator understates
## recall by construction.
##
## Nothing expensive is recomputed. Region membership is rebuilt from the saved
## fit's own null and edge graph (verified identical to the saved region table),
## and qtn_ld_table() -- the costly part -- is built ONCE per cell and reused
## across both surrogate bases and every l_min.
##
## THE INTERACTION TO WATCH. Dedup forgiveness is method-dependent: an arm that
## fragments a true signal into many regions has most of them excluded from the
## FP denominator, while an arm that fragments around a false signal pays in
## full. l_min changes fragmentation directly, so the dedup-loser count is
## reported at every l_min rather than absorbed into precision.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/rescore_with_evaluate_ors.R <scan_dir> [outdir]
## Env: PANEL_DIR (default /Volumes/Nemo/Nemo_sim/analysis_inputs), CORES,
##      LMINS (default 2,3,5,10,20), R2_GRID/D_GRID for the sensitivity check.
## =====================================================================
suppressMessages({ library(data.table); library(LDscnR) })
`%||%` <- function(a, b) if (is.null(a)) b else a

a <- commandArgs(trailingOnly = TRUE)
if (!length(a)) stop("usage: rescore_with_evaluate_ors.R <scan_dir> [outdir]")
SCAN_DIR <- a[1]; OUT <- if (length(a) >= 2) a[2] else SCAN_DIR
PANEL_DIR <- Sys.getenv("PANEL_DIR", "/Volumes/Nemo/Nemo_sim/analysis_inputs")
CORES <- as.integer(Sys.getenv("CORES", "2"))
LMINS <- as.integer(strsplit(Sys.getenv("LMINS", "2,3,5,10,20"), ",")[[1]])
RHO_R2 <- 0.75; RHO_D <- 0.95; DCAP <- 1e5      # score_thresholds(); matches the bundles' stage-2 cap
## Threshold sensitivity. Fixed scoring thresholds are standard in method
## papers; what has to be stable is the ORDERING of the settings compared, not
## the absolute counts. Supplying R2_GRID/D_GRID re-scores the same regions at
## several (r2_match, d_match) pairs so the ordering can be checked. Empty =
## use only the decay-derived values, which is the primary result.
R2_GRID <- { .g <- Sys.getenv("R2_GRID", ""); if (nzchar(.g)) as.numeric(strsplit(.g, ",")[[1]]) else NA_real_ }
D_GRID  <- { .g <- Sys.getenv("D_GRID",  ""); if (nzchar(.g)) as.numeric(strsplit(.g, ",")[[1]]) else NA_real_ }
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

scans <- list.files(SCAN_DIR, pattern = "^scan_.*[.]rds$", full.names = TRUE)
if (!length(scans)) stop("no scan_*.rds in ", SCAN_DIR)
cellof <- function(f) sub("^scan_(V[0-9.]+_c[0-9.]+_env[0-9]+)_.*$", "\\1", basename(f))
cells <- unique(vapply(scans, cellof, character(1)))
cat(sprintf("[1] %d scan(s) across %d cell(s); l_min grid: %s\n",
            length(scans), length(cells), paste(LMINS, collapse = ","))); flush.console()

## span of a member list, computed exactly as .region_table does, so the
## reconstructed regions can be matched to the saved significance table
span_of <- function(r, pos, chr)
  data.table(chr = unname(chr[vapply(r, `[`, character(1), 1L)]),
             lo = vapply(r, function(x) min(pos[x]), numeric(1)),
             hi = vapply(r, function(x) max(pos[x]), numeric(1)),
             size = lengths(r))

out <- list(); k <- 0L
for (cell in cells) {
  pf <- file.path(PANEL_DIR, sprintf("panel_%s.rds", cell))
  if (!file.exists(pf)) { cat(sprintf("  [skip] no panel for %s\n", cell)); next }
  panel <- readRDS(pf); map <- as.data.table(panel$map)
  if (!"true_pos_QTN" %in% names(map)) map <- flag_true_qtns(map)
  n_true <- sum(map$true_pos_QTN %in% TRUE)
  th <- score_thresholds(as.data.table(panel$decay_sum), rho_r2 = RHO_R2, rho_d = RHO_D, dmax_cap = DCAP)
  sc <- scans[vapply(scans, cellof, character(1)) == cell]

  ## union of every marker that any basis clusters at any l_min -- qtn_ld_table
  ## is built once over this and reused throughout the cell
  uni <- unique(unlist(lapply(sc, function(f) { x <- readRDS(f)
    names(x$null$C_obs)[which(x$null$C_obs >= x$params$tau)] }), use.names = FALSE))
  t0 <- Sys.time()
  qtab <- qtn_ld_table(panel$GTs, map, uni, 2e6, cores = CORES)
  cat(sprintf("[2] %s: %d true_pos_QTN of %d QTN | universe %d | qtn_ld_table %.1f min | r2_match %.3f d_match %.0f\n",
              cell, n_true, sum(map$type == "QTN"), length(uni),
              as.numeric(Sys.time() - t0, units = "mins"), th$r2min, th$dmax)); flush.console()

  pos <- stats::setNames(map$Pos, map$marker); chrv <- stats::setNames(as.character(map$Chr), map$marker)
  qtn_pos <- map[true_pos_QTN %in% TRUE, .(Chr = as.character(Chr), Pos)]
  ## A pooled cell is a GENOME of 20 chromosomes, and about half of them carry no
  ## detectable QTN -- as in any real genome. The QTN-free set is therefore NOT
  ## "the Chr2s": across these ten genomes 104 of 200 chromosomes are QTN-free and
  ## 4 of those are Chr1s that happen to carry none. Counting only _Chr2 both
  ## undercounted and implied Chr2 is a designated control rather than simply one
  ## chromosome among many without signal.
  per_chr <- map[, .(nq = sum(true_pos_QTN %in% TRUE)), by = Chr]
  qtn_free_chr <- as.character(per_chr[nq == 0, Chr])
  n_chr <- nrow(per_chr)
  frac_chr_qtn_free <- length(qtn_free_chr) / n_chr   # the null expectation for a call landing there

  for (f in sc) {
    x <- readRDS(f); basis <- x$params$basis; tau <- x$params$tau
    mk <- names(x$null$C_obs)[which(x$null$C_obs >= tau)]
    all_r <- ld_regions(mk, x$edges)
    for (L in LMINS) {
      r <- all_r[lengths(all_r) >= L]
      if (!length(r)) { k <- k + 1L
        out[[k]] <- data.table(cell, basis, l_min = L, n_regions = 0L, n_sig = 0L,
                               TP = 0L, FP = 0L, FN = n_true, extras = 0L,
                               precision = NA_real_, recall = 0, dedup_losers = 0L,
                               contain_TP = 0L, contain_recall = 0,
                               sig_on_qtn_free = 0L, sig_on_qtn_bearing = 0L,
                               n_chr = n_chr, frac_chr_qtn_free = frac_chr_qtn_free,
                               n_true = n_true, lambda_obs = x$null$lambda_obs %||% NA_real_,
                               lambda_surr_med = stats::median(x$null$lambda_surr %||% NA_real_),
                               ## a cell that called NO regions is a primary result -- recall 0,
                               ## precision undefined. Omitting the flag here dropped exactly the
                               ## undetectable cells from the aggregate and biased it upward.
                               r2_match = th$r2min, d_match = th$dmax, primary = TRUE)
        next }
      ## significance from a fresh region scan at this l_min (cheap: re-clusters
      ## the surrogates, but does NOT redo the C-score reduction)
      sr <- ld_region_scan(x$null, x$edges, tau = tau, l_min = L, fdr = x$params$fdr)
      sp <- span_of(r, pos, chrv)
      sig_key <- sr$regions[sig == TRUE, paste(chr, lo, hi)]
      keep <- paste(sp$chr, sp$lo, sp$hi) %in% sig_key
      rs <- r[keep]; sps <- sp[keep]

      ## primary = the decay-derived thresholds; extra settings only if asked for
      grid <- data.table(r2m = th$r2min, dm = th$dmax, primary = TRUE)
      if (!all(is.na(R2_GRID)) || !all(is.na(D_GRID)))
        grid <- rbind(grid, CJ(r2m = if (all(is.na(R2_GRID))) th$r2min else R2_GRID,
                               dm  = if (all(is.na(D_GRID)))  th$dmax  else D_GRID)[, primary := FALSE])
      for (gi in seq_len(nrow(grid))) {
      r2m <- grid$r2m[gi]; dm <- grid$dm[gi]
      ev <- evaluate_ors(rs, map, qtab, r2m, dm)
      dg <- tryCatch(LDscnR:::.diagnose_ors(rs, map, qtab, r2m, dm), error = function(e) NULL)
      dl <- if (!is.null(dg) && nrow(dg)) sum(dg$dropped_by_dedup %in% TRUE) else NA_integer_

      ## positional containment, as the threshold-free secondary
      cTP <- if (nrow(sps)) sum(vapply(seq_len(nrow(sps)), function(i)
        any(qtn_pos$Chr == sps$chr[i] & qtn_pos$Pos >= sps$lo[i] & qtn_pos$Pos <= sps$hi[i]), logical(1))) else 0L
      cRec <- if (n_true) sum(vapply(seq_len(nrow(qtn_pos)), function(j) nrow(sps) > 0 &&
        any(sps$chr == qtn_pos$Chr[j] & sps$lo <= qtn_pos$Pos[j] & sps$hi >= qtn_pos$Pos[j]), logical(1))) / n_true else NA_real_

      k <- k + 1L
      out[[k]] <- data.table(cell, basis, l_min = L, n_regions = length(r), n_sig = length(rs),
        TP = ev$TP, FP = ev$FP, FN = ev$FN, extras = ev$extras,
        precision = ev$Precision, recall = ev$Recall, dedup_losers = dl,
        contain_TP = cTP, contain_recall = cRec,
        sig_on_qtn_free = sum(sps$chr %in% qtn_free_chr),
        sig_on_qtn_bearing = sum(!(sps$chr %in% qtn_free_chr)),
        n_chr = n_chr, frac_chr_qtn_free = frac_chr_qtn_free, n_true = n_true,
        lambda_obs = x$null$lambda_obs %||% NA_real_,
        lambda_surr_med = stats::median(x$null$lambda_surr %||% NA_real_),
        r2_match = r2m, d_match = dm, primary = grid$primary[gi])
      }
    }
    cat(sprintf("    %s / %-9s done\n", cell, basis)); flush.console()
  }
}
res <- rbindlist(out, fill = TRUE)
fwrite(res, file.path(OUT, "rescore_evaluate_ors.csv"))
cat(sprintf("\n[3] wrote rescore_evaluate_ors.csv (%d rows)\n", nrow(res)))

cat("\n=== replicate-averaged over cells (mean +- SE), decay-derived thresholds ===\n")
## recall is averaged over EVERY cell -- a cell that found nothing has recall 0
## and must count. precision is undefined where no region was called (0/0), so it
## is averaged over the cells that called at least one, and both n's are reported.
agg <- res[primary %in% TRUE, .(n_cells = .N, n_prec = sum(!is.na(precision)),
               recall = sprintf("%.3f+-%.3f", mean(recall, na.rm=TRUE), sd(recall, na.rm=TRUE)/sqrt(.N)),
               precision = sprintf("%.3f+-%.3f", mean(precision, na.rm=TRUE),
                                   sd(precision, na.rm=TRUE)/sqrt(sum(!is.na(precision)))),
               TP = round(mean(TP),1), FP = round(mean(FP),1), dedup = round(mean(dedup_losers, na.rm=TRUE),1),
               on_qtn_free = round(mean(sig_on_qtn_free),1),
               pct_on_qtn_free = round(100*sum(sig_on_qtn_free)/pmax(sum(sig_on_qtn_free+sig_on_qtn_bearing),1)),
               pct_expected_by_chance = round(100*mean(frac_chr_qtn_free)),
               contain_recall = round(mean(contain_recall, na.rm=TRUE),3)),
           by = c("basis", "l_min")]
setorderv(agg, c("basis", "l_min")); print(agg, nrow = 100)

if (any(res$primary %in% FALSE)) {
  cat("\n=== threshold sensitivity: is the l_min ORDERING stable? ===\n")
  o <- res[, .(precision = mean(precision, na.rm = TRUE), recall = mean(recall, na.rm = TRUE)),
           by = c("basis", "r2_match", "d_match", "l_min")]
  o[, rank_prec := frank(-precision, ties.method = "min"), by = c("basis", "r2_match", "d_match")]
  rk <- dcast(o, basis + r2_match + d_match ~ l_min, value.var = "rank_prec")
  print(rk, nrow = 100)
  cat("  identical rows => the ranking of l_min by precision does not depend on the thresholds\n")
}
