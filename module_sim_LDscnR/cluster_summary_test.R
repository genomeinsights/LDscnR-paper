## =============================================================================
## cluster_summary_test.R -- four ways to give a stage-2 cluster one p-value.
##
## The module has been using the LD-central pruned representative throughout.
## That choice is defensible but costly: no driving QTN is ever the central
## marker (0 of 21 on the examined panel), because QTN sit in large clusters.
## This compares it against three alternatives on identical clusters, identical
## GRM and identical BH treatment, so the only thing that varies is how a
## cluster's member p-values become one number.
##
##   representative  p of the LD-central marker            (current)
##   eMLG            EMMAX on the cluster's consensus genotype
##   best SNP        min p over members -- ANTI-CONSERVATIVE, no correction,
##                   included as the ceiling rather than as a usable rule
##   Simes           Simes combination over members, which is what the
##                   stickleback analysis uses and is valid under PRDS
##   max ld_w        p of the member with the highest ld_w -- a phenotype-blind
##                   choice of representative, as opposed to the LD-central one
##
## Scored at LOCUS level: clusters flagged, and how many contain a driving QTN.
##
## CROSSED WITH ld_w FILTERING. Every filtering result elsewhere in this module
## used the REPRESENTATIVE p-value, which the comparison below shows is the
## weakest of the summarisation rules. Filtering has therefore only ever been
## evaluated on a lossy summary. Here each summarisation is crossed with
## cluster-ranked ld_w selection at several k, so the two axes are separated.
## Env: SIM_DATA, CELL, TAG, ENV, FILES, ALPHA
## =============================================================================
suppressMessages({library(data.table); library(LDscnR)})
SIM  <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
OUT  <- Sys.getenv("OUT", "module_sim_LDscnR/results/filter_then_test")
CELL <- Sys.getenv("CELL", "V0.5_c1"); TAG <- Sys.getenv("TAG", "nobgs")
ENV  <- as.integer(Sys.getenv("ENV", "3"))
FILES <- as.integer(strsplit(Sys.getenv("FILES", "1,2,3,4,5,6,7,8,9,10"), ",")[[1]])
ALPHA <- as.numeric(Sys.getenv("ALPHA", "0.05"))
## The cluster partition is built from GENOTYPES and never sees a p-value, so
## swapping the association method changes exactly one column. The eMLG arm is
## the exception -- it refits the association on consensus genotypes -- and is
## skipped unless the engine is emmax.
ENG   <- Sys.getenv("ENGINE", "emmax")
PCOL  <- if (ENG == "emmax") "emx_p" else "lfmm_p"
KS    <- as.integer(strsplit(Sys.getenv("KS", "1000,5000,20000,50000"), ",")[[1]])
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

simes <- function(p) { p <- sort(p[is.finite(p)]); n <- length(p)
  if (!n) return(NA_real_); min(n * p / seq_len(n)) }

per_file <- function(i) {
  f <- sprintf("%s/adapt_%s_chr%d_%s_env%d.rds", SIM, TAG, i, CELL, ENV)
  if (!file.exists(f)) return(NULL)
  x <- readRDS(f); m <- flag_true_qtns(as.data.table(x$map))
  pr <- ld_prune_and_eMLG(GTs = x$GTs, stage1 = x$complexity_reduction$stage1,
          LD_decay = x$LD_decay, ld_w_col = "ld_w_095", ld_w_threshold = 0.025,
          score_threshold = 0.80, min_r2_rho = 0.5, distance_threshold = 1e5,
          compute_unflagged_eMLG = TRUE, cores = 1)
  g  <- as.data.table(pr$groups)
  ms <- rbindlist(lapply(seq_len(nrow(g)), function(k)
          data.table(marker = g$members[[k]], CL = paste0(i, "_", g$group_id[k]))))
  mm <- merge(m, ms, by = "marker", all.x = TRUE)[!is.na(CL)]

  ## per-cluster summaries from the stored per-marker p-values
  su <- mm[, .(p_rep  = emx_p[marker %in% pr$pruned][1],
               p_maxldw = get(PCOL)[which.max(ld_w_095)],
               p_best = min(get(PCOL), na.rm = TRUE),
               p_simes= simes(emx_p),
               n      = .N,
               ld_w   = median(ld_w_095, na.rm = TRUE),
               has_qtn= any(true_pos_QTN %in% TRUE),
               Chr = Chr[1], Pos = median(Pos)), by = CL]

  ## eMLG: EMMAX on the consensus genotypes, same GRM as the bundle used
  E <- pr$eMLG
  p_emlg <- rep(NA_real_, nrow(su)); names(p_emlg) <- su$CL
  if (!is.null(E) && ncol(E) > 0) {
    cn <- colnames(E)
    if (!is.null(cn)) {
      keep <- paste0(i, "_", cn) %in% su$CL
      Ek <- E[, keep, drop = FALSE]
      res <- tryCatch({
        ## emmax() returns list(F, pval, Rsq) -- NOT $p. An earlier version
        ## checked for $p, fell through to as.numeric() on a list, and failed.
        pp <- emmax(Y = as.numeric(x$env$env), X = Ek, K = x$GRM, cores = 1)
        stopifnot(is.list(pp), "pval" %in% names(pp))
        as.numeric(pp$pval)
      }, error = function(e) { message("    emmax on eMLG failed (", basename(f), "): ", conditionMessage(e)); NULL })
      if (!is.null(res) && length(res) == ncol(Ek))
        p_emlg[paste0(i, "_", colnames(Ek))] <- res
    }
  }
  su[, p_emlg := p_emlg[CL]]
  su[, file := i][]
}
su <- rbindlist(lapply(FILES, per_file), fill = TRUE)
cat(sprintf("\n  %d clusters, %d contain a driving QTN; eMLG p available for %d (%.0f%%)\n",
            nrow(su), sum(su$has_qtn), sum(is.finite(su$p_emlg)), 100*mean(is.finite(su$p_emlg))))

score <- function(col, keep = NULL, kk = NA_integer_) {
  p <- su[[col]]; ok <- is.finite(p)
  if (!is.null(keep)) ok <- ok & keep
  q <- rep(NA_real_, length(p)); q[ok] <- p.adjust(p[ok], "BH")
  sig <- which(!is.na(q) & q < ALPHA)
  data.table(summary = col, k = kk, tested = sum(ok), flagged = length(sig),
             with_qtn = sum(su$has_qtn[sig]),
             precision = if (length(sig)) mean(su$has_qtn[sig]) else NA_real_,
             recall = sum(su$has_qtn[sig]) / sum(su$has_qtn),
             cutoff = if (length(sig)) max(p[sig]) else NA_real_)
}
SUMS <- c("p_rep","p_emlg","p_best","p_simes","p_maxldw")
res <- rbindlist(lapply(SUMS, score))
res[, PR := precision * recall]
cat("\n=== A. no filtering, all clusters ===\n"); print(res)

## cluster-ranked ld_w selection, crossed with summarisation
ord <- order(-su$ld_w)
cx <- rbindlist(lapply(KS, function(kk) {
  if (kk >= nrow(su)) return(NULL)
  keep <- rep(FALSE, nrow(su)); keep[head(ord, kk)] <- TRUE
  rbindlist(lapply(SUMS, function(cl) score(cl, keep, kk)))
}))
cx[, PR := precision * recall]
cat("\n=== B. crossed with cluster-ranked ld_w filtering ===\n")
print(dcast(cx, k ~ summary, value.var = "with_qtn")[order(k)])
cat("  (cells are QTN-bearing clusters recovered, out of", sum(su$has_qtn), ")\n")
cat("\n  precision:\n"); print(dcast(cx, k ~ summary, value.var = "precision")[order(k)])
cat("\n  clusters flagged:\n"); print(dcast(cx, k ~ summary, value.var = "flagged")[order(k)])
res <- rbind(res, cx, fill = TRUE)
fwrite(cbind(cell = CELL, tag = TAG, env = ENV, res),
       file.path(OUT, sprintf("cluster_summary_test_%s_%s_env%d.csv", CELL, TAG, ENV)))
cat("\n  p_best has NO within-cluster correction and is anti-conservative;\n")
cat("  it is the ceiling, not a usable rule. p_simes is the valid version.\n")
