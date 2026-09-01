## =============================================================================
## module_sim_LDscnR / proxy_T1_local_decay.R
##
## T1 of the pre-registered proxy test: does the LOCAL LD-decay rate track the
## LOCAL recombination rate? PK's argument for making distance_threshold adapt
## to the recombination landscape rests on this, and the sims are the right
## instrument because the true rate is a simulation INPUT, not a pedigree
## estimate -- 11,527 distinct rec_rate values per chromosome against the
## panel's 838 map bins at a median of 489 kb.
##
## NO NEW FIT IS NEEDED. LD_decay$by_chr[[chr]]$decay already holds a per-window
## a/b/c with start/end, which is a(x) as built. Ground truth per window is
## (cM[end] - cM[start]) / span_Mb from the map's own cM column.
##
## SIGN: a is a DECAY RATE, so fast decay (high a) should accompany high
## recombination. The proxy passes T1 if the within-chromosome Spearman is
## clearly positive in a clear majority of chromosomes.
##
## RESULT: T1 PASSES DECISIVELY. Median within-chromosome Spearman +0.757,
## positive in 20 of 20 chromosomes (sign p = 1.9e-06), pooled +0.740. Stronger
## than the panel's +0.605, which is what an exact truth rather than a 489 kb
## map bin should buy.
##
## AND THE TRANSFER FUNCTION FAILS FOR A REASON THAT REVERSES THE OBVIOUS FIX.
## d_from_rho(a, 0.95) over these windows runs 111 kb to 4.7 Mb on the finite
## ones, and 17 of 363 windows have a == 0 EXACTLY, giving an infinite distance.
## The natural reading is that those are failed fits and the cure is to shrink
## a_local toward the chromosome fit where the window is poorly determined.
## THAT READING IS WRONG. The degenerate windows have a median true rate of
## 0.083 cM/Mb against 3.310 elsewhere (Wilcoxon p = 1.1e-09) while their fit
## quality is essentially unchanged (contrast 0.948 against 0.965) -- they are
## COLD BLOCKS being reported correctly. Recombination really is ~zero there, so
## LD really does not decay, so the appropriate merge distance really is very
## large. Shrinking them toward the chromosome mean would destroy precisely the
## signal that motivates a local threshold.
##
## So the fix belongs in the transfer function after all, not in the estimate:
## d = rho/(a(1-rho)) diverges as a -> 0 and a -> 0 is a CORRECT answer. What
## the cap should be dimensioned by is the size of the feature to be bridged --
## measured here at q90 550 kb and max 3.7 Mb for zero-recombination blocks --
## not by anything derived from decay. The local rate says WHERE to widen; the
## block-size distribution says HOW MUCH. Neither is sufficient alone.
##
## T1 IS NOT T2. A quantity can track truth and still not improve the partition
## -- cluster geometry is 30x enriched for containing QTN and no restriction
## built on it beats using every cluster. T2 (does a local threshold beat the
## flat one on QTN precision AND recall) is the test that decides, and T2b runs
## it inside versus outside cold blocks, where a flat threshold should fail if
## it fails anywhere.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/proxy_T1_local_decay.R
## Env: SIM_DATA, CELL, TAG, ENV, FILES, OUT
## =============================================================================
suppressMessages({library(data.table); library(LDscnR)})
SIM   <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
OUT   <- Sys.getenv("OUT", "module_sim_LDscnR/results/operating_points")
CELL  <- Sys.getenv("CELL", "V0.5_c1")
TAG   <- Sys.getenv("TAG", "nobgs")
ENV   <- as.integer(Sys.getenv("ENV", "1"))
FILES <- as.integer(strsplit(Sys.getenv("FILES", "1,2,3,4,5,6,7,8,9,10"), ",")[[1]])
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

W <- rbindlist(lapply(FILES, function(i) {
  f <- sprintf("%s/adapt_%s_chr%d_%s_env%d.rds", SIM, TAG, i, CELL, ENV)
  if (!file.exists(f)) return(NULL)
  x <- readRDS(f); m <- as.data.table(x$map)
  rbindlist(lapply(names(x$LD_decay$by_chr), function(ch) {
    d  <- as.data.table(x$LD_decay$by_chr[[ch]]$decay)
    mm <- m[Chr == ch][order(Pos)]
    ## Truth per window from the map's own cM column. Windows overlap (overlap
    ## 0.5), which inflates neither the correlation's sign nor its direction but
    ## does make the per-chromosome n optimistic -- hence the sign test across
    ## chromosomes rather than a pooled p-value.
    tru <- vapply(seq_len(nrow(d)), function(k) {
      s <- mm[Pos >= d$start[k] & Pos <= d$end[k]]
      if (nrow(s) < 2) return(NA_real_)
      (max(s$cM) - min(s$cM)) / ((d$end[k] - d$start[k]) / 1e6) }, numeric(1))
    data.table(file = i, Chr = ch, w = seq_len(nrow(d)), a = d$a, c = d$c,
               contrast = d$contrast, regime = d$regime, cMMb = tru,
               span = d$end - d$start) })) }))
W <- W[is.finite(a) & is.finite(cMMb)]
fwrite(W, file.path(OUT, "proxy_T1_windows.csv"))

R <- W[, .(n_win = .N, rho = suppressWarnings(stats::cor(a, cMMb, method = "spearman"))),
       by = .(file, Chr)][is.finite(rho)]
fwrite(R, file.path(OUT, "proxy_T1_by_chr.csv"))
cat(sprintf("T1: %d chromosomes, %d windows, median %d windows/chr\n",
            nrow(R), nrow(W), as.integer(median(R$n_win))))
cat(sprintf("  Spearman(a_local, true cM/Mb) within chromosome:\n"))
cat(sprintf("    median %+.3f | IQR %+.3f to %+.3f | range %+.3f to %+.3f\n",
    median(R$rho), quantile(R$rho, .25), quantile(R$rho, .75), min(R$rho), max(R$rho)))
cat(sprintf("    positive in %d of %d chromosomes, sign p = %.5f\n",
    sum(R$rho > 0), nrow(R), stats::binom.test(sum(R$rho > 0), nrow(R))$p.value))
cat(sprintf("\n  pooled across all windows: %+.3f\n",
    stats::cor(W$a, W$cMMb, method = "spearman")))
## Fit quality per window, for the shrinkage transfer function: where the window
## fit is poorly determined, a_local should be pulled toward the chromosome fit.
cat(sprintf("\n  window fit: contrast median %.3f (q10 %.3f) | regimes: %s\n",
    median(W$contrast, na.rm = TRUE), quantile(W$contrast, .1, na.rm = TRUE),
    paste(sprintf("%s %d", names(table(W$regime)), table(W$regime)), collapse = ", ")))
cat(sprintf("  a_local spread %.0f-fold -> d_from_rho(a_local, 0.95) %.0f kb to %.0f kb\n",
    max(W$a) / min(W$a), min(d_from_rho(W$a, 0.95)) / 1e3, max(d_from_rho(W$a, 0.95)) / 1e3))
