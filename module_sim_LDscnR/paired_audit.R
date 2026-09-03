## =============================================================================
## paired_audit.R -- re-score every alpha-vs-C comparison in this module as a
## PAIRED test with a sign count, not a mean of paired differences.
##
## WHY THIS EXISTS. The comparison scripts already differenced within row
## (C_minus_alpha := PR_AUC_C - PR_AUC_alpha), which is the right shape. What
## they did NOT do is report how many arms actually favoured each method. A
## mean of paired differences can be carried by a minority of arms while the
## sign count sits at chance, and that is exactly what happened: the claim
## "C wins in 3 of 4 cells" rested on per-cell means of +0.001 to +0.009 that
## a sign test shows to be null.
##
## Reads only saved CSVs -- no bundle access, runs in seconds.
## =============================================================================
suppressMessages(library(data.table))
R <- Sys.getenv("RESULTS", "module_sim_LDscnR/results")

sgn <- function(x, lab) {
  x <- x[is.finite(x)]
  nz <- x[x != 0]
  if (!length(nz)) { cat(sprintf("  %-22s n=%3d  all ties\n", lab, length(x))); return(invisible(NULL)) }
  p <- binom.test(sum(nz > 0), length(nz))$p.value
  w <- suppressWarnings(wilcox.test(x)$p.value)
  cat(sprintf("  %-22s n=%3d  mean_d=%+.4f  median_d=%+.4f  C>a %d/%d (%.0f%%)  sign p=%.3g  wilcox p=%.3g\n",
              lab, length(x), mean(x), median(x), sum(nz > 0), length(nz),
              100 * sum(nz > 0) / length(nz), p, w))
}

## -- 1. the deployable question: does tau_C = 0.05 beat alpha = 0.05? --------
cat("=== tau_C_sweep, l_min=5: paired C-alpha by tau ===\n")
d <- fread(file.path(R, "bgs5_headtohead/tau_C_sweep.csv"))
for (t in sort(unique(d$tau))) {
  if (t %in% c(0.05, 0.1, 0.2, 0.3, 0.5, 0.7)) sgn(d[l_min == 5 & tau == t, C_PR - a_PR], sprintf("tau=%.2f", t))
}

## -- 2. the withdrawn claim: "C wins in 3 of 4 cells at every rho_ld" -------
cat("\n=== rho_ld_sweep, l_min=2: pooled by rho_ld ===\n")
r <- fread(file.path(R, "bgs5_headtohead/rho_ld_sweep.csv"))
for (v in sort(unique(r$rho_ld))) sgn(r[l_min == 2 & rho_ld == v, C_PR - a_PR], sprintf("rho_ld=%.2f", v))
cat("\n=== rho_ld_sweep, l_min=2, rho_ld=0.75: BY CELL (the withdrawn claim) ===\n")
for (cl in sort(unique(r$cell))) sgn(r[l_min == 2 & rho_ld == 0.75 & cell == cl, C_PR - a_PR], cl)

## -- 3. head-to-head: PR recomputed from stored counts ----------------------
cat("\n=== bgs5_headtohead, l_min=5: paired C-alpha by cell ===\n")
h <- fread(file.path(R, "bgs5_headtohead/bgs5_headtohead.csv"))
h[, `:=`(aPR = (alpha_tp / pmax(alpha_tp + alpha_fp, 1)) * (alpha_tp / pmax(n_qtn, 1)),
         cPR = (C_tp     / pmax(C_tp + C_fp, 1))         * (C_tp     / pmax(n_qtn, 1)))]
for (cl in sort(unique(h$cell))) sgn(h[l_min == 5 & cell == cl, cPR - aPR], cl)
sgn(h[l_min == 5, cPR - aPR], "POOLED")

## -- 4. integrated PR-AUC: the one place C holds up -------------------------
cat("\n=== PR-AUC (all_bgs5): paired C - alpha, pooled ===\n")
v <- fread(file.path(R, "all_bgs5/bgs_vs_nobgs_prauc.csv"))
sgn(v[, PR_AUC_C - PR_AUC_alpha], "POOLED")
cat("\n  NOTE: all_bgs5 is the DCAP=5e5 generation. The dcap=1e5 rerun writes to\n")
cat("  results/all_bgs5_dcap1e5 and should be re-audited here when it lands.\n")
