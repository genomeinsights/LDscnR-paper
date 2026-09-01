## =============================================================================
## module_sim_LDscnR / cap_cost.R
##
## What does a cap on the derived distance COST in proxy signal? I argued to 2c
## that the cap should be dimensioned by the size of the feature to be bridged
## (zero-recombination blocks: q90 550 kb, max 3.7 Mb). This measures what that
## choice throws away, and the answer is that it is far from free.
##
## THE TWO EXCLUSION RULES ARE NOT THE SAME SET. I scored T1 excluding a == 0
## windows (7.0% at w = 20); 2c scored the panel excluding windows above a cap.
## On bgs5 a 2.1 Mb cap excludes 32%, and a 550 kb cap -- my own block q90 --
## excludes 65%. Comparing my a == 0 result against their cap result was not the
## same operation, so this runs theirs.
##
##   w    all    excl a==0   cap 2.1 Mb   cap 550 kb   (% windows kept at 2.1 Mb)
##   5   0.690     0.690       0.660        0.750           57.9
##  10   0.758     0.758       0.639        0.521           67.5
##  20   0.830     0.782       0.639        0.427           67.7
##  50   0.829     0.742       0.556        0.243           66.9
##
## THE CONTRAST WITH THE PANEL SURVIVES THE MATCHED OPERATION AND STRENGTHENS.
## Capping LOWERS the correlation here at every window count and RAISES it on
## the panel (0.633 -> 0.767). Divergent windows are cold blocks correctly
## identified on bgs5 and noise on the panel -- the same two-mechanism
## conclusion the near-zero fractions gave, now from an independent route.
##
## AND IT UNDERMINES THE CAP I PROPOSED. Capping at the block-size q90 costs
## 0.830 -> 0.427 at w = 20. The windows above any plausible cap carry much of
## the proxy's signal, because they ARE the cold blocks -- the features a local
## threshold exists to handle. So the cap is not a harmless safety rail: it is
## most destructive exactly where locality matters most, which is the same shape
## as the shrinkage error and I made it twice.
##
## WHAT THIS DOES AND DOES NOT SHOW. Excluding a window measures how much signal
## lives above the cap. It is NOT the operation a cap performs, which is
## d <- min(d_derived, cap) with the window retained. T2 has to test that, not
## this. The honest statement is that the uncapped quantity is physically
## meaningless (4,961 Mb on a 30 Mb chromosome) while the capped one discards
## the features of interest, and nothing measured so far resolves that tension.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/cap_cost.R
## Env: OUT
## =============================================================================
suppressMessages({library(data.table); library(LDscnR)})
OUT <- Sys.getenv("OUT", "module_sim_LDscnR/results/operating_points")
W <- fread(file.path(OUT, "proxy_T1_by_window_raw.csv"))
W[, d95 := ifelse(a > 0, d_from_rho(a, 0.95), Inf)]
sc <- function(d) {
  r <- d[, .(rho = suppressWarnings(stats::cor(a, cMMb, method = "spearman"))),
         by = .(w, file, Chr)][is.finite(rho)]
  if (!nrow(r)) NA_real_ else median(r$rho)
}
S <- rbindlist(lapply(sort(unique(W$w)), function(ww) {
  d <- W[w == ww]
  data.table(w = ww, n_win = nrow(d),
             pct_a_zero     = 100 * mean(d$a == 0),
             pct_over_2.1Mb = 100 * mean(d$d95 > 2.1e6),
             pct_over_550kb = 100 * mean(d$d95 > 5.5e5),
             rho_all        = sc(d),
             rho_excl_a0    = sc(d[a > 0]),
             rho_cap_2.1Mb  = sc(d[d95 <= 2.1e6]),
             rho_cap_550kb  = sc(d[d95 <= 5.5e5])) }))
fwrite(S, file.path(OUT, "cap_cost.csv"))
print(S)
