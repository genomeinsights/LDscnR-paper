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
## EXCLUSION IS NOT CAPPING, AND IT OVERSTATES THE COST BY ~5x. A cap is
## d <- min(d_derived, cap) with the window RETAINED, so the windows above it are
## tied at the ceiling rather than removed. Measuring the real operation (sweep
## below, w = 20, higher is better):
##
##   cap      250 kb  550 kb  1 Mb   2.1 Mb  3.7 Mb  10 Mb  30 Mb   none
##   at ceil   89.9%   65.1%  45.7%   32.3%   26.2%  17.5%  10.8%   0.0%
##   rho       0.381   0.656  0.755   0.793   0.803  0.826  0.829   0.830
##
## Capping at the block q90 costs 0.830 -> 0.656, not the 0.830 -> 0.427 that
## exclusion implied. So my retraction of the cap was itself premature, and it
## was premature on evidence I had already labelled in the same commit as not
## measuring the operation in question.
##
## THE COST IS MONOTONE IN THE CAP, SO THERE IS NO OPTIMUM AND THE RULE IS
## "GENEROUS". Every cap degrades the correlation and tighter caps degrade it
## more -- no cap wins on this criterion. But a 10 Mb ceiling costs 0.004 while
## bounding a quantity that otherwise reaches 4,961 Mb on a 30 Mb chromosome.
## THE CAP BELONGS WHERE THE QUANTITY STOPS BEING PHYSICALLY MEANINGFUL, NOT
## WHERE THE FEATURES STOP. Dimensioning it from the block-size distribution --
## my original proposal -- puts it an order of magnitude too tight and costs
## 0.174 to buy nothing the generous cap does not already buy.
##
## On the panel the same operation is INERT (2c: moves the correlation by at
## most 0.005 at any cap, against 0.13 for exclusion), because their divergent
## windows are 0.6-1.8% of the total rather than a third. So the two datasets
## disagree about whether capping MATTERS, not about what to cap.
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

## The sweep behind the header. Note the baseline: "none" retains a == 0 windows
## at d = Inf, which Spearman ranks at the top -- so it must be compared against
## capped series that also retain them, not against the a > 0 subset above.
caps <- c(2.5e5, 5.5e5, 1e6, 2.1e6, 3.7e6, 1e7, 3e7, Inf)
SW <- rbindlist(lapply(sort(unique(W$w)), function(ww) rbindlist(lapply(caps, function(cp) {
  d <- copy(W[w == ww]); d[, dc := pmin(d95, cp)]
  r <- d[, .(rho = suppressWarnings(stats::cor(dc, cMMb, method = "spearman"))),
         by = .(file, Chr)][is.finite(rho)]
  data.table(w = ww, cap_kb = cp / 1e3, pct_at_ceiling = 100 * mean(d$d95 > cp),
             rho = -median(r$rho)) }))))
fwrite(SW, file.path(OUT, "cap_sweep.csv"))
cat("\ncap sweep, min(d, cap) with windows retained (rho sign-flipped, higher better):\n")
print(dcast(SW, cap_kb ~ w, value.var = "rho"))
