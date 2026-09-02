## =============================================================================
## module_sim_LDscnR / kinship_rank_diagnostic.R
##
## WHY sigma_e^2 COLLAPSES, resolved against the panel session's contradicting
## result rather than left as "two datasets, two explanations".
##
## The disagreement: on bgs5, random population-constant covariates are NOT
## degenerate (vg_share 0.138, range 0.000-0.646); on the stickleback panel all
## 50 draws are (vg_share 1.000). The panel session attributed theirs to sampling
## geometry -- 35 populations, median 2 fish each, so the kinship encodes
## population identity. BUT bgs5 IS ALSO MEDIAN 2 PER POPULATION (80 populations,
## 160 individuals), so per-population replication cannot be what differs, and by
## a raw span argument bgs5 should be MORE degenerate (80/160 against 35/117),
## not less. It is the reverse.
##
## THE DISCRIMINATING QUANTITY IS THE KINSHIP'S EFFECTIVE RANK RELATIVE TO THE
## NUMBER OF POPULATIONS, which is set by how much WITHIN-population genetic
## variation there is:
##
##   bgs5: effective rank 122.0 of 160, against 80 populations   (ratio 1.53)
##         R^2 of off-diagonal entries on population-pair identity  0.708
##         within-population spread as % of between-population     170%
##
## A GRM that carries substantial within-population variation has dimensions no
## population-constant vector can reach, so such a covariate does not span it and
## REML retains a residual. A GRM at high F_ST collapses toward population
## identity -- effective rank approaching the population count -- and then ANY
## population-constant covariate spans it.
##
## So the portable statement is neither "the environment is the structure axis"
## (true on bgs5, false on the panel: their |cor(ecotype, PC1)| is 0.143 against
## our 0.52-0.86) nor "the covariate is population-constant" (true on the panel,
## false on bgs5). It is that SIGMA_E^2 COLLAPSES WHEN THE COVARIATE LIES IN THE
## SPAN OF A LOW-RANK KINSHIP, and the two datasets get there by different routes
## -- one because the covariate is the leading structure axis, the other because
## the kinship has few dimensions to lie outside of.
##
## Both routes are invisible in the symptom, which is what makes the sigma_e^2
## pre-flight check worth having: it fires either way without needing the cause.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/kinship_rank_diagnostic.R
## =============================================================================
suppressMessages({library(data.table); library(LDscnR)})
REMLE <- utils::getFromNamespace("emma.REMLE", "LDscnR")
EIGR  <- utils::getFromNamespace("emma.eigen.R.wo.Z", "LDscnR")
SIM   <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
CELLS <- strsplit(Sys.getenv("CELLS", "V0.5_c1,V0.5_c2,V1_c1.5,V2_c1"), ",")[[1]]

for (cl in CELLS) {
  x <- readRDS(sprintf("%s/adapt_nobgs_chr1_%s_env1.rds", SIM, cl))
  K <- x$GRM; n <- nrow(K); pop <- as.character(x$env$pop); np <- length(unique(pop))
  ev <- pmax(eigen(K, symmetric = TRUE)$values, 0)
  eff <- sum(ev)^2 / sum(ev^2)
  oi <- which(upper.tri(K), arr.ind = TRUE); v <- K[upper.tri(K)]
  same <- pop[oi[, 1]] == pop[oi[, 2]]
  pp <- paste(pmin(pop[oi[,1]], pop[oi[,2]]), pmax(pop[oi[,1]], pop[oi[,2]]))
  Xo <- matrix(1, n, 1); Kn <- K / mean(diag(K)); eg <- EIGR(Kn, Xo)
  f <- function(y) { r <- REMLE(y, Xo, Kn, eig.R = eg); r$vg / (r$vg + r$ve) }
  set.seed(9); pl <- unique(pop)
  pc <- replicate(50, { vv <- stats::setNames(stats::rnorm(length(pl)), pl); as.numeric(vv[pop]) })
  cat(sprintf("%-9s n=%d pops=%d (median %g/pop) | eff.rank %.1f (ratio to pops %.2f)\n",
              cl, n, np, stats::median(table(pop)), eff, eff / np))
  cat(sprintf("          R2(GRM ~ pop-pair) %.3f | within/between spread %.0f%% | random pop-const vg_share %.3f\n",
              summary(stats::lm(v ~ factor(pp)))$r.squared,
              100 * stats::sd(v[same]) / stats::sd(v[!same]), mean(apply(pc, 2, f))))
}
