## =============================================================================
## module_sim_LDscnR / env_structure_alignment.R
##
## WHY sigma_e^2 = 0, and why no sampling change fixes it.
##
## The obvious explanation -- the covariate is population-constant and 80
## populations nearly span a 160-individual kinship -- IS WRONG, and the test is
## one line: fit REML to RANDOM population-constant vectors.
##
##   observed environment            vg_share 1.000
##   random population-constant      vg_share 0.181  (range 0.000 - 0.534)
##   random individual-level         vg_share 0.176  (range 0.000 - 0.708)
##
## A random population-constant covariate is not degenerate. Only the actual
## environment is. So it is not the sampling design, and changing
## individuals-per-patch or patch count cannot help.
##
## THE CAUSE IS ISOLATION BY ADAPTATION. The environment is not merely
## correlated with structure; it IS the leading axis of it, because the
## populations are locally adapted to it and genome-wide allele frequencies
## track it:
##
##   cell        |cor(env, PC1)|   r2 from PC1   from PC1-5   from PC1-20
##   V0.5_c1          0.706           0.499        0.870        0.886
##   V0.5_c2          0.857           0.734        0.896        0.943
##   V1_c1.5          0.522           0.272        0.884        0.925
##   V2_c1            0.796           0.633        0.832        0.874
##
## The environment lies almost entirely inside the leading ~20 dimensions of
## genetic structure. REML has nothing left to call residual. This is the
## fundamental confound of landscape genomics rather than a simulation defect --
## which is an argument for simulating it, but not at a strength that makes the
## mixed model degenerate.
##
## CONSEQUENCE FOR THE DESIGN QUESTION. Signal cannot be raised by sampling
## without raising power, which defeats the purpose of a marginal-signal
## benchmark. It can only be raised by REDUCING THE COLLINEARITY between the
## environment and the dominant structure axis, which is a demographic change:
## an environment that varies on an axis other than the main demographic one, or
## a mosaic rather than a smooth cline, so migration mixes environmentally
## distinct neighbours and structure stops tracking the environment.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/env_structure_alignment.R
## =============================================================================
suppressMessages({library(data.table); library(LDscnR)})
REMLE <- utils::getFromNamespace("emma.REMLE", "LDscnR")
EIGR  <- utils::getFromNamespace("emma.eigen.R.wo.Z", "LDscnR")
SIM   <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
CELLS <- strsplit(Sys.getenv("CELLS", "V0.5_c1,V0.5_c2,V1_c1.5,V2_c1"), ",")[[1]]

for (cl in CELLS) {
  x <- readRDS(sprintf("%s/adapt_nobgs_chr1_%s_env1.rds", SIM, cl))
  K <- x$GRM / mean(diag(x$GRM)); n <- nrow(K); Xo <- matrix(1, n, 1)
  eg <- EIGR(K, Xo); env <- as.numeric(x$env$env); pop <- x$env$pop
  vs <- function(y) { r <- REMLE(y, Xo, K, eig.R = eg); r$vg / (r$vg + r$ve) }
  set.seed(5); pl <- unique(pop)
  pc0 <- replicate(20, { v <- stats::setNames(stats::rnorm(length(pl)), pl)
                         as.numeric(v[as.character(pop)]) })
  ev <- eigen(x$GRM, symmetric = TRUE)
  cat(sprintf("%-9s vg_share: env %.3f | random pop-constant %.3f | random individual %.3f\n",
              cl, vs(env), mean(apply(pc0, 2, vs)),
              mean(apply(replicate(20, stats::rnorm(n)), 2, vs))))
  cat(sprintf("          |cor(env,PC1)| %.3f | env r2 from PC1 %.3f, PC1-5 %.3f, PC1-20 %.3f\n",
              abs(stats::cor(env, ev$vectors[, 1])),
              summary(stats::lm(env ~ ev$vectors[, 1]))$r.squared,
              summary(stats::lm(env ~ ev$vectors[, 1:5]))$r.squared,
              summary(stats::lm(env ~ ev$vectors[, 1:20]))$r.squared))
}
