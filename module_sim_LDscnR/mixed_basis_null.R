## =============================================================================
## module_sim_LDscnR / mixed_basis_null.R
##
## PK: try a mixed genetic + spatial basis. K(lambda) = lambda*Kg + (1-lambda)*Ks,
## both scaled to mean diagonal 1, surrogates drawn from vg*K(lambda) + ve*I with
## REML components. lambda is FITTED by profiling the REML likelihood, so the
## observed phenotype chooses the mixture rather than the analyst.
##
## THE MIXTURE CALIBRATES THE TAIL, AND THE CALIBRATING WEIGHT IS STABLE.
## % of neutral units below 1e-4, surrogate / observed:
##
##   panel                observed   l=0     l=0.02  l=0.05  l=0.25   REML picks
##   V0.5_c2 nobgs env1    0.0576    5.48x   2.30x   1.31x   0.38x    0.02
##   V0.5_c2 bgs   env2    0.1594    3.44x   1.20x   0.69x   0.15x    0.02
##   V1_c1.5 nobgs env1    0.0251    3.62x   1.83x   1.00x   0.40x    0.00
##
## lambda ~ 0.05 lands within 1.5x of truth on all three, against 2-5x too light
## for pure genetic and 3.4-5.5x too heavy for pure spatial. REML independently
## picks 0.00-0.02, the right region and slightly too spatial. That is a
## non-circular selection rule landing near the calibrating value.
##
## AND IT IS ALL UNUSABLE, BECAUSE THE MIXTURE FAILS THE REGION-LOCKING GATE.
## Overlap between surrogate discoveries and OBSERVED discoveries, against chance:
##
##   lambda   1.00     0.25     0.05             0.00
##   V0.5_c2   0x       67x      157x             68x
##   V0.5_c1    -         -     2218x (46.7%)      -
##
## Only PURE GENETIC (lambda = 1) passes -- 0% overlap, 16 distinct units from 16
## hits -- and pure genetic is the basis that is 2-5x too light.
##
## THE DILEMMA, WHICH IS THE ACTUAL FINDING. A surrogate has power against the
## real confounding only if it RESEMBLES the observed covariate; but resembling
## it makes the surrogate a REPLICATE rather than a null, and it rediscovers the
## same units. Here the environment is a spatial gradient, so spatial surrogates
## are alternative spatial gradients and hit the same clusters -- on V0.5_c1,
## 46.7% of all surrogate discoveries are units the observed scan also found.
## Genetic surrogates are unlike the environment and so discover nothing. No
## mixture escapes this: validity and power move in opposite directions along
## lambda. That is a mechanism for the failure mode PK had already hit.
##
## ONE CAVEAT THAT MAY MAKE THIS A WORST CASE. REML puts vg/(vg+ve) = 1.00 at
## every lambda: the simulated environment has NO residual variance, being a
## deterministic spatial function. A real environment carries measurement noise
## and idiosyncratic variation, so a surrogate can resemble its structure without
## being a replicate of it. The dilemma may be sharper here than on real data,
## and testing that needs an environment with a residual component.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/mixed_basis_null.R
## =============================================================================
NSIM <- 40
lams <- c(0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25)
run <- function(CELL, TAG, ENV) {
  P <- list(); yy <- NULL; GRM <- NULL; XY <- NULL; CT <- list()
  for (i in 1:10) {
    x <- readRDS(sprintf("%s/adapt_%s_chr%d_%s_env%d.rds", SIM, TAG, i, CELL, ENV))
    m <- as.data.table(x$map)
    pr <- ld_prune_and_eMLG(GTs=x$GTs, stage1=x$complexity_reduction$stage1, LD_decay=x$LD_decay,
          ld_w_col="ld_w_095", ld_w_threshold=0.025, score_threshold=0.80, min_r2_rho=0.5,
          distance_threshold=1e5, compute_unflagged_eMLG=TRUE, cores=1)
    E <- pr$eMLG; g <- as.data.table(pr$groups); gk <- g[g$group_id %in% colnames(E)]
    CT[[length(CT)+1]] <- merge(data.table(CL=paste0(i,"_",gk$group_id), Chr=as.character(gk$Chr)),
                                unique(m[,.(Chr=as.character(Chr), chr_type)]), by="Chr")
    colnames(E) <- paste0(i,"_",colnames(E))
    if (is.null(yy)) { yy <- as.numeric(x$env$env); GRM <- x$GRM; XY <- as.matrix(x$env[,.(x,y)]) }
    P[[length(P)+1]] <- list(prep=emmax_setup(E, x$GRM), cl=colnames(E))
  }
  CLall <- unlist(lapply(P,`[[`,"cl")); ntrl <- CLall %in% rbindlist(CT)[chr_type=="ntrl", CL]
  scan1 <- function(v) unlist(lapply(P, function(z) as.numeric(emmax_fast(z$prep, v))))
  n <- length(yy); Xo <- matrix(1,n,1)
  nrm <- function(M) M/mean(diag(M)); Kg <- nrm(GRM)
  Dm <- as.matrix(dist(XY)); l <- median(Dm[lower.tri(Dm)]); Ks <- nrm(exp(-0.5*(Dm/l)^2))
  tf <- function(p) 100*mean(p < 1e-4, na.rm=TRUE)
  ot <- tf(scan1(yy)[ntrl])
  ll <- sapply(lams, function(L){K <- L*Kg+(1-L)*Ks; REMLE(yy,Xo,K,eig.R=EIGR(K,Xo))$REML})
  rt <- sapply(lams, function(L) {
    K <- L*Kg+(1-L)*Ks; re <- REMLE(yy,Xo,K,eig.R=EIGR(K,Xo))
    Sig <- re$vg*K + re$ve*diag(n); e <- eigen(Sig, symmetric=TRUE)
    A <- e$vectors %*% diag(sqrt(pmax(e$values,0)), n); set.seed(11)
    NE <- A %*% matrix(rnorm(n*NSIM), n, NSIM)
    mean(sapply(seq_len(NSIM), function(b) tf(scan1(NE[,b])[ntrl]))) })
  data.table(cell=CELL, tag=TAG, env=ENV, obs=ot, lam=lams, surr=rt, reml=ll)
}
R <- rbindlist(lapply(list(c("V0.5_c2","nobgs",1), c("V0.5_c2","bgs",2), c("V1_c1.5","nobgs",1)),
                      function(z) run(z[1], z[2], as.integer(z[3]))))
for (k in unique(paste(R$cell,R$tag,R$env))) {
  d <- R[paste(cell,tag,env)==k]
  cat(sprintf("\n%s | observed tail %.4f | REML picks lambda %.2f\n", k, d$obs[1], d$lam[which.max(d$reml)]))
  cat(sprintf("  %s\n", paste(sprintf("l=%.2f:%.4f(%.2fx)", d$lam, d$surr, d$surr/d$obs), collapse="  ")))
  i <- which.min(abs(log(d$surr/d$obs)))
  cat(sprintf("  CALIBRATING lambda = %.2f (ratio %.2fx)\n", d$lam[i], d$surr[i]/d$obs[1]))
}
