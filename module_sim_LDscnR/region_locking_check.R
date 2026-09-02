## =============================================================================
## module_sim_LDscnR / region_locking_check.R
##
## REGION-LOCKING CHECK (PK's failure mode, relayed by the panel session): does
## the null discover the SAME units as the observed data? Counting E[V] cannot
## detect this -- only identity can. If surrogate discoveries concentrate on the
## units the observed scan finds, the statistic is tracking cluster structure
## rather than the phenotype, and the null is not a null.
##
## RESULT: THE SPATIAL BASIS FAILS AND ITS NUMBERS ARE WITHDRAWN.
##
##   basis     surrogate hits   distinct units   also an OBSERVED discovery
##   spatial   304 / 52 draws        88                  10.5%
##   vc          8 /  6 draws         8                   0.0%
##   pop         2 /  2 draws         2                   0.0%
##
## Chance overlap is 0.021%, so the spatial basis is 500x enriched for
## rediscovering the observed discoveries, one unit is hit in 27 of 100 draws,
## and five units carry 26% of all hits. Its high E[V] -- 65 to 160 over 40
## panels, which had been read as the basis being appropriately stringent -- is
## REGION LOCKING, not stringency. The estimated FDP of 1.4-3.2 computed from it
## is meaningless and is withdrawn.
##
## WHAT REMAINS IS WORSE THAN A CORRECTION. The two bases that pass the identity
## check are both far too weak: E[V] ~ 0.03, an estimated FDP of ~0.001 against a
## neutral-chromosome FDP of ~0.70. So NO SURROGATE BASIS TESTED SO FAR PRODUCES
## A CALIBRATED FDR -- one is invalid and two are anti-conservative by two orders
## of magnitude. The permutation-null direction is not yet working.
##
## A likely cause worth testing next: the surrogates are residualised against the
## observed phenotype (following structured_null()), which makes them orthogonal
## to it rather than a draw from the null distribution of phenotypes. The real
## null has a heavier tail than these surrogates produce -- 0.185% of neutral
## units below p = 0.001 against the 0.100% a uniform null gives.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/region_locking_check.R
## =============================================================================
suppressMessages({library(data.table); library(LDscnR)})
emma.REMLE <- utils::getFromNamespace("emma.REMLE", "LDscnR")
SIM <- "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5"
CELL <- "V0.5_c1"; TAG <- "nobgs"; ENV <- 1; NSIM <- 100
source_draw <- function(GRM, pop, y, nsim, seed, basis, coords, prep) {
  n <- length(pop)
  if (basis == "spatial") {
    Dm <- as.matrix(stats::dist(coords)); l <- stats::median(Dm[lower.tri(Dm)])
    e <- eigen(exp(-0.5*(Dm/l)^2), symmetric=TRUE)
  } else if (basis == "vc") {
    re <- emma.REMLE(y, prep$Xo, prep$Kn, eig.R=prep$eigR)
    e <- eigen(re$vg*prep$Kn + re$ve*diag(n), symmetric=TRUE)
  } else {
    pl <- sort(unique(pop)); np <- length(pl); idx <- lapply(pl, function(p) which(pop==p))
    Om <- matrix(0,np,np)
    for (a in seq_len(np)) for (b in a:np) { v <- mean(GRM[idx[[a]],idx[[b]]]); Om[a,b]<-v; Om[b,a]<-v }
    ee <- eigen(Om, symmetric=TRUE); L <- ee$vectors %*% diag(sqrt(pmax(ee$values,0)), np)
    set.seed(seed); Z <- L %*% matrix(rnorm(np*nsim), np, nsim)
    out <- matrix(NA_real_, n, nsim)
    for (a in seq_len(np)) out[idx[[a]],] <- rep(Z[a,], each=length(idx[[a]]))
    return(apply(out, 2, function(v) as.numeric(resid(lm(v ~ y)))))
  }
  L <- e$vectors %*% diag(sqrt(pmax(e$values,0)), n)
  set.seed(seed); out <- L %*% matrix(rnorm(n*nsim), n, nsim)
  apply(out, 2, function(v) as.numeric(resid(lm(v ~ y))))
}
P <- list(); yy <- NULL; GRM <- NULL; POP <- NULL; XY <- NULL
for (i in 1:10) {
  x <- readRDS(sprintf("%s/adapt_%s_chr%d_%s_env%d.rds", SIM, TAG, i, CELL, ENV))
  pr <- ld_prune_and_eMLG(GTs=x$GTs, stage1=x$complexity_reduction$stage1, LD_decay=x$LD_decay,
        ld_w_col="ld_w_095", ld_w_threshold=0.025, score_threshold=0.80, min_r2_rho=0.5,
        distance_threshold=1e5, compute_unflagged_eMLG=TRUE, cores=1)
  E <- pr$eMLG; colnames(E) <- paste0(i,"_",colnames(E))
  if (is.null(yy)) { yy <- as.numeric(x$env$env); GRM <- x$GRM; POP <- x$env$pop
                     XY <- as.matrix(x$env[, .(x,y)]) }
  P[[length(P)+1]] <- list(prep=emmax_setup(E, x$GRM), cl=colnames(E))
}
CLall <- unlist(lapply(P, `[[`, "cl"))
scan1 <- function(v) unlist(lapply(P, function(z) as.numeric(emmax_fast(z$prep, v))))
p_obs <- scan1(yy); obs_sig <- CLall[p.adjust(p_obs,"BH") < 0.05]
for (bs in c("spatial","vc","pop")) {
  NE <- source_draw(GRM, POP, yy, NSIM, 7L, bs, XY, P[[1]]$prep)
  sets <- lapply(seq_len(NSIM), function(b) CLall[p.adjust(scan1(NE[,b]),"BH") < 0.05])
  allhits <- unlist(sets); nb <- sum(lengths(sets) > 0)
  ov <- if (length(allhits)) mean(allhits %in% obs_sig) else NA_real_
  cat(sprintf("%-8s obs=%d | surrogate hits: %d total across %d/%d draws, %d DISTINCT units\n",
      bs, length(obs_sig), length(allhits), nb, NSIM, uniqueN(allhits)))
  cat(sprintf("         %% of surrogate hits that are ALSO observed discoveries: %.1f%% (chance %.4f%%)\n",
      100*ov, 100*length(obs_sig)/length(CLall)))
  if (length(allhits)) { tb <- sort(table(allhits), decreasing=TRUE)
    cat(sprintf("         most-repeated unit hit in %d of %d draws | top-5 share %.1f%%\n",
        as.integer(tb[1]), NSIM, 100*sum(head(tb,5))/length(allhits))) }
}
