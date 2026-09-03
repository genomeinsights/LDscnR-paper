## Does ORTHOGONALISING the surrogate against y prevent region locking?
## All locking numbers reported so far used the orthogonal (residualised) form.
suppressMessages({library(data.table); library(LDscnR)})
REMLE <- utils::getFromNamespace("emma.REMLE","LDscnR")
SIM <- "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5"; NSIM <- 80
run <- function(CELL, TAG, ENV) {
  P <- list(); yy <- NULL; GRM <- NULL; XY <- NULL; POP <- NULL
  for (i in 1:10) {
    x <- readRDS(sprintf("%s/adapt_%s_chr%d_%s_env%d.rds", SIM, TAG, i, CELL, ENV))
    pr <- ld_prune_and_eMLG(GTs=x$GTs, stage1=x$complexity_reduction$stage1, LD_decay=x$LD_decay,
          ld_w_col="ld_w_095", ld_w_threshold=0.025, score_threshold=0.80, min_r2_rho=0.5,
          distance_threshold=1e5, compute_unflagged_eMLG=TRUE, cores=1)
    E <- pr$eMLG; colnames(E) <- paste0(i,"_",colnames(E))
    if (is.null(yy)) { yy <- as.numeric(x$env$env); GRM <- x$GRM
                       XY <- as.matrix(x$env[,.(x,y)]); POP <- x$env$pop }
    P[[length(P)+1]] <- list(prep=emmax_setup(E, x$GRM), cl=colnames(E)) }
  CLall <- unlist(lapply(P,`[[`,"cl"))
  scan1 <- function(v) unlist(lapply(P, function(z) as.numeric(emmax_fast(z$prep, v))))
  n <- length(yy); obs <- CLall[p.adjust(scan1(yy),"BH") < 0.05]; ch <- length(obs)/length(CLall)
  nrm <- function(M) M/mean(diag(M))
  Dm <- as.matrix(dist(XY)); l <- median(Dm[lower.tri(Dm)])
  bases <- list(spatial=nrm(exp(-0.5*(Dm/l)^2)), pop=nrm(GRM))
  cat(sprintf("\n%s %s env%d | obs=%d | chance %.4f%%\n", CELL, TAG, ENV, length(obs), 100*ch))
  for (bn in names(bases)) for (rs in c(TRUE, FALSE)) {
    e <- eigen(bases[[bn]], symmetric=TRUE)
    A <- e$vectors %*% diag(sqrt(pmax(e$values,0)), n); set.seed(41)
    NE <- A %*% matrix(rnorm(n*NSIM), n, NSIM)
    if (rs) NE <- apply(NE, 2, function(v) as.numeric(resid(lm(v ~ yy))))
    hits <- unlist(lapply(seq_len(NSIM), function(b) CLall[p.adjust(scan1(NE[,b]),"BH") < 0.05]))
    cat(sprintf("  %-8s orthogonal=%-5s  hits %5d, %4d distinct | overlap %.2f%% = %.0fx chance\n",
        bn, rs, length(hits), uniqueN(hits),
        100*mean(hits %in% obs), if (length(hits)) mean(hits %in% obs)/ch else NA)) }
}
run("V0.5_c1","nobgs",1); run("V0.5_c2","nobgs",1)
