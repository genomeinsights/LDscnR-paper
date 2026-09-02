suppressMessages({library(data.table); library(LDscnR)})
REMLE <- utils::getFromNamespace("emma.REMLE","LDscnR"); EIGR <- utils::getFromNamespace("emma.eigen.R.wo.Z","LDscnR")
SIM <- "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5"; NSIM <- 100
run <- function(CELL, TAG, ENV, lams) {
  P <- list(); yy <- NULL; GRM <- NULL; XY <- NULL
  for (i in 1:10) {
    x <- readRDS(sprintf("%s/adapt_%s_chr%d_%s_env%d.rds", SIM, TAG, i, CELL, ENV))
    pr <- ld_prune_and_eMLG(GTs=x$GTs, stage1=x$complexity_reduction$stage1, LD_decay=x$LD_decay,
          ld_w_col="ld_w_095", ld_w_threshold=0.025, score_threshold=0.80, min_r2_rho=0.5,
          distance_threshold=1e5, compute_unflagged_eMLG=TRUE, cores=1)
    E <- pr$eMLG; colnames(E) <- paste0(i,"_",colnames(E))
    if (is.null(yy)) { yy <- as.numeric(x$env$env); GRM <- x$GRM; XY <- as.matrix(x$env[,.(x,y)]) }
    P[[length(P)+1]] <- list(prep=emmax_setup(E, x$GRM), cl=colnames(E)) }
  CLall <- unlist(lapply(P,`[[`,"cl"))
  scan1 <- function(v) unlist(lapply(P, function(z) as.numeric(emmax_fast(z$prep, v))))
  n <- length(yy); Xo <- matrix(1,n,1); nrm <- function(M) M/mean(diag(M))
  Kg <- nrm(GRM); Dm <- as.matrix(dist(XY)); l <- median(Dm[lower.tri(Dm)]); Ks <- nrm(exp(-0.5*(Dm/l)^2))
  obs <- CLall[p.adjust(scan1(yy),"BH") < 0.05]
  for (L in lams) {
    K <- L*Kg+(1-L)*Ks; re <- REMLE(yy,Xo,K,eig.R=EIGR(K,Xo))
    Sig <- re$vg*K + re$ve*diag(n); e <- eigen(Sig, symmetric=TRUE)
    A <- e$vectors %*% diag(sqrt(pmax(e$values,0)), n); set.seed(23)
    NE <- A %*% matrix(rnorm(n*NSIM), n, NSIM)
    sets <- lapply(seq_len(NSIM), function(b) CLall[p.adjust(scan1(NE[,b]),"BH") < 0.05])
    h <- unlist(sets); ch <- length(obs)/length(CLall)
    cat(sprintf("  %s l=%.2f | obs=%d | hits %d over %d draws, %d distinct | overlap %.2f%% vs chance %.4f%% = %.0fx",
        CELL, L, length(obs), length(h), sum(lengths(sets)>0), uniqueN(h),
        100*mean(h %in% obs), 100*ch, if (length(h)) mean(h %in% obs)/ch else NA))
    if (length(h)) { tb <- sort(table(h), decreasing=TRUE)
      cat(sprintf(" | top unit %d/%d draws\n", as.integer(tb[1]), NSIM)) } else cat("\n") }
}
cat("REGION-LOCKING GATE for the mixed basis (spatial component was locked at 500x)\n")
run("V0.5_c2","nobgs",1, c(0, 0.05, 0.25, 1))
run("V0.5_c1","nobgs",1, c(0.05))
