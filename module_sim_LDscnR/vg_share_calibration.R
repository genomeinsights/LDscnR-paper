## =============================================================================
## module_sim_LDscnR / vg_share_calibration.R
##
## Is the MVN null's miscalibration PREDICTABLE from the REML variance split
## alone? A user can read sigma_e^2 before running any null, so if the split
## predicts the direction of the error, that is more useful than a null which is
## calibrated only in a regime neither of our datasets occupies.
##
## REPLICATED OVER 9 PANELS, AND IT CORRECTS A SINGLE-PANEL NUMBER I HAD ALREADY
## CIRCULATED. Median ratio of surrogate tail to observed tail on neutral
## chromosomes (1.0 = calibrated, below 1 = anti-conservative):
##
##   h2     vg_share   median ratio   panels anti-conservative
##   1.0      1.00         0.37              8 of 9
##   0.9      0.95         0.67              7 of 9
##   0.7      0.87         1.25              2 of 9
##   0.5      0.73         1.42              2 of 9
##   0.3      0.56         1.11              3 of 9
##
## THE DIRECTION REPLICATES AND THE MAGNITUDE I QUOTED DID NOT. I reported 0.17x
## (5.9x too light) from V0.5_c2 env1. Over nine panels the median is 0.37x, and
## the per-panel range at h2 = 1 runs 0.05x to 1.19x. V0.5_c2 is the extreme
## cell, and the ratio-of-pooled-tails (5.45x) is dominated by it. The honest
## statement is ANTI-CONSERVATIVE IN 8 OF 9 PANELS, MEDIAN 2.7x TOO LIGHT.
##
## TWO STATISTICAL TRAPS THIS RUN WALKED INTO, both worth keeping. The MEAN of
## per-panel ratios is 1.5e6, because one panel at h2 = 0.9 has an observed tail
## of exactly 0 and the ratio diverges -- the median is the only usable summary.
## And the ratio of means (5.45x) and the median of ratios (0.37x) point in
## OPPOSITE directions relative to 1, because the panels are heterogeneous; a
## single-number summary of this quantity has to say which it is.
##
## WHAT IS AND IS NOT ESTABLISHED. Established: when sigma_e^2 is near zero the
## null is anti-conservative, consistently enough to license the one-sided
## reading. NOT established: that we have a calibrated null. Calibration appears
## only at h2 <= 0.7, and BOTH datasets sit at h2 = 1 -- the sims because the
## environment is a deterministic spatial function, the stickleback panel because
## a population-constant covariate lies entirely in the kinship span. The
## pathological regime is the one our data occupy.
suppressMessages({library(data.table); library(LDscnR); library(parallel)})
REMLE <- utils::getFromNamespace("emma.REMLE","LDscnR"); EIGR <- utils::getFromNamespace("emma.eigen.R.wo.Z","LDscnR")
SIM <- "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5"; NSIM <- 40
H2 <- c(1, 0.9, 0.7, 0.5, 0.3)
one <- function(CELL, TAG, ENV) {
  P <- list(); y0 <- NULL; GRM <- NULL; CT <- list()
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
    if (is.null(y0)) { y0 <- as.numeric(x$env$env); GRM <- x$GRM }
    P[[length(P)+1]] <- list(prep=emmax_setup(E, x$GRM), cl=colnames(E)) }
  CLall <- unlist(lapply(P,`[[`,"cl")); ntrl <- CLall %in% rbindlist(CT)[chr_type=="ntrl", CL]
  scan1 <- function(v) unlist(lapply(P, function(z) as.numeric(emmax_fast(z$prep, v))))
  tf <- function(p) 100*mean(p < 1e-4, na.rm=TRUE)
  n <- length(y0); Xo <- matrix(1,n,1); K <- GRM/mean(diag(GRM))
  rbindlist(lapply(H2, function(h2) {
    set.seed(97); e <- rnorm(n); e <- e/sd(e)*sd(y0)*sqrt((1-h2)/h2)
    yy <- if (h2 >= 1) y0 else y0 + e
    ot <- tf(scan1(yy)[ntrl])
    re <- REMLE(yy, Xo, K, eig.R=EIGR(K, Xo))
    Sig <- re$vg*K + re$ve*diag(n); ee <- eigen(Sig, symmetric=TRUE)
    A <- ee$vectors %*% diag(sqrt(pmax(ee$values,0)), n); set.seed(98)
    NE <- A %*% matrix(rnorm(n*NSIM), n, NSIM)
    st <- mean(sapply(seq_len(NSIM), function(b) tf(scan1(NE[,b])[ntrl])))
    data.table(cell=CELL, tag=TAG, env=ENV, h2=h2,
               vg_share=re$vg/(re$vg+re$ve), obs=ot, surr=st, ratio=st/max(ot,1e-9)) })) }
g <- CJ(cell=c("V0.5_c1","V0.5_c2","V2_c1"), tag="nobgs", env=1:3, sorted=FALSE)
R <- rbindlist(mclapply(seq_len(nrow(g)), function(k)
       tryCatch(one(g$cell[k], g$tag[k], g$env[k]), error=function(e) NULL), mc.cores=5))
fwrite(R, "module_sim_LDscnR/results/structure_null/vg_share_calibration.csv")
cat("== does the REML variance split predict the null's miscalibration? (9 panels)\n\n")
print(R[, .(panels=.N, vg_share=round(mean(vg_share),2), obs_tail=round(mean(obs),4),
            surr_tail=round(mean(surr),4), ratio=round(mean(ratio),2),
            se=round(sd(ratio)/sqrt(.N),2)), by=h2][order(-h2)])
cat(sprintf("\nSpearman(vg_share, ratio) = %.3f over %d panel-h2 cells\n",
            cor(R$vg_share, R$ratio, method="spearman"), nrow(R)))
cat(sprintf("at vg_share > 0.99: ratio %.2f (n=%d) | at vg_share < 0.7: ratio %.2f (n=%d)\n",
    mean(R[vg_share>0.99]$ratio), sum(R$vg_share>0.99),
    mean(R[vg_share<0.7]$ratio), sum(R$vg_share<0.7)))
