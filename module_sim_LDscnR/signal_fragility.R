## =============================================================================
## module_sim_LDscnR / signal_fragility.R
##
## Why does 10% covariate noise remove every discovery? Two candidate causes:
## the SIGNAL dies, or BH's THRESHOLD tightens because most of the small
## p-values were confounding and noise removes them. Answered by tracking the
## p-values of QTN-tagging units directly, independent of BH.
##
##   h2     BH discoveries   BH threshold   best TRUE-unit p   TRUE units p<1e-4
##   1.00        16            1.06e-05        3.36e-08              5
##   0.99         1            3.21e-07        3.51e-06              4
##   0.95         0              --            8.45e-05              2
##   0.90         0              --            2.23e-04              0
##   0.70         0              --            1.54e-03              0
##
## THE SIGNAL DIES, AND ABSURDLY FAST: ONE PERCENT covariate noise costs TWO
## ORDERS OF MAGNITUDE in the best true p-value. That cannot be a variance
## effect -- the QTN have r2 = 0.39-0.51 against the environment, far above the
## ~0.095 detectable at p < 1e-4 with n = 160.
##
## IT IS THE DEGENERATE MIXED-MODEL FIT. When the covariate lies in the span of
## the kinship (sigma_e^2 = 0), EMMAX whitens against a near-degenerate residual
## and the surviving p-values are extreme for that reason rather than because the
## association is strong. Give the model 1% of unstructured variance and the
## whitening changes completely. SO THE DETECTABLE SIGNAL AT h2 = 1 IS AN
## ARTEFACT OF THE SAME DEGENERACY THAT MAKES THE NULL ANTI-CONSERVATIVE -- one
## cause, and it produces both the discoveries and the miscalibration.
##
## Note also that the median TRUE-unit p is 0.17-0.26 at EVERY h2, and only 5 of
## 76 QTN-tagging units reach p < 1e-4 even at h2 = 1. Most of the truth was
## never detectable; recall of 0.20-0.25 across the grid is that fact.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/signal_fragility.R
## =============================================================================
suppressMessages({library(data.table); library(LDscnR)})
SIM <- "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5"
CELL <- "V0.5_c2"; TAG <- "nobgs"; ENV <- 1
P <- list(); y0 <- NULL; GRM <- NULL; lnk <- list(); CT <- list()
for (i in 1:10) {
  x <- readRDS(sprintf("%s/adapt_%s_chr%d_%s_env%d.rds", SIM, TAG, i, CELL, ENV))
  m <- flag_true_qtns(as.data.table(x$map))
  pr <- ld_prune_and_eMLG(GTs=x$GTs, stage1=x$complexity_reduction$stage1, LD_decay=x$LD_decay,
        ld_w_col="ld_w_095", ld_w_threshold=0.025, score_threshold=0.80, min_r2_rho=0.5,
        distance_threshold=1e5, compute_unflagged_eMLG=TRUE, cores=1)
  E <- pr$eMLG; g <- as.data.table(pr$groups)
  ms <- data.table(marker = unlist(g$members, use.names = FALSE),
            CL = paste0(i,"_",rep.int(g$group_id, lengths(g$members))))
  mm <- merge(m, ms, by="marker", all.x=TRUE)[!is.na(CL)]
  th <- score_thresholds(as.data.table(x$LD_decay$decay_sum), rho_r2=0.75, rho_d=0.95, dmax_cap=1e5)
  drv <- mm[true_pos_QTN %in% TRUE]
  if (nrow(drv)) lnk[[length(lnk)+1]] <- rbindlist(lapply(seq_len(nrow(drv)), function(j) {
    near <- mm[as.character(Chr)==as.character(drv$Chr[j]) & abs(Pos-drv$Pos[j]) < th$dmax]
    if (!nrow(near)) return(NULL)
    r2 <- suppressWarnings(cor(x$GTs[,drv$marker[j]], x$GTs[,near$marker], use="pairwise.complete.obs")^2)
    d <- data.table(CL=near$CL, r2=as.numeric(r2))[is.finite(r2) & r2>=th$r2min]
    if (!nrow(d)) return(NULL); d[, .(CL=paste0(i,"_",sub("^[0-9]+_","",CL[1]))), by=CL][, .(CL)] }))
  CT[[length(CT)+1]] <- data.table(CL=paste0(i,"_",g$group_id[g$group_id %in% colnames(E)]))
  colnames(E) <- paste0(i,"_",colnames(E))
  if (is.null(y0)) { y0 <- as.numeric(x$env$env); GRM <- x$GRM }
  P[[length(P)+1]] <- list(prep=emmax_setup(E, x$GRM), cl=colnames(E))
}
CLall <- unlist(lapply(P,`[[`,"cl")); true <- CLall %in% unique(rbindlist(lnk)$CL)
scan1 <- function(v) unlist(lapply(P, function(z) as.numeric(emmax_fast(z$prep, v))))
n <- length(y0)
cat(sprintf("%s | %d units, %d tag a QTN\n\n", CELL, length(CLall), sum(true)))
for (h2 in c(1, 0.99, 0.95, 0.9, 0.7)) {
  set.seed(97); e <- rnorm(n); e <- e/sd(e)*sd(y0)*sqrt((1-h2)/h2)
  yy <- if (h2>=1) y0 else y0+e
  p <- scan1(yy); q <- p.adjust(p, "BH")
  cat(sprintf("h2=%.2f | BH discoveries %2d | BH p-threshold %.2e | best TRUE-unit p %.2e | median TRUE-unit p %.3f | TRUE units p<1e-4: %d\n",
      h2, sum(q<0.05), max(p[q<0.05], -Inf), min(p[true], na.rm=TRUE),
      median(p[true], na.rm=TRUE), sum(p[true] < 1e-4, na.rm=TRUE)))
}
