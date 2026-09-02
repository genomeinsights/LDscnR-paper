## =============================================================================
## module_sim_LDscnR / null_tail_calibration.R
##
## PK: try the null WITHOUT residualising, and use the ground truth to say which
## nulls are anti-conservative and which conservative.
##
## The ground truth is exact and needs no convention: neutral-chromosome units
## hold no QTN, so their p-values ARE the null distribution under the real
## phenotype. A valid surrogate must reproduce it. Percentages below 1e-4, where
## BH actually operates (a uniform null gives 0.0100):
##
##                          V0.5_c1        V0.5_c2
##   OBSERVED (truth)        0.0131         0.0576     <- 1.3x and 5.8x uniform
##   vc      resid=TRUE      0.0086         0.0123
##   vc      resid=FALSE     0.0077         0.0137
##   pop     resid=TRUE      0.0102         0.0283
##   pop     resid=FALSE     0.0107         0.0298
##   spatial resid=TRUE      0.0182         0.1569
##   spatial resid=FALSE     0.0283         0.3846
##
## THE CLASSIFICATION PK ASKED FOR, on the cell that needs a null (V0.5_c2):
##   ANTI-CONSERVATIVE   vc (4.2-4.7x too light), pop (1.9-2.0x too light)
##   CONSERVATIVE        spatial (2.7x, and 6.7x without residualisation) -- and
##                       invalid anyway, being region-locked at 500x chance.
##
## DROPPING RESIDUALISATION IS DIRECTIONALLY RIGHT AND QUANTITATIVELY INSUFFICIENT.
## It moves every basis toward a heavier tail, which is the direction needed:
## +11% for vc, +5% for pop, +145% for spatial. For the two VALID bases that
## leaves them still 2-4x too light, so it is not the fix.
##
## THE WORST FEATURE OF THIS RESULT IS THE PATTERN, NOT THE MAGNITUDE. On
## V0.5_c1 the observed tail is 0.0131 and vc/pop give 0.0077-0.0107 -- close,
## and that cell's BH is already well calibrated (its neutral tail sits BELOW
## uniform). On V0.5_c2, where BH fails, the surrogates miss by 2-5x. THE NULLS
## AGREE WITH THE DATA EXACTLY WHERE THEY ARE NOT NEEDED AND FAIL WHERE THEY ARE.
##
## The truth sits BETWEEN pop and spatial on the failing cell, which suggests the
## confounding carries both a genetic and a spatial component and no single-basis
## surrogate spans it -- the sim environment is a spatial gradient acting on a
## spatially structured population.
##
## A DIAGNOSTIC-SAMPLING LESSON worth keeping: this was first run on V0.5_c1,
## whose observed tail is BELOW uniform, so there was nothing to reproduce and
## every basis looked fine. Choosing a diagnostic panel without first checking it
## exhibits the problem wasted the run.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/null_tail_calibration.R    (env: CELL, TAG, ENV, NSIM)
## =============================================================================
suppressMessages({library(data.table); library(LDscnR)})
emma.REMLE <- utils::getFromNamespace("emma.REMLE", "LDscnR")
SIM <- "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5"
CELL <- Sys.getenv("CELL", "V0.5_c2"); TAG <- Sys.getenv("TAG", "nobgs")
ENV  <- as.integer(Sys.getenv("ENV", "1")); NSIM <- as.integer(Sys.getenv("NSIM", "60"))
P <- list(); yy <- NULL; GRM <- NULL; POP <- NULL; XY <- NULL; CT <- list()
for (i in 1:10) {
  x <- readRDS(sprintf("%s/adapt_%s_chr%d_%s_env%d.rds", SIM, TAG, i, CELL, ENV))
  m <- as.data.table(x$map)
  pr <- ld_prune_and_eMLG(GTs=x$GTs, stage1=x$complexity_reduction$stage1, LD_decay=x$LD_decay,
        ld_w_col="ld_w_095", ld_w_threshold=0.025, score_threshold=0.80, min_r2_rho=0.5,
        distance_threshold=1e5, compute_unflagged_eMLG=TRUE, cores=1)
  E <- pr$eMLG; g <- as.data.table(pr$groups)
  ct <- unique(m[, .(Chr=as.character(Chr), chr_type)])
  gk <- g[g$group_id %in% colnames(E)]
  cl <- data.table(CL=paste0(i,"_",gk$group_id), Chr=as.character(gk$Chr))
  CT[[length(CT)+1]] <- merge(cl, ct, by="Chr", all.x=TRUE)
  colnames(E) <- paste0(i,"_",colnames(E))
  if (is.null(yy)) { yy <- as.numeric(x$env$env); GRM <- x$GRM; POP <- x$env$pop
                     XY <- as.matrix(x$env[, .(x,y)]) }
  P[[length(P)+1]] <- list(prep=emmax_setup(E, x$GRM), cl=colnames(E))
}
CLall <- unlist(lapply(P, `[[`, "cl")); CTa <- rbindlist(CT)
ntrl <- CLall %in% CTa[chr_type=="ntrl", CL]
scan1 <- function(v) unlist(lapply(P, function(z) as.numeric(emmax_fast(z$prep, v))))
tail_of <- function(p) c(pct01=100*mean(p<0.01, na.rm=TRUE), pct001=100*mean(p<0.001, na.rm=TRUE),
                         pct0001=100*mean(p<1e-4, na.rm=TRUE))
po <- scan1(yy)
cat(sprintf("OBSERVED, neutral units (n=%d): %%<0.01 %.3f | %%<0.001 %.3f | %%<1e-4 %.4f\n",
            sum(ntrl), tail_of(po[ntrl])[1], tail_of(po[ntrl])[2], tail_of(po[ntrl])[3]))
cat("  (a uniform null gives 1.000 / 0.100 / 0.0100)\n\n")
mk <- function(basis, resid) {
  n <- length(POP)
  if (basis=="vc") { re <- emma.REMLE(yy, P[[1]]$prep$Xo, P[[1]]$prep$Kn, eig.R=P[[1]]$prep$eigR)
    e <- eigen(re$vg*P[[1]]$prep$Kn + re$ve*diag(n), symmetric=TRUE)
  } else if (basis=="spatial") { Dm <- as.matrix(dist(XY)); l <- median(Dm[lower.tri(Dm)])
    e <- eigen(exp(-0.5*(Dm/l)^2), symmetric=TRUE)
  } else { pl <- sort(unique(POP)); np <- length(pl); idx <- lapply(pl, function(q) which(POP==q))
    Om <- matrix(0,np,np)
    for (a in seq_len(np)) for (b in a:np) { v <- mean(GRM[idx[[a]],idx[[b]]]); Om[a,b]<-v; Om[b,a]<-v }
    ee <- eigen(Om, symmetric=TRUE); L <- ee$vectors %*% diag(sqrt(pmax(ee$values,0)), np)
    set.seed(11); Z <- L %*% matrix(rnorm(np*NSIM), np, NSIM)
    o <- matrix(NA_real_, n, NSIM); for (a in seq_len(np)) o[idx[[a]],] <- rep(Z[a,], each=length(idx[[a]]))
    return(if (resid) apply(o,2,function(v) as.numeric(resid(lm(v ~ yy)))) else o) }
  L <- e$vectors %*% diag(sqrt(pmax(e$values,0)), n); set.seed(11)
  o <- L %*% matrix(rnorm(n*NSIM), n, NSIM)
  if (resid) apply(o,2,function(v) as.numeric(resid(lm(v ~ yy)))) else o
}
for (bs in c("vc","pop","spatial")) for (rs in c(TRUE,FALSE)) {
  NE <- mk(bs, rs)
  tt <- rowMeans(sapply(seq_len(NSIM), function(b) tail_of(scan1(NE[,b])[ntrl])))
  cat(sprintf("%-8s resid=%-5s  %%<0.01 %.3f | %%<0.001 %.3f | %%<1e-4 %.4f\n",
              bs, rs, tt[1], tt[2], tt[3]))
}
