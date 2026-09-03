## =============================================================================
## module_sim_LDscnR / span_matched_precision.R
##
## I warned the panel session that EcoPeak overlap is span-sensitive -- a wider
## region satisfies it more easily -- so their equal-peak-rate claim needed a
## span-matched check. They then applied the same warning to MY assembly rule on
## their data, where it produces 452 kb median spans and a peak rate rising
## 39.4 -> 50.0 -> 63.6% in lockstep with span. So the warning is owed back:
## tagging counts a region true if ANY member is in LD with a QTN, and a wider
## region has more members and more chances.
##
## IT DOES NOT HOLD FOR THIS METRIC, and the shape is not what six panels showed.
## Precision by span quintile, 72 panels (8 of 80 gave no discoveries in one or
## both arms):
##
##   quintile   median span   precision      six-panel version
##   Q1            ~0            0.155            0.284
##   Q2            93 kb         0.311            0.357
##   Q3           854 kb         0.284            0.296
##   Q4           5.7 Mb         0.222            0.103
##   Q5          14.7 Mb         0.149            0.132
##
## AN INVERTED U. Both extremes are about half as precise as the middle, and the
## NARROWEST quintile -- effectively single-marker regions -- is as weak as the
## widest. The six-panel table showed only the upper fall and was written up as
## evidence that WIDE regions are the problem; the lower tail was simply too
## thinly populated to see. The finding is about EXTREMES, not width.
##
## Tagging is still not inflated by span: precision FALLS in the widest quintile,
## the opposite of an annotation-overlap rate, which a wide region satisfies by
## covering more sequence. That difference is why the two datasets disagree about
## the consensus statistic -- the panel's metric rewards width and this one
## penalises it.
##
## AND THE STAGE-2 ADVANTAGE STRENGTHENS UNDER MATCHING: five of five quintiles
## (0.173/0.142, 0.342/0.287, 0.313/0.263, 0.230/0.216, 0.161/0.141), and
## stage-2 regions are slightly NARROWER in median (788 against 946 kb) where the
## six-panel version had them wider. Overall precision 0.244 against 0.210. So
## the advantage rests on a reversal rather than a near-tie.
##
suppressMessages({library(data.table); library(LDscnR)})
SIM   <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
CELLS <- strsplit(Sys.getenv("CELLS", "V0.5_c1,V0.5_c2"), ",")[[1]]
TAGS  <- strsplit(Sys.getenv("TAGS", "nobgs"), ",")[[1]]
ENVS  <- as.integer(strsplit(Sys.getenv("ENVS", "1,2,3"), ",")[[1]])
OUTF  <- Sys.getenv("OUT", "")
simes <- function(p){p<-sort(p[is.finite(p)]);n<-length(p);if(!n) NA_real_ else min(n*p/seq_len(n))}
allr <- list()
for (CELL in CELLS) for (TAG in TAGS) for (ENV in ENVS) {
  D <- list(); LK <- list()
  for (i in 1:10) {
    f <- sprintf("%s/adapt_%s_chr%d_%s_env%d.rds", SIM, TAG, i, CELL, ENV)
    if (!file.exists(f)) next
    x <- readRDS(f); m <- flag_true_qtns(as.data.table(x$map))
    s1 <- as.data.table(x$complexity_reduction$stage1$map_snp)[, .(marker, CL1=CL_id)]
    pr <- ld_prune_and_eMLG(GTs=x$GTs, stage1=x$complexity_reduction$stage1, LD_decay=x$LD_decay,
      ld_w_col="ld_w_095", ld_w_threshold=0.025, score_threshold=0.80, min_r2_rho=0.5,
      distance_threshold=1e5, compute_unflagged_eMLG=TRUE, cores=1)
    s2 <- ld_group_map(pr, prefix=i)[, .(marker, CL2=group_id)]
    mm <- merge(merge(m, s1, by="marker", all.x=TRUE), s2, by="marker", all.x=TRUE)
    mm <- mm[!is.na(CL1)&!is.na(CL2)][, CL1:=paste0(i,"_",CL1)]
    th <- score_thresholds(as.data.table(x$LD_decay$decay_sum), rho_r2=0.75, rho_d=0.95, dmax_cap=1e5)
    drv <- mm[true_pos_QTN %in% TRUE]
    if (nrow(drv)) LK[[length(LK)+1]] <- rbindlist(lapply(seq_len(nrow(drv)), function(j) {
      near <- mm[as.character(Chr)==as.character(drv$Chr[j]) & abs(Pos-drv$Pos[j]) < th$dmax]
      if (!nrow(near)) return(NULL)
      r2 <- suppressWarnings(cor(x$GTs[,drv$marker[j]], x$GTs[,near$marker], use="pairwise.complete.obs")^2)
      d <- data.table(CL2=near$CL2, r2=as.numeric(r2))[is.finite(r2) & r2>=th$r2min]
      if (!nrow(d)) NULL else unique(d[, .(CL2)]) }))
    D[[length(D)+1]] <- mm[, .(marker, CL1, CL2, Pos, p=emx_p)]
  }
  DD <- rbindlist(D); tagged <- unique(rbindlist(LK)$CL2); m12 <- unique(DD[, .(CL1,CL2)])
  bh <- function(p){ok<-is.finite(p);q<-rep(NA_real_,length(p));q[ok]<-p.adjust(p[ok],"BH");which(!is.na(q)&q<0.05)}
  spans <- DD[, .(span=max(Pos)-min(Pos), n=.N), by=CL2]
  U2 <- DD[, .(p=simes(p)), by=CL2]; r2set <- U2$CL2[bh(U2$p)]
  U1 <- DD[, .(p=simes(p)), by=CL1]; s1s <- U1$CL1[bh(U1$p)]
  r1set <- unique(m12[CL1 %in% s1s]$CL2)
  for (nm in c("stage 2","stage 1")) {
    rs <- if (nm=="stage 2") r2set else r1set
    if (!length(rs)) next
    sp <- spans[CL2 %in% rs]
    allr[[length(allr)+1]] <- data.table(cell=CELL, tag=TAG, env=ENV, arm=nm, CL2=sp$CL2,
      span=sp$span, n=sp$n, tp=sp$CL2 %in% tagged)
  }
}
R <- rbindlist(allr)
if (nzchar(OUTF)) fwrite(R, OUTF)
cat("== region spans and precision by arm\n")
print(R[, .(regions=.N, median_span_kb=round(median(span)/1e3,1), median_markers=median(n),
            precision=round(mean(tp),3)), by=arm])
cat("\n== is TAGGING span-sensitive? precision by span decile, pooled\n")
R[, sb := cut(span, quantile(span, 0:5/5), include.lowest=TRUE, labels=c("Q1","Q2","Q3","Q4","Q5"))]
print(R[, .(regions=.N, median_span_kb=round(median(span)/1e3,1), precision=round(mean(tp),3)), by=sb][order(sb)])
cat("\n== SPAN-MATCHED: precision within each span quintile, by arm\n")
print(dcast(R[, .(prec=round(mean(tp),3), n=.N), by=.(sb,arm)], sb ~ arm, value.var=c("prec","n")))
