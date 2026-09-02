## =============================================================================
## module_sim_LDscnR / ldw_dcap_coupling.R
##
## ld_w_threshold and distance_threshold are not independent axes. dcap splits a
## run at a gap between consecutive FLAGGED stage-1 clusters, and ld_w_threshold
## is what decides which clusters are flagged -- so moving ld_w moves the dcap
## operating point without dcap changing at all.
##
##   ld_w    flagged clusters   gaps > 1e5   OPERATING POINT   degenerate chr
##   0.0125       294,267          2,133         0.73%              5
##   0.0250       103,825          1,598         1.54%              0
##   0.0500        18,401          2,668        14.62%              6
##   0.1000         2,310            593        27.03%             21
##
## 37-fold movement from ld_w alone, monotone. The panel shows 2,000-fold in the
## same direction, so this is a package property rather than a bgs5 one.
##
## DEGENERACY AT BOTH ENDS. A chromosome with no gap above dcap becomes a SINGLE
## RUN and stage 2 cannot preserve any structure on it. That happens at the low
## end because flagged clusters are packed too densely to leave a large gap, and
## at the high end because too few survive -- at ld_w = 0.10, 44 of 160
## chromosomes have fewer than two flagged clusters and drop out entirely. The
## canonical 0.025 is the only level tested with no degenerate chromosome, which
## is luck rather than design.
##
## So the two parameters have a JOINT ADMISSIBLE REGION, dataset-specific, and a
## factorial grid over them has cells that should not be run rather than cells
## that merely correlate.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/ldw_dcap_coupling.R
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(parallel)})
SIM <- "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5"
THR <- c(0.0125, 0.025, 0.05, 0.10)
jobs <- CJ(CELL=c("V0.5_c1","V0.5_c2","V1_c1.5","V2_c1"), TAG=c("nobgs","bgs"), i=1:10, sorted=FALSE)
one <- function(k) {
  r <- jobs[k]
  f <- sprintf("%s/adapt_%s_chr%d_%s_env1.rds", SIM, r$TAG, r$i, r$CELL)
  if (!file.exists(f)) return(NULL)
  x <- readRDS(f); ms <- as.data.table(x$complexity_reduction$stage1$map_snp)
  n_all <- uniqueN(ms$CL_id)
  rbindlist(lapply(THR, function(t) {
    ids <- ms[ld_w_095 > t, unique(CL_id)]
    ext <- ms[CL_id %in% ids, .(Chr=Chr[1], pmin=min(Pos), pmax=max(Pos)), by=CL_id]
    rbindlist(lapply(unique(ext$Chr), function(ch) {
      e <- ext[Chr==ch][order(pmin)]; if (nrow(e) < 2) return(NULL)
      g <- e$pmin[-1] - e$pmax[-nrow(e)]; g <- g[is.finite(g)]
      data.table(CELL=r$CELL, TAG=r$TAG, i=r$i, Chr=ch, ldw=t, n_all=n_all,
                 n_flag=nrow(e), n_gap=length(g), over=sum(g>1e5), mx=max(g)) })) })) }
G <- rbindlist(Filter(Negate(is.null), mclapply(seq_len(nrow(jobs)), one, mc.cores=8)))
fwrite(G, "module_sim_LDscnR/results/operating_points/ldw_dcap_coupling.csv")
S <- G[, .(flagged=sum(n_flag),
           splits=sum(over), op_point=round(100*sum(over)/sum(n_gap),4),
           max_gap_kb=round(max(mx)/1e3), degenerate_chr=sum(over==0), n_chr=.N), by=ldw][order(ldw)]
print(S)
cat(sprintf("\noperating point moves %.0f-fold across ld_w levels (panel: 2000-fold)\n",
            max(S$op_point)/min(S$op_point)))
cat("\ndegenerate chromosomes (zero gaps > 1e5, so the chromosome is ONE run) by cell:\n")
print(dcast(G[, .(deg=sum(over==0), n=.N), by=.(CELL,ldw)], CELL ~ ldw, value.var="deg"))
