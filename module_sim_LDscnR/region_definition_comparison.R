## =============================================================================
## module_sim_LDscnR / region_definition_comparison.R
##
## The panel session asked how a REGION is defined here, because region counts
## appear in both halves of the manuscript and the definitions had never been
## reconciled -- a crossing quantity we both missed.
##
## THIS SIDE, and it differs by arm:
##   stage-1 and single-SNP arms  a region is a distinct STAGE-2 GROUP containing
##                                at least one significant tested unit, assembled
##                                STRICTLY AFTER BH.
##   stage-2 arm                  regions ARE the tested units, so assembly is
##                                necessarily BEFORE testing. That arm is not
##                                test-then-merge at all.
##
## Stage-2 grouping is a HYBRID, not purely LD: flag stage-1 clusters with any
## member ld_w > threshold; within the flagged set split into runs where the gap
## between consecutive stage-1 clusters exceeds 100 kb; merge by LD within a run
## under a dynamic cut. The panel's rule is purely physical, merging significant
## units on a chromosome when the gap is <= 300 kb.
##
## MEASURED ON THE SAME SIGNIFICANT UNITS, the two disagree:
##
##   cell      env   sig units   LD-aware   distance 300 kb   ratio
##   V0.5_c1    1       27          21            13          1.62
##   V0.5_c1    2       22          20            14          1.43
##   V0.5_c1    3       49          32            14          2.29
##   V0.5_c2    1       22          17            19          0.89
##   V0.5_c2    2      249         162           139          1.17
##   V0.5_c2    3      111          83            79          1.05
##
## LD-AWARE GIVES 1.41x MORE REGIONS ON AVERAGE, up to 2.29x, and occasionally
## fewer (0.89) because distance merging can also split where an LD group spans a
## gap over 300 kb. So region counts are NOT comparable between the two halves:
## LD splits what distance fuses, which is the behaviour PK argued for, since two
## causal variants can be adjacent and separable by LD while distance cannot see
## the difference.
##
## IT COMPOUNDS WITH THE SPAN PROBLEM. More and narrower regions lower a
## per-region peak-overlap rate for the SAME underlying discoveries, so
## region-level percentages are not comparable across definitions either. Every
## region-level number needs its definition stated beside it, and the
## significant-UNIT count reported alongside, since that one is definition-free.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/region_definition_comparison.R
## =============================================================================
suppressMessages({library(data.table); library(LDscnR)})
SIM <- "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5"; GAP <- 3e5
simes <- function(p){p<-sort(p[is.finite(p)]);n<-length(p);if(!n) NA_real_ else min(n*p/seq_len(n))}
out <- list()
for (CELL in c("V0.5_c1","V0.5_c2")) for (ENV in 1:3) {
  D <- list()
  for (i in 1:10) {
    f <- sprintf("%s/adapt_nobgs_chr%d_%s_env%d.rds", SIM, i, CELL, ENV)
    if (!file.exists(f)) next
    x <- readRDS(f); m <- flag_true_qtns(as.data.table(x$map))
    s1 <- as.data.table(x$complexity_reduction$stage1$map_snp)[, .(marker, CL1=CL_id)]
    pr <- ld_prune_and_eMLG(GTs=x$GTs, stage1=x$complexity_reduction$stage1, LD_decay=x$LD_decay,
      ld_w_col="ld_w_095", ld_w_threshold=0.025, score_threshold=0.80, min_r2_rho=0.5,
      distance_threshold=1e5, compute_unflagged_eMLG=TRUE, cores=1)
    s2 <- ld_group_map(pr, prefix=i)[, .(marker, CL2=group_id)]
    mm <- merge(merge(m, s1, by="marker", all.x=TRUE), s2, by="marker", all.x=TRUE)
    D[[length(D)+1]] <- mm[!is.na(CL1)&!is.na(CL2)][, `:=`(CL1=paste0(i,"_",CL1),
      chrkey=paste0(i,"_",Chr))][, .(marker, CL1, CL2, chrkey, Pos, p=emx_p)]
  }
  DD <- rbindlist(D)
  U1 <- DD[, .(p=simes(p)), by=CL1]
  ok <- is.finite(U1$p); q <- rep(NA_real_,nrow(U1)); q[ok] <- p.adjust(U1$p[ok],"BH")
  sig <- U1$CL1[which(!is.na(q) & q<0.05)]
  if (!length(sig)) next
  S <- DD[CL1 %in% sig, .(chrkey=chrkey[1], lo=min(Pos), hi=max(Pos), CL2=CL2[1]), by=CL1]
  ## (a) LD-aware: distinct stage-2 groups
  n_ld <- uniqueN(S$CL2)
  ## (b) distance-only: merge units on a chromosome when the gap <= GAP
  setorder(S, chrkey, lo)
  S[, gap := lo - shift(hi), by=chrkey]
  S[, newr := is.na(gap) | gap > GAP]
  n_dist <- S[, sum(newr)]
  out[[length(out)+1]] <- data.table(cell=CELL, env=ENV, sig_units=length(sig),
    regions_LD=n_ld, regions_dist300kb=n_dist)
}
R <- rbindlist(out)
R[, ratio := regions_LD/pmax(regions_dist300kb,1)]
cat("== assembling the SAME significant stage-1 clusters two ways\n\n")
print(R)
cat(sprintf("\n  mean LD-aware %.1f vs distance-300kb %.1f regions  => LD gives %.2fx more\n",
    mean(R$regions_LD), mean(R$regions_dist300kb), mean(R$ratio)))
