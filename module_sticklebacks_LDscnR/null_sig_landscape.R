## =====================================================================
## module_sticklebacks_LDscnR / null_sig_landscape.R
##
## Second-tier C-score for regions. As the marker C-score is the fraction of
## (rho, q*) cells in which a MARKER is a valid outlier candidate, the region score
## here is the fraction of (tau_C, l_min) cells in which a REGION is BH-significant
## (location-matched empirical p, region-level FDR q_R < 0.05) against the global
## permutation null. Integrating over the grid dissolves the "pick the best cell"
## multiple-testing trap: a region significant across many cells is credible; one
## significant in a single cell is not.
##
## Outputs:
##   * a per-cell n_sig heatmap (# regions with q_R<0.05 at each tau_C x l_min), with
##     the a-priori operating point (x) and the post-hoc max-significant cell (o);
##   * a ranked table of consensus loci by their second-tier stability score
##     (fraction of the 25 x 7 = 175 grid cells in which the locus is significant).
##
## Run from the LDscnR-paper root:
##   Rscript module_sticklebacks_LDscnR/null_sig_landscape.R
## Writes figures/null_sig_landscape_popperm.png, figures/region_stability2_manhattan.png,
##        results/null_sig_landscape.csv, results/region_stability2.csv
## =====================================================================
suppressMessages({ library(data.table); library(LDscnR); library(ggplot2) })

BND <- "module_sticklebacks_LDscnR/data/3sp_LDscnR_data.rds"
NULLF <- "module_sticklebacks_LDscnR/results/null_popperm_3sp.rds"
RES <- "module_sticklebacks_LDscnR/results"; OUTFIG <- "module_sticklebacks_LDscnR/figures"
TAUS <- seq(0.02, 0.5, by = 0.02); LMINS <- c(1,2,3,5,10,15,20); NCELL <- length(TAUS)*length(LMINS)
RHO_LD <- 0.60; DCAP <- 1e5; BCAP <- 100L; FDR <- 0.05; OP_TAU <- 0.05; OP_LMIN <- 3L; GAP <- 1e4
ZISSOU <- c("#3B9AB2","#78B7C5","#EBCC2A","#E1AF00","#F21A00")

d <- readRDS(BND); map <- as.data.table(d$map)
null <- readRDS(NULLF); B <- min(BCAP, length(null$C_surr)); surrs <- null$C_surr[seq_len(B)]; C_obs <- null$C_obs
edges <- ld_edges(null$universe, d$GTs, map[, .(marker, Chr, Pos)],
                  as.data.table(d$LD_decay$decay_sum), rho_ld = RHO_LD, dcap = DCAP)
mpos <- stats::setNames(map$Pos, map$marker); mchr <- stats::setNames(as.integer(gsub("Chr","",map$Chr)), map$marker)
cat(sprintf("[1] pop-perm B=%d ; edges over universe=%d ; %d grid cells\n", B, length(null$universe), NCELL)); flush.console()

cluster_at <- function(C, tau) { mk <- names(C)[C >= tau]
  if (!length(mk)) return(data.table(size=integer(), chr=integer(), lo=numeric(), hi=numeric(), score=numeric()))
  r <- ld_regions(mk, edges)
  rbindlist(lapply(r, function(x) data.table(size=length(x), chr=unname(mchr[x[1]]), lo=min(mpos[x]), hi=max(mpos[x]), score=sum(C[x])))) }
emp_q <- function(O, Sl) { ep <- vapply(seq_len(nrow(O)), function(i) {
    best <- vapply(Sl, function(S) { if (!nrow(S)) return(0); h <- S[chr==O$chr[i] & lo<=O$hi[i] & hi>=O$lo[i]]; if (!nrow(h)) 0 else max(h$score) }, numeric(1))
    (1 + sum(best >= O$score[i])) / (1 + B) }, numeric(1)); stats::p.adjust(ep, "BH") }

sig_acc <- list(); counts <- list()
for (tau in TAUS) {
  O_all <- cluster_at(C_obs, tau); S_all <- lapply(surrs, function(C) cluster_at(C, tau))
  for (lm in LMINS) {
    O <- O_all[size >= lm]
    if (!nrow(O)) { counts[[length(counts)+1L]] <- data.table(tau=tau, lmin=lm, n_obs=0L, n_sig=0L); next }
    O[, q_R := emp_q(O, lapply(S_all, function(S) S[size >= lm]))]
    sig <- O[q_R < FDR]
    counts[[length(counts)+1L]] <- data.table(tau=tau, lmin=lm, n_obs=nrow(O), n_sig=nrow(sig))
    if (nrow(sig)) sig_acc[[length(sig_acc)+1L]] <- sig[, .(tau, lmin=lm, chr, lo, hi, score, q_R)]
  }
  cat(sprintf("   tau=%.2f done\n", tau)); flush.console()
}
res <- rbindlist(counts); fwrite(res, file.path(RES, "null_sig_landscape.csv"))
best <- res[which.max(n_sig)]; op <- res[tau==OP_TAU & lmin==OP_LMIN]
NUSE <- res[n_sig > 0, .N]                 # usable parameter space: cells yielding >=1 significant region
cat(sprintf("[1b] usable cells (>=1 significant region) = %d of %d\n", NUSE, NCELL)); flush.console()

## ---- second-tier score: consensus loci across significant cells --------------------
sig_all <- rbindlist(sig_acc)
setorder(sig_all, chr, lo)
sig_all[, cl := { o<-order(lo); s<-lo[o]; e<-hi[o]; g<-c(TRUE, s[-1] > cummax(e)[-length(e)] + GAP); cumsum(g)[order(o)] }, by = chr]
loci <- sig_all[, .(lo=min(lo), hi=max(hi), n_cells_sig=uniqueN(paste(tau,lmin)),
                    best_q=min(q_R), max_score=round(max(score),2),
                    tau_range=sprintf("%.2f-%.2f", min(tau), max(tau)),
                    lmin_max=max(lmin)), by=chr]
loci[, stability := round(n_cells_sig / NUSE, 3)][, locus := sprintf("Chr%d:%.2f-%.2f", chr, lo/1e6, hi/1e6)]
loci[chr==1 & lo<21.93e6 & hi>21.40e6, locus := paste0(locus, " (inversion)")]
loci[chr==4 & lo<12.83e6 & hi>12.79e6, locus := paste0(locus, " (Eda)")]
setorder(loci, -stability, best_q)
fwrite(loci, file.path(RES, "region_stability2.csv")); fwrite(sig_all, file.path(RES, "sig_all_cells.csv"))

## ---- per-region C2: anchor on the structure-null operating point (tau=0.05, l_min=3) --
## The anchor supplies clean, LD-supported region boundaries (out of the small-cluster
## corner); C2 then annotates each with its cross-grid robustness. C2 is insensitive to
## the exact anchor cell -- a real locus scores high whether drawn here or a neighbour.
anchor <- cluster_at(C_obs, OP_TAU)[size >= OP_LMIN]
anchor[, c2 := 0]
for (i in seq_len(nrow(anchor))) { ov <- sig_all[chr==anchor$chr[i] & lo<=anchor$hi[i] & hi>=anchor$lo[i]]
  set(anchor, i, "c2", uniqueN(paste(ov$tau, ov$lmin)) / NUSE) }
anchor[, locus := sprintf("Chr%d:%.2f-%.2f", chr, lo/1e6, hi/1e6)]
anchor[chr==1 & lo<21.93e6 & hi>21.40e6, locus := paste0(locus, " (inversion)")]
anchor[chr==4 & lo<12.83e6 & hi>12.79e6, locus := paste0(locus, " (Eda)")]
setorder(anchor, -c2, -score)
fwrite(anchor, file.path(RES, "region_c2_anchored.csv"))

cat(sprintf("\n=== per-cell: max cell tau=%.2f l_min=%d -> %d/%d significant (grid steps 0.02; op tau=0.05 clustered off-grid for the anchor) ===\n",
            best$tau, best$lmin, best$n_sig, best$n_obs))
cat(sprintf("\n=== PER-REGION C2: %d operating-point regions (tau=%.2f, l_min=%d), ranked by grid-stability ===\n",
            nrow(anchor), OP_TAU, OP_LMIN))
print(anchor[, .(locus, size, score=round(score,2), c2=round(c2,3))], nrow=100)

## ---- figures ----------------------------------------------------------------------
g1 <- ggplot(res, aes(factor(tau), factor(lmin, levels=LMINS), fill=n_sig)) +
  geom_tile(color="white", linewidth=0.2) +
  geom_point(data=data.table(tau=OP_TAU,lmin=OP_LMIN), inherit.aes=FALSE, aes(factor(tau),factor(lmin,levels=LMINS)), shape=4, size=3, stroke=1.2) +
  geom_point(data=best, inherit.aes=FALSE, aes(factor(tau),factor(lmin,levels=LMINS)), shape=1, size=4, stroke=1.2) +
  geom_text(aes(label=ifelse(n_sig>0,n_sig,"")), size=2.3, color="grey15") +
  scale_fill_gradientn(colours=ZISSOU, name="regions\nq_R<0.05") +
  labs(x=expression(tau[C]), y=expression(l[min]),
       title="# regions surviving BH FDR (q_R<0.05) vs global permutation null, per grid point",
       subtitle=sprintf("x = a-priori operating point (%.2f, %d); o = post-hoc max cell (%.2f, %d) -- post-hoc max is grid-selection biased", OP_TAU, OP_LMIN, best$tau, best$lmin)) +
  theme_minimal(base_size=10) + theme(panel.grid=element_blank(), axis.text.x=element_text(angle=90,vjust=0.5,size=6), plot.subtitle=element_text(size=8))
ggsave(file.path(OUTFIG,"null_sig_landscape_popperm.png"), g1, width=13, height=5, dpi=170)

## Manhattan coloured by the second-tier stability score
mh <- copy(map[, .(marker, Chr, Pos, C=C_obs[map$marker])]); mh[is.na(C), C:=0]; mh[, chr:=as.integer(gsub("Chr","",Chr))]
mh[, c2 := NA_real_]
for (i in seq_len(nrow(anchor))) mh[chr==anchor$chr[i] & Pos>=anchor$lo[i] & Pos<=anchor$hi[i], c2 := anchor$c2[i]]
chr_lev <- paste0("Chr", sort(unique(mh$chr))); mh[, Chr:=factor(Chr, levels=chr_lev)]
g2 <- ggplot(mh, aes(Pos, C)) +
  geom_point(data=mh[is.na(c2)], color="grey80", size=0.4, alpha=0.5) +
  geom_point(data=mh[!is.na(c2)], aes(color=c2), size=1.8) +
  facet_wrap(~Chr, nrow=1, scales="free_x") +
  scale_color_gradientn(colours=ZISSOU, name=expression(C^(2)), limits=c(0,1)) +
  labs(x="Genomic position", y="C-score",
       title=expression("3sp operating-point regions coloured by the second-tier C-score "*C^(2)*" (fraction of usable "*tau[C]*" x "*l[min]*" cells BH-significant vs the permutation null)")) +
  theme_minimal(base_size=10) + theme(panel.grid=element_blank(), strip.text=element_text(face="bold",size=6), axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.spacing.x=unit(0.05,"lines"))
ggsave(file.path(OUTFIG,"region_c2_manhattan.png"), g2, width=18, height=4.5, dpi=170)
cat("\n[2] wrote null_sig_landscape_popperm.png + region_c2_manhattan.png\n")
