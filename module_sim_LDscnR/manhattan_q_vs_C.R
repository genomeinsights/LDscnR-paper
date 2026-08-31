## =====================================================================
## module_sim_LDscnR / manhattan_q_vs_C.R
##
## What each rule actually selects, on one genome: -log10(q) with the BH 0.05
## line beside the C-score with the tau_C = 0.05 line, same markers, same axis.
## Plus the called regions coloured by how many SNPs they contain.
##
## Continuous cumulative genome axis with alternating chromosome bands, not 20
## facets -- at 20 facets the points are too small to read.
##
## ld_w colour is CAPPED at 0.5 (values above show at the top colour): the
## distribution is long-tailed and an uncapped scale puts almost every marker
## in the bottom colour bin.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/manhattan_q_vs_C.R
## Env: SIM_INPUTS, NULL_DIR, OUT, ENV (default 1)
## =====================================================================
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(scales); library(LDscnR)})
A <- "/Volumes/Nemo/Nemo_sim/analysis_inputs"
D <- "module_sim_LDscnR/results/nulls_V2_c1"; E <- 1
panel <- readRDS(file.path(A, sprintf("panel_V2_c1_env%d.rds", E)))
map <- flag_true_qtns(as.data.table(panel$map))
lv <- unique(map$Chr)
lv <- lv[order(as.integer(sub("^R([0-9]+)_.*","\\1",lv)), sub("^R[0-9]+_","",lv))]
map[, Chr := factor(Chr, levels = lv)]; setorder(map, Chr, Pos)
off <- map[, .(len = max(Pos)), by = Chr][, start := cumsum(shift(len, fill = 0))]
off[, `:=`(mid = start + len/2, band = seq_len(.N) %% 2 == 1)]
map <- merge(map, off[, .(Chr, start)], by = "Chr", sort = FALSE)
map[, gx := (start + Pos)/1e6]
bands <- off[band == TRUE][, `:=`(xmin = start/1e6, xmax = (start+len)/1e6)]

d <- rbindlist(lapply(c("emmax","lfmm"), function(eng) {
  f <- file.path(A, sprintf("pvals_V2_c1_env%d_%s_%s_B100.rds", E, eng,
                            if (eng=="emmax") "genetic" else "env_orth"))
  p <- readRDS(f)$p_obs
  C <- ld_cscore(p, panel$ld_ws, alpha=0.05, rho=colnames(panel$ld_ws), qstar=seq(0,.95,by=.05))
  Cv <- rep(0, nrow(map)); names(Cv) <- map$marker; Cv[names(C)] <- C
  data.table(engine=eng, Chr=map$Chr, Pos=map$Pos, gx=map$gx,
             true_QTN=map$true_pos_QTN %in% TRUE,
             mlq=-log10(pmax(p.adjust(p,"BH"),1e-12)), C=Cv,
             ld_w=map$ld_w_095, qtnLD=map$max_LD_with_QTN) }))

base <- function(dd, yv, thr, ylab, ttl) list(
  geom_rect(data=bands, inherit.aes=FALSE, aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf),
            fill="grey93", alpha=.5),
  geom_hline(yintercept=thr, colour="#9E4630", linewidth=.55, linetype="22"),
  facet_wrap(~engine, ncol=1, scales="free_y"),
  scale_x_continuous(breaks=off$mid/1e6, labels=sub("^R","",as.character(off$Chr)), expand=c(.005,0)),
  labs(y=ylab, x=NULL, title=ttl),
  theme_bw(base_size=12),
  theme(strip.background=element_blank(), panel.grid=element_blank(),
        axis.text.x=element_text(size=7, angle=90, vjust=.5), plot.title=element_text(size=13)))

## ---- (1) ld_w capped at 0.5 -------------------------------------------------
mk1 <- function(yv, thr, ylab, ttl) { dd <- copy(d); setorder(dd, ld_w)
  ggplot(dd, aes(gx, .data[[yv]], colour=ld_w)) + base(dd,yv,thr,ylab,ttl) +
    geom_point(size=1.05, alpha=.8) +
    geom_point(data=dd[true_QTN==TRUE], shape=3, colour="black", size=3.4, stroke=1) +
    scale_colour_viridis_c(name=expression(ld[w]), limits=c(0,0.5), oob=squish) }
ggsave("module_sim_LDscnR/figures/manhattan_q_vs_C_env1_by_ldw.png",
  (mk1("mlq", -log10(0.05), expression(-log[10]~q~"(BH)"), "BH-adjusted p-values -- dashed line q = 0.05") /
   mk1("C", 0.05, "C-score", "C-score -- dashed line tau_C = 0.05")) +
  plot_annotation(title=sprintf("V2_c1 env%d -- coloured by local LD", E),
    subtitle="ld_w colour CAPPED at 0.5; values above are shown at the top colour. Crosses = detectable QTN.",
    theme=theme(plot.subtitle=element_text(size=10, colour="grey35"))),
  width=20, height=12, dpi=200)
cat("  wrote ..._by_ldw.png (capped at 0.5)\n")

## ---- (2) called regions, coloured by how many SNPs they contain -------------
reg <- rbindlist(lapply(c("genetic"), function(b) {
  r <- fread(file.path(D, sprintf("regions_V2_c1_env%d_emmax_%s_B100.csv", E, b)))
  r[, basis := b] }))
reg <- merge(reg, off[, .(Chr, start)], by.x="chr", by.y="Chr", sort=FALSE)
reg[, `:=`(x1=(start+lo)/1e6, x2=(start+hi)/1e6, xm=(start+(lo+hi)/2)/1e6)]
## markers belonging to a called region, tagged with that region's SNP count
mem <- rbindlist(lapply(seq_len(nrow(reg)), function(i)
  d[engine=="emmax" & Chr==reg$chr[i] & Pos>=reg$lo[i] & Pos<=reg$hi[i]][, size := reg$size[i]]))
pA <- ggplot(d[engine=="emmax"], aes(gx, C)) +
  geom_rect(data=bands, inherit.aes=FALSE, aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf),
            fill="grey93", alpha=.5) +
  geom_point(size=.8, colour="grey78") +
  geom_point(data=mem, aes(colour=size), size=1.5) +
  geom_hline(yintercept=0.05, colour="#9E4630", linewidth=.55, linetype="22") +
  geom_point(data=d[engine=="emmax" & true_QTN==TRUE], shape=3, colour="black", size=3.4, stroke=1) +
  scale_colour_viridis_c(name="SNPs in\nregion", trans="log10", option="C") +
  scale_x_continuous(breaks=off$mid/1e6, labels=sub("^R","",as.character(off$Chr)), expand=c(.005,0)) +
  labs(y="C-score", x=NULL, title="Markers belonging to a called outlier region, coloured by region size") +
  theme_bw(base_size=12) + theme(panel.grid=element_blank(),
    axis.text.x=element_text(size=7, angle=90, vjust=.5), plot.title=element_text(size=13))
pB <- ggplot(reg, aes(x=xm, y=s_R, colour=size)) +
  geom_rect(data=bands, inherit.aes=FALSE, aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf),
            fill="grey93", alpha=.5) +
  geom_segment(aes(x=x1, xend=x2, y=s_R, yend=s_R), linewidth=2.4) +
  geom_point(data=reg[has_qtn==TRUE], shape=3, colour="black", size=3.4, stroke=1) +
  scale_colour_viridis_c(name="SNPs in\nregion", trans="log10", option="C") +
  scale_y_log10() +
  scale_x_continuous(breaks=off$mid/1e6, labels=sub("^R","",as.character(off$Chr)), expand=c(.005,0)) +
  labs(y=expression(s[R]~"(summed C-mass, log)"), x=NULL,
       title="The called regions themselves -- span on x, s_R on y, cross = contains a QTN") +
  theme_bw(base_size=12) + theme(panel.grid=element_blank(),
    axis.text.x=element_text(size=7, angle=90, vjust=.5), plot.title=element_text(size=13))
ggsave("module_sim_LDscnR/figures/regions_by_size_env1.png",
  (pA / pB) + plot_annotation(
    title=sprintf("V2_c1 env%d, EMMAX -- outlier regions by SNP count (tau_C 0.05, l_min 3)", E),
    subtitle=sprintf("%d regions, %d contain a detectable QTN. Grey = markers in no region.",
                     nrow(reg), sum(reg$has_qtn %in% TRUE)),
    theme=theme(plot.subtitle=element_text(size=10, colour="grey35"))),
  width=20, height=11, dpi=200)
cat(sprintf("  wrote regions_by_size_env1.png | %d regions, sizes %d-%d\n",
            nrow(reg), min(reg$size), max(reg$size)))
