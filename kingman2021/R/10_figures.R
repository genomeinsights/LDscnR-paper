## Figures for the kingman2021 report.  Style follows the project default:
## no grey strip background, genomic axes in Mb.
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
P   <- path.expand("~/gitlab/LDscnR-paper/kingman2021")
EP  <- path.expand("~/gitlab/LD-scaling-genome-scans/empirical_data/kingman2021/ecopeaks")
TR  <- path.expand("~/gitlab/LD-scaling-genome-scans/empirical_data/kingman2021/tracks")
FIG <- file.path(P, "figures"); dir.create(FIG, showWarnings = FALSE)
ROMAN <- c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV",
           "XVI","XVII","XVIII","XIX","XX","XXI")
th <- theme_minimal(base_size = 10) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 7),
        panel.grid.minor = element_blank(), legend.position = "top",
        panel.spacing.x = unit(0.06, "lines"))

## ---------- Fig 1: Kingman c155 signal, published peaks, and the 3sp regions ----------
d <- fread(file.path(EP, "c155_global.snp_p.tsv.gz"), select = c("Chr","Pos","p"))
d[, chr_num := match(sub("^chr","",Chr), ROMAN)]
d <- d[!is.na(chr_num)]
set.seed(1)
plotd <- rbind(d[p < 0.01], d[p >= 0.01][sample(.N, min(.N, 300000))])
plotd[, Chr := factor(paste0("chr", ROMAN[chr_num]), levels = paste0("chr", ROMAN))]
pk <- fread(file.path(TR, "gasAcu1-4.c155.specific.50kb.final.peaks.bed"), header = FALSE,
            select = 1:3, col.names = c("Chr","start","end"))
pk[, Chr := factor(Chr, levels = levels(plotd$Chr))]; pk <- pk[!is.na(Chr)]
## Global-specific peaks have a median width of 21 kb, which is sub-pixel on a 30 Mb
## chromosome. Pad the drawn rectangle (not the data) to a minimum visible width.
PAD <- 2e5
pk[, `:=`(dstart = pmax(0, (start+end)/2 - PAD/2), dend = (start+end)/2 + PAD/2)]
pk[end-start > PAD, `:=`(dstart = start, dend = end)]
rg <- fread(file.path(P,"data","liftover","lfmm_g14.bed"), header = FALSE,
            col.names = c("Chr","start","end","name"))
rg[, Chr := factor(Chr, levels = levels(plotd$Chr))]; rg <- rg[!is.na(Chr)]
rg[, `:=`(dstart = pmax(0, (start+end)/2 - PAD/2), dend = (start+end)/2 + PAD/2)]
rg[end-start > PAD, `:=`(dstart = start, dend = end)]
g1 <- ggplot(plotd, aes(Pos/1e6, -log10(p))) +
  geom_rect(data = rg, inherit.aes = FALSE, aes(xmin = dstart/1e6, xmax = dend/1e6, ymin = -Inf, ymax = Inf),
            fill = "#3182bd", alpha = 0.30) +
  geom_rect(data = pk, inherit.aes = FALSE, aes(xmin = dstart/1e6, xmax = dend/1e6, ymin = -Inf, ymax = Inf),
            fill = "firebrick", alpha = 0.55) +
  geom_point(colour = "grey30", size = 0.25, alpha = 0.35) +
  facet_wrap(~ Chr, nrow = 1, scales = "free_x") + th +
  labs(x = "position (Mb)", y = expression(-log[10](italic(p))~"  Kingman c155 marine vs freshwater"),
       title = "Kingman global (c155) SNP-based ecotype signal",
       subtitle = "red = published Global-specific EcoPeaks; blue = 3sp LD-aware outlier regions. Both padded to a 200 kb minimum drawn width for visibility.")
ggsave(file.path(FIG,"fig1_c155_manhattan_with_regions.png"), g1, width = 17, height = 4.2, dpi = 200)

## ---------- Fig 2: validation of the re-implemented test ------------------------------
val <- rbindlist(list(
  data.table(cohort="c155 global", fdr=c(1e-2,1e-3,1e-4,1e-5), n_peaks=c(286,94,46,21),
             recovered=c(36,36,34,20), n_pub=39),
  data.table(cohort="c150 N.E. Pacific", fdr=1e-4, n_peaks=274, recovered=195, n_pub=209)))
val[, pct := 100*recovered/n_pub]
g2a <- ggplot(val, aes(factor(fdr), pct, fill = cohort)) +
  geom_col(position = position_dodge(preserve="single"), width=.65) +
  geom_hline(yintercept = 90, linetype = 2, linewidth = .4) +
  scale_fill_manual(values = c("c155 global"="#2c7fb8","c150 N.E. Pacific"="#d95f0e")) +
  th + labs(x = "BH FDR threshold", y = "% of published specific\nEcoPeaks recovered",
            title = "A  Re-implemented SNP test recovers the published peaks", fill = NULL)
g2b <- ggplot(val[cohort=="c155 global"], aes(n_peaks, pct)) +
  geom_line(colour="grey50") + geom_point(size=2.4, colour="#2c7fb8") +
  geom_vline(xintercept = 39, linetype = 2, colour = "firebrick") +
  annotate("text", x = 60, y = 55, label = "published:\n39 peaks", colour="firebrick", size=3, hjust=0) +
  th + labs(x = "number of peaks I call", y = "% published recovered",
            title = "B  Sensitivity / set-size trade-off (c155)")
ggsave(file.path(FIG,"fig2_validation.png"), g2a + g2b + plot_layout(widths=c(1.3,1)),
       width = 9.5, height = 3.6, dpi = 200)

## ---------- Fig 3: the geographic gradient --------------------------------------------
enr <- fread(file.path(P,"data","enrichment_summary.csv"))
enr[, lab := c(c151_nEur="c151 N. Europe\n(matches 3sp)", c155_global="c155 Global\n(partial match)",
               c150_pacNW="c150 N.E. Pacific\n(different basin)")[cohort]]
enr[, lab := factor(lab, levels = c("c151 N. Europe\n(matches 3sp)","c155 Global\n(partial match)",
                                    "c150 N.E. Pacific\n(different basin)"))]
g3 <- ggplot(enr[top==0.01], aes(lab, fold, fill = lab)) +
  geom_col(width = .6) + geom_hline(yintercept = 1, linetype = 2) +
  geom_text(aes(label = sprintf("%.2fx\np = %.4f", fold, pval)), vjust = -0.25, size = 3.1) +
  scale_fill_manual(values = c("#238b45","#2c7fb8","#d95f0e"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0,.22))) + th +
  labs(x = NULL, y = "fold enrichment of the top 1% of\nKingman p inside 3sp outlier regions",
       title = "Agreement tracks cohort geography",
       subtitle = "top 1% of Kingman p, inside 3sp outlier regions vs a within-chromosome shuffle")
ggsave(file.path(FIG,"fig3_geographic_gradient.png"), g3, width = 7.0, height = 4.2, dpi = 200)

## ---------- Fig 4: the shared chrI locus ----------------------------------------------
W <- list(chr = "chrI", lo = 25.6e6, hi = 27.1e6)
dz <- d[Chr == W$chr & Pos >= W$lo & Pos <= W$hi]
c151 <- fread(file.path(EP,"c151_nEur.snp_p.tsv.gz"), select = c("Chr","Pos","p"))[
  Chr == W$chr & Pos >= W$lo & Pos <= W$hi]
z <- rbind(dz[, .(Pos, p, panel = "Kingman c155 (global)")],
           c151[, .(Pos, p, panel = "Kingman c151 (N. Europe)")])
z[, panel := factor(panel, levels = c("Kingman c155 (global)","Kingman c151 (N. Europe)"))]
pk1 <- pk[Chr == W$chr & end > W$lo & start < W$hi]
rg1 <- rg[Chr == W$chr & end > W$lo & start < W$hi]
g4 <- ggplot(z, aes(Pos/1e6, -log10(p))) +
  geom_rect(data = rg1, inherit.aes = FALSE, aes(xmin=start/1e6, xmax=end/1e6, ymin=-Inf, ymax=Inf),
            fill = "#3182bd", alpha = .22) +
  geom_rect(data = pk1, inherit.aes = FALSE, aes(xmin=start/1e6, xmax=end/1e6, ymin=-Inf, ymax=Inf),
            fill = "firebrick", alpha = .30) +
  geom_point(size = .5, alpha = .5, colour = "grey20") +
  facet_wrap(~ panel, ncol = 1, scales = "free_y") + th +
  labs(x = "chrI position (Mb, gasAcu1-4)", y = expression(-log[10](italic(p))),
       title = "The one locus every analysis agrees on: chrI ~26.1-26.5 Mb",
       subtitle = "blue = 3sp LD-aware outlier region 3; red = published Global-specific EcoPeak\ngenes: tns1, igfbp5, igfbp2a, stk11ip, cx44.2, inha, spega, ctdsp1")
ggsave(file.path(FIG,"fig4_chrI_locus.png"), g4, width = 8, height = 5, dpi = 200)
cat("wrote 4 figures to", FIG, "\n")
