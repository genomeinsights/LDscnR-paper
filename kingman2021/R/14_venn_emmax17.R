## Venn of a 3sp region set's overlap with the Kingman *specific* EcoPeak sets
## (Global-specific vs Pacific-specific), from R/12 (EMMAX) or R/15 (LFMM).
## Usage: Rscript R/14_venn_emmax17.R [tag]      tag = emmax17 (default) | lfmm
## Reads data/overlap_detail_<tag>.csv + data/overlap_summary_<tag>.csv;
## writes figures/fig5_venn_<tag>.png. Small zones are labelled by locus; large
## zones show counts (with the Chr1 inversion / Eda flagged in the intersection).
suppressMessages(library(data.table))
P   <- path.expand("~/gitlab/LDscnR-paper/kingman2021")
a   <- commandArgs(trailingOnly = TRUE); TAG <- if (length(a) >= 1) a[1] else "emmax17"
MLAB <- if (TAG == "lfmm") "LFMM" else "EMMAX"
det <- fread(file.path(P, "data", sprintf("overlap_detail_%s.csv", TAG)))
sm  <- fread(file.path(P, "data", sprintf("overlap_summary_%s.csv", TAG)))
n_total <- sm$n_regions[1]
gp <- function(s) round(sm[set == s, fold_region], 1); pp <- function(s) signif(sm[set == s, p_region], 2)

det[, lab := sprintf("%s:%.2f", Chr, reg_start/1e6)]
det[, notable := (Chr == "Chr1" & reg_start > 21.3e6 & reg_start < 22.0e6) |
                 (Chr == "Chr4" & reg_start > 12.7e6 & reg_start < 12.9e6)]
glob <- unique(det[set == "GlobalSpec", region]); paci <- unique(det[set == "PacSpec", region])
both <- intersect(glob, paci); gonly <- setdiff(glob, paci); ponly <- setdiff(paci, glob)
n_neither <- n_total - length(union(glob, paci))
labof <- function(regs) unique(det[region %in% regs, lab])
nice <- function(v) { v <- sub("Chr1:21.50", "Chr1 inversion", v, fixed = TRUE)
  sub("Chr4:12.81", "Chr4/Eda", v, fixed = TRUE) }
zonetext <- function(regs, cap = 5L) { L <- nice(labof(regs))       # big count number shown separately
  if (length(L) <= cap) paste(L, collapse = "\n")
  else { hot <- nice(labof(intersect(regs, det[notable == TRUE, region])))
    if (length(hot)) paste0("(incl. ", paste(hot, collapse = ", "), ")") else "" } }

circ <- function(x, y, r, col, bord) { a <- seq(0, 2*pi, length.out = 300)
  polygon(x + r*cos(a), y + r*sin(a), col = col, border = bord, lwd = 2.5) }

png(file.path(P, "figures", sprintf("fig5_venn_%s.png", TAG)), width = 2000, height = 1500, res = 220)
par(mar = c(3.6, 1, 3.2, 1))
plot.new(); plot.window(xlim = c(-2.5, 2.5), ylim = c(-1.9, 1.7), asp = 1)
BLU <- "#1F78B4"; ORA <- "#E08214"
circ(-0.62, 0, 1.2, adjustcolor(BLU, 0.16), BLU); circ(0.62, 0, 1.2, adjustcolor(ORA, 0.18), ORA)
text(-1.55, 1.42, "Kingman\nGlobal-specific", col = BLU, font = 2, cex = 1.02)
text( 1.55, 1.42, "Kingman\nPacific-specific", col = ORA, font = 2, cex = 1.02)
text(-1.16, 0.34, zonetext(gonly), cex = 0.82, col = "grey15")
text( 1.16, 0.34, zonetext(ponly), cex = 0.82, col = "grey15")
text( 0.00, 0.40, zonetext(both),  cex = 0.84, font = 2, col = "grey5")
text(-1.16, -0.72, length(gonly), cex = 1.5, col = BLU, font = 2)
text( 1.16, -0.72, length(ponly), cex = 1.5, col = ORA, font = 2)
text( 0.00, -0.34, length(both),  cex = 1.5, col = "grey20", font = 2)
text(0, -1.66, sprintf("%d of %d regions overlap no specific EcoPeak", n_neither, n_total),
     cex = 0.92, col = "grey35")
title(main = sprintf("3sp %s outlier regions (%d, l_min=3)  vs  Kingman EcoPeaks", MLAB, n_total),
      cex.main = 1.05)
mtext(sprintf("Region-level overlap, gasAcu1, rotation null B=2000:  Global-specific %.1fx, p=%s   ·   Pacific-specific %.1fx, p=%s",
              gp("c155.specific"), format(pp("c155.specific")), gp("c150.specific"), format(pp("c150.specific"))),
      side = 1, line = 2.0, cex = 0.80, col = "grey25")
dev.off()
cat(sprintf("[%s] both=%d Global-only=%d Pacific-only=%d neither=%d -> figures/fig5_venn_%s.png\n",
            TAG, length(both), length(gonly), length(ponly), n_neither, TAG))
