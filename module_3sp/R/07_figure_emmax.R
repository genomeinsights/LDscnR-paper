## =============================================================================
## module_3sp/R/07_figure_emmax.R
##
## THE MAIN FIGURE, faithfully ported from LD-pruning-paper's
## R_3sp_blocks/07_figures.R (the design PK chose as the manuscript template),
## sourced from module_3sp's real, current bundle rather than the retired one.
##
## THE Y AXIS IS -log10(q) FOR THE SINGLE MARKER, with the FDR line drawn. The
## claim is not "clustering finds more"; it is "clustering finds regions
## single-marker testing cannot see AT THE SAME FDR", and only a q axis shows
## that -- a p axis would let a reader assume the difference is a threshold
## choice.
##
## EcoPeak status is drawn PER REGION below the axis (triangles), not as
## background shading: shading says peaks exist somewhere, the markers say
## WHICH of our regions map and which do not.
##
## This REPLACES the two-panel EMMAX/LFMM figure at the manuscript's figure4
## slot (PK). That figure remains available as a supplementary/engine-
## comparison figure from 05_manhattan.R; it is not deleted, just not the one
## delivered to figure4 by 06_manuscript.R going forward.
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(ggplot2); library(ggrastr)})
devtools::load_all("~/gitlab/LDscnR")
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_3sp"), "R", "00_config.R"))
STAGE <- "07_figure_emmax"
say("=== %s ===\n\n", STAGE)

b  <- readRDS(file.path(PATHS$out, "02_bundle", "bundle.rds"))
sc <- readRDS(file.path(PATHS$out, "03_scan", "scan.rds"))
map <- b$map
EDA <- list(Chr = "Chr4", from = 12.70e6, to = 12.90e6)

## single-marker scan for the y-axis point cloud (q, not p -- the whole point of the figure)
pm <- emmax_fast(emmax_setup(b$GTs, b$GRM), b$eco)
qm <- p.adjust(pm, "BH")
say("[1] single-marker BH: %d of %s markers reach q <= %.2f (max -log10 q = %.2f)\n",
    sum(qm <= ALPHA), format(length(qm), big.mark=","), ALPHA, max(-log10(qm)))

## regions and EcoPeak status: REUSED from 03_scan.R, not recomputed
R <- copy(sc$consensus$test$regions); R[, rid := .I]
onp <- sc$consensus$rotation$per_region
setkey(onp, chr, from, to); setkey(R, Chr, from, to)
R[, on_peak := onp[.(Chr, from, to)]$on_peak]
eda_hit <- R[Chr == EDA$Chr & to >= EDA$from & from <= EDA$to]
say("[2] %d significant clusters -> %d regions, %d on an EcoPeak\n",
    sum(sc$consensus$test$units$significant), nrow(R), sum(R$on_peak))
say("[3] Eda (%s %.2f-%.2f Mb): %s\n", EDA$Chr, EDA$from/1e6, EDA$to/1e6,
    if (nrow(eda_hit)) sprintf("recovered, region %.2f-%.2f Mb", eda_hit$from[1]/1e6, eda_hit$to[1]/1e6)
    else "not recovered")

PAL <- LDscnR:::default_cluster_colours(); R[, col := NA_character_]; cur <- 0L
for (ch in unique(R$Chr)) { used <- character(0)
  for (i in R[Chr == ch, which = TRUE]) {
    repeat { cur <- cur %% length(PAL) + 1L; if (!(PAL[cur] %in% used)) break }
    R$col[i] <- PAL[cur]; used <- c(used, PAL[cur]) } }
stopifnot(R[, .(ok = uniqueN(col) == .N), by = Chr][, all(ok)])

M <- data.table(Chr = map$Chr, pos = map$Pos, y = -log10(pmax(qm, .Machine$double.xmin)))
setkey(R, Chr, from, to)
M[, rid := foverlaps(M[, .(Chr, from = pos, to = pos)], R,
      by.x = c("Chr","from","to"), type = "within", mult = "first", nomatch = NA)$rid]
OFF <- M[, .(mx = max(pos)), by = Chr][order(Chr)][, off := cumsum(c(0, head(mx,-1)) + 2e6)][]
M   <- merge(M, OFF[, .(Chr, off)], by = "Chr")[, gx := pos + off]
CV  <- setNames(R$col, as.character(R$rid))
chr_lab <- as.integer(gsub("Chr", "", OFF$Chr))
edax <- mean(c(EDA$from, EDA$to)) + OFF[Chr == EDA$Chr, off]
RB <- merge(R, OFF[, .(Chr, off)], by = "Chr")
RB[, gx := (from + to)/2 + off][, status := fifelse(on_peak, "maps to an EcoPeak", "does not")]
chr_n <- as.integer(gsub("Chr", "", M$Chr))
YR <- -0.30

p <- ggplot() +
  rasterise(geom_point(data = M[is.na(rid)], aes(gx, y, colour = factor(chr_n[is.na(M$rid)] %% 2)),
            size = 0.3, alpha = 0.4), dpi = 200) +
  scale_colour_manual(values = c("0"="grey80","1"="grey62"), guide="none") +
  ggnewscale::new_scale_colour() +
  rasterise(geom_point(data = M[!is.na(rid)], aes(gx, y, colour = factor(rid)), size = 1.1), dpi=200) +
  scale_colour_manual(values = CV, guide = "none") +
  geom_hline(yintercept = -log10(ALPHA), linetype = "dashed", linewidth = 0.45, colour = "grey20") +
  annotate("text", x = 0, y = -log10(ALPHA) + 0.09, hjust = -0.02, size = 3.1, colour = "grey20",
           label = sprintf("single-marker FDR %.2f -- %d of %s markers reach it",
                           ALPHA, sum(qm <= ALPHA), format(length(qm), big.mark = ","))) +
  geom_hline(yintercept = -0.10, linewidth = 0.25, colour = "grey75") +
  geom_point(data = RB, aes(gx, YR, shape = status, fill = status), size = 2.1, stroke = 0.4) +
  scale_shape_manual(values = c("maps to an EcoPeak" = 24, "does not" = 25), name = NULL) +
  scale_fill_manual(values = c("maps to an EcoPeak" = "grey15", "does not" = "white"), name = NULL) +
  annotate("segment", x = edax, xend = edax, y = 1.98, yend = 1.83,
           arrow = arrow(length = unit(0.16,"cm")), linewidth = 0.4) +
  annotate("text", x = edax, y = 2.02, label = "italic(Eda)", parse = TRUE, size = 3.4) +
  scale_x_continuous(breaks = OFF$off + OFF$mx/2, labels = chr_lab, expand = c(0.005, 0)) +
  scale_y_continuous(limits = c(-0.52, 2.05), breaks = seq(0, 2, 0.5)) +
  labs(x = "chromosome", y = expression(-log[10](q)*"  (single marker)"),
       title = sprintf("LD-complexity reduction recovers %d regions where single-marker testing recovers almost nothing",
                       nrow(R))) +
  theme_bw(12) +
  theme(panel.grid = element_blank(), legend.position = c(0.055, 0.80),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8),
        plot.title = element_text(size = 12))

OUT_DIR <- file.path(PATHS$out, STAGE); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_PDF <- file.path(OUT_DIR, "emmax_consensus_manhattan.pdf")
ggsave(OUT_PDF, p, width = 15, height = 4.8, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "emmax_consensus_manhattan.png"), p, width = 15, height = 4.8, dpi = 220)
say("\n[4] wrote %s\n", OUT_PDF)
write_receipt(STAGE, inputs = file.path(PATHS$out, "03_scan", "_receipt.rds"),
             params = list(), outputs = OUT_PDF)
say("    receipt: %s\n", receipt_path(STAGE))
