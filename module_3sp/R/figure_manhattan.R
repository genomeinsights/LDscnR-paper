## =============================================================================
## module_3sp/R/08_figure_both.R
##
## PK: both engines, EACH using 07_figure_emmax.R's template (single-marker
## -log10(q), FDR line, cluster-level regions coloured, Eda arrow) -- stacked,
## not the two-scale-per-panel design 05_manhattan.R used for the
## engine-comparison figure. Revised for manuscript use (PK, second pass):
##   - no overall plot title/subtitle -- that text belongs in the LaTeX figure
##     legend (caption), not baked into the image; a suggested caption is
##     printed by this script rather than guessed at silently
##   - no "single-marker FDR..." text above the dashed line -- it did not
##     render legibly at figure scale; the dashed line itself is unchanged
##   - EcoPeak status as a RUG (short tick per region, black = maps to an
##     EcoPeak, grey = does not), replacing the triangle points + shape/fill
##     legend
##   - y-scale FREE per panel, anchored by the FDR line at the same data value
##     in both rather than a shared numeric limit
##   - a small hand-built legend for the rug colours, EMMAX panel only, top
##     right, with extra headroom reserved above EMMAX's own data range so it
##     does not overlap points; LFMM keeps its data range tight
##
## EACH PANEL IS SELF-CONTAINED: its own single-marker q, its own significant
## regions, its own EcoPeak status, its own colour rotation. Colours are NOT
## shared across panels -- there is no joint-region step here, unlike
## 05_manhattan.R, because PK asked for "the above template" applied twice,
## not the joint-colouring scheme. If a shared-colour version is wanted later,
## that is 05_manhattan.R's job, not this one's.
##
## ONE COMBINED COLOUR SCALE PER PANEL for the point cloud, not two stacked via
## ggnewscale::new_scale_colour() -- the two-scale version produced a real
## ggplot2 aesthetics-length error under patchwork's panel combination (see
## 05_manhattan.R's header). The EcoPeak rug uses FIXED colours (two separate
## geom_segment layers, not a mapped/scaled aesthetic) for the same reason: a
## second scale on top of the point cloud's would reintroduce that risk.
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(ggplot2)
  library(patchwork); library(ggrastr)})
devtools::load_all("~/gitlab/LDscnR")
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_3sp"), "R", "00_config.R"))
STAGE <- "08_figure_both"
say("=== %s ===\n\n", STAGE)

b  <- readRDS(file.path(PATHS$out, "02_bundle", "bundle.rds"))
sc <- readRDS(file.path(PATHS$out, "03_EMMAX", "scan.rds"))
lf <- readRDS(file.path(PATHS$out, "04_lfmm", "lfmm.rds"))
map <- b$map
EDA <- list(Chr = "Chr4", from = 12.70e6, to = 12.90e6)

## per-marker q for each engine's point cloud
pm_emmax <- emmax_fast(emmax_setup(b$GTs, b$GRM), b$eco)
Q <- list(EMMAX = p.adjust(pm_emmax, "BH"), LFMM = p.adjust(lf$lfmm_p, "BH"))

## FREE y-scale per panel (PK): each panel ranges over its own data, so neither engine's
## points get clipped (LFMM reaches -log10(q) = 3.76, well past EMMAX's ~1.85) and neither
## is squashed to share the other's range. The FDR dashed line, drawn at the same DATA
## value (-log10(0.05)) in both panels, is what anchors them to a common reference despite
## the differing axis maxima -- not a shared numeric limit.

## regions + EcoPeak status per engine, REUSED from 03_EMMAX.R/04_lfmm.R
REG <- list(EMMAX = copy(sc$consensus$test$regions), LFMM = copy(lf$test$regions))
ROT <- list(EMMAX = sc$consensus$rotation, LFMM = lf$rotation)
for (eng in names(REG)) {
  onp <- ROT[[eng]]$per_region; setkey(onp, chr, from, to)
  setkey(REG[[eng]], Chr, from, to)
  REG[[eng]][, on_peak := onp[.(Chr, from, to)]$on_peak]
  REG[[eng]][, rid := .I]
}
say("[1] EMMAX: %d regions, %d on EcoPeak (fold %.2fx, p=%.4f)\n", nrow(REG$EMMAX),
    sum(REG$EMMAX$on_peak), ROT$EMMAX$fold, ROT$EMMAX$p)
say("[2] LFMM:  %d regions, %d on EcoPeak (fold %.2fx, p=%.4f)\n", nrow(REG$LFMM),
    sum(REG$LFMM$on_peak), ROT$LFMM$fold, ROT$LFMM$p)

OFF <- map[, .(mx = max(Pos)), by = Chr][order(Chr)][, off := cumsum(c(0, head(mx,-1)) + 2e6)][]
chr_lab <- as.integer(gsub("Chr", "", OFF$Chr))
BG <- c("grey80", "grey62")

## known marine-freshwater inversions (Jones et al. 2012; cited by Roberts Kingman et al.
## 2021, i.e. this repo's own kingman2021/) -- chrI/1 (442 kb), chrXI/11 (412 kb), chrXXI/21
## (1,700 kb). PK checked our own regions against Jones et al. 2012's figure by eye: none of
## our chrXI/chrXXI regions coincide with the published inversions there, so those two are
## dropped from the figure entirely (stated in text instead, not annotated on the plot).
## Chr1 is the exception -- our largest EMMAX region (21.49-21.94 Mb) lands almost exactly on
## the informal "Chr1:21.50" reference kingman2021/R/14_venn_emmax17.R hardcodes for this
## inversion -- so it gets an Eda-style arrow annotation (a shaded band read as a line at
## genome-wide scale; an arrow is legible at this width the way the Eda one already is).
INV1 <- list(Chr = "Chr1", pos = 21.5e6, label = "chr1 inversion")

## top_expand: extra headroom ABOVE the data max, as a fraction of the data range. EMMAX
## gets more (room for the rug legend, drawn only on that panel); LFMM stays tight.
## show_legend: draw the black/grey EcoPeak-rug legend in this panel's own top-right
## corner. A SECOND ggplot colour SCALE for the rug (so its own legend could be built
## automatically) is deliberately avoided -- see the header note on why a second
## colour scale broke under patchwork here before. The rug is two FIXED-colour
## geom_segment layers instead (no scale, so no conflict), and its legend is built by
## hand from two short annotate() swatches plus text, positioned from the panel's own
## data range so it sits in the corner regardless of that range's actual size.
mk_panel <- function(eng, lab, top_expand = 0.08, show_legend = FALSE) {
  R <- REG[[eng]]; qv <- Q[[eng]]
  PAL <- LDscnR:::default_cluster_colours(); R[, col := NA_character_]; cur <- 0L
  for (ch in unique(R$Chr)) { used <- character(0)
    for (i in R[Chr == ch, which = TRUE]) {
      repeat { cur <- cur %% length(PAL) + 1L; if (!(PAL[cur] %in% used)) break }
      R$col[i] <- PAL[cur]; used <- c(used, PAL[cur]) } }
  stopifnot(R[, .(ok = uniqueN(col) == .N), by = Chr][, all(ok)])

  M <- data.table(Chr = map$Chr, pos = map$Pos, y = -log10(pmax(qv, .Machine$double.xmin)))
  setkey(R, Chr, from, to)
  M[, rid := foverlaps(M[, .(Chr, from = pos, to = pos)], R,
        by.x = c("Chr","from","to"), type = "within", mult = "first", nomatch = NA)$rid]
  M <- merge(M, OFF[, .(Chr, off)], by = "Chr")[, gx := pos + off]
  chr_n <- as.integer(gsub("Chr", "", M$Chr))
  CV <- setNames(R$col, as.character(R$rid))
  CV_all <- c(CV, setNames(BG, c("..bg0", "..bg1")))
  M[, cc := ifelse(!is.na(rid), as.character(rid), ifelse(chr_n %% 2 == 0, "..bg0", "..bg1"))]
  M[, ord := !is.na(rid)]; setorder(M, ord)

  y_max <- max(M$y)
  x_rng <- range(M$gx)
  edax  <- mean(c(EDA$from, EDA$to)) + OFF[Chr == EDA$Chr, off]
  invx  <- INV1$pos + OFF[Chr == INV1$Chr, off]

  ## rug: one tick per region at its genomic midpoint, on TWO ROWS (PK: overlapping ticks in
  ## one row were hard to tell apart) -- on_peak == TRUE above, FALSE below, so the two
  ## colours never compete for the same y position even when regions share an x position.
  ##
  ## Sized as a fraction of the panel's FULL RENDERED SPAN, not of y_max alone (PK: rug
  ## height differed visibly between panels). BOTTOM_MULT here must match the bottom
  ## expand() fraction below -- the two panels have different top_expand (EMMAX reserves
  ## headroom for its legend, LFMM doesn't), so a size defined as a fraction of y_max alone
  ## is implicitly a SMALLER fraction of the actual rendered height on whichever panel has
  ## the larger top_expand, which is exactly the mismatch PK saw.
  BOTTOM_MULT <- 0.14
  full_span <- y_max * (1 + top_expand + BOTTOM_MULT)
  RB <- merge(R, OFF[, .(Chr, off)], by = "Chr")[, gx := (from + to)/2 + off]
  rug_h       <- 0.018 * full_span
  rug_y_true  <- -0.028 * full_span
  rug_y_false <- -0.066 * full_span

  ## legend swatch geometry, corner-anchored from this panel's own data range so it
  ## does not need to know the range in advance
  leg_x0 <- x_rng[2] - 0.14 * diff(x_rng); leg_x1 <- leg_x0 + 0.018 * diff(x_rng)
  leg_y1 <- y_max * (1 + top_expand * 0.80); leg_y2 <- y_max * (1 + top_expand * 0.55)

  p <- ggplot(M) +
    rasterise(geom_point(aes(gx, y, colour = cc, size = ord, alpha = ord)), dpi = 200) +
    scale_colour_manual(values = CV_all, guide = "none") +
    scale_size_manual(values = c("TRUE" = 1.1, "FALSE" = 0.3), guide = "none") +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4), guide = "none") +
    geom_hline(yintercept = -log10(ALPHA), linetype = "dashed", linewidth = 0.45, colour = "grey20") +
    geom_segment(data = RB[on_peak == TRUE],
                 aes(x = gx, xend = gx, y = rug_y_true - rug_h, yend = rug_y_true + rug_h),
                 colour = "black", linewidth = 0.6) +
    geom_segment(data = RB[on_peak == FALSE],
                 aes(x = gx, xend = gx, y = rug_y_false - rug_h, yend = rug_y_false + rug_h),
                 colour = "grey70", linewidth = 0.6) +
    annotate("segment", x = edax, xend = edax, y = y_max * (1 + top_expand*0.15),
             yend = y_max * 0.98, arrow = arrow(length = unit(0.16,"cm")), linewidth = 0.4) +
    annotate("text", x = edax, y = y_max * (1 + top_expand*0.28),
             label = "italic(Eda)", parse = TRUE, size = 3.4) +
    annotate("segment", x = invx, xend = invx, y = y_max * (1 + top_expand*0.15),
             yend = y_max * 0.98, arrow = arrow(length = unit(0.16,"cm")), linewidth = 0.4) +
    annotate("text", x = invx, y = y_max * (1 + top_expand*0.28),
             label = INV1$label, size = 3.0) +
    scale_x_continuous(breaks = OFF$off + OFF$mx/2, labels = chr_lab, expand = c(0.005, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0.14, top_expand))) +
    labs(x = NULL, y = expression(-log[10](q)*"  (single marker)"), title = lab) +
    theme_bw(12) +
    theme(panel.grid = element_blank(), plot.title = element_text(size = 12))

  if (show_legend) p <- p +
    annotate("segment", x = leg_x0, xend = leg_x1, y = leg_y1, yend = leg_y1,
             colour = "black", linewidth = 0.6) +
    annotate("text", x = leg_x1, y = leg_y1, label = "  maps to an EcoPeak",
             hjust = 0, vjust = 0.5, size = 3.0, colour = "grey20") +
    annotate("segment", x = leg_x0, xend = leg_x1, y = leg_y2, yend = leg_y2,
             colour = "grey70", linewidth = 0.6) +
    annotate("text", x = leg_x1, y = leg_y2, label = "  does not",
             hjust = 0, vjust = 0.5, size = 3.0, colour = "grey20")
  p
}

pA <- mk_panel("EMMAX", sprintf("A  EMMAX, consensus test -- %d regions (%d on EcoPeak)",
                                nrow(REG$EMMAX), sum(REG$EMMAX$on_peak)),
               top_expand = 0.32, show_legend = TRUE)
pB <- mk_panel("LFMM", sprintf("B  LFMM, Simes -- %d regions (%d on EcoPeak)",
                               nrow(REG$LFMM), sum(REG$LFMM$on_peak)),
               top_expand = 0.10, show_legend = FALSE) + labs(x = "chromosome")
FIG <- pA / pB

OUT_DIR <- file.path(PATHS$out, STAGE); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_PDF <- file.path(OUT_DIR, "both_engines_manhattan.pdf")
ggsave(OUT_PDF, FIG, width = 15, height = 8.4, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "both_engines_manhattan.png"), FIG, width = 15, height = 8.4, dpi = 200)
say("\n[3] wrote %s\n", OUT_PDF)
say("\n[4] suggested LaTeX figure legend (the title text removed from the image):\n")
say('    LD-complexity reduction recovers regions where single-marker testing\n')
say('    recovers almost nothing, for both association engines. Each panel plots\n')
say('    single-marker $-\\log_{10}(q)$ against the dashed FDR line at $\\alpha$ =\n')
say('    %.2f; points inside a significant cluster are coloured by region, one\n', ALPHA)
say('    colour never repeated within a chromosome. The rug below each axis marks\n')
say('    every region, black where it maps to a published EcoPeak and grey where it\n')
say('    does not (legend, panel A). A: EMMAX, consensus test, %d regions, %d on an\n',
    nrow(REG$EMMAX), sum(REG$EMMAX$on_peak))
say('    EcoPeak. B: LFMM, Simes (its p-values are precomputed, so no consensus\n')
say('    test is available for this engine), %d regions, %d on an EcoPeak.\n',
    nrow(REG$LFMM), sum(REG$LFMM$on_peak))
say('    The chr1 inversion arrow marks one of three chromosomal inversions long\n')
say('    known to distinguish marine and freshwater ecotypes (Jones et al. 2012;\n')
say('    chrI, chrXI, chrXXI -- also recovered as LD-clusters 6, 22 and 29 by Fang\n')
say('    et al. 2020). Our largest EMMAX region (chr1:21.49-21.94 Mb) coincides\n')
say('    closely with it. Text, not figure: none of our chrXI or chrXXI regions\n')
say('    coincided with the published locations of the other two inversions\n')
say('    (checked visually against Jones et al. 2012); we did not source precise\n')
say('    breakpoints for any of the three (Jones et al. 2012 Supplementary Table\n')
say('    10 was not accessible to us) so this comparison is approximate.\n')
write_receipt(STAGE, inputs = c(file.path(PATHS$out, "03_EMMAX", "_receipt.rds"),
                                file.path(PATHS$out, "04_lfmm", "_receipt.rds")),
             params = list(), outputs = OUT_PDF)
say("    receipt: %s\n", receipt_path(STAGE))
