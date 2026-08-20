## =====================================================================
## module_C2 / R/12_manhattan.R    [Question 5]
##
## Two companion Manhattan figures on an identical chromosome layout and a
## comparable colour scale:
##   A  anchored : markers inside the (0.05, 3) REFERENCE-POINT regions coloured by
##                 that region's detection support across the null-admissible grid;
##   B  anchor-free : every marker coloured by its own marker-level support S_m.
## Other markers light grey in both. Labelled and unlabelled versions of each.
##
## The colour is DETECTION support. It is not called a second-tier significance
## score, because significance is exactly redundant with detection here (0
## detected-but-not-significant cells in all 19 configurations; see 10_support.log).
##
## Rscript module_C2/R/12_manhattan.R
## =====================================================================
source("module_C2/R/00_helpers.R")
suppressMessages({ library(ggplot2); library(ggrepel) })
core <- readRDS(file.path(C2$CACHE, "grid_core.rds"))
KEEP <- readRDS(file.path(C2$CACHE, "admissible_cells.rds"))
SM   <- readRDS(file.path(C2$CACHE, "marker_support.rds"))
SUP  <- readRDS(file.path(C2$CACHE, "support.rds"))
anc  <- c2_anchors(core)
map  <- as.data.table(core$D$map)
C_obs <- core$D$C_obs

chr_lev <- paste0("Chr", sort(unique(as.integer(gsub("Chr", "", map$Chr)))))
BASE <- map[, .(marker, Chr = factor(Chr, levels = chr_lev), Pos)]
BASE[, C := ifelse(marker %in% names(C_obs), C_obs[marker], 0)][is.na(C), C := 0]

th <- theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 7),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.spacing.x = unit(0.05, "lines"),
        plot.subtitle = element_text(size = 7, colour = "grey30", lineheight = 1.15))

## a shared colour scale so the two panels are directly comparable
LIMS <- c(0, 0.65)

manh <- function(dt, lab, title, sub, labels = NULL, lims = LIMS) {
  g <- ggplot(dt, aes(Pos, C)) +
    geom_point(data = dt[is.na(sup)], colour = "grey82", size = 0.45, alpha = 0.55) +
    geom_point(data = dt[!is.na(sup)], aes(colour = sup), size = 1.6) +
    facet_wrap(~ Chr, nrow = 1, scales = "free_x") +
    scale_colour_gradientn(colours = C2$ZISSOU, name = lab, limits = lims,
                           oob = scales::squish) +
    labs(x = "Genomic position", y = "C-score", title = title,
         subtitle = paste(strwrap(sub, width = 190), collapse = "\n")) + th
  if (!is.null(labels) && nrow(labels))
    g <- g + geom_text_repel(data = labels, aes(Pos, C, label = tag), inherit.aes = FALSE,
                             size = 2.1, min.segment.length = 0, segment.size = 0.2,
                             box.padding = 0.4, max.overlaps = 30, colour = "grey15")
  g
}

for (GN in c("P20", "P10")) {
  ncell <- nrow(KEEP[[GN]])
  eps <- if (GN == "P20") 0.20 else 0.10

  ## ---- A. anchored: reference-point regions coloured by their support ----
  D <- SUP[[paste0(GN, ".rec50")]][, .(sup = sum(detected) / ncell), by = label]
  dtA <- copy(BASE)[, sup := NA_real_]
  for (i in seq_len(nrow(D))) {
    mk <- anc$mk[[D$label[i]]]
    dtA[marker %in% mk, sup := D$sup[i]]
  }
  labA <- rbindlist(lapply(c("*inv", "*Eda"), function(p) {
    lb <- grep(p, D$label, fixed = TRUE, value = TRUE); if (!length(lb)) return(NULL)
    mk <- anc$mk[[lb]]; d <- dtA[marker %in% mk][which.max(C)]
    d[, tag := if (p == "*inv") "Chr1 inversion" else "Eda"][, .(Chr, Pos, C, tag)] }))
  ## the strongest-supported region also labelled
  top <- D[which.max(sup)]
  dtop <- dtA[marker %in% anc$mk[[top$label]]][which.max(C)]
  labA <- rbind(labA, dtop[, .(Chr, Pos, C,
    tag = sprintf("top support %.2f", top$sup))])

  subA <- sprintf("colour = fraction of the %d null-admissible cells (Pr0(N>0) <= %.2f) in which the region is recovered; grey = all other markers. The reference point supplies the displayed boundaries only -- it does not define the admissible grid or the evidence.",
                  ncell, eps)
  gA0 <- manh(dtA, "support", "Reference-point regions coloured by support across the null-admissible operating grid", subA)
  gA1 <- manh(dtA, "support", "Reference-point regions coloured by support across the null-admissible operating grid", subA, labA)
  ggsave(sprintf("%s/fig13_manhattan_reference_support_%s.png", C2$FIG, GN), gA0, width = 18, height = 4.5, dpi = 170)
  ggsave(sprintf("%s/fig13_manhattan_reference_support_%s_labelled.png", C2$FIG, GN), gA1, width = 18, height = 4.5, dpi = 170)
  ## region support maxes far below marker support, so the shared scale flattens this
  ## panel; an own-scale copy is emitted for reading the anchored structure itself.
  ownl <- c(0, max(D$sup))
  subAo <- sprintf("%s SCALE DIFFERS from the anchor-free panel (0-%.2f here vs 0-%.2f there): a region must be recovered >=50%% intact, a marker need only fall in ANY region.",
                   subA, ownl[2], LIMS[2])
  gA2 <- manh(dtA, "support", "Reference-point regions coloured by support across the null-admissible operating grid", subAo, labA, lims = ownl)
  ggsave(sprintf("%s/fig13_manhattan_reference_support_%s_ownscale.png", C2$FIG, GN), gA2, width = 18, height = 4.5, dpi = 170)

  ## ---- B. anchor-free: every marker by its own S_m ----------------------
  s <- SM[[GN]]
  dtB <- copy(BASE)[, sup := NA_real_]
  dtB[marker %in% names(s), sup := s[marker]]
  labB <- dtB[!is.na(sup)][order(-sup)][1:3][, tag := sprintf("S=%.2f", sup)][, .(Chr, Pos, C, tag)]
  subB <- sprintf("colour = S_m, the fraction of the same %d admissible cells in which the marker lies in ANY empirical region. No anchor set is used. Identical layout and colour scale to the panel above; S_m is operating-grid support, not independent significance (Spearman with the marker C-score = %.2f).",
                  ncell, stats::cor(s, C_obs[names(s)], method = "spearman"))
  gB0 <- manh(dtB, expression(S[m]), "Anchor-free marker-level support across the null-admissible operating grid", subB)
  gB1 <- manh(dtB, expression(S[m]), "Anchor-free marker-level support across the null-admissible operating grid", subB, labB)
  ggsave(sprintf("%s/fig14_manhattan_marker_support_%s.png", C2$FIG, GN), gB0, width = 18, height = 4.5, dpi = 170)
  ggsave(sprintf("%s/fig14_manhattan_marker_support_%s_labelled.png", C2$FIG, GN), gB1, width = 18, height = 4.5, dpi = 170)
  c2_msg("[%s] wrote 4 figures (anchored + anchor-free, labelled + plain); %d cells, max region support %.3f, max S_m %.3f\n",
         GN, ncell, max(D$sup), max(s))
}
c2_msg("[done] figures in %s\n", C2$FIG)
