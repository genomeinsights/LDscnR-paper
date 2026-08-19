## =====================================================================
## module_C2 / R/06_make_figures.R
## Figures 1-5 for the operating-grid stability exploration.
## Rscript module_C2/R/06_make_figures.R
## =====================================================================
source("module_C2/R/00_helpers.R")
suppressMessages({ library(ggplot2) })
th <- theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 7),
        plot.subtitle = element_text(size = 7, colour = "grey30"))

sg <- readRDS(file.path(C2$CACHE, "grid_scored.rds"))
ag <- readRDS(file.path(C2$CACHE, "anchor_grid.rds")); B <- ag$B
g  <- ag$G$rec50
NCELL <- length(C2$TAUS) * length(C2$LMINS)

## ---- Fig 1: significant-region counts across the grid ---------------
## The primary operating point tau_C = 0.05 is NOT on the prototype's grid
## (seq(0.02, 0.50, 0.02) steps over it), so its column is computed separately
## from the augmented core and drawn set apart, to make that explicit.
res <- copy(sg$resF)[, off_grid := FALSE]
core <- readRDS(file.path(C2$CACHE, "grid_core.rds"))
ct <- core$by_tau[["0.05"]]
if (!is.null(ct)) {
  op <- rbindlist(lapply(C2$LMINS, function(lm) {
    keep <- which(ct$obs$size >= lm)
    if (!length(keep)) return(data.table(tau = 0.05, lmin = lm, n_obs = 0L, n_sig = 0L))
    O <- ct$obs[keep]; pq <- c2_emp_pq(O, ct$surr[b <= B & size >= lm], B)
    data.table(tau = 0.05, lmin = lm, n_obs = nrow(O), n_sig = sum(pq$q_R < C2$FDR)) }))
  res <- rbind(res, op[, off_grid := TRUE])
}
NUSE <- sg$resF[n_sig > 0, .N]
tlev <- c(sort(unique(sg$resF$tau)), 0.05)
res[, taul := factor(tau, levels = tlev)]
f1 <- ggplot(res, aes(taul, factor(lmin, levels = C2$LMINS), fill = n_sig)) +
  geom_tile(aes(colour = off_grid), linewidth = 0.4) +
  geom_text(aes(label = ifelse(n_sig > 0, n_sig, "")), size = 2.1, colour = "grey15") +
  scale_colour_manual(values = c(`FALSE` = "white", `TRUE` = "grey20"), guide = "none") +
  scale_fill_gradientn(colours = C2$ZISSOU, name = expression(regions~q[R]<0.05)) +
  annotate("text", x = length(tlev), y = length(C2$LMINS) + 0.75, label = "operating\npoint", size = 2.2) +
  coord_cartesian(clip = "off") +
  labs(x = expression(tau[C]), y = expression(l[min]),
       title = sprintf("Significant regions per operating-grid cell (B = %d, location-matched empirical p, BH within cell)", B),
       subtitle = sprintf("|U| = %d usable of |G| = %d cells (availability %.2f); every cell producing a region produces a SIGNIFICANT one. The primary operating point tau_C = 0.05 (outlined, far right) is NOT on the grid seq(0.02, 0.50, 0.02).",
                          NUSE, length(C2$TAUS) * length(C2$LMINS), NUSE / (length(C2$TAUS) * length(C2$LMINS)))) +
  th + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
             plot.margin = margin(5, 20, 5, 5))
ggsave(file.path(C2$FIG, "fig1_grid_nsig.png"), f1, width = 12, height = 4.2, dpi = 170)

## ---- Fig 2: every anchor locus across every cell --------------------
ord <- g[, sum(sig), by = label][order(-V1)]$label
g[, state := fifelse(sig, "significant", fifelse(detected, "detected only", "absent"))]
g[, label := factor(label, levels = rev(ord))]
f2 <- ggplot(g, aes(factor(tau), label, fill = state)) +
  geom_tile(colour = "white", linewidth = 0.15) +
  facet_wrap(~ factor(lmin, levels = C2$LMINS, labels = paste0("l_min = ", C2$LMINS)), nrow = 1) +
  scale_fill_manual(values = c(absent = "grey92", `detected only` = "#78B7C5",
                               significant = "#F21A00"), name = NULL) +
  labs(x = expression(tau[C]), y = NULL,
       title = "Each prespecified anchor locus across the operating grid (marker-membership matching, recover >= 0.5)",
       subtitle = sprintf("'detected only' (called but not location-matched significant) occurs in %d of %d cells: the significance layer never fires independently of detection",
                          sum(g$detected & !g$sig), nrow(g))) +
  th + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 4),
             axis.text.y = element_text(size = 5.5), legend.position = "top")
ggsave(file.path(C2$FIG, "fig2_anchor_grid_tiles.png"), f2, width = 15, height = 4.5, dpi = 170)

## ---- Fig 3: denominator variants ------------------------------------
V <- fread(file.path(C2$RES, "C2_variant_comparison.csv"))
Vl <- melt(V[, .(label, S_U, S_G, D, Q_over_D)], id.vars = "label",
           variable.name = "variant", value.name = "score")
Vl[, label := factor(label, levels = V[order(S_G)]$label)]
Vl[, variant := factor(variant, levels = c("S_U", "S_G", "D", "Q_over_D"),
     labels = c("A: S^U (usable cells)", "B: S^G (all cells)",
                "D: detection D", "D: Q/D (null confirmation)"))]
f3 <- ggplot(Vl, aes(score, label, fill = variant)) +
  geom_col(position = "dodge", width = 0.8) +
  scale_fill_manual(values = c("#E1AF00", "#F21A00", "#78B7C5", "#3B9AB2"), name = NULL) +
  labs(x = "score", y = NULL, title = "Denominator variants on the fixed anchor loci",
       subtitle = sprintf("S^U is exactly S^G x |G|/|U| = %.2f x S^G -- a rescaling, not a re-ranking; on this dataset D == Q so Q/D == 1 for every detected locus",
                          NCELL / max(sg$resF[n_sig > 0, .N], 1L))) +
  th + theme(axis.text.y = element_text(size = 6), legend.position = "top")
ggsave(file.path(C2$FIG, "fig3_denominator_variants.png"), f3, width = 9, height = 5.5, dpi = 170)

## ---- Fig 4: detection vs significance stability ---------------------
f4 <- ggplot(V, aes(D, Q)) +
  geom_abline(slope = 1, linetype = 2, colour = "grey70") +
  geom_point(aes(size = Q_over_D, colour = Q_over_D)) +
  ggrepel::geom_text_repel(aes(label = sub("_Chr", " Chr", sub("Mb.*", "", label))),
                           size = 2.1, max.overlaps = 30, segment.size = 0.2) +
  scale_colour_gradientn(colours = C2$ZISSOU, name = "Q/D", limits = c(0, 1)) +
  scale_size_continuous(range = c(1, 4), guide = "none") +
  labs(x = "detection stability  D  (called as a region)",
       y = "significance stability  Q  (called AND location-matched significant)",
       title = "Detection vs significance stability: the two are not interchangeable",
       subtitle = "points on the diagonal are always-significant-when-found; points far below are region-calling-robust but null-fragile") +
  th
ggsave(file.path(C2$FIG, "fig4_detection_vs_significance.png"), f4, width = 7, height = 5.5, dpi = 170)

## ---- Fig 5: grid sensitivity ----------------------------------------
GS <- fread(file.path(C2$RES, "grid_sensitivity.csv"))
gcols <- setdiff(names(GS), "label")
RK <- copy(GS); for (n in gcols) RK[[n]] <- frank(-GS[[n]], ties.method = "min")
RL <- melt(RK, id.vars = "label", variable.name = "grid", value.name = "rank")
RL[, label := factor(label, levels = GS[order(G1_original)]$label)]
f5 <- ggplot(RL, aes(grid, rank, group = label, colour = label)) +
  geom_line(alpha = 0.7) + geom_point(size = 1.4) +
  scale_y_reverse(breaks = seq(1, nrow(GS), 2)) +
  labs(x = NULL, y = "rank (1 = most stable)",
       title = "Anchor-locus ranking under six prespecified grids",
       subtitle = "flat bundle = grid-robust ranking. A column collapsing to rank 1 is fully TIED: under G5 (tau_C > 0.26) every locus scores 0, and under G2 (coarse) 13 of 17 do -- the ranking is not merely re-ordered, it is destroyed") +
  th + theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1, size = 6))
ggsave(file.path(C2$FIG, "fig5_grid_sensitivity_ranks.png"), f5, width = 8, height = 5, dpi = 170)

c2_msg("[fig] wrote 5 figures to %s\n", C2$FIG)
