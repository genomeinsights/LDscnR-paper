## Figures for the LD-clustering result: where the flagged regions go, and why
## the single-SNP scan's extra output is mostly noise.
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
R   <- "module_sim_LDscnR/results/filter_then_test"
OUT <- "module_sim_LDscnR/figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
d <- fread(file.path(R, "snp_vs_cluster_dedup_allpanels.csv"))
d[, `:=`(raw_PR = raw_prec*raw_rec, dedup_PR = dedup_prec*dedup_rec)]
lab <- c(snp = "single SNP", rep = "representative", simes = "Simes", emlg = "eMLG")
d[, an := factor(lab[analysis], levels = rev(lab))]
pal <- c("single SNP" = "#B0392B", "representative" = "#7A6A1F",
         "Simes" = "#1F6F8B", "eMLG" = "#2E7156")

## A -- where flagged clusters go: contains a QTN / satellite / tags nothing
comp <- melt(d[, .(an, contain, tag_not_contain, neither)], id.vars = "an",
             variable.name = "kind", value.name = "n")
comp[, kind := factor(c(contain = "contains a QTN", tag_not_contain = "satellite (tags a QTN)",
                        neither = "tags nothing")[as.character(kind)],
                      levels = c("contains a QTN","satellite (tags a QTN)","tags nothing"))]
cs <- comp[, .(n = median(n)), by = .(an, kind)]
pA <- ggplot(cs, aes(n, an, fill = kind)) +
  geom_col(width = .62) +
  scale_fill_manual(values = c("contains a QTN" = "#1F6F8B",
                               "satellite (tags a QTN)" = "#C9A227",
                               "tags nothing" = "grey78"), name = NULL) +
  labs(x = "flagged stage-2 clusters (median over 80 panels)", y = NULL,
       title = "A  Most of the single-SNP scan's extra output tags no QTN at all",
       subtitle = "Satellites -- real LD neighbours of a QTN -- are a minor share of the errors in every analysis.") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.subtitle = element_text(colour = "grey30", size = 7.4))

## B -- PR by cell, paired
pb <- d[, .(PR = median(dedup_PR, na.rm = TRUE)), by = .(cell, an)]
pb[, an := factor(as.character(an), levels = as.character(lab))]
pB <- ggplot(pb, aes(PR, cell, colour = an)) +
  geom_line(aes(group = cell), colour = "grey80", linewidth = 2.2, alpha = .5) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = pal, name = NULL, limits = as.character(lab)) +
  labs(x = "PR after satellite removal (median over 20 panels per cell)", y = NULL,
       title = "B  Decisive where there is signal, indistinguishable where there is not",
       subtitle = "Simes beats the single-SNP scan in 19 of 19 panels in V0.5_c1. In the other three cells the paired tests are NULL (26 W / 22 L, all p > 0.4), not adverse.") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.subtitle = element_text(colour = "grey30", size = 7.4))

## C -- what satellite removal does
sr <- melt(d[, .(an, raw = raw_prec, dedup = dedup_prec)], id.vars = "an",
           variable.name = "scoring", value.name = "prec")
sr[, scoring := factor(c(raw = "satellites counted as errors", dedup = "satellites removed")[as.character(scoring)],
                       levels = c("satellites counted as errors","satellites removed"))]
pC <- ggplot(sr[is.finite(prec)], aes(prec, an, fill = scoring)) +
  geom_boxplot(outlier.size = .5, linewidth = .3, width = .6) +
  scale_fill_manual(values = c("grey82", "#BBD3DE"), name = NULL) +
  labs(x = "precision", y = NULL,
       title = "C  Removing satellites lifts precision but changes no ordering",
       subtitle = "One region per QTN, the evaluate_ors convention. Gains are ~0.03 (single SNP) to ~0.13 (cluster-level).") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.subtitle = element_text(colour = "grey30", size = 7.4))

ggsave(file.path(OUT, "clustering_gain.png"), pA / pB / pC, width = 9.5, height = 10, dpi = 190)
cat(sprintf("  written: %s\n", file.path(OUT, "clustering_gain.png")))
cat("\n  medians used:\n"); print(cs)
