## module_sim/R_figures/fig_precision_recall_scatter.R
##
## Precision-vs-recall scatter, medians per (tag, arm) -- mirrors the
## reference figure's Panel A. Circles = nobgs, triangles = bgs, coloured by
## method.
##
## [!] FIXED 2026-09-06 (external audit, item 8): the five methods used to be
## joined by a line per tag. Dropped -- ARM_LEVELS is a reduction-degree
## ORDERING (clustered methods first, single-SNP last), not a continuous
## trajectory a system moves along, so a connecting line implied a
## relationship (interpolation, a shared axis) the five points don't have.
##
## RESTRICTED TO HIGH-DISPERSAL CELLS, same scope as fig_fbeta.R (PK, reading
## fig_pooled_pr.R's Precision*Recall row: "for everything but high
## dispersal the results are useless regardless, so no point in comparing
## there").
suppressMessages({library(data.table); library(ggplot2)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim"), "R", "00_config.R"))
say("=== fig_precision_recall_scatter ===\n\n")

POOL_PATH <- file.path(PATHS$out, "05_pool", "pooled_pr.rds")
if (!file.exists(POOL_PATH)) stop("R/05_pool.R has not produced: ", POOL_PATH)
pooled <- readRDS(POOL_PATH)$pooled

HIGH_DISP_CELLS <- c("V0.5_c1", "V1_c1", "V2_c1")
pr <- pooled[cell %in% HIGH_DISP_CELLS, .(tag, cell, arm, Precision, Recall)]
say("[0] %d (tag,cell,arm) panels, cells restricted to high dispersal: %s\n",
    nrow(pr), paste(HIGH_DISP_CELLS, collapse = ", "))

ARM_LEVELS <- c("emmax_consensus", "emmax_simes", "lfmm_simes", "emmax_snp", "lfmm_snp")
## [!] FIXED 2026-09-06 (external audit, item 8): was "EMMAX representative" --
## the arm uses the consensus dosage across the cluster's markers, not a
## single representative SNP.
ARM_LABELS <- c(emmax_consensus = "EMMAX consensus", emmax_simes = "EMMAX Simes",
                lfmm_simes = "LFMM Simes", emmax_snp = "EMMAX single-SNP", lfmm_snp = "LFMM single-SNP")
ARM_COLOURS <- c(emmax_consensus = "#1565C0", emmax_simes = "#26A69A",
                 lfmm_simes = "#7B1FA2", emmax_snp = "#F9A825", lfmm_snp = "#C0392B")

med <- pr[, .(Precision = median(Precision), Recall = median(Recall)), by = .(tag, arm)]
say("[1] median over %d high-dispersal cells, by (tag, arm)\n", length(HIGH_DISP_CELLS))
med[, arm := factor(arm, levels = ARM_LEVELS)]
setorder(med, tag, arm)
med[, arm_label := factor(ARM_LABELS[as.character(arm)], levels = ARM_LABELS[ARM_LEVELS])]
med[, tag := factor(tag, levels = c("nobgs", "bgs"))]
names(ARM_COLOURS) <- ARM_LABELS[names(ARM_COLOURS)]

p <- ggplot(med, aes(Recall, Precision, colour = arm_label)) +
  geom_point(aes(shape = tag), size = 3.5) +
  scale_colour_manual(values = ARM_COLOURS, name = "method") +
  scale_shape_manual(values = c(nobgs = 16, bgs = 17), name = NULL) +
  labs(x = "recall", y = "precision",
      title = "Clustering trades recall for precision",
      subtitle = sprintf("Medians over %d high-dispersal cells per (tag, method). Circles = nobgs, triangles = bgs.",
                          length(HIGH_DISP_CELLS))) +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  theme_bw(11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top", legend.box = "vertical")

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "fig_precision_recall_scatter.pdf")
ggsave(OUT, p, width = 8, height = 7, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 8, height = 7, dpi = 200)
say("\n[2] wrote %s (+ .png)\n", OUT)
