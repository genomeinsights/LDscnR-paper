## module_sim/R_figures/fig_fp_neutral_chr.R
##
## Visualize R/05_pool.R's neutral-chromosome (Chr2) FP analysis. PK:
## "analyse the neutral chromosomes separately (only FPs of course)." Chr2
## never carries a QTN (R/04_score.R), so every significant region there is
## an assumption-free false positive -- an empirical false-discovery count
## under the real genotypes/kinship/spatial structure, no permutation
## needed. Two panels: (a) mean FP count per environment, by cell and arm;
## (b) the same, broken down by cluster size (mirrors fig_fp_by_size, but a
## raw count here since the FP proportion on Chr2 is trivially 1 everywhere).
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim"), "R", "00_config.R"))
say("=== fig_fp_neutral_chr ===\n\n")

POOL_PATH <- file.path(PATHS$out, "05_pool", "pooled_pr.rds")
if (!file.exists(POOL_PATH)) stop("R/05_pool.R has not produced: ", POOL_PATH)
pool <- readRDS(POOL_PATH)
fp_neutral <- pool$fp_neutral
fp_neutral_size <- pool$fp_neutral_size

## Same cell labels/ordering and arm colours as fig_pooled_pr.R / fig_fp_by_size.R.
SI_MAP   <- c(V0.5 = "strong", V1 = "medium", V2 = "weak")
DISP_MAP <- c(c1 = "high", c1.5 = "medium", c2 = "low")
fp_neutral[, c("V_raw", "c_raw") := tstrsplit(cell, "_", fixed = TRUE)]
fp_neutral[, SI   := factor(SI_MAP[V_raw], levels = c("strong", "medium", "weak"))]
fp_neutral[, disp := factor(DISP_MAP[c_raw], levels = c("high", "medium", "low"))]
fp_neutral[, cell_label := sprintf("%s SI, %s disp", SI, disp)]
cell_order <- unique(fp_neutral[, .(cell, cell_label, disp, SI)])
setorder(cell_order, disp, SI)
fp_neutral[, cell_label := factor(cell_label, levels = cell_order$cell_label)]

ARM_LEVELS <- c("emmax_consensus", "emmax_simes", "lfmm_simes",
                "emmax_snp", "emmax_snp_clustered", "lfmm_snp", "lfmm_snp_clustered")
ARM_COLOURS <- c(emmax_consensus = "#1565C0", emmax_simes = "#26A69A", lfmm_simes = "#7B1FA2",
                 emmax_snp = "#F9A825", emmax_snp_clustered = "#FFCC80",
                 lfmm_snp = "#C0392B", lfmm_snp_clustered = "#EF9A9A")
fp_neutral[, arm := factor(arm, levels = ARM_LEVELS)]
fp_neutral[, tag := factor(tag, levels = c("nobgs", "bgs"))]
fp_neutral_size[, arm := factor(arm, levels = ARM_LEVELS)]
fp_neutral_size[, tag := factor(tag, levels = c("nobgs", "bgs"))]
SIZE_LABELS <- c("1", "2", "3", "4-5", "6-10", "11-20", "21-50", "50+")
fp_neutral_size[, size_bin := factor(size_bin, levels = SIZE_LABELS)]

say("[1] %d (tag,cell,arm) points, mean Chr2 FP count range %.3f-%.2f\n",
    nrow(fp_neutral), min(fp_neutral$mean_FP), max(fp_neutral$mean_FP))

p1 <- ggplot(fp_neutral, aes(cell_label, mean_FP, colour = arm, group = arm)) +
  geom_line(alpha = 0.5, linewidth = 0.4, position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(ymin = pmax(0, mean_FP - FP_SE), ymax = mean_FP + FP_SE),
                  position = position_dodge(width = 0.4), size = 0.3) +
  facet_wrap(~tag, nrow = 1) +
  scale_colour_manual(values = ARM_COLOURS, name = "method") +
  labs(x = NULL, y = "mean FP count per environment (Chr2 only)",
      title = "Neutral chromosome: false discoveries with zero QTN present") +
  theme_bw(11) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 7, angle = 40, hjust = 1), legend.position = "top")

p2 <- ggplot(fp_neutral_size, aes(size_bin, mean_FP, colour = arm, group = arm)) +
  geom_line(alpha = 0.5, linewidth = 0.4, position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(ymin = pmax(0, mean_FP - FP_SE), ymax = mean_FP + FP_SE),
                  position = position_dodge(width = 0.4), size = 0.3) +
  facet_wrap(~tag, nrow = 1) +
  scale_colour_manual(values = ARM_COLOURS, name = "method", guide = "none") +
  scale_y_log10() +
  labs(x = "stage-1 cluster size (n markers)", y = "mean FP count per environment (log scale)",
      title = "Same, by cluster size (pooled across cells)",
      subtitle = "all Chr2 significant regions are FP by construction (no QTN on Chr2) -- no permutation, an assumption-free empirical false-discovery count") +
  theme_bw(11) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8))

p <- p1 / p2 + plot_layout(heights = c(1, 1))

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "fig_fp_neutral_chr.pdf")
ggsave(OUT, p, width = 10, height = 10, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 10, height = 10, dpi = 200)
say("\n[2] wrote %s (+ .png)\n", OUT)
