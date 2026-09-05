## module_sim/R_figures/fig_fp_by_size.R
##
## Visualize R/05_pool.R's fp_by_size: pooled (sum first, then divide) FP
## proportion among significant stage-1 clusters/units, mean +- SE across the
## 10 environments (the replicate axis), as a function of cluster size. PK:
## "what is the proportion of FP as a function of cluster size" / "the main
## comparison is between methods, not between bgs" -- coloured by arm
## (method), faceted by tag.
suppressMessages({library(data.table); library(ggplot2)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim"), "R", "00_config.R"))
say("=== fig_fp_by_size ===\n\n")

POOL_PATH <- file.path(PATHS$out, "05_pool", "pooled_pr.rds")
if (!file.exists(POOL_PATH)) stop("R/05_pool.R has not produced: ", POOL_PATH)
fp <- readRDS(POOL_PATH)$fp_by_size

## [!] FIXED 2026-09-06: "1" added -- single-SNP arms now produce genuine
## size-1 regions (R/04_score.R), and R/05_pool.R's bins now start at 0 so
## n_loci==1 lands in its own bin instead of falling into cut()'s NA gap.
SIZE_LABELS <- c("1", "2", "3", "4-5", "6-10", "11-20", "21-50", "50+")
ARM_LEVELS <- c("emmax_consensus", "emmax_simes", "lfmm_simes", "emmax_snp", "lfmm_snp")
ARM_COLOURS <- c(emmax_consensus = "#1565C0", emmax_simes = "#26A69A",
                 lfmm_simes = "#7B1FA2", emmax_snp = "#F9A825", lfmm_snp = "#C0392B")
fp[, size_bin := factor(size_bin, levels = SIZE_LABELS)]
fp[, arm := factor(arm, levels = ARM_LEVELS)]
fp[, tag := factor(tag, levels = c("nobgs", "bgs"))]

say("[1] %d (tag,arm,size_bin) points, n_sig range %d-%d\n", nrow(fp), min(fp$n_sig), max(fp$n_sig))

p <- ggplot(fp, aes(size_bin, FP_proportion, colour = arm, group = arm)) +
  geom_line(alpha = 0.5, linewidth = 0.4, position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(ymin = pmax(0, FP_proportion - FP_proportion_SE),
                     ymax = pmin(1, FP_proportion + FP_proportion_SE)),
                  position = position_dodge(width = 0.4), size = 0.3) +
  facet_wrap(~tag, nrow = 1) +
  scale_colour_manual(values = ARM_COLOURS, name = "method") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "stage-1 cluster size (n markers)", y = "FP proportion among significant clusters",
      title = "False-positive proportion falls with cluster size",
      subtitle = "mean +/- SE across the 10 environments (the replicate axis); pooled counts (sum TP/FP first, then divide) within each environment") +
  theme_bw(11) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8), legend.position = "top")

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "fig_fp_by_size.pdf")
ggsave(OUT, p, width = 10, height = 5, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 10, height = 5, dpi = 200)
say("\n[2] wrote %s (+ .png)\n", OUT)
