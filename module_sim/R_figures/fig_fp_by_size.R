## module_sim/R_figures/fig_fp_by_size.R
##
## Visualize R/05_pool.R's fp_by_size: pooled (sum first, then divide) FP
## proportion among significant stage-1 clusters/units, as a function of
## cluster size, faceted by arm, coloured by tag (bgs/nobgs). PK: "what is
## the proportion of FP as a function of cluster size."
suppressMessages({library(data.table); library(ggplot2)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim"), "R", "00_config.R"))
say("=== fig_fp_by_size ===\n\n")

POOL_PATH <- file.path(PATHS$out, "05_pool", "pooled_pr.rds")
if (!file.exists(POOL_PATH)) stop("R/05_pool.R has not produced: ", POOL_PATH)
fp <- readRDS(POOL_PATH)$fp_by_size

SIZE_LABELS <- c("2", "3", "4-5", "6-10", "11-20", "21-50", "50+")
fp[, size_bin := factor(size_bin, levels = SIZE_LABELS)]
fp[, arm := factor(arm, levels = c("emmax_consensus", "emmax_simes", "lfmm_simes", "emmax_snp", "lfmm_snp"))]
fp[, tag := factor(tag, levels = c("nobgs", "bgs"))]

say("[1] %d (tag,arm,size_bin) points, n_sig range %d-%d\n", nrow(fp), min(fp$n_sig), max(fp$n_sig))

p <- ggplot(fp, aes(size_bin, FP_proportion, colour = tag, group = tag)) +
  geom_line(alpha = 0.5, linewidth = 0.4) +
  geom_point(aes(size = n_sig)) +
  facet_wrap(~arm, nrow = 1) +
  scale_colour_manual(values = c(nobgs = "#1565C0", bgs = "#C0392B")) +
  scale_size_continuous(name = "n significant\nclusters", range = c(1.5, 5)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "stage-1 cluster size (n markers)", y = "FP proportion among significant clusters",
      title = "False-positive proportion falls with cluster size",
      subtitle = "pooled across all reps x envs (sum TP/FP first, then divide); point size = number of significant clusters in that bin") +
  theme_bw(11) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8), legend.position = "top")

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "fig_fp_by_size.pdf")
ggsave(OUT, p, width = 13, height = 4.5, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 13, height = 4.5, dpi = 200)
say("\n[2] wrote %s (+ .png)\n", OUT)
