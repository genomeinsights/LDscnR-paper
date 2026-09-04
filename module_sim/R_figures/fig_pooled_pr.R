## module_sim/R_figures/fig_pooled_pr.R
##
## Visualize R/05_pool.R's replicate-averaged Precision/Recall (mean +- SE)
## across cells, faceted by engine/arm, coloured by tag (bgs/nobgs).
suppressMessages({library(data.table); library(ggplot2)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim"), "R", "00_config.R"))
say("=== fig_pooled_pr ===\n\n")

POOL_PATH <- file.path(PATHS$out, "05_pool", "pooled_pr.rds")
if (!file.exists(POOL_PATH)) stop("R/05_pool.R has not produced: ", POOL_PATH)
pooled <- readRDS(POOL_PATH)$pooled

## Cell order = selection strength (V0.5 strongest -> V2 weakest, per PK: V
## inversely proportional to selection intensity) then dispersal kernel c.
CELL_ORDER <- c("V0.5_c1", "V0.5_c1.5", "V0.5_c2", "V1_c1", "V1_c1.5", "V2_c1", "V2_c1.5")
pooled[, cell := factor(cell, levels = CELL_ORDER)]
pooled[, arm  := factor(arm, levels = c("emmax_consensus", "emmax_simes", "lfmm_consensus", "lfmm_simes"))]
pooled[, tag  := factor(tag, levels = c("nobgs", "bgs"))]

## long format: one row per (tag, cell, arm, metric), metric in {Precision, Recall},
## with its own SE -- lets both share one facet_grid(metric ~ arm) rather than two
## separately-assembled plots.
long <- rbindlist(list(
  pooled[, .(tag, cell, arm, metric = "Recall",    value = Recall,    SE = Recall_SE)],
  pooled[, .(tag, cell, arm, metric = "Precision", value = Precision, SE = Precision_SE)]
))
long[, metric := factor(metric, levels = c("Recall", "Precision"))]

say("[1] %d (tag,cell,arm,metric) points ; %d Precision NAs (zero-discovery replicates all NA -- see 05_pool.R)\n",
    nrow(long), sum(is.na(long$value) & long$metric == "Precision"))

p <- ggplot(long, aes(cell, value, colour = tag, group = tag)) +
  geom_line(alpha = 0.5, linewidth = 0.4, position = position_dodge(width = 0.3)) +
  geom_pointrange(aes(ymin = pmax(0, value - SE), ymax = pmin(1, value + SE)),
                  position = position_dodge(width = 0.3), size = 0.35, fatten = 2) +
  facet_grid(metric ~ arm) +
  scale_colour_manual(values = c(nobgs = "#1565C0", bgs = "#C0392B")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = NULL,
      title = "Pooled TP/FP scoring: replicate-averaged Precision/Recall by cell (mean +/- SE, N=10 reps)",
      subtitle = "cells ordered by selection strength (V0.5 strongest) then dispersal kernel c") +
  theme_bw(11) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "fig_pooled_pr.pdf")
ggsave(OUT, p, width = 10, height = 5.5, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 10, height = 5.5, dpi = 200)
say("\n[2] wrote %s (+ .png)\n", OUT)
