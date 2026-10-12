## final_analysis/R_figures/figureS_simulation_truth_sensitivity.R
##
## Supplementary figure (CLAUDE_REANALYSIS_INSTRUCTIONS.md, "Final outputs"
## item 8: "only if it adds information beyond its table"). Reads
## results/simulation_truth_sensitivity.tsv (R/09_truth_sensitivity.R):
## grand-pooled precision/recall for the 5 primary Stage-2 region methods
## across the prespecified detectable-QTN additive-variance-share
## threshold grid (0.01-0.20, includes the primary 5% value). The value
## this adds over the table: the method RANKING is visibly stable across
## the whole grid -- no crossing lines -- which is the actual claim this
## sensitivity check is meant to support.
suppressMessages({library(data.table); library(ggplot2)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R_figures", "palette.R"))
say("=== figureS_simulation_truth_sensitivity ===\n\n")

ts <- fread("results/simulation_truth_sensitivity.tsv")
METHOD_LEVELS <- c("emmax_snp_region", "emmax_simes_region", "emmax_consensus_region", "lfmm_snp_region", "lfmm_simes_region")
METHOD_LABELS <- c("emmax_snp_region" = "EMMAX marker-wise", "emmax_simes_region" = "EMMAX Simes",
                   "emmax_consensus_region" = "EMMAX consensus", "lfmm_snp_region" = "LFMM marker-wise",
                   "lfmm_simes_region" = "LFMM Simes")
ts[, method := factor(method, levels = METHOD_LEVELS, labels = METHOD_LABELS[METHOD_LEVELS])]

long <- rbindlist(list(
  ts[, .(method, va_share_threshold, metric = "Precision", value = precision, lo = precision_ci_lo, hi = precision_ci_hi)],
  ts[, .(method, va_share_threshold, metric = "Recall", value = recall, lo = recall_ci_lo, hi = recall_ci_hi)]
))
long[, metric := factor(metric, levels = c("Precision", "Recall"))]

METHOD_COLOURS <- c("EMMAX marker-wise" = COL_MARKER, "EMMAX Simes" = COL_SIMES, "EMMAX consensus" = COL_CONSENSUS,
                    "LFMM marker-wise" = COL_MARKER, "LFMM Simes" = COL_SIMES)
METHOD_LINETYPES <- c("EMMAX marker-wise" = "solid", "EMMAX Simes" = "solid", "EMMAX consensus" = "solid",
                      "LFMM marker-wise" = "22", "LFMM Simes" = "22")
METHOD_SHAPES <- c("EMMAX marker-wise" = 16, "EMMAX Simes" = 16, "EMMAX consensus" = 16,
                   "LFMM marker-wise" = 1, "LFMM Simes" = 1)

p <- ggplot(long, aes(va_share_threshold, value, colour = method, linetype = method, shape = method)) +
  geom_vline(xintercept = VA_SHARE_DETECTABLE, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = method), alpha = 0.12, colour = NA, linetype = 0) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
  facet_wrap(~metric, scales = "free_y") +
  scale_colour_manual(values = METHOD_COLOURS, name = NULL) +
  scale_fill_manual(values = METHOD_COLOURS, guide = "none") +
  scale_linetype_manual(values = METHOD_LINETYPES, name = NULL) +
  scale_shape_manual(values = METHOD_SHAPES, name = NULL) +
  scale_x_continuous(breaks = VA_SHARE_SENSITIVITY_GRID) +
  labs(x = "Detectable-QTN Va-share threshold (dashed = primary, 5%)",
       y = "Grand-pooled estimate (95% map-cluster bootstrap CI)",
       title = "Truth-threshold sensitivity: method ranking is stable across the grid") +
  theme_bw(11) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(), legend.position = "top")

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "figureS_simulation_truth_sensitivity.pdf")
ggsave(OUT, p, width = 9, height = 5, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 9, height = 5, dpi = 200)
say("wrote %s (+ .png)\n", OUT)
