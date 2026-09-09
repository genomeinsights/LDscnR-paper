## final_analysis/R_figures/figureS_simulation_absolute_performance.R
##
## Supplementary figure (CLAUDE_REANALYSIS_INSTRUCTIONS.md, "Final outputs"
## item 6): absolute precision and recall for every primary method and
## cell -- the four EMMAX-arm methods, all seven cells, both tags.
suppressMessages({library(data.table); library(ggplot2)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))
say("=== figureS_simulation_absolute_performance ===\n\n")

pf <- fread("results/simulation_performance.tsv")
CELL_LEVELS <- c("V0.5_c1", "V1_c1", "V2_c1", "V0.5_c1.5", "V1_c1.5", "V2_c1.5", "V0.5_c2")
METHOD_LEVELS <- c("emmax_snp", "emmax_snp_nonsingleton", "emmax_simes", "emmax_consensus")
METHOD_LABELS <- c("emmax_snp" = "Unrestricted marker-wise", "emmax_snp_nonsingleton" = "Marker-wise, singletons excluded",
                   "emmax_simes" = "Stage-1 Simes", "emmax_consensus" = "Stage-1 consensus")
pf[, cell := factor(cell, levels = CELL_LEVELS)]
pf[, tag := factor(tag, levels = c("nobgs", "bgs"), labels = c("no BGS", "BGS"))]
pf[, method := factor(method, levels = METHOD_LEVELS, labels = METHOD_LABELS[METHOD_LEVELS])]

long <- rbindlist(list(
  pf[, .(tag, cell, method, metric = "Precision", value = precision, lo = precision_ci_lo, hi = precision_ci_hi)],
  pf[, .(tag, cell, method, metric = "Recall", value = recall, lo = recall_ci_lo, hi = recall_ci_hi)]
))
long[, metric := factor(metric, levels = c("Precision", "Recall"))]

METHOD_COLOURS <- c("Unrestricted marker-wise" = "#26A69A", "Marker-wise, singletons excluded" = "#FFCC80",
                    "Stage-1 Simes" = "#7B1FA2", "Stage-1 consensus" = "#1565C0")

p <- ggplot(long, aes(cell, value, colour = method)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.6), size = 0.3, fatten = 2) +
  facet_grid(metric ~ tag) +
  scale_colour_manual(values = METHOD_COLOURS, name = NULL) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL, y = "Pooled estimate (95% map-cluster bootstrap CI)",
       title = "Absolute precision and recall, every primary method and cell") +
  theme_bw(11) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top",
        legend.text = element_text(size = 8))

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "figureS_simulation_absolute_performance.pdf")
ggsave(OUT, p, width = 8, height = 6, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 8, height = 6, dpi = 200)
say("wrote %s (+ .png)\n", OUT)
