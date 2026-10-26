## final_analysis/R_figures/figureS_simulation_absolute_performance.R
##
## Supplementary figure (CLAUDE_REANALYSIS_INSTRUCTIONS.md, "Final outputs"
## item 6): absolute STAGE-2 REGION precision and recall for every primary
## method and cell -- the three EMMAX-arm region methods, all seven cells,
## both tags. emmax_snp_region is the unrestricted marker-wise comparator
## made comparable via post hoc Stage-2 assembly of every discovered
## (>=1 BH-significant marker) phenotype-blind cluster, including
## singletons -- not a plain marker-wise count.
suppressMessages({library(data.table); library(ggplot2)})
source(file.path(file.path(Sys.getenv("LDSCNR_PAPER_ROOT", "."), "module_sim_3sp53/final_analysis"), "R", "00_config.R"))
source(file.path(file.path(Sys.getenv("LDSCNR_PAPER_ROOT", "."), "module_sim_3sp53/final_analysis"), "R_figures", "palette.R"))
say("=== figureS_simulation_absolute_performance ===\n\n")

pf <- fread("results/simulation_performance.tsv")
CELL_LEVELS <- c("V0.5_c1", "V1_c1", "V2_c1", "V0.5_c1.5", "V1_c1.5", "V2_c1.5", "V0.5_c2")
METHOD_LEVELS <- c("emmax_snp_region", "emmax_simes_region", "emmax_consensus_region")
METHOD_LABELS <- c("emmax_snp_region" = "Marker-wise", "emmax_simes_region" = "Simes", "emmax_consensus_region" = "Consensus")
pf[, cell := factor(cell, levels = CELL_LEVELS, labels = CELL_LABELS[CELL_LEVELS])]
pf[, tag := factor(tag, levels = c("nobgs", "bgs"), labels = c("no BGS", "BGS"))]
pf[, method := factor(method, levels = METHOD_LEVELS, labels = METHOD_LABELS[METHOD_LEVELS])]

long <- rbindlist(list(
  pf[, .(tag, cell, method, metric = "Precision", value = precision, lo = precision_ci_lo, hi = precision_ci_hi)],
  pf[, .(tag, cell, method, metric = "Recall", value = recall, lo = recall_ci_lo, hi = recall_ci_hi)]
))
long[, metric := factor(metric, levels = c("Precision", "Recall"))]

METHOD_COLOURS <- APPROACH_COLOURS[c("Marker-wise", "Simes", "Consensus")]

p <- ggplot(long, aes(cell, value, colour = method)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.6), size = 0.3, fatten = 2) +
  facet_grid(metric ~ tag) +
  scale_colour_manual(values = METHOD_COLOURS, name = NULL) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL, y = "Pooled estimate (95% CI)") +
  theme_bw(11) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "top",
        legend.text = element_text(size = 8))

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "figureS_simulation_absolute_performance.pdf")
ggsave(OUT, p, width = 8, height = 6, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 8, height = 6, dpi = 200)
say("wrote %s (+ .png)\n", OUT)
