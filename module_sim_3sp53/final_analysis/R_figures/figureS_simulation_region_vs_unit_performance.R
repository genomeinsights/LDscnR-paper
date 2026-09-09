## final_analysis/R_figures/figureS_simulation_region_vs_unit_performance.R
##
## Supplementary figure, NOT named in the instructions' output list: absolute
## precision/recall for the three Stage-1-unit-scored restricted methods
## (emmax_simes, emmax_consensus, lfmm_simes) paired against their new
## Stage-2-region-scored counterparts (05_score_truth.R/06_summarise.R,
## 2026-09-10 -- "score at the Stage-2 region level too", PK), all 7 cells,
## both tags. Same visual style as figureS_simulation_absolute_performance.R.
suppressMessages({library(data.table); library(ggplot2)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))
say("=== figureS_simulation_region_vs_unit_performance ===\n\n")

## simulation_performance.tsv/_lfmm.tsv are now the PRIMARY region-level
## tables (2026-09-10 update) -- unit-level diagnostics moved to
## simulation_performance_diagnostic.tsv.
pd <- fread("results/simulation_performance_diagnostic.tsv")
pr <- fread("results/simulation_performance_region.tsv")
pf <- rbindlist(list(
  pd[method %in% c("emmax_simes", "emmax_consensus", "lfmm_simes")],
  pr[method %in% c("emmax_simes_region", "emmax_consensus_region", "lfmm_simes_region")]
), use.names = TRUE, fill = TRUE)

CELL_LEVELS <- c("V0.5_c1", "V1_c1", "V2_c1", "V0.5_c1.5", "V1_c1.5", "V2_c1.5", "V0.5_c2")
METHOD_LEVELS <- c("emmax_simes", "emmax_simes_region", "emmax_consensus", "emmax_consensus_region",
                   "lfmm_simes", "lfmm_simes_region")
METHOD_LABELS <- c("emmax_simes" = "EMMAX Simes (unit)", "emmax_simes_region" = "EMMAX Simes (region)",
                   "emmax_consensus" = "EMMAX consensus (unit)", "emmax_consensus_region" = "EMMAX consensus (region)",
                   "lfmm_simes" = "LFMM Simes (unit)", "lfmm_simes_region" = "LFMM Simes (region)")
pf[, cell := factor(cell, levels = CELL_LEVELS)]
pf[, tag := factor(tag, levels = c("nobgs", "bgs"), labels = c("no BGS", "BGS"))]
pf[, method := factor(method, levels = METHOD_LEVELS, labels = METHOD_LABELS[METHOD_LEVELS])]

long <- rbindlist(list(
  pf[, .(tag, cell, method, metric = "Precision", value = precision, lo = precision_ci_lo, hi = precision_ci_hi)],
  pf[, .(tag, cell, method, metric = "Recall", value = recall, lo = recall_ci_lo, hi = recall_ci_hi)]
))
long[, metric := factor(metric, levels = c("Precision", "Recall"))]
long[, engine_stat := sub(" \\((unit|region)\\)$", "", method)]
long[, granularity := ifelse(grepl("\\(region\\)$", method), "region", "unit")]

METHOD_COLOURS <- c("EMMAX Simes" = "#7B1FA2", "EMMAX consensus" = "#1565C0", "LFMM Simes" = "#2E7D32")

p <- ggplot(long, aes(cell, value, colour = engine_stat, shape = granularity, alpha = granularity)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.6), size = 0.3, fatten = 2) +
  facet_grid(metric ~ tag) +
  scale_colour_manual(values = METHOD_COLOURS, name = NULL) +
  scale_shape_manual(values = c(unit = 16, region = 17), name = "Hypothesis") +
  scale_alpha_manual(values = c(unit = 0.55, region = 1), guide = "none") +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL, y = "Pooled estimate (95% map-cluster bootstrap CI)",
       title = "Absolute precision and recall: Stage-1 unit vs Stage-2 region scoring",
       subtitle = "Filled triangle = one Stage-2 assembled region is one hypothesis; faded circle = one Stage-1 unit is one hypothesis (same method, same combo, different scoring granularity)") +
  theme_bw(11) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top",
        legend.text = element_text(size = 8), plot.subtitle = element_text(size = 8, colour = "grey30"))

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "figureS_simulation_region_vs_unit_performance.pdf")
ggsave(OUT, p, width = 9, height = 6.5, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 9, height = 6.5, dpi = 200)
say("wrote %s (+ .png)\n", OUT)
