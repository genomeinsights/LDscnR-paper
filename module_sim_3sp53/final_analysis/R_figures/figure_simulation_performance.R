## final_analysis/R_figures/figure_simulation_performance.R
##
## The single main simulation figure (CLAUDE_REANALYSIS_INSTRUCTIONS.md,
## "Final outputs" item 5): two aligned panels showing the PAIRED change in
## precision and recall from unrestricted marker-wise EMMAX (emmax_snp) for
## Stage-1 Simes and consensus, across all seven cells, with BGS treatment
## distinguishable but not visually dominant -- method is colour (the
## primary comparison), tag is small facet columns, not a competing colour
## dimension. Shows ALL SEVEN designed cells, not a high-dispersal subset,
## per instructions: "If some regimes fail, that is part of the method's
## operating range and belongs in the result."
##
## Reads results/simulation_method_contrasts.tsv (R/06_summarise.R's paired
## bootstrap contrasts vs emmax_snp -- already the primary map-cluster CI,
## not the across-environment SE the instructions reject).
suppressMessages({library(data.table); library(ggplot2)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))
say("=== figure_simulation_performance ===\n\n")

cn <- fread("results/simulation_method_contrasts.tsv")
cn <- cn[method %in% c("emmax_simes", "emmax_consensus")]
CELL_LEVELS <- c("V0.5_c1", "V1_c1", "V2_c1", "V0.5_c1.5", "V1_c1.5", "V2_c1.5", "V0.5_c2")
cn[, cell := factor(cell, levels = CELL_LEVELS)]
cn[, tag := factor(tag, levels = c("nobgs", "bgs"), labels = c("no BGS", "BGS"))]
cn[, method := factor(method, levels = c("emmax_simes", "emmax_consensus"),
                      labels = c("Stage-1 Simes", "Stage-1 consensus"))]

long <- rbindlist(list(
  cn[, .(tag, cell, method, metric = "Precision", diff = diff_precision, lo = diff_precision_ci_lo, hi = diff_precision_ci_hi)],
  cn[, .(tag, cell, method, metric = "Recall", diff = diff_recall, lo = diff_recall_ci_lo, hi = diff_recall_ci_hi)]
))
long[, metric := factor(metric, levels = c("Precision", "Recall"))]

p <- ggplot(long, aes(cell, diff, colour = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.5), size = 0.35, fatten = 2.2) +
  facet_grid(metric ~ tag, scales = "free_y") +
  scale_colour_manual(values = c("Stage-1 Simes" = "#7B1FA2", "Stage-1 consensus" = "#1565C0"), name = NULL) +
  labs(x = NULL, y = "Change from unrestricted marker-wise EMMAX (95% map-cluster bootstrap CI)",
       title = "Paired change in precision and recall from phenotype-blind Stage-1 clustering") +
  theme_bw(11) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top")

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "figure_simulation_performance.pdf")
ggsave(OUT, p, width = 8, height = 6, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 8, height = 6, dpi = 200)
say("wrote %s (+ .png)\n", OUT)
