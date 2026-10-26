## final_analysis/R_figures/figure_simulation_performance.R
##
## The single main simulation figure (CLAUDE_REANALYSIS_INSTRUCTIONS.md,
## "Final outputs" item 5): two aligned panels showing the PAIRED change in
## STAGE-2 REGION precision and recall from the unrestricted marker-wise
## comparator (emmax_snp_region -- every phenotype-blind Stage-1 cluster,
## including singletons, discovered when it contains >=1 BH-significant
## marker, then assembled the SAME way as Simes/consensus, so it is a
## comparable reported-call unit, not "every significant SNP") for
## Stage-1 Simes and consensus, across all seven cells, with BGS treatment
## distinguishable but not visually dominant -- method is colour (the
## primary comparison), tag is small facet columns, not a competing colour
## dimension. Shows ALL SEVEN designed cells, not a high-dispersal subset,
## per instructions: "If some regimes fail, that is part of the method's
## operating range and belongs in the result."
##
## Reads results/simulation_method_contrasts.tsv (R/06_summarise.R's paired
## bootstrap contrasts vs emmax_snp_region -- already the primary map-cluster
## CI, not the across-environment SE the instructions reject).
## 2026-10-12 (PK review): a pooled primary (c = 1) column added on the left --
## the estimate quoted in the text -- readable cell labels, the shared
## Simes-blue / consensus-orange palette, and no in-figure title (it was
## clipped; the caption carries it).
suppressMessages({library(data.table); library(ggplot2)})
source(file.path(file.path(Sys.getenv("LDSCNR_PAPER_ROOT", "."), "module_sim_3sp53/final_analysis"), "R", "00_config.R"))
source(file.path(file.path(Sys.getenv("LDSCNR_PAPER_ROOT", "."), "module_sim_3sp53/final_analysis"), "R_figures", "palette.R"))
say("=== figure_simulation_performance ===\n\n")

POOLED_X <- "All primary runs\n(high gene flow)"
cn <- fread("results/simulation_method_contrasts.tsv")
cn <- cn[method %in% c("emmax_simes_region", "emmax_consensus_region")]
cn[, panel := fifelse(tag == "nobgs", "No BGS", "BGS")]
cn[, x := CELL_LABELS[cell]]
pc <- fread("results/simulation_method_contrasts_c1.tsv")
pc <- pc[method_b == "emmax_snp_region" & method_a %in% c("emmax_simes_region", "emmax_consensus_region")]
pc[, `:=`(method = method_a, panel = "Pooled", x = POOLED_X)]
keep <- c("panel", "x", "method", "diff_precision", "diff_precision_ci_lo", "diff_precision_ci_hi",
          "diff_recall", "diff_recall_ci_lo", "diff_recall_ci_hi")
cn <- rbind(pc[, ..keep], cn[, ..keep])
cn[, panel := factor(panel, levels = c("Pooled", "No BGS", "BGS"))]
cn[, x := factor(x, levels = c(POOLED_X, CELL_LABELS[CELL_ORDER]))]
cn[, method := factor(method_approach(method), levels = c("Simes", "Consensus"))]

long <- rbindlist(list(
  cn[, .(panel, x, method, metric = "Precision", diff = diff_precision, lo = diff_precision_ci_lo, hi = diff_precision_ci_hi)],
  cn[, .(panel, x, method, metric = "Recall", diff = diff_recall, lo = diff_recall_ci_lo, hi = diff_recall_ci_hi)]
))
long[, metric := factor(metric, levels = c("Precision", "Recall"))]

p <- ggplot(long, aes(x, diff, colour = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.5), size = 0.35, fatten = 2.2) +
  facet_grid(metric ~ panel, scales = "free", space = "free_x") +
  scale_colour_manual(values = APPROACH_COLOURS[c("Simes", "Consensus")], name = NULL) +
  labs(x = NULL, y = "Change relative to marker-wise EMMAX (95% CI)") +
  theme_bw(11) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "top")

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "figure_simulation_performance.pdf")
ggsave(OUT, p, width = 9, height = 6.5, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 9, height = 6.5, dpi = 200)
say("wrote %s (+ .png)\n", OUT)
