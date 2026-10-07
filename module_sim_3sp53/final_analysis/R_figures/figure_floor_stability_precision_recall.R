## =============================================================================
## final_analysis/R_figures/figure_floor_stability_precision_recall.R
##
## EXPLORATORY (PK, 2026-09-27) -- precision-recall trade-off across the
## cumulative stability tiers (>=1/3 = all canonical regions, >=2/3, =3/3),
## per method, read purely from results/floor_stability_summary.tsv (no
## model rerun). The existing floor=2 canonical reference (results/
## floor_stability_reference_floor2.tsv) is plotted as a separate, distinctly
## shaped point -- changing the canonical floor and imposing a stability
## filter are two different interventions (PK's explicit framing) and must
## not look like one continuous path on this figure.
## =============================================================================
suppressMessages({library(data.table); library(ggplot2)})
MODULE_ROOT <- path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis")

summ <- fread(file.path(MODULE_ROOT, "results", "floor_stability_summary.tsv"))
tier_dt <- summ[call_set %in% c("all_ge1", "ge2", "eq3")]
tier_dt[, call_set := factor(call_set, levels = c("all_ge1", "ge2", "eq3"),
                             labels = c("all canonical (>=1/3)", ">=2/3 settings", "3/3 settings"))]

ref_file <- file.path(MODULE_ROOT, "results", "floor_stability_reference_floor2.tsv")
ref <- if (file.exists(ref_file)) fread(ref_file) else NULL

p <- ggplot(tier_dt, aes(recall, precision, colour = method, group = method)) +
  geom_errorbar(aes(ymin = precision_ci_lo, ymax = precision_ci_hi), width = 0, alpha = 0.5) +
  geom_errorbar(aes(xmin = recall_ci_lo, xmax = recall_ci_hi), width = 0, alpha = 0.5, orientation = "y") +
  geom_path(linewidth = 0.6) +
  geom_point(aes(shape = call_set), size = 3) +
  { if (!is.null(ref)) geom_point(data = ref, aes(recall, precision),
                                 inherit.aes = FALSE, colour = "grey30", shape = 4, size = 3) } +
  theme_bw(11) + theme(strip.background = element_blank()) +
  labs(x = "Recall (unique detectable QTN recovered / total detectable QTN)",
      y = "Precision (TP regions / all regions)",
      shape = "Stability tier (cumulative)", colour = "Method",
      title = "Floor-stability trade-off: precision vs. recall across cumulative tiers",
      subtitle = "Grey X (if present) = existing floor=2 canonical reference -- a separate intervention, not a fourth tier point",
      caption = "600 primary c=1 simulations (V0.5_c1/V1_c1/V2_c1, both BGS treatments, 10 maps x 10 environmental continuations).\nErrorbars are 95% map-cluster-bootstrap intervals (clustered by rep). EXPLORATORY -- not yet in the manuscript.")

ggsave(file.path(MODULE_ROOT, "figures", "figure_floor_stability_precision_recall.pdf"), p, width = 7.5, height = 5.5)
cat("Wrote figures/figure_floor_stability_precision_recall.pdf\n")
