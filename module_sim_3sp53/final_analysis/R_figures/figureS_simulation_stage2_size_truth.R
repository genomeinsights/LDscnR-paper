## final_analysis/R_figures/figureS_simulation_stage2_size_truth.R
##
## Analysis 2 figure: FP proportion by number of constituent markers, FP
## proportion by number of constituent Stage-1 units, and the opportunity-
## null comparison (does the observed FP-by-size gradient exceed what pure
## genomic-coverage opportunity predicts). Reads R/14's three output TSVs;
## no association model or Stage-2 assembly is rerun here.
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))
say("=== figureS_simulation_stage2_size_truth ===\n\n")

desc <- fread("results/simulation_stage2_size_truth.tsv")
METHOD_LABELS <- c(emmax_snp_region = "EMMAX unrestricted", emmax_simes_region = "EMMAX Simes",
                   emmax_consensus_region = "EMMAX consensus", lfmm_snp_region = "LFMM unrestricted",
                   lfmm_simes_region = "LFMM Simes")
METHOD_COLOURS <- c("EMMAX unrestricted" = "#26A69A", "EMMAX Simes" = "#7B1FA2",
                    "EMMAX consensus" = "#1565C0", "LFMM unrestricted" = "#EF6C00", "LFMM Simes" = "#2E7D32")
desc[, method := factor(METHOD_LABELS[method], levels = METHOD_LABELS)]
desc[, tag := factor(tag, levels = c("nobgs", "bgs"), labels = c("no BGS", "BGS"))]

MARKER_BIN_LABELS <- c("1", "2", "3-5", "6-10", "11-20", "21-50", "51+")
UNIT_BIN_LABELS <- c("1", "2", "3-4", "5-9", "10+")

p_markers <- ggplot(desc[size_measure == "n_markers"][, bin := factor(bin, levels = MARKER_BIN_LABELS)],
                    aes(bin, fp_prop, colour = method, group = method)) +
  geom_pointrange(aes(ymin = fp_prop_ci_lo, ymax = fp_prop_ci_hi), position = position_dodge(width = 0.5),
                  size = 0.3, fatten = 2) +
  geom_line(position = position_dodge(width = 0.5), alpha = 0.5) +
  facet_wrap(~tag) +
  scale_colour_manual(values = METHOD_COLOURS, name = NULL) +
  labs(x = "Constituent markers", y = "False-positive proportion", title = "By number of constituent markers") +
  theme_bw(10) + theme(strip.background = element_blank(), panel.grid.minor = element_blank(), legend.position = "top")

p_units <- ggplot(desc[size_measure == "n_units"][, bin := factor(bin, levels = UNIT_BIN_LABELS)],
                  aes(bin, fp_prop, colour = method, group = method)) +
  geom_pointrange(aes(ymin = fp_prop_ci_lo, ymax = fp_prop_ci_hi), position = position_dodge(width = 0.5),
                  size = 0.3, fatten = 2) +
  geom_line(position = position_dodge(width = 0.5), alpha = 0.5) +
  facet_wrap(~tag) +
  scale_colour_manual(values = METHOD_COLOURS, name = NULL, guide = "none") +
  labs(x = "Constituent Stage-1 units", y = "False-positive proportion", title = "By number of constituent Stage-1 units") +
  theme_bw(10) + theme(strip.background = element_blank(), panel.grid.minor = element_blank())

## opportunity-null comparison, if available
opp_file <- "results/simulation_stage2_opportunity_null.tsv"
if (file.exists(opp_file)) {
  opp <- fread(opp_file)
  opp[, marker_bin := factor(marker_bin, levels = MARKER_BIN_LABELS)]
  opp_long <- rbindlist(list(
    opp[, .(marker_bin, series = "Observed TP rate", rate = observed_tp_rate)],
    opp[, .(marker_bin, series = "Opportunity-null hit rate", rate = opportunity_hit_rate)]))
  p_opp <- ggplot(opp_long, aes(marker_bin, rate, colour = series, group = series)) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
    scale_colour_manual(values = c("Observed TP rate" = "#1B9E77", "Opportunity-null hit rate" = "grey50"), name = NULL) +
    labs(x = "Constituent markers", y = "Rate",
         title = "Observed truth rate vs. chromosome-preserving random-window opportunity null") +
    theme_bw(10) + theme(strip.background = element_blank(), panel.grid.minor = element_blank(), legend.position = "top")
} else {
  p_opp <- ggplot() + annotate("text", x = 0, y = 0, label = "Opportunity-null comparison not available", size = 4) +
    theme_void()
}

p <- (p_markers | p_units) / p_opp + plot_layout(heights = c(1, 1))
FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "figureS_simulation_stage2_size_truth.pdf")
ggsave(OUT, p, width = 10, height = 8.5, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 10, height = 8.5, dpi = 200)
say("wrote %s (+ .png)\n", OUT)
