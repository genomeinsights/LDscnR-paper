## final_analysis/R_figures/figureS_simulation_stage2_size_truth.R
##
## Analysis 2 figure: FP proportion by cluster size (markers in the region's
## significant Stage-1 units), FP proportion by number of significant Stage-1 units, and the opportunity-
## null comparison (does the observed FP-by-size gradient exceed what pure
## genomic-coverage opportunity predicts). Reads R/14's three output TSVs;
## no association model or Stage-2 assembly is rerun here.
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R_figures", "palette.R"))
say("=== figureS_simulation_stage2_size_truth ===\n\n")

desc <- fread("results/simulation_stage2_size_truth.tsv")
METHOD_LABELS <- c(emmax_snp_region = "EMMAX marker-wise", emmax_simes_region = "EMMAX Simes",
                   emmax_consensus_region = "EMMAX consensus", lfmm_snp_region = "LFMM marker-wise",
                   lfmm_simes_region = "LFMM Simes")
METHOD_COLOURS <- c("EMMAX marker-wise" = COL_MARKER, "EMMAX Simes" = COL_SIMES, "EMMAX consensus" = COL_CONSENSUS,
                    "LFMM marker-wise" = COL_MARKER, "LFMM Simes" = COL_SIMES)
METHOD_LINETYPES <- c("EMMAX marker-wise" = "solid", "EMMAX Simes" = "solid", "EMMAX consensus" = "solid",
                      "LFMM marker-wise" = "22", "LFMM Simes" = "22")
METHOD_SHAPES <- c("EMMAX marker-wise" = 16, "EMMAX Simes" = 16, "EMMAX consensus" = 16,
                   "LFMM marker-wise" = 1, "LFMM Simes" = 1)
desc[, method := factor(METHOD_LABELS[method], levels = METHOD_LABELS)]
desc[, tag := factor(tag, levels = c("nobgs", "bgs"), labels = c("no BGS", "BGS"))]

MARKER_BIN_LABELS <- c("1", "2", "3-5", "6-10", "11-20", "21-50", "51+")
UNIT_BIN_LABELS <- c("1", "2", "3-4", "5-9", "10+")
## print bin ranges with an en dash ("3–5") and tilt the marker-bin labels so they do not overlap
en_dash <- scale_x_discrete(labels = function(x) gsub("-", "\u2013", x, fixed = TRUE))
tilt_x <- theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_markers <- ggplot(desc[size_measure == "n_markers"][, bin := factor(bin, levels = MARKER_BIN_LABELS)],
                    aes(bin, fp_prop, colour = method, linetype = method, shape = method, group = method)) +
  geom_pointrange(aes(ymin = fp_prop_ci_lo, ymax = fp_prop_ci_hi), position = position_dodge(width = 0.5),
                  size = 0.3, fatten = 2) +
  geom_line(position = position_dodge(width = 0.5), alpha = 0.5) +
  facet_wrap(~tag) +
  scale_colour_manual(values = METHOD_COLOURS, name = NULL) +
  scale_linetype_manual(values = METHOD_LINETYPES, name = NULL) +
  scale_shape_manual(values = METHOD_SHAPES, name = NULL) +
  labs(x = "Cluster size (markers)", y = "False-discovery proportion", title = "By cluster size", tag = "A") +
  en_dash + theme_bw(11) + tilt_x + theme(plot.tag = element_text(face = "bold", size = 12), strip.background = element_blank(), panel.grid.minor = element_blank(), legend.position = "top")

p_units <- ggplot(desc[size_measure == "n_units"][, bin := factor(bin, levels = UNIT_BIN_LABELS)],
                  aes(bin, fp_prop, colour = method, linetype = method, shape = method, group = method)) +
  geom_pointrange(aes(ymin = fp_prop_ci_lo, ymax = fp_prop_ci_hi), position = position_dodge(width = 0.5),
                  size = 0.3, fatten = 2) +
  geom_line(position = position_dodge(width = 0.5), alpha = 0.5) +
  facet_wrap(~tag) +
  scale_colour_manual(values = METHOD_COLOURS, name = NULL, guide = "none") +
  scale_linetype_manual(values = METHOD_LINETYPES, guide = "none") +
  scale_shape_manual(values = METHOD_SHAPES, guide = "none") +
  labs(x = "Significant Stage-1 units in region", y = "False-discovery proportion", title = "By number of significant Stage-1 units", tag = "B") +
  en_dash + theme_bw(11) + theme(plot.tag = element_text(face = "bold", size = 12), strip.background = element_blank(), panel.grid.minor = element_blank())

## opportunity-null comparison, if available
opp_file <- "results/simulation_stage2_opportunity_null.tsv"
if (file.exists(opp_file)) {
  opp <- fread(opp_file)
  opp[, marker_bin := factor(marker_bin, levels = MARKER_BIN_LABELS)]
  opp_long <- rbindlist(list(
    opp[, .(marker_bin, series = "Observed", rate = observed_tp_rate)],
    opp[, .(marker_bin, series = "Expected after relocation", rate = opportunity_hit_rate)]))
  p_opp <- ggplot(opp_long, aes(marker_bin, rate, colour = series, group = series)) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
    scale_colour_manual(values = c("Observed" = "#1B9E77", "Expected after relocation" = "grey50"), name = NULL) +
    labs(x = "Cluster size (markers)", y = "Proportion linked to a detectable QTN",
         title = "Observed and expected overlap with detectable QTN", tag = "C") +
    en_dash + theme_bw(11) + theme(plot.tag = element_text(face = "bold", size = 12), strip.background = element_blank(), panel.grid.minor = element_blank(), legend.position = "top")
  ## D: observed / expected ratio, which makes the proportionally faster rise of
  ## the expected rate (C) visible directly
  p_ratio <- ggplot(opp, aes(marker_bin, excess_ratio, group = 1)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    geom_line(linewidth = 0.7, colour = "#1B9E77") + geom_point(size = 1.6, colour = "#1B9E77") +
    scale_y_log10(limits = c(1, NA)) +
    labs(x = "Cluster size (markers)", y = "Observed / expected (log scale)",
         title = "Enrichment over genomic coverage", tag = "D") +
    en_dash + theme_bw(11) + theme(plot.tag = element_text(face = "bold", size = 12), strip.background = element_blank(), panel.grid.minor = element_blank())
} else {
  p_opp <- ggplot() + annotate("text", x = 0, y = 0, label = "Opportunity-null comparison not available", size = 4) +
    theme_void()
  p_ratio <- ggplot() + theme_void()
}

top <- (p_markers | p_units) + plot_layout(guides = "collect") & theme(legend.position = "top")
p <- wrap_elements(top) / (p_opp | p_ratio) + plot_layout(heights = c(1.45, 1))
FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "figureS_simulation_stage2_size_truth.pdf")
ggsave(OUT, p, width = 10, height = 10, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 10, height = 10, dpi = 200)
say("wrote %s (+ .png)\n", OUT)
