## final_analysis/R_figures/figure_simulation_null_truth_calibration.R
##
## Analysis 1 figure: does the null-to-observed Stage-2-region ratio predict
## the known false-discovery proportion? Reads R/13's
## results/simulation_null_truth_calibration.tsv (map/burn-in-level rows) --
## no null replicates are rerun here. Produces a manuscript-ready version
## (no title, per the spec) and an exploratory "_labelled" version with a
## descriptive title/subtitle and per-point annotation of the balanced
## design used (see ADDITIONAL_ANALYSES_AUDIT.md for why the full grid was
## not run). The raw plotted points are also exported as TSV.
suppressMessages({library(data.table); library(ggplot2)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))
say("=== figure_simulation_null_truth_calibration ===\n\n")

ml <- fread("results/simulation_null_truth_calibration.tsv")
METHOD_LABELS <- c(emmax_consensus = "EMMAX consensus", emmax_simes = "EMMAX Simes")
METHOD_COLOURS <- c("EMMAX consensus" = "#1565C0", "EMMAX Simes" = "#7B1FA2")
SCHEME_LABELS <- c(group = "Group permutation", mvn = "Kinship-matched (MVN)", spatial = "Spatial MVN")

ml[, method := factor(METHOD_LABELS[method], levels = METHOD_LABELS)]
ml[, scheme := factor(SCHEME_LABELS[scheme], levels = SCHEME_LABELS)]
ml[, above_one := !is.na(R_null_obs) & R_null_obs > 1]

## export the exact plotted data (including NA/zero-observed rows, per spec)
fwrite(ml, "results/figure_simulation_null_truth_calibration_data.tsv", sep = "\t")

plot_dt <- ml[!is.na(R_null_obs) & !is.na(FDP_truth)]
n_above_one <- sum(plot_dt$above_one)
n_zero_obs <- sum(is.na(ml$R_null_obs))

.base_plot <- function(with_title) {
  g <- ggplot(plot_dt, aes(FDP_truth, R_null_obs, colour = method, shape = above_one)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_point(size = 1.9, alpha = 0.8) +
    facet_wrap(~scheme) +
    scale_colour_manual(values = METHOD_COLOURS, name = NULL) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    coord_cartesian(xlim = c(0, 1)) +
    labs(x = "Known false-discovery proportion (truth)",
         y = "Null-to-observed Stage-2 region ratio (raw, uncapped)") +
    theme_bw(11) +
    theme(strip.background = element_blank(), panel.grid.minor = element_blank(), legend.position = "top")
  if (with_title) {
    g <- g + labs(title = "Does the null-calibrated discovery ratio predict the known false-discovery proportion?",
                  subtitle = sprintf("Balanced design: %d map/burn-in groups per scheme x method; %d points above the 1:1 line (triangles); %d zero-observed groups excluded (NA ratio, not floored)",
                                     nrow(plot_dt) / (nlevels(ml$scheme) * nlevels(ml$method)), n_above_one, n_zero_obs))
  }
  g
}

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

p_manuscript <- .base_plot(with_title = FALSE)
OUT <- file.path(FIG_DIR, "figure_simulation_null_truth_calibration.pdf")
ggsave(OUT, p_manuscript, width = 9, height = 3.6, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p_manuscript, width = 9, height = 3.6, dpi = 200)

p_labelled <- .base_plot(with_title = TRUE)
OUT_L <- file.path(FIG_DIR, "figure_simulation_null_truth_calibration_labelled.pdf")
ggsave(OUT_L, p_labelled, width = 9, height = 4.4, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT_L), p_labelled, width = 9, height = 4.4, dpi = 200)

say("wrote %s (+ .png) and %s (+ .png)\n", OUT, OUT_L)
say("wrote results/figure_simulation_null_truth_calibration_data.tsv\n")
