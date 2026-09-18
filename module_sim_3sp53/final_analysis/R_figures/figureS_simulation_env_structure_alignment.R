## final_analysis/R_figures/figureS_simulation_env_structure_alignment.R
##
## PK's environment-genetic-structure alignment hypothesis (2026-09-18, see
## env-structure-alignment-hypothesis.md memory): x = population-level R^2 of
## env ~ top 5 GRM eigenvectors (per (tag,cell,rep,env), NOT pooled across
## envs). Panel A: known FP proportion among reported regions. Panel B:
## spatial-null-to-observed region ratio, log2 scale, hline at 1. Colour =
## method (consensus/Simes); facet_grid(V ~ c); BGS treatment shown by
## shape, secondary to colour; point size ~ number of reported regions.
##
## Two outputs, matching this module's manuscript/exploratory convention:
## figureS_simulation_env_structure_alignment.pdf (no title, no in-plot
## caption -- the caption text belongs in the manuscript's own LaTeX
## \caption{}, see the companion _caption.txt this script also writes) and
## the _labelled.pdf twin (title + in-plot caption, for internal review).
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))
say("=== figureS_simulation_env_structure_alignment ===\n\n")

gp <- fread("results/simulation_env_structure_alignment.tsv")
gp[, V := sub("_c.*$", "", cell)]
gp[, c_val := sub("^V[0-9.]+_", "", cell)]
V_LEVELS <- c("V0.5 (strong selection)", "V1 (medium selection)", "V2 (weak selection)")
C_LEVELS <- c("c1 (high gene flow)", "c1.5 (medium gene flow)", "c2 (low gene flow)")
gp[, V := factor(V, levels = c("V0.5", "V1", "V2"), labels = V_LEVELS)]
gp[, c_val := factor(c_val, levels = c("c1", "c1.5", "c2"), labels = C_LEVELS)]
gp[, method_label := ifelse(method == "emmax_consensus", "Consensus", "Simes")]
gp[, tag_label := ifelse(tag == "bgs", "BGS", "No BGS")]

METHOD_COLOURS <- c(Consensus = "#1565C0", Simes = "#7B1FA2")

## "not simulated" labels for the two (V,c) combinations this design doesn't
## include (V1_c2, V2_c2 -- this module is 7 of 9 possible (V,c) cells, see
## nemo-3sp53-canonical-production memory).
MISSING_PANELS <- data.table(V = factor(c("V1 (medium selection)", "V2 (weak selection)"), levels = V_LEVELS),
                             c_val = factor("c2 (low gene flow)", levels = C_LEVELS))
x_mid <- mean(range(gp$r2_axes))

## =============================================================================
## [!] Two corrections from PK's third review of this figure:
## 1. Confidence ribbons removed entirely (both panels). They were built
## from ordinary model/OLS standard errors, which treat every row (an
## environmental continuation, or a BGS/no-BGS pair) as an independent
## observation -- directly contradicting this figure's own caption, which
## says each cell has only 5 independent map/burn-in histories. A ribbon
## built from the actual map-cluster bootstrap would be the correct fix,
## but reusing that machinery here (per-x-grid-point bootstrap prediction
## intervals) is a bigger undertaking than this figure needs; omitting the
## ribbon and showing only the fitted line is honest and sufficient, per
## PK's own suggested resolution.
## 2. Panel A's fitted lines now come from ONE glm per cell,
## cbind(n_FP, n_regions - n_FP) ~ r2_axes + method_label -- a single
## shared alignment slope with a per-method intercept offset, predicted
## separately for each method -- matching EXACTLY the within-cell model
## saved in simulation_env_structure_alignment_models.tsv (r2_axes +
## method_f). The previous version fit a fully separate slope+intercept
## per (cell, method), which does not correspond to any reported model.
## =============================================================================
.fp_curve_cell <- function(d) {
  if (nrow(d) < 10 || uniqueN(d$method_label) < 2) return(NULL)
  fit <- tryCatch(stats::glm(cbind(n_FP, n_regions - n_FP) ~ r2_axes + method_label, data = d, family = stats::binomial()),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  xg <- seq(min(d$r2_axes), max(d$r2_axes), length.out = 40)
  newdat <- CJ(r2_axes = xg, method_label = unique(d$method_label))
  newdat[, fit := predict(fit, newdata = newdat, type = "response")]
  newdat
}
curves_a <- gp[, { cv <- .fp_curve_cell(.SD); if (is.null(cv)) NULL else cv }, by = .(cell, V, c_val)]

## plain ASCII only in axis text -- an en-dash here previously rendered as
## "..." in the PNG device on the mini (same class of issue PK flagged for
## the R-squared symbol in an earlier draft of this label).
X_LAB <- "Environment-structure alignment\n(environmental variation explained by five main axes of relatedness)"

pA <- ggplot(gp, aes(r2_axes, fp_prop, colour = method_label)) +
  geom_point(aes(shape = tag_label, size = n_regions), alpha = 0.3) +
  geom_line(data = curves_a, aes(x = r2_axes, y = fit, colour = method_label), linewidth = 0.8) +
  geom_text(data = MISSING_PANELS, aes(x = x_mid, y = 0.5, label = "not simulated"), inherit.aes = FALSE,
           colour = "grey50", fontface = "italic", size = 3.2) +
  facet_grid(V ~ c_val) +
  scale_colour_manual(values = METHOD_COLOURS, name = NULL) +
  scale_shape_manual(values = c(`No BGS` = 16, BGS = 1), name = NULL) +
  scale_size_area(max_size = 3, name = "Reported\nregions") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "False-positive proportion among reported regions") +
  theme_bw(10) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(), legend.position = "top")

## Panel B: spatial-null / observed ratio vs. alignment, log2 scale, hline at
## 1, intuitive ratio ticks. Same shared-slope-per-cell correction as panel
## A: ONE lm per cell, log2_ratio ~ r2_axes + method_label, predicted
## separately per method -- matching the reported within-cell model
## (r2_axes + method_f) exactly, not a per-method independent OLS fit.
gp_b <- gp[ratio_null_obs > 0]
RATIO_BREAKS <- c(0.0625, 0.25, 1, 4, 16)
.ratio_curve_cell <- function(d) {
  if (nrow(d) < 10 || uniqueN(d$method_label) < 2) return(NULL)
  fit <- tryCatch(stats::lm(log2(ratio_null_obs) ~ r2_axes + method_label, data = d), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  xg <- seq(min(d$r2_axes), max(d$r2_axes), length.out = 40)
  newdat <- CJ(r2_axes = xg, method_label = unique(d$method_label))
  newdat[, fit := 2^predict(fit, newdata = newdat)]
  newdat
}
curves_b <- gp_b[, { cv <- .ratio_curve_cell(.SD); if (is.null(cv)) NULL else cv }, by = .(cell, V, c_val)]

pB <- ggplot(gp_b, aes(r2_axes, ratio_null_obs, colour = method_label, shape = tag_label)) +
  geom_hline(yintercept = 1, linetype = "22", colour = "grey40") +
  geom_point(aes(size = n_regions), alpha = 0.3) +
  geom_line(data = curves_b, aes(x = r2_axes, y = fit, colour = method_label, shape = NULL), linewidth = 0.8) +
  geom_text(data = MISSING_PANELS, aes(x = x_mid, y = 2.8, label = "not simulated"), inherit.aes = FALSE,
           colour = "grey50", fontface = "italic", size = 3.2) +
  facet_grid(V ~ c_val) +
  scale_y_continuous(trans = "log2", breaks = RATIO_BREAKS, labels = as.character(RATIO_BREAKS)) +
  scale_colour_manual(values = METHOD_COLOURS, name = NULL) +
  scale_shape_manual(values = c(`No BGS` = 16, BGS = 1), name = NULL) +
  scale_size_area(max_size = 3, guide = "none") +
  labs(x = X_LAB, y = "Spatial null : observed region ratio") +
  theme_bw(10) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")

n_zero_null <- sum(gp$ratio_null_obs == 0, na.rm = TRUE)
CAPTION <- paste0("Each within-cell analysis contains only 5 independent map/burn-in histories (reps); the many\n",
                  "points shown per cell are environmental continuations and paired BGS/no-BGS analyses of those\n",
                  "same 5 histories, not independent replicates. Fitted lines have no confidence ribbon: ordinary\n",
                  "model standard errors would treat those non-independent points as independent, which this\n",
                  "design does not support. In both panels, the two methods' lines within a cell share one\n",
                  "alignment slope (only their intercept differs), matching the reported within-cell model exactly.\n",
                  sprintf("%d observation%s with zero expected null regions cannot be shown on the logarithmic\n", n_zero_null, if (n_zero_null == 1) "" else "s"),
                  "lower-panel axis and are omitted from panel B only.")
CAPTION_THEME <- theme(plot.caption = element_text(hjust = 0, size = 8.5, lineheight = 1.2, margin = margin(t = 8)),
                       plot.margin = margin(t = 5.5, r = 5.5, b = 12, l = 5.5))

## manuscript version: no title, no in-plot caption -- that text belongs in
## the manuscript's own LaTeX \caption{}, written out separately below.
p <- (pA / pB)
p_labelled <- (pA / pB) +
  plot_annotation(title = "Environment-genetic-structure alignment vs. known FP proportion and spatial-null calibration",
                  caption = CAPTION, theme = CAPTION_THEME)

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "figureS_simulation_env_structure_alignment.pdf")
ggsave(OUT, p, width = 9, height = 10.2, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 9, height = 10.2, dpi = 200)
OUT_L <- sub("\\.pdf$", "_labelled.pdf", OUT)
ggsave(OUT_L, p_labelled, width = 9, height = 11.4, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT_L), p_labelled, width = 9, height = 11.4, dpi = 200)

CAPTION_TXT <- sub("\\.pdf$", "_caption.txt", OUT)
writeLines(gsub("\n", " ", CAPTION), CAPTION_TXT)
say("wrote %s (+ .png), %s (+ .png), and %s (LaTeX-ready caption text)\n", OUT, OUT_L, CAPTION_TXT)
