## final_analysis/R_figures/figureS_simulation_env_structure_alignment.R
##
## PK's environment-genetic-structure alignment hypothesis (2026-09-18, see
## env-structure-alignment-hypothesis.md memory): x = population-level R^2 of
## env ~ top 5 GRM eigenvectors (per (tag,cell,rep,env), NOT pooled across
## envs). Panel A: known FP proportion. Panel B: spatial-null-to-observed
## region ratio, log2 scale, hline at 1. Raw points light + a linear trend;
## colour = method (consensus/Simes); facet_grid(V ~ c); BGS treatment shown
## by shape, secondary to colour.
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))
say("=== figureS_simulation_env_structure_alignment ===\n\n")

gp <- fread("results/simulation_env_structure_alignment.tsv")
gp[, V := sub("_c.*$", "", cell)]
gp[, c_val := sub("^V[0-9.]+_", "", cell)]
gp[, V := factor(V, levels = c("V0.5", "V1", "V2"), labels = c("V0.5 (strong sel.)", "V1 (medium sel.)", "V2 (weak sel.)"))]
gp[, c_val := factor(c_val, levels = c("c1", "c1.5", "c2"), labels = c("c1 (high gene flow)", "c1.5 (medium gene flow)", "c2 (low gene flow)"))]
gp[, method_label := ifelse(method == "emmax_consensus", "Consensus", "Simes")]
gp[, tag_label := ifelse(tag == "bgs", "BGS", "No BGS")]

METHOD_COLOURS <- c(Consensus = "#1565C0", Simes = "#7B1FA2")

## Panel A: known FP proportion vs. alignment
pA <- ggplot(gp, aes(r2_axes, fp_prop, colour = method_label, shape = tag_label)) +
  geom_point(alpha = 0.25, size = 1.3) +
  geom_smooth(aes(group = method_label), method = "lm", se = TRUE, linewidth = 0.8) +
  facet_grid(V ~ c_val) +
  scale_colour_manual(values = METHOD_COLOURS, name = NULL) +
  scale_shape_manual(values = c(`No BGS` = 16, BGS = 1), name = NULL) +
  labs(x = NULL, y = "Known FP proportion") +
  theme_bw(10) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(), legend.position = "top")

## Panel B: spatial-null / observed ratio vs. alignment, log2 scale, hline at 1
gp_b <- gp[ratio_null_obs > 0]
pB <- ggplot(gp_b, aes(r2_axes, ratio_null_obs, colour = method_label, shape = tag_label)) +
  geom_hline(yintercept = 1, linetype = "22", colour = "grey40") +
  geom_point(alpha = 0.25, size = 1.3) +
  geom_smooth(aes(group = method_label), method = "lm", se = TRUE, linewidth = 0.8) +
  facet_grid(V ~ c_val) +
  scale_y_continuous(trans = "log2") +
  scale_colour_manual(values = METHOD_COLOURS, name = NULL) +
  scale_shape_manual(values = c(`No BGS` = 16, BGS = 1), name = NULL) +
  labs(x = "Environment-structure alignment (pop-level R², top 5 GRM axes)",
       y = "Spatial null : observed region ratio") +
  theme_bw(10) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")

p <- (pA / pB) +
  plot_annotation(title = "Environment-genetic-structure alignment vs. known FP proportion and spatial-null calibration")

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "figureS_simulation_env_structure_alignment.pdf")
ggsave(OUT, p, width = 9, height = 10, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 9, height = 10, dpi = 200)
say("wrote %s (+ .png)\n", OUT)
