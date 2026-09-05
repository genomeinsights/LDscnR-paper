## module_sim_bgs2/R_audit/fig_fp_by_size_fixed.R
##
## Audit response (audit_simulations.tex, item 3): corrected FP-by-cluster-size
## figure, NOBGS ONLY, built from fp_by_size_fixed.R's matched-denominator
## recomputation (module_sim_bgs2/results/fp_by_size_fixed_nobgs.rds). Two
## panels: (a) pooled across the 7 cells, one line per method -- the
## corrected version of the original archived figure; (b) the SAME pooled
## quantity broken out by cell, to show whether the cluster-size trend is
## stable within cells or partly an artifact of which cells contribute how
## many clusters at each size (the audit's second point re this figure).
suppressMessages({library(data.table); library(ggplot2)})

IN <- path.expand("~/gitlab/LDscnR-paper/module_sim_bgs2/results/fp_by_size_fixed_nobgs.rds")
d <- readRDS(IN)
fp <- d$fp_by_size
fp_cell <- d$fp_by_size_cell

SIZE_LABELS <- c("2", "3", "4-5", "6-10", "11-20", "21-50", "50+")
ARM_LEVELS <- c("emmax_consensus", "emmax_simes", "lfmm_simes", "emmax_snp", "lfmm_snp")
ARM_COLOURS <- c(emmax_consensus = "#1565C0", emmax_simes = "#26A69A",
                 lfmm_simes = "#7B1FA2", emmax_snp = "#F9A825", lfmm_snp = "#C0392B")
CELL_LEVELS <- c("V0.5_c1", "V0.5_c1.5", "V0.5_c2", "V1_c1", "V1_c1.5", "V2_c1", "V2_c1.5")

fp[, size_bin := factor(size_bin, levels = SIZE_LABELS)]
fp[, arm := factor(arm, levels = ARM_LEVELS)]
fp_cell[, size_bin := factor(size_bin, levels = SIZE_LABELS)]
fp_cell[, arm := factor(arm, levels = ARM_LEVELS)]
fp_cell[, cell := factor(cell, levels = CELL_LEVELS)]

cat(sprintf("[1] pooled: %d (arm,size_bin) points, n_sig range %d-%d\n", nrow(fp), min(fp$n_sig), max(fp$n_sig)))
cat(sprintf("[2] by-cell: %d (cell,arm,size_bin) points, n_sig range %d-%d\n", nrow(fp_cell), min(fp_cell$n_sig), max(fp_cell$n_sig)))

FIG_DIR <- path.expand("~/gitlab/LDscnR-paper/module_sim_bgs2/figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

p_pooled <- ggplot(fp, aes(size_bin, FP_proportion, colour = arm, group = arm)) +
  geom_line(alpha = 0.5, linewidth = 0.4, position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(ymin = pmax(0, FP_proportion - FP_proportion_SE),
                     ymax = pmin(1, FP_proportion + FP_proportion_SE)),
                  position = position_dodge(width = 0.4), size = 0.3) +
  scale_colour_manual(values = ARM_COLOURS, name = "method") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "stage-1 cluster size (n markers)", y = "FP proportion among significant clusters",
      title = "False-positive proportion falls with cluster size (nobgs, pooled over cells)",
      subtitle = "matched denominator TP/(TP+FP) at point estimate and SE; mean +/- SE across 10 environments") +
  theme_bw(11) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8), legend.position = "top")

OUT1 <- file.path(FIG_DIR, "fig_fp_by_size_nobgs_fixed.pdf")
ggsave(OUT1, p_pooled, width = 8, height = 5.5, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT1), p_pooled, width = 8, height = 5.5, dpi = 200)
cat(sprintf("[3] wrote %s (+ .png)\n", OUT1))

p_cell <- ggplot(fp_cell, aes(size_bin, FP_proportion, colour = arm, group = arm)) +
  geom_line(alpha = 0.5, linewidth = 0.35, position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(ymin = pmax(0, FP_proportion - FP_proportion_SE),
                     ymax = pmin(1, FP_proportion + FP_proportion_SE)),
                  position = position_dodge(width = 0.4), size = 0.18, fatten = 1.3) +
  facet_wrap(~cell, nrow = 2) +
  scale_colour_manual(values = ARM_COLOURS, name = "method") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "stage-1 cluster size (n markers)", y = "FP proportion among significant clusters",
      title = "Same quantity, stratified by cell",
      subtitle = "checks whether the pooled trend reflects a within-cell size effect or shifting cell composition") +
  theme_bw(10) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 6), legend.position = "top")

OUT2 <- file.path(FIG_DIR, "fig_fp_by_size_nobgs_bycell.pdf")
ggsave(OUT2, p_cell, width = 11, height = 7, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT2), p_cell, width = 11, height = 7, dpi = 200)
cat(sprintf("[4] wrote %s (+ .png)\n", OUT2))
