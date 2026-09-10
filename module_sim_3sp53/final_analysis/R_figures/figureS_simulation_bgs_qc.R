## final_analysis/R_figures/figureS_simulation_bgs_qc.R
##
## The one compact BGS validation figure (CLAUDE_REANALYSIS_INSTRUCTIONS.md,
## "One compact BGS validation" + "Final outputs" item 7). Reads
## results/simulation_bgs_validation.tsv (R/08_bgs_validation.R).
##
## Two panels, same x-axis (recombination-rate bin):
##  - LEFT (the validated BGS signature): pre-MAF-filter retained-marker
##    ratio (bgs/nobgs) -- BGS retains ~17% fewer segregating markers
##    overall, MOST at low recombination (linked selection's footprint is
##    largest there) -- the classic background-selection signature.
##  - RIGHT (the trap this confirms): post-MAF-filter mean heterozygosity
##    difference (bgs - nobgs) is POSITIVE, not negative -- among markers
##    that SURVIVE the MAF>0.10 cut, BGS looks more diverse, backwards
##    from panel 1's aggregate effect. MAF filtering removes most of what
##    BGS actually did (drove rare variants to loss) before He ever sees
##    it; the mild positive residual among survivors is consistent with
##    associative overdominance. Do not read panel 2 alone as "BGS
##    validation failed" -- panel 1 is the validated result.
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))
say("=== figureS_simulation_bgs_qc ===\n\n")

bv <- fread("results/simulation_bgs_validation.tsv")
BIN_LEVELS <- c("overall", "zero", "low", "medium", "high")
bv[, bin := factor(bin, levels = BIN_LEVELS)]

pre <- bv[filt == "pre_maf"]
post <- bv[filt == "post_maf"]

p1 <- ggplot(pre[bin != "overall"], aes(bin, ratio_n_markers)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_point(size = 2.2, colour = "#B71C1C") +
  geom_segment(aes(y = 1, yend = ratio_n_markers, x = bin, xend = bin), colour = "#B71C1C", linewidth = 0.6) +
  labs(x = "Recombination-rate bin", y = "Retained markers, bgs / nobgs\n(pre-MAF-filter)",
       title = "Segregating-marker retention", subtitle = "BGS's validated effect: fewer polymorphic sites, most at low recombination") +
  theme_bw(10) + theme(panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 7.5, colour = "grey30"))

p2 <- ggplot(post[bin != "overall"], aes(bin, mean_diff_he)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_pointrange(aes(ymin = diff_ci_lo, ymax = diff_ci_hi), size = 0.4, colour = "#1565C0") +
  labs(x = "Recombination-rate bin", y = "Mean He difference, bgs - nobgs\n(post-MAF-filter survivors)",
       title = "Heterozygosity among MAF-filter survivors", subtitle = "The trap: positive, backwards from panel 1's effect") +
  theme_bw(10) + theme(panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 7.5, colour = "grey30"))

p <- p1 + p2 + plot_annotation(
  title = "BGS quality control: reduced polymorphism (correctly measured) vs the post-filter trap",
  subtitle = sprintf("Grand-pooled (700 map/tag/env pairs, rep-cluster bootstrap CI). Overall retention ratio = %.3f [markers: nobgs=%s, bgs=%s]",
                     pre[bin == "overall", ratio_n_markers], format(round(pre[bin=="overall",n_markers_nobgs]), big.mark=","), format(round(pre[bin=="overall",n_markers_bgs]), big.mark=",")),
  theme = theme(plot.title = element_text(size = 11, face = "bold"), plot.subtitle = element_text(size = 8, colour = "grey30")))

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "figureS_simulation_bgs_qc.pdf")
ggsave(OUT, p, width = 9, height = 4.5, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 9, height = 4.5, dpi = 200)
say("wrote %s (+ .png)\n", OUT)
