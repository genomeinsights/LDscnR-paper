## final_analysis/R_figures/figureS_simulation_manhattan_example.R
##
## Illustrative supplementary figure: one example (tag, cell, env), all 10
## reps concatenated into a 20-"chromosome" genome (see
## manhattan_example_data.R for why env=3/nobgs/V0.5_c1 was picked), EMMAX
## on top / LFMM below, every marker coloured by which Stage-2 assembled
## outlier region (if any) it physically falls in -- LDscnR's own
## default_cluster_colours(), the same palette/style ld_manhattan() itself
## uses. QTN are a "+", coloured the SAME as every other marker so a QTN's
## own region membership (or lack of one) is directly visible.
##
## Stage 2 is illustrative only here (as everywhere in this module) -- it
## is NOT what R/06_summarise.R's precision/recall are computed from; see
## figureS_simulation_manhattan_tpfp.R for the actually-scored Stage-1
## units, coloured by TP/FP.
suppressMessages({library(data.table); library(ggplot2); library(LDscnR)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))
source(file.path(PATHS$module, "R_figures", "manhattan_example_data.R"))
say("=== figureS_simulation_manhattan_example ===\n\n")

g <- build_example_genome(compute_stage2 = TRUE)
dt <- g$dt

dt[, region := NA_character_]
## joined on (marker, rep), NOT marker alone -- each rep is an independent
## simulated genome, so "Chr1:12345"-style marker names collide across
## reps (same physical-position discretisation, different underlying
## individuals); marker alone would silently cross-contaminate rows.
for (r in names(g$rep_data)) {
  rd <- g$rep_data[[r]]
  rn <- as.integer(r)
  reg_em <- assign_stage2_region(rd$map, rd$test_em$regions, rn)
  reg_lf <- assign_stage2_region(rd$map, rd$test_lf$regions, rn)
  lut <- data.table(marker = rd$map$marker, rep = rn, EMMAX = reg_em, LFMM = reg_lf)
  dt[lut, region := fifelse(engine == "EMMAX", i.EMMAX, i.LFMM), on = c("marker", "rep")]
}
dt[, region := factor(region)]

n_regions <- sum(vapply(g$rep_data, function(rd) nrow(rd$test_em$regions) + nrow(rd$test_lf$regions), numeric(1)))
say("total Stage-2 regions across %d reps: %d\n", length(g$reps), n_regions)

region_levels <- levels(dt$region)
region_pal <- stats::setNames(rep(LDscnR:::default_cluster_colours(), length.out = length(region_levels)), region_levels)

d_bg_dot <- dt[is_qtn == FALSE & is.na(region)]
d_fg_dot <- dt[is_qtn == FALSE & !is.na(region)]
d_bg_qtn <- dt[is_qtn == TRUE & is.na(region)]
d_fg_qtn <- dt[is_qtn == TRUE & !is.na(region)]

ALPHA_LINE <- ALPHA

p <- ggplot(dt) +
  geom_rect(data = g$band_rects, aes(xmin = xmin - 1e5, xmax = xmax + 1e5, ymin = -Inf, ymax = Inf),
            fill = "grey90", alpha = 0.5, inherit.aes = FALSE) +
  geom_hline(yintercept = -log10(ALPHA_LINE), linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_point(data = d_bg_dot, aes(pos_cum, -log10(q)), colour = "grey78", size = 0.9, alpha = 0.85) +
  geom_point(data = d_fg_dot, aes(pos_cum, -log10(q), colour = region), size = 1.6, alpha = 0.95) +
  geom_point(data = d_bg_qtn, aes(pos_cum, -log10(q)), shape = 3, size = 3, stroke = 1, colour = "grey78") +
  geom_point(data = d_fg_qtn, aes(pos_cum, -log10(q), colour = region), shape = 3, size = 3, stroke = 1) +
  facet_grid(engine ~ ., scales = "free_y") +
  scale_colour_manual(values = region_pal, na.value = "grey78", name = "Stage-2 outlier region",
                       guide = guide_legend(override.aes = list(size = 3, shape = 16))) +
  scale_x_continuous(breaks = g$chr_mid$mid, labels = g$chr_mid$global_chr) +
  labs(x = "Chromosome (odd = real, grey/even = near-neutral; 2 per rep x 10 reps)", y = expression(-log[10](q)),
       title = sprintf("%s / %s / env %d -- 10 reps as a 20-chromosome illustrative genome", g$tag, g$cell, g$envn),
       subtitle = "+ = QTN, coloured the same as every other marker (Stage-2 outlier region membership, illustrative only); grey78 = outside any call") +
  theme_bw(10) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        axis.text.x = element_text(size = 7),
        legend.key.size = unit(0.35, "cm"), legend.text = element_text(size = 7))

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "figureS_simulation_manhattan_example.pdf")
ggsave(OUT, p, width = 14, height = 7, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 14, height = 7, dpi = 220)
say("wrote %s (+ .png)\n", OUT)
