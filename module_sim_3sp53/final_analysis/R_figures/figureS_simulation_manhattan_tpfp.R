## final_analysis/R_figures/figureS_simulation_manhattan_tpfp.R
##
## Illustrative supplementary figure: same example genome as
## figureS_simulation_manhattan_example.R (10 reps of nobgs/V0.5_c1/env=3
## concatenated into a 20-"chromosome" genome), EMMAX on top / LFMM below,
## but now coloured by TP/FP status of the SIGNIFICANT STAGE-1 UNIT
## (emmax_simes/lfmm_simes) a marker belongs to -- i.e. exactly the
## hypothesis-level classification R/05_score_truth.R counts towards
## R/06_summarise.R's pooled precision/recall (manhattan_example_data.R's
## score_units_tpfp() reproduces that logic verbatim: same
## flag_true_qtns()/qtn_ld_table()/score_thresholds()/
## LDscnR:::.ld_outlier_units() primitives, same PARAMS). This is
## deliberately NOT the Stage-2 assembled regions from
## figureS_simulation_manhattan_example.R -- Stage 2 is not what is scored
## (CLAUDE_REANALYSIS_INSTRUCTIONS.md: "not needed for the primary
## simulation score").
suppressMessages({library(data.table); library(ggplot2); library(LDscnR)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))
source(file.path(PATHS$module, "R_figures", "manhattan_example_data.R"))
say("=== figureS_simulation_manhattan_tpfp ===\n\n")

g <- build_example_genome(compute_stage2 = FALSE)
dt <- g$dt

## per-rep TP/FP status tables, rbound across all 10 reps, THEN joined onto
## dt on (marker, rep) -- not marker alone, since marker names ("Chr1:pos")
## collide across reps (independent genomes sharing the same physical-
## position discretisation); a bare marker join would cross-contaminate
## status labels between reps.
status_em <- list(); status_lf <- list()
for (r in names(g$rep_data)) {
  rd <- g$rep_data[[r]]
  rn <- as.integer(r)
  tpfp <- score_units_tpfp(rd)
  if (nrow(tpfp$em)) status_em[[r]] <- tpfp$em[, .(marker, rep = rn, status)]
  if (nrow(tpfp$lf)) status_lf[[r]] <- tpfp$lf[, .(marker, rep = rn, status)]
}
status_em <- if (length(status_em)) rbindlist(status_em) else data.table(marker = character(), rep = integer(), status = character())
status_lf <- if (length(status_lf)) rbindlist(status_lf) else data.table(marker = character(), rep = integer(), status = character())

n_tp <- sum(status_em$status == "TP") + sum(status_lf$status == "TP")
n_fp <- sum(status_em$status == "FP") + sum(status_lf$status == "FP")
say("EMMAX: %d TP, %d FP significant units' markers across %d reps\n",
    sum(status_em$status == "TP"), sum(status_em$status == "FP"), length(g$reps))
say("LFMM:  %d TP, %d FP significant units' markers across %d reps\n",
    sum(status_lf$status == "TP"), sum(status_lf$status == "FP"), length(g$reps))

dt[, status := NA_character_]
dt[engine == "EMMAX"][status_em, status := i.status, on = c("marker", "rep")] -> em_rows
dt[engine == "LFMM"][status_lf, status := i.status, on = c("marker", "rep")] -> lf_rows
dt <- rbindlist(list(em_rows, lf_rows))
dt[, engine := factor(engine, levels = c("EMMAX", "LFMM"))]
dt[, status := factor(status, levels = c("TP", "FP"))]

status_pal <- c(TP = "#1B9E77", FP = "#D95F02")

d_bg_dot <- dt[is_qtn == FALSE & is.na(status)]
d_fg_dot <- dt[is_qtn == FALSE & !is.na(status)]
d_bg_qtn <- dt[is_qtn == TRUE & is.na(status)]
d_fg_qtn <- dt[is_qtn == TRUE & !is.na(status)]

ALPHA_LINE <- ALPHA

p <- ggplot(dt) +
  geom_rect(data = g$band_rects, aes(xmin = xmin - 1e5, xmax = xmax + 1e5, ymin = -Inf, ymax = Inf),
            fill = "grey90", alpha = 0.5, inherit.aes = FALSE) +
  geom_hline(yintercept = -log10(ALPHA_LINE), linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_point(data = d_bg_dot, aes(pos_cum, -log10(q)), colour = "grey78", size = 0.9, alpha = 0.85) +
  geom_point(data = d_fg_dot, aes(pos_cum, -log10(q), colour = status), size = 1.8, alpha = 0.95) +
  geom_point(data = d_bg_qtn, aes(pos_cum, -log10(q)), shape = 3, size = 3, stroke = 1, colour = "grey78") +
  geom_point(data = d_fg_qtn, aes(pos_cum, -log10(q), colour = status), shape = 3, size = 3, stroke = 1) +
  facet_grid(engine ~ ., scales = "free_y") +
  scale_colour_manual(values = status_pal, na.value = "grey78", name = "Significant Stage-1 unit",
                       guide = guide_legend(override.aes = list(size = 3, shape = 16))) +
  scale_x_continuous(breaks = g$chr_mid$mid, labels = g$chr_mid$global_chr) +
  labs(x = "Chromosome (odd = real, grey/even = near-neutral; 2 per rep x 10 reps)", y = expression(-log[10](q)),
       title = sprintf("%s / %s / env %d -- 10 reps as a 20-chromosome illustrative genome", g$tag, g$cell, g$envn),
       subtitle = sprintf("+ = QTN; coloured markers = members of a significant emmax_simes/lfmm_simes unit, TP (%d) vs FP (%d), exactly as scored for precision/recall; grey78 = not in a significant unit", n_tp, n_fp)) +
  theme_bw(10) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        axis.text.x = element_text(size = 7))

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "figureS_simulation_manhattan_tpfp.pdf")
ggsave(OUT, p, width = 14, height = 7, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 14, height = 7, dpi = 220)
say("wrote %s (+ .png)\n", OUT)
