## final_analysis/R_figures/figureS_simulation_manhattan_tpfp_stage2.R
##
## Third variant of the illustrative TP/FP Manhattan: same example genome,
## same colour scheme, but scored at the STAGE-2 ASSEMBLED-REGION level --
## a whole region is ONE hypothesis, so every marker inside a region shares
## its status (TP if ANY member marker links to a detectable QTN, else
## FP). manhattan_example_data.R's score_stage2_tpfp().
##
## Motivation (PK): scoring at the Stage-1-UNIT level
## (figureS_simulation_manhattan_tpfp.R) can show both TP and FP markers
## inside what is visually a single Stage-2 peak, when Stage-2 has merged
## a truth-linked unit with an adjacent non-linked one -- an artefact of
## Stage-1 granularity that makes the region look like it "creates FPs
## that shouldn't be there." Scoring at the region level removes that
## artefact for this illustration.
##
## Still illustrative only -- R/06_summarise.R's actual pooled
## precision/recall remain hypothesis-level on Stage-1 units (or plain
## markers for the unrestricted arm), NOT Stage-2 regions, per
## CLAUDE_REANALYSIS_INSTRUCTIONS.md's explicit rejection of the
## dedup-neutral region-level metric for the primary simulation score.
suppressMessages({library(data.table); library(ggplot2); library(LDscnR)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))
source(file.path(PATHS$module, "R_figures", "manhattan_example_data.R"))
say("=== figureS_simulation_manhattan_tpfp_stage2 ===\n\n")

g <- build_example_genome(compute_stage2 = TRUE)
dt <- g$dt

status_em <- list(); status_lf <- list()
for (r in names(g$rep_data)) {
  rd <- g$rep_data[[r]]
  rn <- as.integer(r)
  tpfp <- score_stage2_tpfp(rd)
  if (nrow(tpfp$em)) status_em[[r]] <- tpfp$em[, .(marker, rep = rn, status)]
  if (nrow(tpfp$lf)) status_lf[[r]] <- tpfp$lf[, .(marker, rep = rn, status)]
}
status_em <- if (length(status_em)) rbindlist(status_em) else data.table(marker = character(), rep = integer(), status = character())
status_lf <- if (length(status_lf)) rbindlist(status_lf) else data.table(marker = character(), rep = integer(), status = character())

n_tp <- sum(status_em$status == "TP") + sum(status_lf$status == "TP")
n_fp <- sum(status_em$status == "FP") + sum(status_lf$status == "FP")
say("EMMAX (Stage-2 regions): %d TP, %d FP region-member markers across %d reps\n",
    sum(status_em$status == "TP"), sum(status_em$status == "FP"), length(g$reps))
say("LFMM  (Stage-2 regions): %d TP, %d FP region-member markers across %d reps\n",
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
  scale_colour_manual(values = status_pal, na.value = "grey78", name = "Stage-2 region\n(whole-region status)",
                       guide = guide_legend(override.aes = list(size = 3, shape = 16))) +
  scale_x_continuous(breaks = g$chr_mid$mid, labels = g$chr_mid$global_chr) +
  labs(x = "Chromosome (odd = real, grey/even = near-neutral; 2 per rep x 10 reps)", y = expression(-log[10](q)),
       title = sprintf("%s / %s / env %d -- Stage-2 assembled regions, scored as ONE hypothesis each", g$tag, g$cell, g$envn),
       subtitle = sprintf("+ = QTN; every marker in a region shares that region's status, TP (%d) vs FP (%d) -- compare to figureS_simulation_manhattan_tpfp.R's unit-level mixed colouring within one peak; grey78 = no called region", n_tp, n_fp)) +
  theme_bw(10) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        axis.text.x = element_text(size = 7))

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "figureS_simulation_manhattan_tpfp_stage2.pdf")
ggsave(OUT, p, width = 14, height = 7, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 14, height = 7, dpi = 220)
say("wrote %s (+ .png)\n", OUT)
