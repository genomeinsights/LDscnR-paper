## module_sim_3sp53/R_figures/fig_random_removal_control.R
##
## Visualize R/14_random_removal_control.R: size-based vs matched-
## cardinality-random realised_fdr by floor. PK, 2026-09-09: "if we
## randomly remove the same number of outlier clusters (but independently
## of size), would we also see a decline in FPs, and is the decline in
## fig_structured_null_sizesweep relative to that?" then "run it on
## emmax_simes too. Also add SE's if possible."
##
## Full grid: all 7 cells x both tags (bgs, nobgs) x both arms
## (emmax_consensus, emmax_simes), full 10x10 rep x env each, group null
## only. Faceted by arm x tag (2x2) -- the two comparisons that matter:
## bgs vs nobgs (despite very different ground-truth size-truth strength,
## they behave the same -- rules out a recombination/BGS confound) and
## consensus vs simes (does a genuinely multi-comparison combining rule
## behave differently from consensus_dosage's single lower-noise
## variable -- see README.md's 2026-09-09 mechanism discussion).
##
## SE = sd/sqrt(n) across the per-combo (rep x env) values contributing
## to each mean, n = combos with n_obs_size>0 at that floor (shown by
## point size, same as before) -- NOT the SE of the N_RANDOM=100 draws
## within a combo (already averaged away inside fdr_random_conditional
## before this script ever sees it).
suppressMessages({library(data.table); library(ggplot2)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53"), "R", "00_config.R"))
say("=== fig_random_removal_control ===\n\n")

FILES <- Sys.glob(file.path(PATHS$out, "14_random_removal_control", "rrc_*.rds"))
if (!length(FILES)) stop("R/14_random_removal_control.R has not produced any out/14_random_removal_control/rrc_*.rds")
all <- rbindlist(lapply(FILES, readRDS), fill = TRUE)
## old consensus-only files predate the `arm` column
if (!"arm" %in% names(all)) all[, arm := NA_character_]
all[is.na(arm), arm := "emmax_consensus"]
all[, fdr_size := n_surr_size / max(n_obs_size, 1), by = seq_len(nrow(all))]
say("[1] %d rows, %d cells x %d tags x %d arms\n", nrow(all), uniqueN(all$cell), uniqueN(all$tag), uniqueN(all$arm))

.se <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) return(NA_real_); stats::sd(x) / sqrt(length(x)) }

by_tag_arm <- all[n_obs_size > 0, .(
    fdr_size = mean(fdr_size), fdr_size_se = .se(fdr_size),
    fdr_random = mean(fdr_random_conditional, na.rm = TRUE), fdr_random_se = .se(fdr_random_conditional),
    n_combo = .N
  ), by = .(arm, tag, size_floor)]

plot_long <- melt(by_tag_arm, id.vars = c("arm", "tag", "size_floor", "n_combo"),
                  measure.vars = list(mean = c("fdr_size", "fdr_random"), se = c("fdr_size_se", "fdr_random_se")))
plot_long[, variant := factor(variable, levels = 1:2, labels = c("size-based (>= floor)", "random, same cardinality"))]
plot_long[, arm := factor(arm, levels = c("emmax_consensus", "emmax_simes"))]
plot_long[, tag := factor(tag, levels = c("bgs", "nobgs"))]

p <- ggplot(plot_long, aes(size_floor, mean, colour = variant, fill = variant)) +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(aes(size = n_combo)) +
  facet_grid(arm ~ tag) +
  scale_x_continuous(breaks = c(2, 3, 5, 10, 20, 50)) +
  scale_colour_manual(values = c("size-based (>= floor)" = "#C62828", "random, same cardinality" = "#1565C0"), name = NULL) +
  scale_fill_manual(values = c("size-based (>= floor)" = "#C62828", "random, same cardinality" = "#1565C0"), guide = "none") +
  scale_size_continuous(name = "n combos\n(n_obs>0)", range = c(1, 4)) +
  labs(x = "minimum stage-1 unit size (markers) / matched random-subset cardinality",
       y = "group-null realised FDR (mean +/- SE across rep x env combos)",
       title = "Size-based restriction vs. a matched-cardinality random control",
       subtitle = "Full 7-cell x 2-tag x 10x10 rep-env grid. Identical at floor=2 by construction (the whole candidate pool).\nBeyond that, random is equal or BETTER at every floor -- in every arm x tag combination checked.") +
  theme_bw(11) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(), legend.position = "top")

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "fig_random_removal_control.pdf")
ggsave(OUT, p, width = 10, height = 8, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 10, height = 8, dpi = 200)
say("\n[2] wrote %s (+ .png)\n", OUT)
