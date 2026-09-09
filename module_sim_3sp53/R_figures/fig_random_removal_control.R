## module_sim_3sp53/R_figures/fig_random_removal_control.R
##
## Visualize R/14_random_removal_control.R: size-based vs matched-
## cardinality-random realised_fdr by floor. PK, 2026-09-09: "if we
## randomly remove the same number of outlier clusters (but independently
## of size), would we also see a decline in FPs, and is the decline in
## fig_structured_null_sizesweep relative to that?"
##
## Full grid as of this run: all 7 cells x both tags (bgs, nobgs), full
## 10x10 rep x env each, emmax_consensus / group null only. Faceted by
## TAG (pooled across cells within each) rather than by cell -- the
## bgs-vs-nobgs comparison is the interesting one: despite nobgs having a
## MUCH stronger ground-truth size-truth relationship than bgs (fig_fp_by_
## size: 55% vs 80% FP at size 50+), the two tags show essentially
## IDENTICAL null-vs-size behaviour here, which is what retracted the
## original recombination/BGS-confound explanation for this effect -- see
## README.md.
suppressMessages({library(data.table); library(ggplot2)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53"), "R", "00_config.R"))
say("=== fig_random_removal_control ===\n\n")

FILES <- Sys.glob(file.path(PATHS$out, "14_random_removal_control", "rrc_*.rds"))
if (!length(FILES)) stop("R/14_random_removal_control.R has not produced any out/14_random_removal_control/rrc_*.rds")
all <- rbindlist(lapply(FILES, readRDS), fill = TRUE)
all[, fdr_size := n_surr_size / max(n_obs_size, 1), by = seq_len(nrow(all))]
say("[1] %d rows, %d cells x %d tags\n", nrow(all), uniqueN(all$cell), uniqueN(all$tag))

by_tag <- all[n_obs_size > 0, .(fdr_size = mean(fdr_size), fdr_random = mean(fdr_random_conditional, na.rm = TRUE),
                                 n_combo = .N), by = .(tag, size_floor)]
by_tag[, panel := tag]

pooled <- all[n_obs_size > 0, .(fdr_size = mean(fdr_size), fdr_random = mean(fdr_random_conditional, na.rm = TRUE),
                                 n_combo = .N), by = size_floor]
pooled[, panel := "pooled (both tags)"]

plot_dt <- rbindlist(list(pooled, by_tag[, .(panel, size_floor, fdr_size, fdr_random, n_combo)]), use.names = TRUE)
plot_dt[, panel := factor(panel, levels = c("pooled (both tags)", "bgs", "nobgs"))]
plot_long <- melt(plot_dt, id.vars = c("panel", "size_floor", "n_combo"),
                  measure.vars = c("fdr_size", "fdr_random"),
                  variable.name = "variant", value.name = "realised_fdr")
plot_long[, variant := factor(variant, levels = c("fdr_size", "fdr_random"),
                              labels = c("size-based (>= floor)", "random, same cardinality"))]

p <- ggplot(plot_long, aes(size_floor, realised_fdr, colour = variant)) +
  geom_line(linewidth = 0.7) +
  geom_point(aes(size = n_combo)) +
  facet_wrap(~panel, nrow = 1) +
  scale_x_continuous(breaks = c(2, 3, 5, 10, 20, 50)) +
  scale_colour_manual(values = c("size-based (>= floor)" = "#C62828", "random, same cardinality" = "#1565C0"),
                      name = NULL) +
  scale_size_continuous(name = "n combos\n(n_obs>0)", range = c(1, 4)) +
  labs(x = "minimum stage-1 unit size (markers) / matched random-subset cardinality",
       y = "group-null realised FDR",
       title = "Size-based restriction vs. a matched-cardinality random control",
       subtitle = "emmax_consensus, group null, full 7-cell x 2-tag x 10x10 rep-env grid. Identical at floor=2 by construction\n(the whole candidate pool). Beyond that, random is equal or BETTER at every floor -- in bgs AND nobgs, despite\nnobgs having a much stronger ground-truth size-truth relationship (fig_fp_by_size). Size shows no advantage.") +
  theme_bw(11) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(), legend.position = "top")

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "fig_random_removal_control.pdf")
ggsave(OUT, p, width = 10, height = 5, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 10, height = 5, dpi = 200)
say("\n[2] wrote %s (+ .png)\n", OUT)
