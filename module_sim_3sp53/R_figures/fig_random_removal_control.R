## module_sim_3sp53/R_figures/fig_random_removal_control.R
##
## Visualize R/14_random_removal_control.R: size-based vs matched-
## cardinality-random realised_fdr by floor. PK, 2026-09-09: "if we
## randomly remove the same number of outlier clusters (but independently
## of size), would we also see a decline in FPs, and is the decline in
## fig_structured_null_sizesweep relative to that?"
suppressMessages({library(data.table); library(ggplot2)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53"), "R", "00_config.R"))
say("=== fig_random_removal_control ===\n\n")

FILES <- Sys.glob(file.path(PATHS$out, "14_random_removal_control", "rrc_*.rds"))
if (!length(FILES)) stop("R/14_random_removal_control.R has not produced any out/14_random_removal_control/rrc_*.rds")
all <- rbindlist(lapply(FILES, readRDS))
all[, fdr_size := n_surr_size / max(n_obs_size, 1), by = seq_len(nrow(all))]
say("[1] %d rows, cells: %s\n", nrow(all), paste(sort(unique(all$cell)), collapse = ", "))

## pooled (across the cells actually run) -- primary panel
pooled <- all[n_obs_size > 0, .(fdr_size = mean(fdr_size), fdr_random = mean(fdr_random_conditional, na.rm = TRUE),
                                 n_combo = .N), by = size_floor]
pooled_long <- melt(pooled, id.vars = c("size_floor", "n_combo"),
                    measure.vars = c("fdr_size", "fdr_random"),
                    variable.name = "variant", value.name = "realised_fdr")
pooled_long[, variant := factor(variant, levels = c("fdr_size", "fdr_random"),
                                labels = c("size-based (>= floor)", "random, same cardinality"))]
pooled_long[, cell := "pooled (all cells run)"]

## per-cell panels -- so the pooled line isn't hiding cell-to-cell variation
percell <- all[n_obs_size > 0, .(fdr_size = mean(fdr_size), fdr_random = mean(fdr_random_conditional, na.rm = TRUE),
                                  n_combo = .N), by = .(cell, size_floor)]
percell_long <- melt(percell, id.vars = c("cell", "size_floor", "n_combo"),
                     measure.vars = c("fdr_size", "fdr_random"),
                     variable.name = "variant", value.name = "realised_fdr")
percell_long[, variant := factor(variant, levels = c("fdr_size", "fdr_random"),
                                 labels = c("size-based (>= floor)", "random, same cardinality"))]

plot_dt <- rbindlist(list(pooled_long, percell_long), use.names = TRUE)
plot_dt[, cell := factor(cell, levels = c("pooled (all cells run)", sort(unique(percell$cell))))]

p <- ggplot(plot_dt, aes(size_floor, realised_fdr, colour = variant)) +
  geom_line(linewidth = 0.7) +
  geom_point(aes(size = n_combo)) +
  facet_wrap(~cell, nrow = 1) +
  scale_x_continuous(breaks = c(2, 3, 5, 10, 20, 50)) +
  scale_colour_manual(values = c("size-based (>= floor)" = "#C62828", "random, same cardinality" = "#1565C0"),
                      name = NULL) +
  scale_size_continuous(name = "n combos\n(n_obs>0)", range = c(1, 4)) +
  labs(x = "minimum stage-1 unit size (markers) / matched random-subset cardinality",
       y = "group-null realised FDR",
       title = "Size-based restriction vs. a matched-cardinality random control",
       subtitle = "bgs tag, emmax_consensus, group null, full 10x10 rep x env grid per cell. Identical at floor=2 by construction\n(the whole candidate pool). Beyond that, random is equal or BETTER at every floor in every cell checked --\nsize-based restriction shows no advantage, and the gap widens with floor.") +
  theme_bw(11) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(), legend.position = "top")

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "fig_random_removal_control.pdf")
ggsave(OUT, p, width = 12, height = 5, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 12, height = 5, dpi = 200)
say("\n[2] wrote %s (+ .png)\n", OUT)
