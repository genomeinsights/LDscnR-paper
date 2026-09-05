## module_sim/R_figures/fig_bgs_windows.R
##
## The n_snp-ratio BGS-vs-recombination estimator from
## /Volumes/Nemo/Nemo_sim/bgs_effect_newmaps/measure_bgs_effect.R, applied to
## this repo's own grid: B_obs = n_snp(bgs)/n_snp(nobgs) per window, by
## recombination-rate quintile (Q1 = lowest), one line per cell (never pooled
## across cells -- see R/11_bgs_windows_pool.R's header).
suppressMessages({library(data.table); library(ggplot2)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim"), "R", "00_config.R"))
say("=== fig_bgs_windows ===\n\n")

POOL_PATH <- file.path(PATHS$out, "11_bgs_windows_pool", "bgs_windows_summary.rds")
if (!file.exists(POOL_PATH)) stop("R/11_bgs_windows_pool.R has not produced: ", POOL_PATH)
x <- readRDS(POOL_PATH)
quintiles <- x$quintiles
summary_by_cell <- x$summary_by_cell

SI_MAP   <- c(V0.5 = "strong", V1 = "medium", V2 = "weak")
DISP_MAP <- c(c1 = "high", c1.5 = "medium", c2 = "low")
quintiles[, c("V_raw", "c_raw") := tstrsplit(cell, "_", fixed = TRUE)]
quintiles[, SI   := factor(SI_MAP[V_raw], levels = c("strong", "medium", "weak"))]
quintiles[, disp := factor(DISP_MAP[c_raw], levels = c("high", "medium", "low"))]
quintiles[, cell_label := sprintf("%s SI, %s disp", SI, disp)]
cell_order <- unique(quintiles[, .(cell, cell_label, disp, SI)])
setorder(cell_order, disp, SI)
say("[0] cell order (disp high->low, then SI high->low): %s\n", paste(cell_order$cell, collapse = ", "))
quintiles[, cell_label := factor(cell_label, levels = cell_order$cell_label)]

say("[1] %d (cell,quintile) points\n", nrow(quintiles))

p <- ggplot(quintiles, aes(q, B, colour = cell_label, group = cell_label)) +
  geom_hline(yintercept = 1, colour = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:5, labels = paste0("Q", 1:5)) +
  scale_colour_brewer(palette = "Dark2", name = NULL) +
  labs(x = "recombination-rate quintile (Q1 = lowest)", y = expression(median~B[obs] == n[snp](bgs) / n[snp](nobgs)),
      title = "BGS reduces segregating-site count more where recombination is low",
      subtitle = "paired per 500kb window, same window/rep/env in both arms; 1 = no BGS effect") +
  theme_bw(11) +
  theme(panel.grid.minor = element_blank(), legend.position = "right")

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "fig_bgs_windows.pdf")
ggsave(OUT, p, width = 9, height = 5.5, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 9, height = 5.5, dpi = 200)
say("\n[2] wrote %s (+ .png)\n", OUT)
