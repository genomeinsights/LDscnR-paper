## module_sim/R_figures/fig_popgen_summary.R
##
## ONE figure summarizing the simulations themselves (not detection
## performance): Fst, mean number of detectable QTN, total Va, and a
## local-adaptation proxy (breeding value ~ environment R^2), by cell,
## coloured by tag -- the tag contrast at each cell IS the BGS-effect
## estimate PK asked for, rather than a fifth derived quantity. Same cell
## order (dispersal high->low, then selection intensity high->low) and
## style as fig_pooled_pr.R, for a directly comparable reading.
suppressMessages({library(data.table); library(ggplot2)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim"), "R", "00_config.R"))
say("=== fig_popgen_summary ===\n\n")

POOL_PATH <- file.path(PATHS$out, "07_popgen_pool", "popgen_summary.rds")
if (!file.exists(POOL_PATH)) stop("R/07_popgen_pool.R has not produced: ", POOL_PATH)
pooled <- readRDS(POOL_PATH)$pooled

## same SI/disp mapping and cell ordering as fig_pooled_pr.R
SI_MAP   <- c(V0.5 = "strong", V1 = "medium", V2 = "weak")
DISP_MAP <- c(c1 = "high", c1.5 = "medium", c2 = "low")
pooled[, c("V_raw", "c_raw") := tstrsplit(cell, "_", fixed = TRUE)]
pooled[, SI   := factor(SI_MAP[V_raw], levels = c("strong", "medium", "weak"))]
pooled[, disp := factor(DISP_MAP[c_raw], levels = c("high", "medium", "low"))]
pooled[, cell_label := sprintf("%s SI, %s disp", SI, disp)]
cell_order <- unique(pooled[, .(cell, cell_label, disp, SI)])
setorder(cell_order, disp, SI)
say("[0] cell order (disp high->low, then SI high->low): %s\n", paste(cell_order$cell, collapse = ", "))
pooled[, cell_label := factor(cell_label, levels = cell_order$cell_label)]
pooled[, tag := factor(tag, levels = c("nobgs", "bgs"))]

long <- rbindlist(list(
  pooled[, .(tag, cell = cell_label, metric = "Fst",                        value = Fst,              SE = Fst_SE)],
  pooled[, .(tag, cell = cell_label, metric = "detectable QTN (n)",         value = n_qtn_detectable, SE = n_qtn_detectable_SE)],
  pooled[, .(tag, cell = cell_label, metric = "total Va",                   value = Va_total,         SE = Va_total_SE)],
  pooled[, .(tag, cell = cell_label, metric = "local adaptation (R2)",      value = local_adapt_r2,   SE = local_adapt_r2_SE)]
))
long[, metric := factor(metric, levels = c("Fst", "detectable QTN (n)", "total Va", "local adaptation (R2)"))]

say("[1] %d (tag,cell,metric) points\n", nrow(long))

p <- ggplot(long, aes(cell, value, colour = tag, group = tag)) +
  geom_line(alpha = 0.5, linewidth = 0.4, position = position_dodge(width = 0.3)) +
  geom_pointrange(aes(ymin = value - SE, ymax = value + SE),
                  position = position_dodge(width = 0.3), size = 0.3) +
  facet_wrap(~metric, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = c(nobgs = "#1565C0", bgs = "#C0392B")) +
  labs(x = NULL, y = NULL,
      title = "Simulation properties by cell: Fst, detectable QTN, Va, local adaptation",
      subtitle = "mean +/- SE across the 10 environments; cells sorted by dispersal (high->low) then selection intensity (high->low); tag contrast = BGS effect") +
  theme_bw(11) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8, angle = 30, hjust = 1), legend.position = "top")

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "fig_popgen_summary.pdf")
ggsave(OUT, p, width = 8, height = 11, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 8, height = 11, dpi = 200)
say("\n[2] wrote %s (+ .png)\n", OUT)
