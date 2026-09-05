## module_sim/R_figures/fig_bgs_effect.R
##
## The BGS effect itself, as its own figure (PK, separate from
## fig_popgen_summary.R's tag-contrast lines): log2(bgs/nobgs), PAIRED at
## (cell,rep,env) before pooling (R/07_popgen_pool.R's bgs_effect table) --
## every (cell,rep,env) has both a bgs and a nobgs run, so this is a genuine
## paired contrast, not two independently-pooled tag summaries subtracted
## afterward. 0 = no effect; negative = bgs lower than nobgs. Same cell
## order as fig_pooled_pr.R/fig_popgen_summary.R.
suppressMessages({library(data.table); library(ggplot2)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim"), "R", "00_config.R"))
say("=== fig_bgs_effect ===\n\n")

POOL_PATH <- file.path(PATHS$out, "07_popgen_pool", "popgen_summary.rds")
if (!file.exists(POOL_PATH)) stop("R/07_popgen_pool.R has not produced: ", POOL_PATH)
bgs_effect <- readRDS(POOL_PATH)$bgs_effect

SI_MAP   <- c(V0.5 = "strong", V1 = "medium", V2 = "weak")
DISP_MAP <- c(c1 = "high", c1.5 = "medium", c2 = "low")
bgs_effect[, c("V_raw", "c_raw") := tstrsplit(cell, "_", fixed = TRUE)]
bgs_effect[, SI   := factor(SI_MAP[V_raw], levels = c("strong", "medium", "weak"))]
bgs_effect[, disp := factor(DISP_MAP[c_raw], levels = c("high", "medium", "low"))]
bgs_effect[, cell_label := sprintf("%s SI, %s disp", SI, disp)]
cell_order <- unique(bgs_effect[, .(cell, cell_label, disp, SI)])
setorder(cell_order, disp, SI)
say("[0] cell order (disp high->low, then SI high->low): %s\n", paste(cell_order$cell, collapse = ", "))
bgs_effect[, cell_label := factor(cell_label, levels = cell_order$cell_label)]

long <- rbindlist(list(
  bgs_effect[, .(cell = cell_label, metric = "Fst",                   value = log2_Fst,              SE = log2_Fst_SE)],
  bgs_effect[, .(cell = cell_label, metric = "detectable QTN (n)",    value = log2_n_qtn_detectable, SE = log2_n_qtn_detectable_SE)],
  bgs_effect[, .(cell = cell_label, metric = "total Va",              value = log2_Va_total,          SE = log2_Va_total_SE)],
  bgs_effect[, .(cell = cell_label, metric = "local adaptation (R2)", value = log2_local_adapt_r2,    SE = log2_local_adapt_r2_SE)]
))
long[, metric := factor(metric, levels = c("Fst", "detectable QTN (n)", "total Va", "local adaptation (R2)"))]

say("[1] %d (cell,metric) points\n", nrow(long))

p <- ggplot(long, aes(cell, value)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_pointrange(aes(ymin = value - SE, ymax = value + SE), colour = "#7B1FA2", size = 0.4) +
  facet_wrap(~metric, scales = "free_y", ncol = 1) +
  labs(x = NULL, y = expression(log[2] * (bgs / nobgs)),
      title = "The effect of background selection (BGS)",
      subtitle = "paired log2(bgs/nobgs) at (cell,rep,env), mean +/- SE over the 10 environments; 0 = no effect, negative = bgs lower; cells sorted by dispersal (high->low) then selection intensity (high->low)") +
  theme_bw(11) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8, angle = 30, hjust = 1))

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "fig_bgs_effect.pdf")
ggsave(OUT, p, width = 8, height = 10, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 8, height = 10, dpi = 200)
say("\n[2] wrote %s (+ .png)\n", OUT)
