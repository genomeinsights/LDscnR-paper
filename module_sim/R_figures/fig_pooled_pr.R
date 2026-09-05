## module_sim/R_figures/fig_pooled_pr.R
##
## Visualize R/05_pool.R's POOLED (sum TP/FP/FN, then one Precision/Recall)
## Precision/Recall by cell, mean +- SE across the 10 ENVIRONMENTS (the true
## replicate axis -- see 05_pool.R's header). PK: "the main comparison is
## between methods, not between bgs" -- coloured by arm (method), faceted by
## tag instead of the other way around.
suppressMessages({library(data.table); library(ggplot2)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim"), "R", "00_config.R"))
say("=== fig_pooled_pr ===\n\n")

POOL_PATH <- file.path(PATHS$out, "05_pool", "pooled_pr.rds")
if (!file.exists(POOL_PATH)) stop("R/05_pool.R has not produced: ", POOL_PATH)
pooled <- readRDS(POOL_PATH)$pooled

## PK: V is selection variance, inversely proportional to selection intensity
## (SI) -- V0.5/V1/V2 = strong/medium/weak SI. c is the dispersal kernel,
## inversely proportional to gene flow (disp) -- c1/c1.5/c2 = high/medium/low
## disp. Cell labels below are built from these, not the raw V/c values.
SI_MAP   <- c(V0.5 = "strong", V1 = "medium", V2 = "weak")
DISP_MAP <- c(c1 = "high", c1.5 = "medium", c2 = "low")
pooled[, c("V_raw", "c_raw") := tstrsplit(cell, "_", fixed = TRUE)]
pooled[, SI   := factor(SI_MAP[V_raw], levels = c("strong", "medium", "weak"))]
pooled[, disp := factor(DISP_MAP[c_raw], levels = c("high", "medium", "low"))]
pooled[, cell_label := sprintf("%s SI, %s disp", SI, disp)]

## Cells sorted by performance (mean Recall across all tag x arm, descending --
## best-performing parameter combination first), per PK's request, rather than
## by SI/disp parameter value.
perf_order <- pooled[, .(mean_recall = mean(Recall)), by = .(cell, cell_label)][order(-mean_recall)]
say("[0] cells ranked by mean Recall (all tag x arm): %s\n",
    paste(sprintf("%s=%.3f", perf_order$cell, perf_order$mean_recall), collapse = ", "))
pooled[, cell_label := factor(cell_label, levels = perf_order$cell_label)]
## NO lfmm_consensus (PK, 2026-09-05: LFMM on complexity-reduced/pooled
## genotypes is not how the method is meant to be used, dropped from
## R/03_scan.R and R/04_score.R). Two single-SNP arms added instead.
ARM_LEVELS <- c("emmax_consensus", "emmax_simes", "lfmm_simes", "emmax_snp", "lfmm_snp")
ARM_COLOURS <- c(emmax_consensus = "#1565C0", emmax_simes = "#26A69A",
                 lfmm_simes = "#7B1FA2", emmax_snp = "#F9A825", lfmm_snp = "#C0392B")
pooled[, arm := factor(arm, levels = ARM_LEVELS)]
pooled[, tag := factor(tag, levels = c("nobgs", "bgs"))]

## long format: one row per (tag, cell, arm, metric), metric in {Precision, Recall}, with its SE.
long <- rbindlist(list(
  pooled[, .(tag, cell = cell_label, arm, metric = "Recall",    value = Recall,    SE = Recall_SE)],
  pooled[, .(tag, cell = cell_label, arm, metric = "Precision", value = Precision, SE = Precision_SE)]
))
long[, metric := factor(metric, levels = c("Recall", "Precision"))]

say("[1] %d (tag,cell,arm,metric) points ; %d Precision NAs (zero-TP-zero-FP cells -- see 05_pool.R)\n",
    nrow(long), sum(is.na(long$value) & long$metric == "Precision"))

p <- ggplot(long, aes(cell, value, colour = arm, group = arm)) +
  geom_line(alpha = 0.5, linewidth = 0.4, position = position_dodge(width = 0.4)) +
  geom_pointrange(aes(ymin = pmax(0, value - SE), ymax = pmin(1, value + SE)),
                  position = position_dodge(width = 0.4), size = 0.3, fatten = 2) +
  facet_grid(metric ~ tag) +
  scale_colour_manual(values = ARM_COLOURS, name = "method") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = NULL,
      title = "Pooled TP/FP scoring: Precision/Recall by cell (all chromosomes pooled -- sum TP/FP/FN, then one ratio)",
      subtitle = "mean +/- SE across the 10 environments (the replicate axis); cells sorted by mean Recall, best-performing first (SI = selection intensity, disp = dispersal/gene flow)") +
  theme_bw(11) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 7, angle = 40, hjust = 1), legend.position = "top")

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "fig_pooled_pr.pdf")
ggsave(OUT, p, width = 10, height = 6.5, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 10, height = 6.5, dpi = 200)
say("\n[2] wrote %s (+ .png)\n", OUT)
