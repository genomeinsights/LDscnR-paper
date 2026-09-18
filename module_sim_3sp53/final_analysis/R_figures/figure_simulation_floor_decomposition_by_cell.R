## final_analysis/R_figures/figure_simulation_floor_decomposition_by_cell.R
##
## Follow-up to figure_simulation_floor_decomposition.R (grand-pooled across
## all 7 cells): does "SNP-based (unrestricted/floor-filtered marker) never
## beats LD-aggregation (Stage-1 Simes/consensus) on Precision x Recall, at a
## given floor" hold within EVERY demographic cell individually, or is the
## pooled result driven by one or two cells? Same pooled-count-then-bootstrap
## construction as the grand-pooled figure, just with `cell` added to the
## grouping (still pooled across BOTH BGS treatments and all envs within a
## cell -- splitting further by tag as well would leave too few reps per
## cell x tag x env stratum for a meaningful bootstrap).
##
## Bootstrap CI is clustered by rep ONLY (not tag x rep): BGS/no-BGS rows for
## the same (cell, rep) are a paired map/burn-in draw and must resample
## together, same principle as R/13's and R/14's bootstrap-pairing fixes
## (PK, 2026-09-18 review, points 3/5).
suppressMessages({library(data.table); library(ggplot2)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "helpers_stage2_truth.R"))
say("=== figure_simulation_floor_decomposition_by_cell ===\n\n")

raw <- readRDS("out_final_v1/12_floor_decomposition_raw.rds")

ARM_LABELS <- c(A_unrestricted_marker = "Unrestricted marker", B_floor_filtered_marker = "Floor-filtered marker",
                C_stage1_simes = "Stage-1 Simes", D_stage1_consensus = "Stage-1 consensus")
ARM_COLOURS <- c("Unrestricted marker" = "#26A69A", "Floor-filtered marker" = "#546E7A",
                 "Stage-1 Simes" = "#7B1FA2", "Stage-1 consensus" = "#1565C0")

## pool BOTH tag values and all envs within (cell, rep) -- i.e. one row per
## (cell, method, arm, floor, rep), keeping rep as the bootstrap-resampling
## unit (paired across tag by construction, since both tag rows are already
## summed into it before the bootstrap ever sees them).
rep_dt <- raw[, .(TP = sum(TP), FP = sum(FP), n_detectable_qtn = sum(n_detectable_qtn),
                  n_recovered = sum(n_recovered), n_tests = sum(n_tests)),
             by = .(cell, method, arm, floor, rep)]
strata <- unique(rep_dt[, .(cell, method, arm, floor)])
rows <- vector("list", nrow(strata))
for (i in seq_len(nrow(strata))) {
  s <- strata[i]
  rd <- rep_dt[cell == s$cell & method == s$method & arm == s$arm & floor == s$floor][match(REPS_ALL, rep)]
  rd[is.na(TP), `:=`(TP = 0L, FP = 0L, n_detectable_qtn = 0L, n_recovered = 0L, n_tests = 0L)]
  bs <- bootstrap_rep_matrix(as.matrix(rd[, .(TP, FP, n_detectable_qtn, n_recovered, n_tests)]),
                             N_BOOTSTRAP, SEEDS[["bootstrap"]])
  prec_b <- bs[, "TP"] / pmax(bs[, "TP"] + bs[, "FP"], 1)
  rec_b <- bs[, "n_recovered"] / pmax(bs[, "n_detectable_qtn"], 1)
  prxrec_b <- prec_b * rec_b
  point_precision <- sum(rd$TP) / pmax(sum(rd$TP) + sum(rd$FP), 1)
  point_recall <- sum(rd$n_recovered) / pmax(sum(rd$n_detectable_qtn), 1)
  rows[[i]] <- data.table(
    cell = s$cell, method = s$method, arm = s$arm, floor = s$floor,
    precision = point_precision, recall = point_recall,
    prec_x_rec = point_precision * point_recall,
    prec_x_rec_ci_lo = ci_quantile(prxrec_b)[1], prec_x_rec_ci_hi = ci_quantile(prxrec_b)[2])
}
gp <- rbindlist(rows)
fwrite(gp, "results/figure_simulation_floor_decomposition_by_cell_data.tsv", sep = "\t")

gp[, engine := ifelse(method == "lfmm", "LFMM", "EMMAX")]
gp <- gp[!(method == "lfmm" & arm == "D_stage1_consensus")]   ## no LFMM consensus arm
gp[, arm_label := factor(ARM_LABELS[arm], levels = ARM_LABELS)]

## Cell naming: V is INVERSELY proportional to selection intensity (V0.5 =
## strong, V1 = medium, V2 = weak) and c is INVERSELY proportional to gene
## flow (c1 = high, c1.5 = medium, c2 = low) (PK, 2026-09-18). Facets sorted
## first by gene flow (high -> low), then by selection intensity (strong ->
## weak) within each gene-flow level; the design is only 7 of the 9 possible
## (V, c) combinations (no V1_c2/V2_c2), so the c2 group has one facet.
SELECTION_LABEL <- c(V0.5 = "Strong selection", V1 = "Medium selection", V2 = "Weak selection")
FLOW_LABEL <- c(c1 = "High gene flow", c1.5 = "Medium gene flow", c2 = "Low gene flow")
CELL_ORDER <- c("V0.5_c1", "V1_c1", "V2_c1", "V0.5_c1.5", "V1_c1.5", "V2_c1.5", "V0.5_c2")
cell_v <- sub("_c.*$", "", CELL_ORDER)
cell_c <- sub("^V[0-9.]+_", "", CELL_ORDER)
CELL_LABELS <- setNames(paste0(SELECTION_LABEL[cell_v], "\n", FLOW_LABEL[cell_c]), CELL_ORDER)
gp[, cell := factor(cell, levels = CELL_ORDER, labels = CELL_LABELS[CELL_ORDER])]

p <- ggplot(gp, aes(floor, prec_x_rec, colour = arm_label, linetype = engine)) +
  geom_ribbon(aes(ymin = prec_x_rec_ci_lo, ymax = prec_x_rec_ci_hi, fill = arm_label), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.2) +
  facet_wrap(~cell, ncol = 4, scales = "free_y") +
  scale_x_continuous(breaks = c(1, 2, 3, 5, 10, 20), trans = "log2") +
  scale_colour_manual(values = ARM_COLOURS, name = NULL) +
  scale_fill_manual(values = ARM_COLOURS, guide = "none") +
  scale_linetype_manual(values = c(EMMAX = "solid", LFMM = "22"), name = NULL) +
  labs(x = "Minimum Stage-1 unit size (floor)", y = "Precision x Recall",
       title = "Precision x Recall across the floor grid, split by demographic cell") +
  theme_bw(11) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(), legend.position = "top")

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "figure_simulation_floor_decomposition_by_cell.pdf")
ggsave(OUT, p, width = 11, height = 7, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 11, height = 7, dpi = 200)
say("wrote %s (+ .png) and results/figure_simulation_floor_decomposition_by_cell_data.tsv\n", OUT)
