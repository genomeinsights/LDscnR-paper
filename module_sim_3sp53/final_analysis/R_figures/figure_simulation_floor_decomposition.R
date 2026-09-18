## final_analysis/R_figures/figure_simulation_floor_decomposition.R
##
## Analysis 3 figure: how much of the precision/recall change across the
## floor grid comes from filtering small Stage-1 units (arms A -> B) versus
## from LD-aggregation itself (arms B -> C/D). Four panels (tests as % of the
## unrestricted-marker baseline, precision, recall, precision x recall)
## across FLOOR_GRID, one line per arm, EMMAX arms A-D and LFMM arms A-C
## shown together (linetype = engine).
## Grand-pooled across all 7 cells x 2 BGS treatments, with map-cluster
## bootstrap CIs computed here from R/12's saved raw per-combo table (the
## per-(tag,cell) CIs in results/simulation_floor_sweep.tsv are a different,
## finer grain than the grand-pooled view this figure needs).
##
## [!] The "Tests" panel plots tests as a PERCENTAGE of that method's own
## unrestricted-marker (arm A) test count at the same floor, not a raw total
## across all 1,400 combos (PK, 2026-09-17: the raw total is dominated by
## the marker-count scale difference between arms and is hard to read as a
## cost).
##
## [!] "Markers retained" panel REMOVED, "Precision x Recall" panel ADDED
## (PK, 2026-09-17, follow-up request). frac_markers_retained stays in the
## underlying TSV for reference (same convention as realised FDP below), just
## without its own panel. Precision x Recall is a single-number summary of
## the trade-off this analysis is about; its CI comes from the SAME
## bootstrap draws as the Precision/Recall panels (prec_b * rec_b per
## replicate), not from naively multiplying the two marginal CIs.
suppressMessages({library(data.table); library(ggplot2)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "helpers_stage2_truth.R"))
say("=== figure_simulation_floor_decomposition ===\n\n")

raw <- readRDS("out_final_v1/12_floor_decomposition_raw.rds")

ARM_LABELS <- c(A_unrestricted_marker = "Unrestricted marker", B_floor_filtered_marker = "Floor-filtered marker",
                C_stage1_simes = "Stage-1 Simes", D_stage1_consensus = "Stage-1 consensus")
ARM_COLOURS <- c("Unrestricted marker" = "#26A69A", "Floor-filtered marker" = "#546E7A",
                 "Stage-1 Simes" = "#7B1FA2", "Stage-1 consensus" = "#1565C0")

## grand-pooled point estimate + map-cluster bootstrap CI, per (method, arm, floor)
rep_dt <- raw[, .(TP = sum(TP), FP = sum(FP), n_detectable_qtn = sum(n_detectable_qtn),
                  n_recovered = sum(n_recovered), n_tests = sum(n_tests)),
             by = .(method, arm, floor, rep)]
strata <- unique(rep_dt[, .(method, arm, floor)])
rows <- vector("list", nrow(strata))
for (i in seq_len(nrow(strata))) {
  s <- strata[i]
  rd <- rep_dt[method == s$method & arm == s$arm & floor == s$floor][match(REPS_ALL, rep)]
  rd[is.na(TP), `:=`(TP = 0L, FP = 0L, n_detectable_qtn = 0L, n_recovered = 0L, n_tests = 0L)]
  bs <- bootstrap_rep_matrix(as.matrix(rd[, .(TP, FP, n_detectable_qtn, n_recovered, n_tests)]),
                             N_BOOTSTRAP, SEEDS[["bootstrap"]])
  prec_b <- bs[, "TP"] / pmax(bs[, "TP"] + bs[, "FP"], 1)
  rec_b <- bs[, "n_recovered"] / pmax(bs[, "n_detectable_qtn"], 1)
  fdp_b <- bs[, "FP"] / pmax(bs[, "TP"] + bs[, "FP"], 1)
  prxrec_b <- prec_b * rec_b
  point_precision <- sum(rd$TP) / pmax(sum(rd$TP) + sum(rd$FP), 1)
  point_recall <- sum(rd$n_recovered) / pmax(sum(rd$n_detectable_qtn), 1)
  rows[[i]] <- data.table(
    method = s$method, arm = s$arm, floor = s$floor,
    n_tests = sum(rd$n_tests), precision = point_precision,
    precision_ci_lo = ci_quantile(prec_b)[1], precision_ci_hi = ci_quantile(prec_b)[2],
    recall = point_recall,
    recall_ci_lo = ci_quantile(rec_b)[1], recall_ci_hi = ci_quantile(rec_b)[2],
    fdp = sum(rd$FP) / pmax(sum(rd$TP) + sum(rd$FP), 1),
    fdp_ci_lo = ci_quantile(fdp_b)[1], fdp_ci_hi = ci_quantile(fdp_b)[2],
    prec_x_rec = point_precision * point_recall,
    prec_x_rec_ci_lo = ci_quantile(prxrec_b)[1], prec_x_rec_ci_hi = ci_quantile(prxrec_b)[2])
}
gp <- rbindlist(rows)

## mean fraction of markers retained (a ratio -- averaged across combos, not
## summed) and tests as a percentage of the SAME method's unrestricted-marker
## (arm A) test count at the SAME floor.
mfr <- raw[, .(mean_frac_markers_retained = mean(frac_markers_retained)), by = .(method, arm, floor)]
gp <- merge(gp, mfr, by = c("method", "arm", "floor"))
tests_a <- gp[arm == "A_unrestricted_marker", .(method, floor, n_tests_a = n_tests)]
gp <- merge(gp, tests_a, by = c("method", "floor"))
gp[, pct_of_unrestricted_tests := 100 * n_tests / n_tests_a]

fwrite(gp, "results/figure_simulation_floor_decomposition_data.tsv", sep = "\t")

gp[, engine := ifelse(method == "lfmm", "LFMM", "EMMAX")]
gp <- gp[!(method == "lfmm" & arm == "D_stage1_consensus")]   ## no LFMM consensus arm
gp[, arm_label := factor(ARM_LABELS[arm], levels = ARM_LABELS)]

## No separate "Realised FDP" panel: with pooled-count ratios, FDP = FP/(TP+FP)
## = 1 - TP/(TP+FP) = 1 - Precision EXACTLY (same TP/FP denominator for both),
## so it would be a mirror-flipped duplicate of the Precision panel, not new
## information (PK, 2026-09-17). realised_fdp/fdp_ci_* stay in the underlying
## TSV (figure_simulation_floor_decomposition_data.tsv, simulation_floor_
## sweep.tsv) since that's this pipeline's established column name, just not
## given its own panel here.
long <- rbindlist(list(
  gp[, .(engine, arm_label, floor, metric = "Tests (% of unrestricted)", value = pct_of_unrestricted_tests, lo = NA_real_, hi = NA_real_)],
  gp[, .(engine, arm_label, floor, metric = "Precision", value = precision, lo = precision_ci_lo, hi = precision_ci_hi)],
  gp[, .(engine, arm_label, floor, metric = "Recall", value = recall, lo = recall_ci_lo, hi = recall_ci_hi)],
  gp[, .(engine, arm_label, floor, metric = "Precision x Recall", value = prec_x_rec, lo = prec_x_rec_ci_lo, hi = prec_x_rec_ci_hi)]
))
long[, metric := factor(metric, levels = c("Tests (% of unrestricted)", "Precision", "Recall", "Precision x Recall"))]

p <- ggplot(long, aes(floor, value, colour = arm_label, linetype = engine)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = arm_label), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.4) +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = c(1, 2, 3, 5, 10, 20), trans = "log2") +
  scale_colour_manual(values = ARM_COLOURS, name = NULL) +
  scale_fill_manual(values = ARM_COLOURS, guide = "none") +
  scale_linetype_manual(values = c(EMMAX = "solid", LFMM = "22"), name = NULL) +
  labs(x = "Minimum Stage-1 unit size (floor)", y = NULL,
       title = "Decomposing filtering (floor) from LD-aggregation (Simes/consensus) across the floor grid") +
  theme_bw(11) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(), legend.position = "top")

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "figure_simulation_floor_decomposition.pdf")
ggsave(OUT, p, width = 9, height = 7, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 9, height = 7, dpi = 200)
say("wrote %s (+ .png) and results/figure_simulation_floor_decomposition_data.tsv\n", OUT)
