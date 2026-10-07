## final_analysis/R_figures/figureS_simulation_floor_decomposition_c1.R
##
## c=1-pooled version of the floor-decomposition figure (PK, 2026-09-19):
## same four panels as the grand-pooled figure_simulation_floor_decomposition.R
## (tests as % of unrestricted marker-wise, precision, recall, precision x
## recall), but restricted to the 600 primary high-gene-flow simulations
## (V0.5_c1, V1_c1, V2_c1) instead of all 7 cells. Reads
## results/simulation_floor_sweep_c1.tsv (already computed by
## R/16_floor_decomposition_c1.R, itself a downstream summary of R/12's
## cached raw table) -- no bootstrap is recomputed here, this script only
## plots numbers that already exist on disk.
##
## Manuscript version only, no title, no in-plot caption -- that text
## belongs in the manuscript's own LaTeX \caption{} and is written out
## verbatim to the companion _caption.txt file (same convention as
## figureS_simulation_env_structure_alignment.R).
##
## The existing figure_simulation_floor_decomposition_by_cell.pdf (per-cell
## breakdown across ALL 7 cells, showing where methods fail at c=1.5/c=2) is
## untouched by this script and remains the complementary view.
suppressMessages({library(data.table); library(ggplot2)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))
say("=== figureS_simulation_floor_decomposition_c1 ===\n\n")

gp <- fread("results/simulation_floor_sweep_c1.tsv")

ARM_LABELS <- c(A_unrestricted_marker = "Unrestricted marker", B_floor_filtered_marker = "Floor-filtered marker",
                C_stage1_simes = "Stage-1 Simes", D_stage1_consensus = "Stage-1 consensus")
ARM_COLOURS <- c("Unrestricted marker" = "#26A69A", "Floor-filtered marker" = "#546E7A",
                 "Stage-1 Simes" = "#7B1FA2", "Stage-1 consensus" = "#1565C0")

## tests as % of the SAME method's unrestricted-marker (arm A) test count at
## the SAME floor -- same convention as the grand-pooled figure.
tests_a <- gp[arm == "A_unrestricted_marker", .(method, floor, n_tests_a = n_tests)]
gp <- merge(gp, tests_a, by = c("method", "floor"))
gp[, pct_of_unrestricted_tests := 100 * n_tests / n_tests_a]

gp[, engine := ifelse(method == "lfmm", "LFMM", "EMMAX")]
gp[, arm_label := factor(ARM_LABELS[arm], levels = ARM_LABELS)]

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
  labs(x = "Minimum Stage-1 unit size (floor)", y = NULL) +
  theme_bw(11) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(), legend.position = "top")

CAPTION <- paste0("Floor decomposition pooled over the 600 primary high-gene-flow (c=1) simulations only ",
                  "(V0.5_c1, V1_c1, V2_c1; both BGS treatments, all ten map/burn-in identities, all ten ",
                  "environmental continuations). Precision, recall and precision x recall shown with 95% ",
                  "percentile map-cluster bootstrap intervals (2,000 replicates, ten map identities resampled ",
                  "as the cluster, all three selection strengths/both BGS treatments/all envs sharing a map ",
                  "identity always resampled together; the same bootstrap draws are used for every method, ",
                  "arm and floor, so contrasts are paired). No confidence band on the tests panel, which is a ",
                  "deterministic percentage of the unrestricted marker-wise test count at the same floor, not ",
                  "a bootstrapped quantity. EMMAX (solid) and LFMM (dashed) shown together; LFMM has no ",
                  "Stage-1 consensus arm. See results/simulation_floor_sweep_c1.tsv and ",
                  "results/simulation_floor_contrasts_c1.tsv for the full point estimates, intervals and paired ",
                  "contrasts underlying this figure, and ADDITIONAL_ANALYSES_AUDIT.md for the interpretation.")

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "figureS_simulation_floor_decomposition_c1.pdf")
ggsave(OUT, p, width = 9, height = 7, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 9, height = 7, dpi = 200)
CAPTION_TXT <- sub("\\.pdf$", "_caption.txt", OUT)
writeLines(CAPTION, CAPTION_TXT)
say("wrote %s (+ .png) and %s (LaTeX-ready caption text)\n", OUT, CAPTION_TXT)
