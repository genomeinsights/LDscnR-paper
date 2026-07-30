######################################################
## Simulated-data OR / PR-AUC driver
##
## Reproduces the two validated test figures:
##   * AUC-PR* trajectories  (Rplot12.png)  -> plot_auc_trajectories()
##   * genome-wide C-score + q panels (Rplot13.png) -> plot_C_score_genome()
##
## Runs on the per-file parsed simulation objects (GTs, map, env,
## LD_decay, ld_ws) produced by Parse_sim_data.R. The OR-clustering /
## C-score engine is paper-specific (not in the LDscnR package); the
## LD-decay and ld_w that it consumes were computed at parse time with
## LDscnR::compute_LD_decay() / LDscnR::compute_ld_w().
##
## Run from the repository root.
######################################################

library(LDscnR)
library(data.table)
library(igraph)
library(parallel)
library(ggplot2)
library(wesanderson)
library(PRROC)
library(patchwork)

source("R/define_ORs_functions.R")        # generic clustering / C-score engine
source("R/Outlier_regions_simulation.R")  # simulation layer (ground truth, PR/AUC)

#----------------------------------------------------------
# 1) Configuration
#----------------------------------------------------------
parsed_folder <- "/Volumes/Nemo/Nemo_sim/parsed_sim_data2"

## One replicate "genome" = one sim/bgs/V/c/env across its chromosomes
## (chr1..chr10). The genome C-score figure facets by chromosome, so it
## expects exactly one file per chromosome -- fix everything but chr here.
sel_pattern <- "^adapt_bgs_chr[0-9]+_V0.5_c2_env2\\.rds$"

p_cols       <- c(EMX = "emx_p", LFMM = "lfmm_p")
r2_grid      <- c(0.9)                                  # interpreted as rho, not raw r^2
th_ldw_grid  <- c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
lmin_grid    <- c(1, 2, 4, 8, 16)
alpha_grid   <- c(0.0001, 0.001, 0.01, 0.05, 0.1)

#----------------------------------------------------------
# 2) Run the per-file pipeline ONCE per file, caching intermediates
#    (scored / map / el / qtn_ld_table / thr) so any downstream
#    C-score analysis is just extraction, no recomputation.
#----------------------------------------------------------
all_results <- run_and_score_all(
  parsed_folder = parsed_folder, p_cols = p_cols,
  th_ldw_grid = th_ldw_grid, r2_grid = r2_grid, lmin_grid = lmin_grid,
  alpha_grid = alpha_grid, r2_grid_scale = "rho",
  return_intermediates = TRUE, pattern = sel_pattern
)

## drop files with no potential outliers (cached as NULL)
dropped <- names(all_results)[vapply(all_results, is.null, logical(1))]
if (length(dropped) > 0)
  message("Dropping ", length(dropped), " file(s) with no potential outliers: ",
          paste(dropped, collapse = ", "))
all_results <- all_results[!vapply(all_results, is.null, logical(1))]

all_scored <- combine_scored_list(all_results)

#----------------------------------------------------------
# 3) Figure 1 (Rplot12): AUC-PR* trajectories -- full LD-filtering
#    vs th_ldw==0 baseline, for EMMAX and LFMM.
#----------------------------------------------------------
auc_results <- list(
  EMMAX   = auc_cummax_PR(all_scored,               p_name = "EMX",  n_perm = 1000, seed = 1),
  LFMM    = auc_cummax_PR(all_scored,               p_name = "LFMM", n_perm = 1000, seed = 1),
  EMMAX_0 = auc_cummax_PR(all_scored[th_ldw == 0],  p_name = "EMX",  n_perm = 1000, seed = 1),
  LFMM_0  = auc_cummax_PR(all_scored[th_ldw == 0],  p_name = "LFMM", n_perm = 1000, seed = 1)
)
fig_auc <- plot_auc_trajectories(auc_results)
fig_auc

#----------------------------------------------------------
# 4) Figure 2 (Rplot13): genome-wide C-score and -log10(q) panels,
#    coloured by LD with the nearest true QTN.
#----------------------------------------------------------
C_genome_EMX  <- compute_C_score_genome_from_results(all_results, p_name = "EMX")
C_genome_LFMM <- compute_C_score_genome_from_results(all_results, p_name = "LFMM")

C_genome_EMX$map[,  log_emx_q  := -log10(p.adjust(emx_p,  "fdr"))]
C_genome_LFMM$map[, log_lfmm_q := -log10(p.adjust(lfmm_p, "fdr"))]

p1 <- plot_C_score_genome(C_genome_EMX,  title = "EMMAX | full LD-filtering (C-score)",
                          color_by = "max_LD_with_QTN", nrow_facets = 1, show_auc_in_strip = FALSE)
p2 <- plot_C_score_genome(C_genome_LFMM, title = "LFMM | full LD-filtering (C-score)",
                          color_by = "max_LD_with_QTN", nrow_facets = 1, show_auc_in_strip = FALSE)
p3 <- plot_C_score_genome(C_genome_EMX,  title = "EMMAX | q", value_col = "log_emx_q",
                          value_label = "-log10(q)", color_by = "max_LD_with_QTN",
                          nrow_facets = 1, show_auc_in_strip = FALSE) +
  geom_hline(yintercept = 1.31, linetype = 2)
p4 <- plot_C_score_genome(C_genome_LFMM, title = "LFMM | q", value_col = "log_lfmm_q",
                          value_label = "-log10(q)", color_by = "max_LD_with_QTN",
                          nrow_facets = 1, show_auc_in_strip = FALSE) +
  geom_hline(yintercept = 1.31, linetype = 2)

fig_genome <- p1 / p2 / p3 / p4
fig_genome
