######################################################
## Simulated-data OR pipeline -- example figures
##
## Companion driver to analyse_sim2.R, on the same engine. Produces:
##   1. PR heatmap over (rho x th_ldw), faceted by l_min x method
##   2. Manhattan genome plots (-log10(q) and ld_w) per chromosome
##   3. genome-wide C-score panels
##
## Input: per-file parsed simulation objects (GTs, map, env, LD_decay,
## ld_ws) from Parse_sim_data.R. Run from the repository root.
######################################################

library(LDscnR)
library(data.table)
library(igraph)
library(parallel)
library(ggplot2)
library(wesanderson)
library(PRROC)
library(patchwork)

source("R/define_ORs_functions.R")
source("R/Outlier_regions_simulation.R")

#----------------------------------------------------------
# Configuration
#----------------------------------------------------------
parsed_folder <- "/Volumes/Nemo/Nemo_sim/parsed_sim_data2"

## one replicate genome (one sim/bgs/V/c/env, chr1..chr10) for the
## per-genome Manhattan / C-score figures
genome_pattern <- "^adapt_bgs_chr[0-9]+_V0.5_c2_env2\\.rds$"

r2_grid_rho <- c(0.25, 0.5, 0.75, 0.9, 0.95)   # interpreted as rho, not raw r^2
th_ldw_grid <- c(0, 0.25, 0.5, 0.75, 0.9, 0.95)
lmin_grid   <- c(1, 2, 4, 8, 16)
alpha_grid  <- c(0.1, 0.05, 0.005)

#----------------------------------------------------------
# 1) PR heatmap across the LD-filtering grid, pooled by V
#----------------------------------------------------------
## span V0.5/V1/V2 (one bgs/c/env, chr1..chr10) so PR can be grouped by V
heatmap_pattern <- "^adapt_bgs_chr[0-9]+_V[0-9.]+_c2_env2\\.rds$"

all_results <- run_and_score_all(
  parsed_folder = parsed_folder, p_cols = c(EMX = "emx_p", LFMM = "lfmm_p"),
  th_ldw_grid = th_ldw_grid, r2_grid = r2_grid_rho, lmin_grid = lmin_grid,
  r2_grid_scale = "rho", bp_th_cluster = 1e6, pattern = heatmap_pattern
)
pooled_by_V <- aggregate_PR(all_results, p_names = c("EMX", "LFMM"), group_vars = "V")

## the rho-scaled clustering threshold column was historically named
## r2_grid_rho and is now r2_th -- support either
r2col <- if ("r2_grid_rho" %in% names(pooled_by_V)) "r2_grid_rho" else "r2_th"

dt <- pooled_by_V[get(r2col) == 0.25,
                  .(Mean_PR = mean(PR)), by = .(rho = factor(rho), th_ldw, l_min, method)]
dt[, th_ldw_f := factor(round(th_ldw, 3), levels = round(th_ldw_grid, 3), ordered = TRUE)]
dt[, rho_f    := factor(rho, levels = unique(rho), ordered = TRUE)]
dt[, l_min    := factor(paste0("l_min=", l_min),
                        levels = paste0("l_min=", lmin_grid), ordered = TRUE)]

fig_heatmap <- ggplot(dt, aes(x = rho_f, y = th_ldw_f, fill = Mean_PR)) +
  geom_tile(color = "white", linewidth = 0.15) +
  facet_grid(l_min ~ method) +
  scale_fill_gradientn(colors = wes_palette("Zissou1", 100, type = "continuous"),
                       name = "Mean PR") +
  labs(x = expression("LD window size relative to LD decay (" * rho * ")"),
       y = "Local LD filtering threshold") +
  theme_minimal(base_size = 18) +
  theme(legend.position = "right", panel.grid = element_blank(),
        strip.text = element_text(face = "bold"), axis.title = element_text(face = "bold"))
fig_heatmap

#----------------------------------------------------------
# 2) Manhattan genome plots (per-genome: one file per chromosome)
#----------------------------------------------------------
files_genome <- list.files(parsed_folder, pattern = genome_pattern, full.names = TRUE)

## NOTE: sel_rho must be a real ld_ws column name. LDscnR::compute_ld_w()
## names its columns "rho_<value>", so use "rho_0.95", not "0.95".
mk_lfmm_manhattan <- function() {
  sg0 <- score_OR_genome(files_genome, p_cols = c(LFMM = "lfmm_p"), p_name = "LFMM",
                         sel_rho = "rho_0.95", sel_th_ldw = 0, sel_r2_th = 0.3,
                         sel_l_min = 1, bp_th_cluster = 1e6)
  a <- plot_OR_manhattan_genome_cached(sg0, value_col = "fdr_q",
         value_label = expression(-log[10](italic(q))), p_col = "lfmm_p")
  b <- plot_OR_manhattan_genome_cached(sg0, value_col = "ld_w",
         value_label = "ld_w | rho=0.95", color_by = "max_LD_with_QTN")
  a / b
}
fig_manhattan_lfmm <- mk_lfmm_manhattan()
fig_manhattan_lfmm

#----------------------------------------------------------
# 3) Genome-wide C-score
#----------------------------------------------------------
fetched <- fetch_C_score_data_all_files(
  file_paths = files_genome, p_cols = c(LFMM = "lfmm_p"),
  th_ldw_grid = th_ldw_grid, r2_grid = r2_grid_rho, lmin_grid = lmin_grid,
  alpha_grid = alpha_grid, r2_grid_scale = "rho"
)
C_score_genome_all <- compute_C_score_genome_from_data(fetched, p_name = "LFMM")

fig_cscore <- plot_C_score_genome(C_score_genome_all,
                value_col = "C_score", value_label = "C-score | LFMM",
                color_by = "max_LD_with_QTN", nrow_facets = 1, show_auc_in_strip = FALSE)
fig_cscore
