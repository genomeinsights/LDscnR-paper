######################################################
## Example: pooling the 10 chr1..chr10 replicates when
## estimating PR and AUC-PR*, for one simulated V/env setting.
##
## Requires outlier_regions_generic.R and outlier_regions_simulation.R
## to be sourced first.
######################################################

source("./R_LDscnR/define_ORs_functions.R")
source("./R_LDscnR/PR_AUC_functions.R")

#----------------------------------------------------------
# 1) Setup
#----------------------------------------------------------
p_cols <- c(EMX = "emx_p", LFMM = "lfmm_p")

th_ldw_grid <- c(0, 0.5, 0.75, 0.9, 0.95)
r2_grid     <- seq(0.6, 0.9, by = 0.1)
lmin_grid   <- c(1, 5, 10, 20)
alpha_grid  <- 0.05   # single value here -- see note at the bottom for a real alpha sweep

#----------------------------------------------------------
# 2) Run the per-file pipeline for the 10 chr1..chr10 replicates
#    sharing one V/env combination (e.g. V2, env1). Each file gets
#    its OWN FDR correction/clustering -- nothing is pooled yet at
#    this stage, this just produces one `scored` table per file.
#----------------------------------------------------------
files_v2_env1 <- list.files("./parsed_sim_data", pattern = "_V2_c1_env1\\.rds$", full.names = TRUE)
stopifnot(length(files_v2_env1) == 10)

r2_grid <- c(0.25,0.5,0.75, 0.9, 0.95)   # now interpreted as rho, not raw r^2
th_ldw_grid <- c(0, 0.25,0.5, 0.75, 0.9, 0.95)
lmin_grid   <- c(1, 2, 4, 8)
alpha_grid <- c(0.1,0.05,0.005)
p_cols <- c(LFMM = "lfmm_p",EMX = "emx_p")

files_v2_env1 <- list.files("./parsed_sim_data", pattern = "_V2_c1_env1\\.rds$", full.names = TRUE)
all_scored <- rbindlist(lapply(files_v2_env1, function(f) {
  message("Processing ", basename(f))
  tryCatch(
    run_and_score_one_sim_file(
      f,p_cols = p_cols, th_ldw_grid = th_ldw_grid, r2_grid = r2_grid, lmin_grid = lmin_grid,alpha_grid = alpha_grid,r2_grid_scale = "rho"
    ),
    error = function(e) { message("  FAILED (", basename(f), "): ", conditionMessage(e)); NULL }
  )
}))

all_scored[,r2_th := r2_grid_rho]

#gc()
## all_scored now has 10 files' worth of rows stacked, one row per
## (file, rho, th_ldw, alpha, r2_th, l_min) combination, each with its
## own TP_<nm>/FP_<nm>/FN_<nm> counts computed independently per file.

#----------------------------------------------------------
# 3) Pooled Precision/Recall/PR: sum TP/FP/FN across the 10 files
#    WITHIN each grid point, THEN compute ratios -- not the other way
#    around. This is what avoids the "many 0's" problem: a single
#    chromosome's grid point very often has TP=0 (few/no true
#    positives findable on that one chromosome at that setting),
#    which mechanically forces PR=0 if computed per file first.
#----------------------------------------------------------
pooled_PR <- aggregate_PR(all_scored, p_names = names(p_cols))

#pooled_PR[,hist(TP)]




## one row per (method, rho, th_ldw, r2_th, l_min, alpha), pooled
## across all 10 files -- inspect a slice, e.g.:
pooled_PR[,.(method, rho, th_ldw, TP, FP, FN, Precision, Recall, PR)]

#----------------------------------------------------------
# 4) AUC-PR*: random-search AUC over the FULL (rho x th_ldw x l_min x
#    r2_th x alpha) grid. auc_cummax_PR() pools TP/FP/FN across
#    whatever rows share the same grid point INTERNALLY -- since
#    all_scored has 10 files' rows per grid point, this happens
#    automatically here, using the exact same pooling logic as
#    aggregate_PR() above (so the two are consistent with each other).
#----------------------------------------------------------

auc_results <- list(
  EMMAX = auc_cummax_PR(all_scored, p_name = "EMX",  n_perm = 1000, seed = 1),
  LFMM  = auc_cummax_PR(all_scored, p_name = "LFMM", n_perm = 1000, seed = 1),
  EMMAX_0 = auc_cummax_PR(all_scored[th_ldw==0], p_name = "EMX",  n_perm = 1000, seed = 1),
  LFMM_0  = auc_cummax_PR(all_scored[th_ldw==0], p_name = "LFMM", n_perm = 1000, seed = 1)
)

plot_auc_trajectories(auc_results)


