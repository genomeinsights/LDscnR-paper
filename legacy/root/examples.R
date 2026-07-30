
#----------------------------------------------------------
# examples
#----------------------------------------------------------

r2_grid_rho <- c(0.25,0.5,0.75, 0.9, 0.95)   # now interpreted as rho, not raw r^2
th_ldw_grid <- c(0, 0.25,0.5, 0.75, 0.9, 0.95)
lmin_grid   <- c(1,2, 4, 8,16)
#r2_grid     <- seq(0.2, 0.9, by = 0.1)

p_cols <- c(EMX = "emx_p",LFMM = "lfmm_p")
all_results_1e6 <- all_results <- run_and_score_all(
  parsed_folder = "./parsed_sim_data", p_cols = p_cols,
  th_ldw_grid = th_ldw_grid, r2_grid = r2_grid_rho, lmin_grid = lmin_grid,
  r2_grid_scale = "rho",bp_th_cluster = 1e6
)

pooled_by_V <- aggregate_PR(all_results, p_names = names(p_cols), group_vars = "V")

pooled_by_V[l_min==8,.(V, rho, th_ldw, r2_grid_rho,TP, FP, FN, Precision, Recall, PR,method)]


ggplot(pooled_by_V[ ],aes(factor(th_ldw),PR,fill=method)) +
  geom_boxplot() +
  theme_bw()

dt <- pooled_by_V[r2_grid_rho==0.25 ,.("Mean_PR"=mean(PR)),by=.(rho=factor(rho),th_ldw,l_min,method)]
#dt <- rbind(tmp1.2)

dt[, th_ldw_f := factor(
  round(th_ldw,3),
  levels = round(th_ldw_grid,3),
  ordered = TRUE
)]

dt[, rho_f := factor(
  rho,
  levels = unique(rho),
  ordered = TRUE
)]
dt[,l_min:=paste0("l_min=",l_min)]
dt[, l_min := factor(
  l_min,
  levels = unique(l_min),
  ordered = TRUE
)]

p1 <- ggplot(dt, aes(
  x = rho_f,
  y = th_ldw_f,
  fill = Mean_PR
)) +
  geom_tile(color = "white", linewidth = 0.15) +
  facet_grid(l_min ~ method) +
  scale_fill_gradientn(
    colors = wes_palette("Zissou1", 100, type = "continuous"),
    name = "Mean PR"
  ) +
  scale_x_discrete(
    breaks = factor(seq(0,0.9,by=0.1)),
    labels = as.character(seq(0,0.9,by=0.1))
  ) +
  # scale_y_discrete(
  #   breaks = round(th_ldw_grid,3)[c(1,4,8,12,16,20)],
  #   labels = paste(as.character(round(th_ldw_grid,2)[c(1,4,8,12,16,20)]*100),"%")
  # ) +
  labs(
    x = expression("LD window size relative to LD decay (" * rho * ")"),
    y = "Local LD filtering threshold"#,
    #title = "Number of outlier regions across LD-filtering parameters"
  ) +
  # scale_x_continuous(
  #   breaks = seq(0.2, 2, by = 0.1)
  # ) +
  # scale_y_discrete(
  #   breaks = th_ldw_grid[c(1, 6, 11, 16, 20)],
  #   labels = scales::number(th_ldw_grid[c(1, 6, 11, 16, 20)], accuracy = 0.01)
  # ) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.margin = margin(1, 1, 1, 1)
  )
p1

#+
#  facet_grid(method~.)
#----------------------------------------------------------
# Manhattan plot
#----------------------------------------------------------


files_v2_env1 <- list.files("./parsed_sim_data", pattern = "_V2_c1_env1\\.rds$", full.names = TRUE)
#map

### LFMM
p_cols <- c(LFMM = "lfmm_p")
# Expensive part -- run once per grid point you care about
scored_genome <- score_OR_genome(
  file_paths = files_v2_env1, p_cols = p_cols, p_name = "LFMM",
  sel_rho = "0.95", sel_th_ldw = 0, sel_r2_th = 0.3, sel_l_min = 1,bp_th_cluster = 1e6
)
p1 <- plot_OR_manhattan_genome_cached(
  scored_genome, value_col = "fdr_q", value_label = expression(-log[10](italic(q))),
  p_col = "lfmm_p"
)
#p1


p2 <- plot_OR_manhattan_genome_cached(scored_genome, value_col = "ld_w", value_label = "ld_w | rho=0.95",color_by = "max_LD_with_QTN")

scored_genome <- score_OR_genome(
  file_paths = files_v2_env1, p_cols = p_cols, p_name = "LFMM",
  sel_rho = "0.95", sel_th_ldw = 0.95, sel_r2_th = 0.3, sel_l_min = 1
)
p3 <- plot_OR_manhattan_genome_cached(
  scored_genome, value_col = "fdr_q", value_label = expression(-log[10](italic(q))),
  p_col = "lfmm_p"
)

p_lfmm <- p1 / p2 / p3

### EMX
p_cols <- c(EMX = "emx_p")
# Expensive part -- run once per grid point you care about
scored_genome <- score_OR_genome(
  file_paths = files_v2_env1, p_cols = p_cols, p_name = "EMX",
  sel_rho = "0.95", sel_th_ldw = 0, sel_r2_th = 0.3, sel_l_min = 1,bp_th_cluster = 1e6
)
p1 <- plot_OR_manhattan_genome_cached(
  scored_genome, value_col = "fdr_q", value_label = expression(-log[10](italic(q))),
  p_col = "emx_p"
)
#p1


p2 <- plot_OR_manhattan_genome_cached(scored_genome, value_col = "ld_w", value_label = "ld_w | rho=0.95",color_by = "max_LD_with_QTN")

scored_genome <- score_OR_genome(
  file_paths = files_v2_env1, p_cols = p_cols, p_name = "EMX",
  sel_rho = "0.95", sel_th_ldw = 0.95, sel_r2_th = 0.3, sel_l_min = 1
)
p3 <- plot_OR_manhattan_genome_cached(
  scored_genome, value_col = "fdr_q", value_label = expression(-log[10](italic(q))),
  p_col = "emx_p"
)

p_emx <- p1 / p2 / p3
p_lfmm | p_emx
#----------------------------------------------------------
# C-scores and AUC
#----------------------------------------------------------
r2_grid <- c(0.25,0.5,0.75, 0.9, 0.95)   # now interpreted as rho, not raw r^2
th_ldw_grid <- c(0, 0.25,0.5, 0.75, 0.9, 0.95)
lmin_grid   <- c(1, 2, 4, 8)
alpha_grid <- c(0.1,0.05,0.005)

p_cols <- c(LFMM = "lfmm_p")

files_v2_env1 <- list.files("./parsed_sim_data", pattern = "_V2_c1_env1\\.rds$", full.names = TRUE)


fetched <- fetch_C_score_data_all_files(
  file_paths = files_v2_env1, p_cols = p_cols,
  th_ldw_grid = th_ldw_grid, r2_grid = r2_grid, lmin_grid = lmin_grid,
  alpha_grid = alpha_grid,r2_grid_scale = "rho"
)

C_score_genome_all      <- compute_C_score_genome_from_data(fetched, p_name = "LFMM")
C_score_genome_strict_a <- compute_C_score_genome_from_data(
  fetched, p_name = "LFMM", filter_fun = function(x) x[alpha <= 0.05]
)

C_score_genome_r06      <- compute_C_score_genome_from_data(
  fetched, p_name = "LFMM", fixed_r2_th = 0.6
)

plot_C_score_genome(C_score_genome_all)                       # C_score as both y and color
plot_C_score_genome(C_score_genome_all, value_col = "lfmm_F", value_label = "LFMM F",
                    color_by = "C_score")                     # F-stat y-axis, C_score coloring


## OR, to try several fixed_r2_th / filters without re-running the
## expensive per-file pipeline each time:


fetched <- fetch_C_score_data_all_files(
  file_paths = files_v2_env1, p_cols = p_cols,
  th_ldw_grid = th_ldw_grid, r2_grid = r2_grid_rho, lmin_grid = lmin_grid,
  alpha_grid = c(0.01, 0.05, 0.1, 0.2),r2_grid_scale = "rho"
)

C_score_genome_all      <- compute_C_score_genome_from_data(fetched, p_name = "LFMM")

C_score_genome_strict_a <- compute_C_score_genome_from_data(
  fetched, p_name = "LFMM", filter_fun = function(x) x[alpha <= 0.05]
)
C_score_genome_r06      <- compute_C_score_genome_from_data(
  fetched, p_name = "LFMM", fixed_r2_th = 0.6
)

plot_C_score_genome(C_score_genome_all)                       # C_score as both y and color
plot_C_score_genome(C_score_genome_all, value_col = "lfmm_F", value_label = "LFMM F",
                    color_by = "C_score")                     # F-stat y-axis, C_score coloring


p_cols <- c(LFMM = "lfmm_p")
C_score_genome_lfmm <- compute_C_score_genome(
  file_paths = files_v2_env1, p_cols = p_cols, p_name = "LFMM",
  th_ldw_grid = th_ldw_grid, r2_grid = r2_grid_rho, lmin_grid = lmin_grid,alpha_grid=alpha_grid
  ## fixed_r2_th = NULL by default -- pools r2_th into the C-score/AUC too
)

C_score_genome_lfmm$map[,log_lfmm_p:=-log10(p.adjust(lfmm_p,"fdr"))]

p1 <- plot_C_score_genome(C_score_genome_lfmm, value_col = "log_lfmm_p", value_label = "-log10(q) | LFMM", color_by = "C_score", nrow_facets = 1,show_auc_in_strip = FALSE) + geom_hline(yintercept = 1.31)                    # F-stat y-axis, C_score coloring
p2 <- plot_C_score_genome(C_score_genome_lfmm, value_col = "C_score", value_label = "C-score | LFMM", color_by = "max_LD_with_QTN", nrow_facets = 1,show_auc_in_strip = FALSE)                    # F-stat y-axis, C_score coloring

p_cols <- c(LFMM = "lfmm_p")
scored_genome <- score_OR_genome(
  file_paths = files_v2_env1, p_cols = p_cols, p_name = "LFMM",
  sel_rho = "0.75", sel_th_ldw = 0, sel_r2_th = 0.3, sel_l_min = 1,bp_th_cluster = 1e6
)

p3 <- plot_OR_manhattan_genome_cached(scored_genome, value_col = "ld_w", value_label = "ld_w | rho=0.75",color_by = "max_LD_with_QTN")

p_lfmm <- p1 / p2 / p3



p_cols <- c(EMX = "emx_p")
C_score_genome_emmax <- compute_C_score_genome(
  file_paths = files_v2_env1, p_cols = p_cols, p_name = "EMX",
  th_ldw_grid = th_ldw_grid, r2_grid = r2_grid_rho, lmin_grid = lmin_grid,alpha_grid=alpha_grid
  ## fixed_r2_th = NULL by default -- pools r2_th into the C-score/AUC too
)
C_score_genome_emmax$map[,log_emx_p:=-log10(p.adjust(emx_p,"fdr"))]
p1 <- plot_C_score_genome(C_score_genome_emmax, value_col = "log_emx_p", value_label = "-log10(q) | EMMAX", color_by = "C_score", nrow_facets = 1,show_auc_in_strip = FALSE) + geom_hline(yintercept = 1.31)                    # F-stat y-axis, C_score coloring
p2 <- plot_C_score_genome(C_score_genome_emmax, value_col = "C_score", value_label = "C-score | EMMAX", color_by = "max_LD_with_QTN", nrow_facets = 1,show_auc_in_strip = FALSE)                    # F-stat y-axis, C_score coloring

scored_genome <- score_OR_genome(
  file_paths = files_v2_env1, p_cols = p_cols, p_name = "EMX",
  sel_rho = "0.75", sel_th_ldw = 0, sel_r2_th = 0.3, sel_l_min = 1,bp_th_cluster = 1e6
)

p3 <- plot_OR_manhattan_genome_cached(scored_genome, value_col = "ld_w", value_label = "ld_w | rho=0.75",color_by = "max_LD_with_QTN")


p_emx <- p1 / p2 / p3

p_lfmm | p_emx

#auc_by_Cscore <- auc_cummax_PR_Cscore_sweep(C_sweep, n_perm = 200, seed = 1)
