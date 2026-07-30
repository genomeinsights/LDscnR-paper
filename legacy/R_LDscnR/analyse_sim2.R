######################################################
## Example: C-score vs raw-alpha comparison across ALL simulated
## data, running the expensive per-file analysis ONLY ONCE per
## file (via run_full_analysis_all_files()), then deriving both
## the raw-alpha baseline and the C-score curve from that single
## pass -- grouped by V so you can compare across selection
## intensities without re-running anything.
##
## Requires outlier_regions_generic.R and outlier_regions_simulation.R
## to be sourced first.
######################################################

library(PRROC)

source("./R_LDscnR/define_ORs_functions.R")
source("./R_LDscnR/PR_AUC_functions.R")

#----------------------------------------------------------
# 1) Setup
#----------------------------------------------------------
p_cols <- c(EMX = "emx_p", LFMM = "lfmm_p")

r2_grid <- c(0.9)   # now interpreted as rho, not raw r^2
th_ldw_grid <- c(0, 0.25,0.5, 0.75, 0.9, 0.95, 0.99)
lmin_grid   <- c(1, 2, 4, 8, 16)
alpha_grid   <- 0.05
alpha_grid   <- c(0.0001,0.001, 0.01, 0.05, 0.1)
C_score_grid <- seq(0, 1, by = 0.05)

sel_r2_th <- 0.9   # FIXED for both curves -- not integrated, not swept
sel_l_min <- 1     # FIXED for both curves

gc()
#----------------------------------------------------------
# 2) Run the expensive per-file pipeline EXACTLY ONCE per file,
#    caching everything (scored/map/el/qtn_ld_table/thr) needed for
#    ANY later C-score analysis -- no C-score grid committed yet.
#----------------------------------------------------------

all_results <- run_and_score_all(
  parsed_folder = "./parsed_sim_data", p_cols = p_cols,
  th_ldw_grid = th_ldw_grid, r2_grid = r2_grid, lmin_grid = lmin_grid,
  alpha_grid = alpha_grid, return_intermediates = TRUE,r2_grid_scale = "rho",pattern = "c1"
)

all_scored <- combine_scored_list(all_results)
pooled_PR <- aggregate_PR(all_scored, p_names = names(p_cols))

ggplot(pooled_PR[,.(mean_PR=mean(PR,na.rm = TRUE)),by=.(r2_th)], aes(r2_th,mean_PR)) +
  geom_line()

#----------------------------------------------------------
# 3) Cheap: stack just the raw grid tables (all_scored) for the
#    alpha baseline. No recomputation -- just extraction + rbind.
#----------------------------------------------------------
all_scored <- combine_scored_list(all_results)

auc_results <- list(
  EMMAX = auc_cummax_PR(all_scored, p_name = "EMX",  n_perm = 1000, seed = 1),
  LFMM  = auc_cummax_PR(all_scored, p_name = "LFMM", n_perm = 1000, seed = 1),
  EMMAX_0 = auc_cummax_PR(all_scored[th_ldw==0], p_name = "EMX",  n_perm = 1000, seed = 1),
  LFMM_0  = auc_cummax_PR(all_scored[th_ldw==0], p_name = "LFMM", n_perm = 1000, seed = 1)
)


plot_auc_trajectories(auc_results)

C_labels_LFMM     <- get_marker_Cscore_labels_from_results(all_results, p_name = "LFMM",filter_fun = function(x) x[l_min>4])
alpha_labels_LFMM <- get_marker_alpha_labels_from_results(all_results, p_col = "lfmm_p")

pr_C     <- auc_pr_PRROC(C_labels_LFMM, curve = TRUE)
pr_alpha <- auc_pr_PRROC(alpha_labels_LFMM, curve = TRUE)

pr_C$auc.davis.goadrich       # preferred over auc.integral -- see note above
pr_alpha$auc.davis.goadrich

## drop any files with zero potential outliers (return_intermediates
## caches NULL for those -- compute_C_score_genome_from_results() needs
## every entry to have a real map to plot, so these can't be included)
dropped <- names(all_results)[vapply(all_results, is.null, logical(1))]
if (length(dropped) > 0) {
  message("Dropping ", length(dropped), " file(s) with no potential outliers: ",
          paste(dropped, collapse = ", "))
}
all_results <- all_results[!vapply(all_results, is.null, logical(1))]

#----------------------------------------------------------
# 3) Build the genome-wide C-score display, cheaply, for each of the
#    4 method x filtering combinations -- no recomputation, just
#    C-score + relabeling from the cache built in step 2.
#----------------------------------------------------------

C_genome_EMX_full   <- compute_C_score_genome_from_results(all_results, p_name = "EMX")
C_genome_LFMM_full  <- compute_C_score_genome_from_results(all_results, p_name = "LFMM")

C_genome_EMX_0  <- compute_C_score_genome_from_results(
  all_results, p_name = "EMX",  filter_fun = function(x) x[th_ldw == 0]
)
C_genome_LFMM_0 <- compute_C_score_genome_from_results(
  all_results, p_name = "LFMM", filter_fun = function(x) x[th_ldw == 0]
)

#----------------------------------------------------------
# 4) Plot -- C-score as both y and color (AUC shown per chromosome
#    in the facet strip labels, since that's a per-file quantity)
#----------------------------------------------------------
C_genome_EMX_full$map[,log_emx_q := -log10(p.adjust(emx_p,"fdr"))]
C_genome_LFMM_full$map[,log_lfmm_q := -log10(p.adjust(lfmm_p,"fdr"))]
C_genome_EMX_0$map[,log_emx_q := -log10(p.adjust(emx_p,"fdr"))]
C_genome_LFMM_0$map[,log_lfmm_q := -log10(p.adjust(lfmm_p,"fdr"))]
p1 <- plot_C_score_genome(C_genome_EMX_full,  title = "EMMAX | full LD-filtering (C-score)",color_by = "max_LD_with_QTN",nrow_facets = 1,show_auc_in_strip = FALSE)
p2 <- plot_C_score_genome(C_genome_LFMM_full, title = "LFMM | full LD-filtering (C-score)",color_by = "max_LD_with_QTN",nrow_facets = 1,show_auc_in_strip = FALSE)
p3 <- plot_C_score_genome(C_genome_EMX_full,  title = "EMMAX | q",color_by = "max_LD_with_QTN",nrow_facets = 1,show_auc_in_strip = FALSE,value_col = "log_emx_q",value_label = "-log10(q)") + geom_hline(yintercept = 1.31,linetype=2)
p4 <- plot_C_score_genome(C_genome_LFMM_full, title = "LFMM | q",color_by = "max_LD_with_QTN",nrow_facets = 1,show_auc_in_strip = FALSE,value_col = "log_lfmm_q",value_label = "-log10(q)") + geom_hline(yintercept = 1.31,linetype=2)
p1 / p2 / p3 / p4

#p3 <- plot_C_score_genome(C_genome_EMX_0,  title = "EMMAX | th_ldw=0 only (C-score)",color_by = "max_LD_with_QTN",nrow_facets = 1,show_auc_in_strip = FALSE)
#p4 <-plot_C_score_genome(C_genome_LFMM_0, title = "LFMM | th_ldw=0 only (C-score)",color_by = "max_LD_with_QTN",nrow_facets = 1,show_auc_in_strip = FALSE)


#plot(pr_C)                    # PRROC objects have a built-in plot method

marker_C_curves <- run_marker_C_score_curve_from_results(all_results, p_name = "LFMM", C_score_grid = C_score_grid)
pooled_C_curve  <- marker_C_curves[, .(TP = sum(TP), FP = sum(FP), FN = sum(FN)), by = C_score_threshold]

marker_alpha_curves <- extract_marker_alpha_baseline_from_results(all_results, p_name = "LFMM")
pooled_alpha_curve  <- marker_alpha_curves[, .(TP = sum(TP), FP = sum(FP), FN = sum(FN)), by = alpha]

pr_curve_auc(pooled_C_curve)$AUC
pr_curve_auc(pooled_alpha_curve)$AUC

C_curve_by_V <- pr_curve_auc(pooled_C_curve)$curve
alpha_by_V <- pr_curve_auc(pooled_alpha_curve)$curve
C_curve_by_V[, group := paste("C-score")]
alpha_by_V[,  group := paste("raw alpha")]

combined <- rbind(
  C_curve_by_V[, .(TP, FP, FN, group)],
  alpha_by_V[,   .(TP, FP, FN,group)],
  fill = TRUE
)

p1 <- plot_pr_curve(combined[V == "V1"], group_col = "group",
                    title = "C-score vs raw alpha (Recall/Precision), V2")

alpha_by_V <- aggregate_PR(all_scored, p_names = names(p_cols), group_vars = "V")[
  th_ldw == 0 & r2_th == sel_r2_th & l_min == sel_l_min
]

alpha_by_V[,plot(alpha,PR,type="l",col=as.numeric(as.factor(method)))]

C_sweep_all <- run_C_score_sweep_from_results(
  all_results, p_names = names(p_cols),
  r2_th_grid = sel_r2_th, C_score_grid = C_score_grid, l_min = sel_l_min,filter_fun = function(x){x[alpha==0.2]}
)

C_sweep_all[,hist(TP)]
C_curve_by_V <- C_sweep_all[, .(TP = sum(TP), FP = sum(FP), FN = sum(FN)), by = .(C_score_threshold, method)]
C_curve_by_V[,hist(TP)]

C_curve_by_V[,PR:=ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_)*ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)]
C_curve_by_V[,plot(C_score_threshold,PR,type="l",col=as.numeric(as.factor(method)))]

#----------------------------------------------------------
# 5) Standard monotonic PR-curve AUC, for both criteria, both
#    methods, EACH V level
#----------------------------------------------------------
V_levels <- unique(all_scored$V)

auc_by_V <- rbindlist(lapply(V_levels, function(v) {
  rbindlist(lapply(names(p_cols), function(nm) {
    C_res     <- pr_curve_auc(C_curve_by_V[method == nm & V == v])
    alpha_res <- pr_curve_auc(alpha_by_V[method == nm & V == v])
    data.table(
      V = v, method = nm,
      AUC_Cscore = C_res$AUC, n_Cscore_thresholds = C_res$n_thresholds,
      AUC_alpha  = alpha_res$AUC, n_alpha_thresholds = alpha_res$n_thresholds
    )
  }))
}))
auc_by_V[order(V, method)]

auc_long <- melt(auc_by_V, id.vars = c("V", "method"),
                 measure.vars = c("AUC_Cscore", "AUC_alpha"),
                 variable.name = "criterion", value.name = "AUC")
auc_long[, criterion := fifelse(criterion == "AUC_Cscore", "C-score", "raw alpha")]

ggplot(auc_long, aes(x = V, y = AUC, color = criterion, linetype = criterion,group=interaction(criterion))) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  facet_wrap(~ method) +
  labs(x = "V (selection variance -- higher = harder detection)",
       y = "AUC (standard PR-curve)", color = NULL, linetype = NULL,
       title = sprintf("C-score vs raw alpha across V (r2_th=%.2f, l_min=%d)", sel_r2_th, sel_l_min)) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")


C_curve_by_V[, group := paste(method, "C-score")]
alpha_by_V[,  group := paste(method, "raw alpha")]

combined <- rbind(
  C_curve_by_V[, .(TP, FP, FN, method, V, group)],
  alpha_by_V[,   .(TP, FP, FN, method, V, group)],
  fill = TRUE
)

p1 <- plot_pr_curve(combined[V == "V2"], group_col = "group",
                    title = "C-score vs raw alpha (Recall/Precision), V2")
#----------------------------------------------------------
# 2) Run the expensive per-file pipeline EXACTLY ONCE per file.
#    Derives both all_scored (raw grid) and C_sweep_all (C-score
#    sweep at sel_r2_th) from the same run_and_score_one_sim_file()
#    call -- no duplicated computation.
#----------------------------------------------------------
full_res <- run_full_analysis_all_files(
  file_paths = all_files, p_cols = p_cols,
  th_ldw_grid = th_ldw_grid, r2_grid = r2_grid, lmin_grid = lmin_grid,
  r2_th_grid_Cscore = sel_r2_th, C_score_grid = C_score_grid,
  l_min_cluster = sel_l_min, alpha_grid = alpha_grid
)
all_scored  <- full_res$all_scored    # raw grid, all files, all V/c/env

pooled_PR <- aggregate_PR(all_scored, p_names = names(p_cols))

ggplot(pooled_PR,aes(factor(r2_th),PR,fill=method)) +
  geom_boxplot()
#all_scored[,r2_th := r2_grid_rho]


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

p1 <- ggplot(dt, aes(
  x = rho_f,
  y = th_ldw_f,
  fill = `Number of outliers`
)) +
  geom_tile(color = "white", linewidth = 0.15) +
  facet_grid(. ~ data) +
  scale_fill_gradientn(
    colors = wes_palette("Zissou1", 100, type = "continuous"),
    name = "Outlier \nregions"
  ) +
  scale_x_discrete(
    breaks = factor(seq(0,0.9,by=0.1)),
    labels = as.character(seq(0,0.9,by=0.1))
  ) +
  scale_y_discrete(
    breaks = round(th_ldw_grid,3)[c(1,4,8,12,16,20)],
    labels = paste(as.character(round(th_ldw_grid,2)[c(1,4,8,12,16,20)]*100),"%")
  ) +
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

auc_results <- list(
  EMMAX = auc_cummax_PR(all_scored, p_name = "EMX",  n_perm = 1000, seed = 1),
  LFMM  = auc_cummax_PR(all_scored, p_name = "LFMM", n_perm = 1000, seed = 1),
  EMMAX_0 = auc_cummax_PR(all_scored[th_ldw==0], p_name = "EMX",  n_perm = 1000, seed = 1),
  LFMM_0  = auc_cummax_PR(all_scored[th_ldw==0], p_name = "LFMM", n_perm = 1000, seed = 1)
)

plot_auc_trajectories(auc_results)


C_sweep_all <- full_res$C_sweep_all   # C-score sweep at sel_r2_th, all files, all V/c/env

#----------------------------------------------------------
# 3) Raw-alpha baseline, pooled WITHIN each V level (across that
#    V's replicate files), same fixed r2_th/l_min, th_ldw = 0.
#    aggregate_PR()'s group_vars pools by V (and method) directly.
#----------------------------------------------------------
alpha_by_V <- aggregate_PR(all_scored, p_names = names(p_cols), group_vars = "V")[
  th_ldw == 0 & r2_th == sel_r2_th & l_min == sel_l_min & rho==0.05
]

#----------------------------------------------------------
# 4) C-score curve, pooled WITHIN each V level (r2_th is already
#    fixed at sel_r2_th, so only V and C_score_threshold remain
#    to group by when pooling TP/FP/FN across that V's files)
#----------------------------------------------------------
# C_curve_by_V <- C_sweep_all[, .(TP = sum(TP), FP = sum(FP), FN = sum(FN)), by = .(C_score_threshold, method, V)]
# dt=C_curve_by_V
# TP_col = "TP"
# FP_col = "FP"
# FN_col = "FN"
#
# TP <- dt[[TP_col]]
# FP <- dt[[FP_col]]
# FN <- dt[[FN_col]]
#
# Precision <- ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_)
# Recall    <- ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)
#
# curve <- data.table(Recall = Recall, Precision = Precision)
# dt_C <- cbind(dt,curve,PR=curve$Recall*curve$Precision)
# dt_C[,stat:="C"]
# dt_C[,th:=C_score_threshold]
#
# alpha_by_V[,stat:="-log10(alpha)"]
# alpha_by_V[,th:=-log10(alpha)]
#
# dt <- rbind(dt_C,alpha_by_V,fill=TRUE)
#
# ggplot(dt, aes(th,PR, col=method)) +
#   geom_line() +
#   facet_wrap(stat~.,scales="free_x") +
#   theme_bw()

C_curve_by_V[, group := paste(method, "C-score")]
alpha_by_V[,  group := paste(method, "raw alpha")]

combined <- rbind(
  C_curve_by_V[, .(TP, FP, FN, method, V, group)],
  alpha_by_V[,   .(TP, FP, FN, method, V, group)],
  fill = TRUE
)

p1 <- plot_pr_curve(combined[V == "V2"], group_col = "group",
              title = "C-score vs raw alpha (Recall/Precision), V2")

## collapse ties in Recall to the highest Precision observed at that Recall
curve <- curve[, .(Precision = max(Precision)), by = Recall]
setorder(curve, Recall)


#----------------------------------------------------------
# 5) Standard monotonic PR-curve AUC, for both criteria, both
#    methods, EACH V level
#----------------------------------------------------------
V_levels <- unique(all_scored$V)

v <- "V2"
nm <- names(p_cols)[1]
auc_by_V <- rbindlist(lapply(V_levels, function(v) {
  rbindlist(lapply(names(p_cols), function(nm) {
    C_res     <- pr_curve_auc(C_curve_by_V[method == nm & V == v])
    alpha_res <- pr_curve_auc(alpha_by_V[method == nm & V == v])
    data.table(
      V = v, method = nm,
      AUC_Cscore = C_res$AUC, n_Cscore_thresholds = C_res$n_thresholds,
      AUC_alpha  = alpha_res$AUC, n_alpha_thresholds = alpha_res$n_thresholds
    )
  }))
}))
#auc_by_V[order(V, method)]

#----------------------------------------------------------
# 6) Plot: AUC (C-score vs raw alpha) as a function of V, one
#    panel/line per method -- the actual "does LD-filtering help
#    more as detection gets harder" comparison, using V directly
#    instead of alpha as the difficulty axis.
#----------------------------------------------------------
auc_long <- melt(auc_by_V, id.vars = c("V", "method"),
                 measure.vars = c("AUC_Cscore", "AUC_alpha"),
                 variable.name = "criterion", value.name = "AUC")
auc_long[, criterion := fifelse(criterion == "AUC_Cscore", "C-score", "raw alpha")]

p2 <- ggplot(auc_long, aes(x = V, y = AUC, color = criterion, linetype = criterion)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  facet_wrap(~ method) +
  labs(x = "V (selection variance -- higher = harder detection)",
       y = "AUC (standard PR-curve)", color = NULL, linetype = NULL,
       title = sprintf("C-score vs raw alpha across V (r2_th=%.2f, l_min=%d)", sel_r2_th, sel_l_min)) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")
#library(patchwork)
p1 | p2
#----------------------------------------------------------
# Note: this single run_full_analysis_all_files() call now covers
# BOTH the C-score-vs-alpha comparison above AND the raw AUC-PR*
# comparison from the earlier example (auc_cummax_PR() on
# all_scored, per V via group_vars, etc.) -- no need to re-run
# the per-file pipeline again for that either.



# parse_raw_data(gz_files,
#                outfolder,
#                min_maf,
#                side=48,
#                keep_inds,
#                cores = 4)
#
# if(FALSE){
#   map[,indx:=.I]
#   par(mfcol=c(2,1))
#   rho = "0.9"
#   map[, q_ldw := ecdf(ld_ws[,rho])(ld_ws[,rho])]
#   map[, ld_w := ld_ws[,rho]]
#   #map[,plot(ld_w)]
#   map[, true_pos := rho_d <=0.99 & ld_rel > 0.25]
#   map[chr_type=="ntrl",true_pos:=FALSE]
#   map[, SNP_class := fifelse(true_pos, "SNPs near causal loci", "Neutral SNPs")]
#
#
#   #map[,plot(MAF,ld_w)]
#   plot(map$indx,-log10(p.adjust(map$lfmm_p,"fdr")),pch=ifelse(map$type=="QTN",3,20),cex=ifelse(map$type=="QTN",3,1),col=ifelse(map$true_pos,"firebrick4","grey"),
#        main="all SNPs | LFMM")
#   abline(h=1.3)
#   qt = 0.8
#   keep <-  map[,which(q_ldw>qt)]
#   #map[keep & type=="QTN"]
#   #hist(-log10(p.adjust(map$lfmm_p[keep],"fdr")))
#   #hist(-log10(p.adjust(map$lfmm_p,"fdr"))[keep])
#
#   #map[,plot(ld_ws,-log10(lfmm_p))]
#
#   plot(map$indx[keep],-log10(p.adjust(map$lfmm_p[keep],"fdr")),
#        pch=ifelse(map[keep]$type=="QTN",3,20),
#        cex=ifelse(map[keep]$type=="QTN",3,1),col=ifelse(map[keep]$true_pos,"firebrick4","grey"),
#        main="ld_w>ld_w[Q(95)] | LFMM",xlim=c(0,nrow(map)))
#   abline(h=1.3)
#
#   plot(map$indx,map$ld_w,pch=ifelse(map$type=="QTN",3,20),cex=ifelse(map$type=="QTN",3,1),col=ifelse(map$true_pos,"firebrick4","grey"))
#
#   plot(map$indx,-log10(p.adjust(map$emx_p,"fdr")),pch=ifelse(map$type=="QTN",3,20),cex=ifelse(map$type=="QTN",3,1),col=ifelse(map$true_pos,"firebrick4","grey"),
#        main="all SNPs | EMMAX")
#   abline(h=1.3)
#
#   plot(map$indx[keep],-log10(p.adjust(map$emx_p[keep],"fdr")),
#        pch=ifelse(map[keep]$type=="QTN",3,20),
#        cex=ifelse(map[keep]$type=="QTN",3,1),col=ifelse(map[keep]$true_pos,"firebrick4","grey"),
#        main="ld_w>ld_w[Q(95)] | EMMAX",xlim=c(0,nrow(map)))
#   abline(h=1.3)
#
#   # map[, true_pos := rho_d <=0.99 & ld_rel > 0.25]
#   # map[chr_type=="ntrl",true_pos:=FALSE]
#   # map[, SNP_class := fifelse(true_pos, "SNPs near causal loci", "Neutral SNPs")]
#
#
#   p1 <- ggplot(map,aes(q_ldw,lfmm_F,col=max_LD_with_QTN)) +
#     geom_point(data=map[chr_type=="ntrl"],alpha=0.5,shape = 21, fill="grey40",size = 3.5, col = "black", stroke = 0.3) +
#     #geom_point(data = map[!which(true_pos)],size = 1.5, col = "grey70", alpha = 0.35)+
#     #geom_smooth(se=FALSE) +
#     geom_point(data = map[which(true_pos)],
#                aes(fill = max_LD_with_QTN),
#                shape = 21, size = 3.5, col = "black", stroke = 0.3)+
#     scale_shape_manual(values = c(20,3),name=NULL) +
#     scale_fill_gradientn(
#       colors = wes_palette("Zissou1", 100, type = "continuous"),name = expression("LD with QTN (" * italic(r)^2 * ")")
#     ) +
#     facet_grid(.~SNP_class)+
#     scale_size_identity() +
#     theme_bw(base_size = 22) +
#     theme(strip.background = element_blank(),
#           legend.background = element_blank(),
#           legend.position = "inside",
#           legend.position.inside = c(0.65,0.7))+
#     ylab("LFMM association statistic") +
#     xlab("Local LD rank quantile")
#   p1
#
#   ggsave(p1,filename="the_premise.png",height = 5,width = 12,units = "in",dpi=300)
# }
#
