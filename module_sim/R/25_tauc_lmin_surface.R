## module_sim/R/25_tauc_lmin_surface.R
## tau_C x l_min POWER SURFACE (sim, truth-based, from the cache 20). Tests the claim:
## some FALSE-POSITIVE single SNPs carry high C (structure noise can be locally
## consistent), so raising l_min prunes FP faster than TP -> more TPs per FP. Outputs:
##  (A) PR frontier (recall vs precision) per l_min, sweeping tau_C -- does higher l_min
##      dominate? (B) TP vs FP as tau_C varies, per l_min -- the "TPs per FP" view.
##  (C) PR-AUC by l_min. (D) at a mid tau_C, FP/TP composition of small vs large clusters
##      -- the direct mechanism test. Replicate-averaged over env1-5, alpha-swept C (sim).
## Run: Rscript module_sim/R/25_tauc_lmin_surface.R [V c]

source("module_sim/R/21_estimate.R")
suppressMessages({ library(ggplot2); library(patchwork) })
a  <- commandArgs(trailingOnly = TRUE)
V  <- if (length(a) >= 1) a[1] else "2"
cc <- if (length(a) >= 2) a[2] else "1"; ENV <- as.character(1:5)
ALPHA_SIM <- c(0.001, 0.01, 0.05, 0.1); TAUC <- seq(0.05, 0.95, by = 0.05); LMINS <- c(1L, 2L, 5L, 10L, 20L)
CA <- lapply(ENV, function(e) load_cache(V, cc, e))

## per (env, method, tau, l_min): TP, FP, recall, precision, PR (dedup-neutral)
grid <- rbindlist(lapply(seq_along(CA), function(i) { ca <- CA[[i]]
  rbindlist(lapply(c("emx_p", "lfmm_p"), function(p) {
    Cv <- cscore_count(ca$map[[p]], ca$LDW, ca$grid$RHO, ca$grid$QSTAR, ALPHA_SIM)
    rbindlist(lapply(TAUC, function(t) { reg <- cluster_from_cache(names(Cv)[Cv >= t], ca$edges)
      rbindlist(lapply(LMINS, function(lm) {
        ev <- evaluate_ORs_qtn(reg[lengths(reg) >= lm], ca$map, ca$qtab, ca$th$r2min, ca$th$dmax)
        data.table(env = i, method = sub("_p$", "", p), tau = t, l_min = lm,
                   TP = ev$TP, FP = ev$FP, recall = ev$Recall,
                   precision = if (is.na(ev$Precision)) NA_real_ else ev$Precision,
                   PR = if (is.na(ev$PR)) 0 else ev$PR) })) })) })) }))
agg <- grid[, .(TP = mean(TP), FP = mean(FP), recall = mean(recall),
                precision = mean(precision, na.rm = TRUE), PR = mean(PR)), by = .(method, tau, l_min)]
agg[, l_min := factor(l_min)]

auc <- agg[, .(PR_AUC = round(pr_auc(recall, precision), 4)), by = .(method, l_min)]
cat("=== PR-AUC by l_min (mean over 5 env, alpha-swept C) ===\n"); print(dcast(auc, method ~ l_min, value.var = "PR_AUC"))

## (D) mechanism test: at a mid tau_C, FP/TP composition of small (<5) vs large (>=10) clusters
midtau <- 0.3
comp <- rbindlist(lapply(seq_along(CA), function(i) { ca <- CA[[i]]
  rbindlist(lapply(c("emx_p", "lfmm_p"), function(p) {
    Cv <- cscore_count(ca$map[[p]], ca$LDW, ca$grid$RHO, ca$grid$QSTAR, ALPHA_SIM)
    reg <- cluster_from_cache(names(Cv)[Cv >= midtau], ca$edges); if (!length(reg)) return(data.table())
    dg <- diagnose_OR_classification(reg, ca$map, ca$qtab, ca$th$r2min, ca$th$dmax)
    dg <- as.data.table(dg); dg[, sz := lengths(reg)]
    dg[, .(method = sub("_p$","",p), env = i, bin = ifelse(sz == 1, "singleton", ifelse(sz < 10, "2-9", ">=10")),
           is_TP = is_TP)] })) }))
if (nrow(comp)) { tab <- comp[, .(regions = .N, TP = sum(is_TP), FP = sum(!is_TP)), by = .(method, bin)]
  tab[, FP_rate := round(FP / regions, 3)]
  cat(sprintf("\n=== cluster TP/FP composition at tau_C=%.2f, by cluster size (pooled 5 env) ===\n", midtau))
  print(tab[order(method, bin)]) }

saveRDS(list(agg = agg, auc = auc, comp = if (nrow(comp)) tab else NULL),
        file.path(dir_data, sprintf("tauc_lmin_surface_V%s_c%s.rds", V, cc)))
pA <- ggplot(agg, aes(recall, precision, color = l_min)) + geom_path(linewidth = 0.8) +
  facet_wrap(~method) + scale_color_viridis_d(end = 0.9) + coord_cartesian(xlim = c(0,0.6), ylim = c(0,1)) +
  labs(title = "(A) PR frontier by l_min (sweep tau_C)", x = "Recall", y = "Precision")
pB <- ggplot(agg, aes(FP, TP, color = l_min)) + geom_path(linewidth = 0.8) + geom_point(size = 0.8) +
  facet_wrap(~method) + scale_color_viridis_d(end = 0.9) +
  labs(title = "(B) TP vs FP by l_min (more TP per FP = up/left)", x = "FP regions", y = "TP regions")
pC <- ggplot(auc, aes(l_min, PR_AUC, fill = method)) + geom_col(position = "dodge") +
  scale_fill_manual(values = c(emx = "#1D3557", lfmm = "#E63946")) +
  labs(title = "(C) PR-AUC by l_min", x = "l_min", y = "PR-AUC")
ggsave(file.path(dir_fig, sprintf("tauc_lmin_surface_V%s_c%s.png", V, cc)),
       (pA / pB / pC) + patchwork::plot_annotation(title = sprintf("V%s_c%s: tau_C x l_min power surface (mean 5 env)", V, cc)) &
         theme_bw(base_size = 10), width = 10, height = 12, dpi = 150)
cat("\nwrote tauc_lmin_surface figure\n")
