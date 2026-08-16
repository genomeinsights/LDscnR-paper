## module_sim/R/26_lmin_truthfree.R
## TRUTH-FREE l_min diagnostic, validated against truth. The structured null estimates
## the FP region count at each (tau_C, l_min); so from the bundle alone:
##   est_TP   = n_obs - n_null      (observed regions minus null-expected)
##   est_prec = est_TP / n_obs      (= 1 - FDR)
## Both need NO ground truth. We classify the SAME observed regions against the sim
## truth (cache qtab) and check est_TP ~ real_TP and est_prec ~ real_precision, and
## whether the null-based objective picks the same optimal l_min as the true PR-AUC.
## Run: Rscript module_sim/R/26_lmin_truthfree.R [V c env]

source("module_sim/R/21_estimate.R")
suppressMessages(library(ggplot2))
a   <- commandArgs(trailingOnly = TRUE)
V   <- if (length(a) >= 1) a[1] else "2"
cc  <- if (length(a) >= 2) a[2] else "1"
env <- if (length(a) >= 3) a[3] else "1"
TAU <- seq(0.05, 0.95, by = 0.05); LMINS <- c(1L, 2L, 5L, 10L, 20L)

bd <- readRDS(file.path(dir_data, sprintf("null_bundle_V%s_c%s_env%s.rds", V, cc, env)))
ca <- load_cache(V, cc, env)
n_reg <- function(Csp, tau, lmin) { mk <- names(Csp)[Csp >= tau]
  if (!length(mk)) return(list(k = 0L, reg = list())); reg <- cluster_from_cache(mk, bd$edges)
  reg <- reg[lengths(reg) >= lmin]; list(k = length(reg), reg = reg) }

res <- rbindlist(lapply(LMINS, function(lm) rbindlist(lapply(TAU, function(t) {
  ob <- n_reg(bd$C_obs, t, lm); n_obs <- ob$k
  n_null <- mean(vapply(bd$C_surr, function(cs) n_reg(cs, t, lm)$k, numeric(1)))
  est_TP <- max(0, n_obs - n_null); est_prec <- if (n_obs > 0) est_TP / n_obs else NA_real_
  ## truth: classify the SAME observed regions against the sim QTN truth
  ev <- evaluate_ORs_qtn(ob$reg, ca$map, ca$qtab, ca$th$r2min, ca$th$dmax)
  data.table(l_min = lm, tau = t, n_obs = n_obs, n_null = round(n_null, 2),
             est_TP = est_TP, est_prec = est_prec,
             real_TP = ev$TP, real_prec = if (is.na(ev$Precision)) NA_real_ else ev$Precision,
             real_recall = ev$Recall) }))))

ok <- res[n_obs > 0 & is.finite(est_prec) & is.finite(real_prec)]
cat(sprintf("=== est (truth-free, from null) vs real (truth), V%s_c%s_env%s ===\n", V, cc, env))
cat(sprintf("cor(est_TP, real_TP)       = %.3f\n", cor(ok$est_TP, ok$real_TP)))
cat(sprintf("cor(est_prec, real_prec)   = %.3f\n", cor(ok$est_prec, ok$real_prec)))

## per-l_min objectives: truth PR-AUC vs null-based yield (max est_TP at est FDR<=0.10)
by_l <- res[, .(
  truth_PR_AUC = round(pr_auc(real_recall, ifelse(is.na(real_prec), 0, real_prec)), 4),
  null_yield_FDR10 = max(c(0, est_TP[ (1 - est_prec) <= 0.10 ])),   # truth-free
  truth_TP_at_bestF = max(real_TP)
), by = l_min]
cat("\n=== optimal l_min: truth PR-AUC vs null-based yield (truth-free) ===\n"); print(by_l)
cat(sprintf("\ntruth-optimal l_min (max PR-AUC) = %d | null-optimal l_min (max yield@FDR<=0.10) = %d\n",
            by_l$l_min[which.max(by_l$truth_PR_AUC)], by_l$l_min[which.max(by_l$null_yield_FDR10)]))

saveRDS(list(res = res, by_l = by_l), file.path(dir_data, sprintf("lmin_truthfree_V%s_c%s_env%s.rds", V, cc, env)))
p <- ggplot(ok, aes(real_TP, est_TP, color = factor(l_min))) + geom_abline(linetype = 2, color = "grey") +
  geom_point(size = 1.6) + scale_color_viridis_d(end = 0.9, name = "l_min") +
  labs(title = sprintf("V%s_c%s_env%s: null-estimated TP vs true TP (truth-free l_min diagnostic)", V, cc, env),
       subtitle = "est_TP = observed - null regions; points on y=x => null predicts TP without truth",
       x = "real TP (from QTN truth)", y = "est TP (observed - null)") + theme_bw(base_size = 11)
ggsave(file.path(dir_fig, sprintf("lmin_truthfree_V%s_c%s_env%s.png", V, cc, env)), p, width = 7, height = 5.5, dpi = 150)
cat("wrote lmin_truthfree figure\n")
