## module_sim/R/26b_lmin_truthfree_aggregate.R
## REPLICATE-AVERAGED truth-free l_min diagnostic (env1 alone flukes to l_min=1; the
## 5-env truth optimum is l_min=2, see 25). For each env, from the null BUNDLE (18d) +
## cache: est_TP = n_obs - n_null and est_FDR = n_null/n_obs are truth-free; the SAME
## observed regions are scored against the QTN truth. Aggregate over env1-5 and ask:
##  (a) does est_TP track real_TP (pooled cor)?
##  (b) does the null-based objective (max est_TP at est_FDR<=0.10) pick the same l_min
##      as the true PR-AUC?
## Run: Rscript module_sim/R/26b_lmin_truthfree_aggregate.R [V c]

source("module_sim/R/21_estimate.R")
suppressMessages(library(ggplot2))
a  <- commandArgs(trailingOnly = TRUE)
V  <- if (length(a) >= 1) a[1] else "2"
cc <- if (length(a) >= 2) a[2] else "1"; ENV <- as.character(1:5)
TAU <- seq(0.05, 0.95, by = 0.05); LMINS <- c(1L, 2L, 5L, 10L)

grid <- rbindlist(lapply(ENV, function(env) {
  bd <- readRDS(file.path(dir_data, sprintf("null_bundle_V%s_c%s_env%s.rds", V, cc, env)))
  ca <- load_cache(V, cc, env)
  n_reg <- function(Csp, tau, lmin) { mk <- names(Csp)[Csp >= tau]
    if (!length(mk)) return(list(k = 0L, reg = list())); reg <- cluster_from_cache(mk, bd$edges)
    reg <- reg[lengths(reg) >= lmin]; list(k = length(reg), reg = reg) }
  rbindlist(lapply(LMINS, function(lm) rbindlist(lapply(TAU, function(t) {
    ob <- n_reg(bd$C_obs, t, lm); n_obs <- ob$k
    n_null <- mean(vapply(bd$C_surr, function(cs) n_reg(cs, t, lm)$k, numeric(1)))
    ev <- evaluate_ORs_qtn(ob$reg, ca$map, ca$qtab, ca$th$r2min, ca$th$dmax)
    data.table(env = env, l_min = lm, tau = t, n_obs = n_obs, n_null = n_null,
               est_TP = max(0, n_obs - n_null), est_FDR = n_null / max(n_obs, 1),
               real_TP = ev$TP, real_prec = if (is.na(ev$Precision)) NA_real_ else ev$Precision,
               real_recall = ev$Recall) }))))
}))

ok <- grid[n_obs > 0 & is.finite(real_prec)]
ok[, est_prec := est_TP / n_obs]
cat(sprintf("=== truth-free (null) vs truth, POOLED over %d env (V%s_c%s) ===\n", length(ENV), V, cc))
cat(sprintf("cor(est_TP, real_TP)     = %.3f\n", cor(ok$est_TP, ok$real_TP)))
cat(sprintf("cor(est_prec, real_prec) = %.3f\n", cor(ok$est_prec, ok$real_prec)))

## per (env,l_min): truth PR-AUC + null-based yield (max est_TP at est_FDR<=0.10); then average over env
se <- function(x){ x <- x[is.finite(x)]; if (length(x) < 2) NA_real_ else sd(x)/sqrt(length(x)) }
perenv <- grid[, .(truth_PR_AUC = pr_auc(real_recall, ifelse(is.na(real_prec), 0, real_prec)),
                   null_yield = max(c(0, est_TP[est_FDR <= 0.10])),
                   real_yield = max(c(0, real_TP[est_FDR <= 0.10]))), by = .(env, l_min)]
byl <- perenv[, .(truth_PR_AUC = mean(truth_PR_AUC), PR_se = se(truth_PR_AUC),
                  null_yield = mean(null_yield), real_yield = mean(real_yield)), by = l_min][order(l_min)]
cat("\n=== per-l_min, mean over env: truth PR-AUC vs null-based yield ===\n")
print(byl[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)])
cat(sprintf("\ntruth-optimal l_min (max PR-AUC) = %d | null-optimal l_min (max yield@estFDR<=0.10) = %d\n",
            byl$l_min[which.max(byl$truth_PR_AUC)], byl$l_min[which.max(byl$null_yield)]))

saveRDS(list(grid = grid, perenv = perenv, byl = byl), file.path(dir_data, sprintf("lmin_truthfree_agg_V%s_c%s.rds", V, cc)))
p <- ggplot(ok, aes(real_TP, est_TP, color = factor(l_min))) + geom_abline(linetype = 2, color = "grey") +
  geom_jitter(width = 0.1, height = 0.1, size = 1.2, alpha = 0.5) + scale_color_viridis_d(end = 0.9, name = "l_min") +
  labs(title = sprintf("V%s_c%s: null-estimated TP vs true TP, pooled 5 env (truth-free l_min diagnostic)", V, cc),
       subtitle = sprintf("cor(est_TP, real_TP) = %.2f; est_TP = observed - null regions", cor(ok$est_TP, ok$real_TP)),
       x = "real TP (QTN truth)", y = "est TP (observed - null)") + theme_bw(base_size = 11)
ggsave(file.path(dir_fig, sprintf("lmin_truthfree_agg_V%s_c%s.png", V, cc)), p, width = 7, height = 5.5, dpi = 150)
cat("wrote lmin_truthfree_agg figure\n")
