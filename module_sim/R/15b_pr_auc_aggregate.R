## module_sim/R/15b_pr_auc_aggregate.R
## Replicate-average the PR-AUC headline over env1-5 (run 15_pr_auc.R for each first).
## Reports mean +/- SE PR-AUC per (method, curve) and draws (a) the mean PR-AUC with
## SE, (b) the per-env PR curves faceted, so the C-score-vs-single-SNP comparison is
## replicate-averaged (single env is unreliable -- see module policy).
## Run from LDscnR-paper/:  Rscript module_sim/R/15b_pr_auc_aggregate.R [V c]
## Output (git-ignored): data/pr_auc_summary_V*.rds, figures/pr_auc_summary_V*.png

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
mod <- "/Users/petrikem/gitlab/LDscnR-paper/module_sim"
dir_data <- file.path(mod, "data"); dir_fig <- file.path(mod, "figures")
a  <- commandArgs(trailingOnly = TRUE)
V  <- if (length(a) >= 1) a[1] else "2"
cc <- if (length(a) >= 2) a[2] else "1"
ENV <- 1:5

f <- sprintf("%s/pr_auc_V%s_c%s_env%d.rds", dir_data, V, cc, ENV)
f <- f[file.exists(f)]
if (!length(f)) stop("no pr_auc_*.rds -- run 15_pr_auc.R across env first")
L <- lapply(seq_along(f), function(i) { x <- readRDS(f[i])
  list(auc = as.data.table(x$auc)[, env := i], curves = as.data.table(x$curves)[, env := i]) })
AUC <- rbindlist(lapply(L, `[[`, "auc")); CUR <- rbindlist(lapply(L, `[[`, "curves"))
AUC[, method := toupper(method)]; CUR[, method := toupper(method)]

se <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) NA_real_ else sd(x) / sqrt(length(x)) }
agg <- AUC[, .(PR_AUC = mean(PR_AUC), se = se(PR_AUC), n = .N), by = .(method, curve, l_min)]
cat(sprintf("=== PR-AUC mean over %d env (V%s_c%s), by l_min ===\n", length(f), V, cc))
print(dcast(agg, method + l_min ~ curve, value.var = "PR_AUC")[, lapply(.SD, function(x) if (is.numeric(x)) round(x,3) else x)])

lt <- c("1" = "dotted", "2" = "solid", "4" = "longdash", "8" = "dotdash")
pA <- ggplot(agg, aes(factor(l_min), PR_AUC, color = curve, group = curve)) +
  geom_line() + geom_point(size = 2) +
  geom_errorbar(aes(ymin = PR_AUC - se, ymax = PR_AUC + se), width = 0.15) +
  facet_wrap(~method) + scale_color_manual(values = c("C-score" = "#D62828", "single-SNP" = "#457B9D"), name = NULL) +
  labs(title = "Mean PR-AUC +/- SE vs l_min", x = "l_min", y = "PR-AUC") + theme_bw(base_size = 10)
pB <- ggplot(CUR, aes(recall, precision, color = curve, linetype = factor(l_min),
                      group = interaction(curve, l_min, env))) +
  geom_path(alpha = 0.4, linewidth = 0.35) + facet_wrap(~method) +
  scale_color_manual(values = c("C-score" = "#D62828", "single-SNP" = "#457B9D"), name = NULL) +
  scale_linetype_manual(values = lt, name = "l_min") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "Per-env PR curves (l_min as linetype)", x = "Recall", y = "Precision") + theme_bw(base_size = 10)
ggsave(file.path(dir_fig, sprintf("pr_auc_summary_V%s_c%s.png", V, cc)),
       pA / pB + patchwork::plot_annotation(
         title = sprintf("V%s_c%s: C-score vs single-SNP, standard PR-AUC (mean over %d env)", V, cc, length(f))),
       width = 11, height = 9, dpi = 150)
saveRDS(list(agg = agg, AUC = AUC, CUR = CUR), file.path(dir_data, sprintf("pr_auc_summary_V%s_c%s.rds", V, cc)))
cat("wrote pr_auc_summary figure\n")
