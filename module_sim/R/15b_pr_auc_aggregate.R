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

pA <- ggplot(agg, aes(factor(l_min), PR_AUC, color = curve, group = curve)) +
  geom_line() + geom_point(size = 2) +
  geom_errorbar(aes(ymin = PR_AUC - se, ymax = PR_AUC + se), width = 0.15) +
  facet_wrap(~method) + scale_color_manual(values = c("C-score" = "#D62828", "single-SNP" = "#457B9D"), name = NULL) +
  labs(title = "Mean PR-AUC +/- SE vs l_min", x = "l_min", y = "PR-AUC") + theme_bw(base_size = 10)
## bottom: mean PR curve over env, faceted by l_min; C-score coloured by tau_C.
## l_min just relabels the tau_C axis -> same frontier across facets, tau_C colours
## slide toward lower values as l_min rises (l_min and tau_C are the same FP budget).
cagg <- CUR[, .(recall = mean(recall, na.rm = TRUE), precision = mean(precision, na.rm = TRUE)),
            by = .(method, curve, l_min, knob)]
csB <- cagg[curve == "C-score"]; ssB <- cagg[curve == "single-SNP"]
pB <- ggplot() +
  geom_path(data = ssB, aes(recall, precision), color = "#457B9D", linewidth = 0.5) +
  geom_point(data = ssB, aes(recall, precision), color = "#457B9D", size = 0.9) +
  geom_path(data = csB, aes(recall, precision), color = "grey70", linewidth = 0.4) +
  geom_point(data = csB, aes(recall, precision, color = knob), size = 1.8) +
  scale_color_viridis_c(name = "tau_C") + facet_grid(method ~ l_min, labeller = label_both) +
  coord_cartesian(xlim = c(0, 0.6), ylim = c(0, 1)) +
  labs(title = "Mean PR curve over env, faceted by l_min (C-score coloured by tau_C; blue = single-SNP)",
       x = "Recall", y = "Precision") + theme_bw(base_size = 10)
ggsave(file.path(dir_fig, sprintf("pr_auc_summary_V%s_c%s.png", V, cc)),
       pA / pB + patchwork::plot_annotation(
         title = sprintf("V%s_c%s: C-score vs single-SNP, standard PR-AUC (mean over %d env)", V, cc, length(f))),
       width = 13, height = 10, dpi = 150)
saveRDS(list(agg = agg, AUC = AUC, CUR = CUR), file.path(dir_data, sprintf("pr_auc_summary_V%s_c%s.rds", V, cc)))
cat("wrote pr_auc_summary figure\n")
