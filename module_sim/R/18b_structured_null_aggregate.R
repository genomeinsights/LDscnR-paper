## module_sim/R/18b_structured_null_aggregate.R
## Replicate-average the STRUCTURED-null tau_C calibration over env1-5 (run 18 first
## for each env). Single env repeatedly flukes (module policy) -- the calibration
## tau_C must be a mean over env. Reports mean +/- SE empirical FDR at each tau_C and
## the tau_C hitting FDR 0.05 / 0.10 on the MEAN curve; draws the per-env FDR curves
## (grey) with the mean +/- SE ribbon and the empirical PR-optimum band.
## Run from LDscnR-paper/:  Rscript module_sim/R/18b_structured_null_aggregate.R [V c]
## Output (git-ignored): data/structured_null_summary_V*.rds, figures/structured_null_summary_V*.png

suppressMessages({ library(data.table); library(ggplot2) })
mod <- "/Users/petrikem/gitlab/LDscnR-paper/module_sim"
dir_data <- file.path(mod, "data"); dir_fig <- file.path(mod, "figures")
a  <- commandArgs(trailingOnly = TRUE)
V  <- if (length(a) >= 1) a[1] else "2"
cc <- if (length(a) >= 2) a[2] else "1"
ENV <- 1:5

f <- sprintf("%s/structured_null_V%s_c%s_env%d.rds", dir_data, V, cc, ENV)
f <- f[file.exists(f)]
if (!length(f)) stop("no structured_null_*.rds -- run 18_structured_null.R across env first")
R <- rbindlist(lapply(seq_along(f), function(i) {
  x <- as.data.table(readRDS(f[i])$res); x[, env := i]; x }))

se <- function(x){ x <- x[is.finite(x)]; if (length(x) < 2) NA_real_ else sd(x)/sqrt(length(x)) }
## POOLED FDR (ratio-of-means): sum counts across env then divide. Robust to the
## small per-env denominators (n_obs 10-70) that make mean-of-ratios blow up when a
## single env fluctuates. se here is the per-env FDR spread, shown for honesty only.
agg <- R[, .(FDR = sum(n_null) / max(sum(n_obs), 1), se = se(FDR),
             n_obs = sum(n_obs), n_null = sum(n_null), n = .N), by = tau]
setorder(agg, tau)
agg[, FDR := rev(cummax(rev(FDR)))]   # enforce monotone-decreasing FDR in tau (step-down)
tau_at <- function(q){ ok <- agg[FDR <= q & n_obs > 0, tau]; if (length(ok)) min(ok) else NA_real_ }
t05 <- tau_at(0.05); t10 <- tau_at(0.10)

## per-env tau_C -- expose the replicate spread (env3 high, env4 low: calibration is noisy)
pe <- R[, .(tau05 = { ok <- tau[FDR <= 0.05 & n_obs > 0]; if (length(ok)) min(ok) else NA_real_ },
            tau10 = { ok <- tau[FDR <= 0.10 & n_obs > 0]; if (length(ok)) min(ok) else NA_real_ }), by = env]
cat(sprintf("=== STRUCTURED-null tau_C, POOLED over %d env (V%s_c%s) ===\n", length(f), V, cc))
print(agg[tau %in% seq(0.1, 1, 0.1), .(tau, FDR = round(FDR,3), n_obs, n_null = round(n_null,1))])
cat("\nper-env tau_C (spread):\n"); print(pe)
cat(sprintf("\nPooled tau_C: FDR<=0.10 -> %.2f | FDR<=0.05 -> %.2f  (l_min=2, alpha in {.001-.1}; PR-optimum ~0.35-0.5)\n",
            t10, t05))

p <- ggplot() +
  geom_line(data = R, aes(tau, FDR, group = env), color = "grey75", linewidth = 0.35) +
  geom_line(data = agg, aes(tau, FDR), color = "#D62828", linewidth = 0.9) +
  geom_hline(yintercept = c(0.05, 0.10), linetype = 2, color = "grey40") +
  annotate("rect", xmin = 0.35, xmax = 0.5, ymin = -Inf, ymax = Inf, fill = "#457B9D", alpha = 0.10) +
  annotate("text", x = 0.425, y = 0.92, label = "empirical\nPR-optimum", size = 3, color = "#457B9D") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = sprintf("V%s_c%s: C-score FDR vs tau_C (STRUCTURED null, pooled over %d env)", V, cc, length(f)),
       subtitle = sprintf("grey = per-env; red = pooled (ratio-of-means, step-down); FDR<=0.05 at tau_C=%.2f, <=0.10 at %.2f (l_min=2)", t05, t10),
       x = "tau_C", y = "empirical FDR") + theme_bw(base_size = 11)
ggsave(file.path(dir_fig, sprintf("structured_null_summary_V%s_c%s.png", V, cc)), p, width = 8, height = 5, dpi = 150)
saveRDS(list(agg = agg, per_env = R, tau05 = t05, tau10 = t10), file.path(dir_data, sprintf("structured_null_summary_V%s_c%s.rds", V, cc)))
cat("wrote structured_null_summary figure\n")
