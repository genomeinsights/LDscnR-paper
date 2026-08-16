## module_sim/R/22_tauc_by_method.R
## Does the best tau_C differ between EMMAX and LFMM, and can ONE tau_C serve both?
## Ad-hoc from the cache (20): per method, PR(tau_C) = Precision x Recall (engine
## convention) over a fine tau_C grid, l_min=2, alpha-SWEPT (sim convention),
## replicate-averaged over env1-5. Reports (a) each method's own best tau_C, (b) the
## C-score distribution shift that explains any gap, (c) the loss from forcing a single
## shared tau_C. The principled resolution -- calibrate to a common FDR TARGET via the
## structured null, which maps that target to a METHOD-SPECIFIC tau_C automatically.
## Run from LDscnR-paper/:  Rscript module_sim/R/22_tauc_by_method.R [V c]

source("module_sim/R/21_estimate.R")   # brings load_cache + helpers + _config
suppressMessages(library(ggplot2))
a  <- commandArgs(trailingOnly = TRUE)
V  <- if (length(a) >= 1) a[1] else "2"
cc <- if (length(a) >= 2) a[2] else "1"; ENV <- as.character(1:5)
ALPHA_SIM <- c(0.001, 0.01, 0.05, 0.1)          # sim keeps the alpha sweep
TAUC <- seq(0.05, 0.95, by = 0.025); LMIN <- 2L

pr_at <- function(cache, pcol, tau) {           # PR = Prec x Rec at one tau_C
  Cv <- cscore_count(cache$map[[pcol]], cache$LDW, cache$grid$RHO, cache$grid$QSTAR, ALPHA_SIM)
  vapply(tau, function(t) { mk <- names(Cv)[Cv >= t]
    reg <- if (length(mk)) cluster_from_cache(mk, cache$edges) else list()
    ev <- evaluate_ORs_qtn(reg[lengths(reg) >= LMIN], cache$map, cache$qtab, cache$th$r2min, cache$th$dmax)
    if (is.na(ev$PR)) 0 else ev$PR }, numeric(1)) }

CA <- lapply(ENV, function(e) load_cache(V, cc, e))
R <- rbindlist(lapply(seq_along(CA), function(i)
  rbindlist(lapply(c("emx_p", "lfmm_p"), function(p)
    data.table(env = i, method = sub("_p$", "", p), tau = TAUC, PR = pr_at(CA[[i]], p, TAUC))))))
agg <- R[, .(PR = mean(PR), se = sd(PR)/sqrt(.N)), by = .(method, tau)]

best <- agg[, .SD[which.max(PR)], by = method]
cat("=== best tau_C per method (max mean PR=Prec x Rec, l_min=2, 5 env) ===\n"); print(best[, .(method, tau, PR = round(PR,4))])

## shared tau_C: the single tau maximizing the mean PR across methods; report each
## method's loss vs its OWN best
joint <- agg[, .(mean_PR = mean(PR)), by = tau][which.max(mean_PR)]
tsh <- joint$tau
loss <- agg[tau == tsh][best, on = "method"][, .(method, tau_shared = tsh, PR_shared = round(PR,4),
            tau_own = i.tau, PR_own = round(i.PR,4), loss = round(i.PR - PR, 4))]
cat(sprintf("\n=== shared tau_C = %.3f (max joint PR) : loss vs own-best ===\n", tsh)); print(loss)

## C-score distribution shift (why the bests differ): quantiles of C>0 per method, env1
c1 <- CA[[1]]
q_emx  <- quantile(cscore_count(c1$map$emx_p,  c1$LDW, c1$grid$RHO, c1$grid$QSTAR, ALPHA_SIM), c(.9,.99,.999,1))
q_lfmm <- quantile(cscore_count(c1$map$lfmm_p, c1$LDW, c1$grid$RHO, c1$grid$QSTAR, ALPHA_SIM), c(.9,.99,.999,1))
cat("\n=== C-score upper quantiles (env1): LFMM shifted right => higher tau_C for same operating point ===\n")
print(rbind(EMMAX = round(q_emx,3), LFMM = round(q_lfmm,3)))

saveRDS(list(agg = agg, best = best, shared = loss), file.path(dir_data, sprintf("tauc_by_method_V%s_c%s.rds", V, cc)))
p <- ggplot(agg, aes(tau, PR, color = method)) +
  geom_ribbon(aes(ymin = PR - se, ymax = PR + se, fill = method), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) + geom_vline(data = best, aes(xintercept = tau, color = method), linetype = 2) +
  geom_vline(xintercept = tsh, linetype = 3, color = "grey30") +
  annotate("text", x = tsh, y = 0, label = sprintf(" shared %.2f", tsh), hjust = 0, vjust = 0, size = 3, color = "grey30") +
  scale_color_manual(values = c(emx = "#1D3557", lfmm = "#E63946")) +
  scale_fill_manual(values = c(emx = "#1D3557", lfmm = "#E63946")) +
  labs(title = sprintf("V%s_c%s: PR=Prec x Rec vs tau_C by method (mean+/-SE, 5 env, l_min=2)", V, cc),
       subtitle = "dashed = each method's best tau_C; dotted = shared (max joint PR). Sim: alpha swept.",
       x = "tau_C", y = "PR (Precision x Recall)") + theme_bw(base_size = 11)
ggsave(file.path(dir_fig, sprintf("tauc_by_method_V%s_c%s.png", V, cc)), p, width = 8, height = 5, dpi = 150)
cat("\nwrote tauc_by_method figure\n")
