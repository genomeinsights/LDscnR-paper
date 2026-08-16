## module_sim/R/24_winning_range.R
## HOW ROBUST is the C-score win? tau_C beats single-SNP alpha over a WIDE RANGE, not a
## knife-edge -- that is the point of integrating over tau_C. For each method: C-score
## PR=Prec x Rec across a tau_C grid vs the single-SNP BEST PR (max over its alpha
## sweep), replicate-averaged over env1-5, l_min=2, alpha-swept C (sim). Reports the
## "winning range": the span of tau_C where mean C-score PR >= mean single-SNP best PR.
## Ad-hoc from the cache (20). Run from LDscnR-paper/:
##   Rscript module_sim/R/24_winning_range.R [V c]

source("module_sim/R/21_estimate.R")
suppressMessages(library(ggplot2))
a  <- commandArgs(trailingOnly = TRUE)
V  <- if (length(a) >= 1) a[1] else "2"
cc <- if (length(a) >= 2) a[2] else "1"; ENV <- as.character(1:5)
ALPHA_SIM <- c(0.001, 0.01, 0.05, 0.1); LMIN <- 2L
TAUC   <- seq(0.05, 0.95, by = 0.025)
ALPHA_S <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1)

CA <- lapply(ENV, function(e) load_cache(V, cc, e))
pr_gate <- function(ca, mk) { reg <- if (length(mk)) cluster_from_cache(mk, ca$edges) else list()
  ev <- evaluate_ORs_qtn(reg[lengths(reg) >= LMIN], ca$map, ca$qtab, ca$th$r2min, ca$th$dmax)
  if (is.na(ev$PR)) 0 else ev$PR }

R <- rbindlist(lapply(seq_along(CA), function(i) { ca <- CA[[i]]
  rbindlist(lapply(c("emx_p", "lfmm_p"), function(p) {
    Cv <- cscore_count(ca$map[[p]], ca$LDW, ca$grid$RHO, ca$grid$QSTAR, ALPHA_SIM)
    cscore <- data.table(method = sub("_p$", "", p), env = i, tau = TAUC,
                         PR = vapply(TAUC, function(t) pr_gate(ca, names(Cv)[Cv >= t]), numeric(1)))
    ss_best <- max(vapply(ALPHA_S, function(al) pr_gate(ca, ca$map$marker[p.adjust(ca$map[[p]], "BH") < al]), numeric(1)))
    cscore[, ss_best := ss_best][] })) }))
se <- function(x){ x <- x[is.finite(x)]; if (length(x) < 2) NA_real_ else sd(x)/sqrt(length(x)) }
agg <- R[, .(PR = mean(PR), PR_se = se(PR), ss_best = mean(ss_best), ss_se = se(ss_best)), by = .(method, tau)]

for (mth in unique(agg$method)) { d <- agg[method == mth]; win <- d[PR >= ss_best, tau]
  ssb <- d$ss_best[1]
  if (length(win)) cat(sprintf("[%s] single-SNP best PR = %.3f | C-score wins for tau_C in [%.3f, %.3f] (%d/%d grid pts, %.0f%% of [0.05,0.95]); peak PR %.3f @ tau=%.3f\n",
        toupper(mth), ssb, min(win), max(win), length(win), nrow(d), 100*length(win)/nrow(d),
        max(d$PR), d$tau[which.max(d$PR)]))
  else cat(sprintf("[%s] single-SNP best PR = %.3f | C-score never exceeds it\n", toupper(mth), ssb)) }

saveRDS(list(agg = agg, per_env = R), file.path(dir_data, sprintf("winning_range_V%s_c%s.rds", V, cc)))
p <- ggplot(agg, aes(tau, PR)) +
  geom_ribbon(aes(ymin = PR - PR_se, ymax = PR + PR_se), fill = "#D62828", alpha = 0.15) +
  geom_line(color = "#D62828", linewidth = 0.9) +
  geom_hline(aes(yintercept = ss_best), linetype = 2, color = "#457B9D") +
  geom_ribbon(data = agg[PR >= ss_best], aes(ymin = ss_best, ymax = PR), fill = "#2A9D8F", alpha = 0.25) +
  facet_wrap(~method) +
  labs(title = sprintf("V%s_c%s: C-score PR vs tau_C (red) vs single-SNP best PR (blue dashed)", V, cc),
       subtitle = "green = winning range where C-score beats single-SNP; mean over 5 env, l_min=2",
       x = "tau_C", y = "PR (Precision x Recall)") + theme_bw(base_size = 11)
ggsave(file.path(dir_fig, sprintf("winning_range_V%s_c%s.png", V, cc)), p, width = 9, height = 5, dpi = 150)
cat("wrote winning_range figure\n")
