## module_sim/R/28_fold_validation.R
## RECONCILE sim <-> poster: does folding clustering+l_min INTO the C-score (poster
## style) beat the per-SNP C (sim style) on the sim, where we have truth?
##   per-SNP C   : frac of (rho,q*) cells where SNP is candidate & BH-FDR<0.05   (clustering/l_min AFTER)
##   folded  C   : frac of (rho,q*) cells where SNP is in a >=l_min cluster of the sig candidates (in-cell)
## Fair test: same alpha=0.05, same FINAL scoring (gate tau_C -> cluster -> >=l_min ->
## evaluate); only the GATING statistic differs. Compute folded C for all l_min in ONE
## pass (cache per-cell clusters). Replicate-averaged PR-AUC over env1-5, both methods.
## Run: Rscript module_sim/R/28_fold_validation.R [V c]

source("module_sim/R/21_estimate.R")
suppressMessages(library(ggplot2))
a  <- commandArgs(trailingOnly = TRUE)
V  <- if (length(a) >= 1) a[1] else "2"
cc <- if (length(a) >= 2) a[2] else "1"; ENV <- as.character(1:5)
ALPHA <- 0.05; QSTAR <- seq(0, 0.95, by = 0.05); TAUC <- seq(0.05, 1, by = 0.05); LMINS <- c(1L, 2L, 5L, 10L)

## folded C: for each cell, cluster the sig candidates once; a SNP counts for l_min if it
## sits in a cluster of size >= l_min. Returns matrix [n_marker x length(LMINS)].
cscore_folded <- function(pv, LDW, edges, rho, qstar, alpha, lmins) {
  idx <- setNames(seq_len(nrow(LDW)), rownames(LDW)); ncell <- length(rho) * length(qstar)
  cnt <- matrix(0, nrow(LDW), length(lmins), dimnames = list(rownames(LDW), as.character(lmins)))
  for (rc in rho) { lw <- LDW[, rc]
    for (q in qstar) { thr <- stats::quantile(lw, q, na.rm = TRUE); cand <- which(lw >= thr)
      if (!length(cand)) next; qv <- stats::p.adjust(pv[cand], "BH"); sig <- cand[qv < alpha]
      if (!length(sig)) next
      reg <- cluster_from_cache(rownames(LDW)[sig], edges); sz <- lengths(reg)
      for (j in seq_along(lmins)) { keep <- unlist(reg[sz >= lmins[j]]); if (length(keep)) cnt[idx[keep], j] <- cnt[idx[keep], j] + 1L } } }
  cnt / ncell
}
prauc <- function(Cv, cache, lmin) { pc <- rbindlist(lapply(TAUC, function(t) {
    mk <- names(Cv)[Cv >= t]; reg <- if (length(mk)) cluster_from_cache(mk, cache$edges) else list()
    ev <- evaluate_ORs_qtn(reg[lengths(reg) >= lmin], cache$map, cache$qtab, cache$th$r2min, cache$th$dmax)
    data.table(recall = ev$Recall, precision = if (is.na(ev$Precision)) NA_real_ else ev$Precision) }))
  pr_auc(pc$recall, pc$precision) }

R <- rbindlist(lapply(ENV, function(env) { ca <- load_cache(V, cc, env)
  rbindlist(lapply(c("emx_p", "lfmm_p"), function(p) {
    Cper <- cscore_count(ca$map[[p]], ca$LDW, ca$grid$RHO, QSTAR, ALPHA)           # per-SNP, alpha=0.05
    Cfold <- cscore_folded(ca$map[[p]], ca$LDW, ca$edges, ca$grid$RHO, QSTAR, ALPHA, LMINS)  # folded, all l_min
    rbindlist(lapply(seq_along(LMINS), function(j) { lm <- LMINS[j]
      data.table(env = env, method = sub("_p$", "", p), l_min = lm,
                 perSNP = prauc(setNames(Cper, rownames(ca$LDW)), ca, lm),
                 folded = prauc(setNames(Cfold[, j], rownames(ca$LDW)), ca, lm)) })) })) }))
se <- function(x){ x <- x[is.finite(x)]; if (length(x) < 2) NA_real_ else sd(x)/sqrt(length(x)) }
agg <- R[, .(perSNP = mean(perSNP), perSNP_se = se(perSNP), folded = mean(folded), folded_se = se(folded),
             delta = mean(folded - perSNP), delta_se = se(folded - perSNP)), by = .(method, l_min)]
cat("=== PR-AUC: folded (poster) vs per-SNP (sim) C, mean over 5 env, alpha=0.05 ===\n")
print(agg[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)])
cat("\n(delta = folded - per-SNP; >0 => folding clustering+l_min into C helps)\n")

saveRDS(list(R = R, agg = agg), file.path(dir_data, sprintf("fold_validation_V%s_c%s.rds", V, cc)))
m <- melt(R, id.vars = c("env","method","l_min"), measure.vars = c("perSNP","folded"), variable.name = "Cdef", value.name = "PR_AUC")
p <- ggplot(m[, .(PR_AUC = mean(PR_AUC), se = se(PR_AUC)), by = .(method, l_min, Cdef)],
            aes(factor(l_min), PR_AUC, color = Cdef, group = Cdef)) +
  geom_line() + geom_point(size = 2) + geom_errorbar(aes(ymin = PR_AUC - se, ymax = PR_AUC + se), width = 0.15) +
  facet_wrap(~method) + scale_color_manual(values = c(perSNP = "#457B9D", folded = "#D62828"), name = "C-score") +
  labs(title = sprintf("V%s_c%s: folded (poster) vs per-SNP (sim) C-score PR-AUC (mean 5 env)", V, cc),
       subtitle = "folded = cluster+l_min inside C; per-SNP = clustering/l_min as post-filter; same alpha=0.05, same final scoring",
       x = "l_min", y = "PR-AUC") + theme_bw(base_size = 11)
ggsave(file.path(dir_fig, sprintf("fold_validation_V%s_c%s.png", V, cc)), p, width = 9, height = 5, dpi = 150)
cat("wrote fold_validation figure\n")
