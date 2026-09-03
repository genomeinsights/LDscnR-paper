## module_sim/R/15_pr_auc.R
## !! REPLICATE-AVERAGE: run across env1-5, aggregate with 15b_pr_auc_aggregate.R. !!
## The principled headline: standard trapezoidal PR-AUC (one ORDERED knob per curve)
## -- retires the random-search AUC-PR* (which existed only for the raw UNORDERED grid).
##
## FIXED (not swept): clustering = decay-relative rho_ld=0.75 / rho_d=0.95 capped at
## 500kb (== the TP-match r2 and distance; self-calibrating per chromosome); l_min=2
## as a post-C filter. SWEPT into the C-score: rho (ld_w window) x q* x alpha.
##   C-score (ours): C = frac of (rho,q*,alpha) cells a SNP is called; sweep tau_C.
##   single-SNP (legacy): sweep its own alpha (FDR).
## Each threshold value: gate -> cluster (fixed) -> drop <l_min regions -> score
## (dedup-neutral) -> one (recall, precision) point. PR-AUC = trapezoid over the sweep.
## Run from LDscnR-paper/:  Rscript module_sim/R/15_pr_auc.R [V c env]
## Output (git-ignored): data/pr_auc_V*.rds, figures/pr_auc_V*.png

source("module_sim/R/_config.R")
suppressMessages(library(ggplot2))
a   <- commandArgs(trailingOnly = TRUE)
V   <- if (length(a) >= 1) a[1] else "2"
cc  <- if (length(a) >= 2) a[2] else "1"
env <- if (length(a) >= 3) a[3] else "1"
QSTAR   <- seq(0, 0.95, by = 0.05)
ALPHA_C <- c(0.001, 0.01, 0.05, 0.1)          # alpha SWEPT on sim (reviewer evidence; motivates fixing alpha empirically -- 15c)
TAUC    <- seq(0.05, 1.0, by = 0.05)          # C-gate sweep (the ordered knob)
ALPHA_S <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1)  # single-SNP alpha sweep
LMIN_SET <- c(1L, 2L, 4L, 8L)                 # post-C region-size filter, shown as linetype
RHO_LD <- 0.75; RHO_D <- 0.95; DCAP <- 1e5

## pool GTs + map + full ld_ws + decay
files <- list.files(SIM_DATA, full.names = TRUE,
                    pattern = sprintf("^adapt_bgs_chr[0-9]+_V%s_c%s_env%s\\.rds$", V, cc, env))
reps <- as.integer(sub(".*_chr([0-9]+)_.*", "\\1", basename(files))); files <- files[order(reps)]
maps <- gts <- ldws <- decs <- vector("list", length(files))
for (i in seq_along(files)) {
  d <- readRDS(files[i]); m <- as.data.table(d$map)
  m[, `:=`(Chr = paste0("R", i, "_", Chr), marker = paste0("R", i, "_", marker))]
  G <- d$GTs; colnames(G) <- paste0("R", i, "_", colnames(G))
  lw <- d$ld_ws; rownames(lw) <- m$marker
  ds <- as.data.table(d$LD_decay$decay_sum); ds[, Chr := paste0("R", i, "_", Chr)]
  maps[[i]] <- m; gts[[i]] <- G; ldws[[i]] <- lw; decs[[i]] <- ds
}
map <- flag_true_positive_QTNs(rbindlist(maps, fill = TRUE))
GTs <- do.call(cbind, gts); LDW <- do.call(rbind, ldws)[map$marker, ]
decay_sum <- rbindlist(decs, fill = TRUE)
th  <- score_thresholds(decay_sum, rho_r2 = 0.75, rho_d = 0.95, dmax_cap = DCAP)
RHO <- colnames(LDW); ncell <- length(RHO) * length(QSTAR) * length(ALPHA_C)
cat(sprintf("V%s_c%s_env%s: %d SNPs, true_pos_QTN=%d | TP-match r2>%.3f dist<%.0fkb | C-cells=%d\n",
    V, cc, env, nrow(map), map[true_pos_QTN==TRUE,.N], th$r2min, th$dmax/1e3, ncell))

clust <- function(mk) cluster_regions(mk, map, GTs, decay_sum = decay_sum,
                                      rho_ld = RHO_LD, rho_d = RHO_D, dcap_max = DCAP)

pr_one <- function(pcol) {
  pv <- map[[pcol]]
  ## 2D+alpha C-score: fraction of (rho,q*,alpha) cells where SNP is candidate & FDR<alpha
  calls <- unlist(lapply(RHO, function(rc) { lw <- LDW[, rc]
    unlist(lapply(QSTAR, function(q) { thr <- stats::quantile(lw, q, na.rm = TRUE); cand <- which(lw >= thr)
      qv <- stats::p.adjust(pv[cand], "BH")
      unlist(lapply(ALPHA_C, function(al) map$marker[cand[qv < al]])) })) }))
  C <- table(calls); Cv <- setNames(as.numeric(C) / ncell, names(C))
  single_sets <- lapply(ALPHA_S, function(al) map$marker[p.adjust(pv, "BH") < al])
  qtab <- precompute_QTN_LD(GTs, map, unique(c(names(Cv), unlist(single_sets))), 2e6, cores = 4)
  score <- function(mk) { reg <- if (length(mk)) clust(mk) else list()   # cluster once, filter per l_min
    rbindlist(lapply(LMIN_SET, function(lm) {
      ev <- evaluate_ORs_qtn(reg[lengths(reg) >= lm], map, qtab, th$r2min, th$dmax)
      data.table(l_min = lm, recall = ev$Recall, precision = ev$Precision) })) }
  pc_C <- rbindlist(lapply(TAUC, function(t) score(names(Cv)[Cv >= t])[, `:=`(curve = "C-score", knob = t)]))
  pc_S <- rbindlist(lapply(seq_along(ALPHA_S), function(i)
    score(single_sets[[i]])[, `:=`(curve = "single-SNP", knob = ALPHA_S[i])]))
  rbind(pc_C, pc_S)[, method := sub("_p$", "", pcol)][]
}

curves <- rbind(pr_one("emx_p"), pr_one("lfmm_p"))
auc <- curves[, .(PR_AUC = round(pr_auc(recall, precision), 4)), by = .(method, curve, l_min)]
cat("\n=== PR-AUC (trapezoidal) by l_min ===\n"); print(dcast(auc, method + l_min ~ curve, value.var = "PR_AUC"))
saveRDS(list(curves = curves, auc = auc, params = list(V = V, c = cc, env = env)),
        file.path(dir_data, sprintf("pr_auc_V%s_c%s_env%s.rds", V, cc, env)))

p <- ggplot(curves, aes(recall, precision, color = curve, linetype = factor(l_min))) +
  geom_path(linewidth = 0.6) + facet_wrap(~method) +
  scale_color_manual(values = c("C-score" = "#D62828", "single-SNP" = "#457B9D"), name = NULL) +
  scale_linetype_manual(values = c("1" = "dotted", "2" = "solid", "4" = "longdash", "8" = "dotdash"),
                        name = "l_min") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = sprintf("V%s_c%s_env%s: PR curves + trapezoidal PR-AUC (l_min as linetype)", V, cc, env),
       subtitle = "C-score: sweep tau_C (rho,q*,alpha folded in) | single-SNP: sweep alpha | clustering fixed rho0.75/0.95<=500kb",
       x = "Recall", y = "Precision") + theme_bw(base_size = 10)
ggsave(file.path(dir_fig, sprintf("pr_auc_V%s_c%s_env%s.png", V, cc, env)), p, width = 11, height = 5, dpi = 150)
cat("wrote pr_auc figure\n")
