## module_sim/R/10_cscore.R
## q*-sweep consistency (C-score) as a q*-robust alternative to a single RMSC q*.
## Instead of picking one ld_w quantile threshold, sweep a grid of them; a SNP's
## C-score = fraction of thresholds at which it is called an outlier (candidate &
## BH-FDR-significant within candidates). SNPs recovered across many thresholds
## (high C) are stable; fragile terminal-spike calls get low C. We then gate on C
## (>= tau_C) to define outlier regions -- no single q* needed.
##
## Two questions this answers:
##   (1) does high C concentrate on QTN-linked SNPs? (C vs max_LD_with_QTN)
##   (2) does C-gating recover the recall that a fragile single q* loses?
##       (C_th sweep vs single-q* ld_w vs single-SNP, dedup-neutral scoring)
##
## Run from LDscnR-paper/:  Rscript module_sim/R/10_cscore.R [V c env] [method]
##   method = emx_p (default) or lfmm_p
## Output (git-ignored): data/cscore_V*_c*_env*_<m>.rds, figures/cscore_*.png

source("module_sim/R/_config.R")
suppressMessages(library(ggplot2))
a    <- commandArgs(trailingOnly = TRUE)
V    <- if (length(a) >= 1) a[1] else "2"
cc   <- if (length(a) >= 2) a[2] else "1"
env  <- if (length(a) >= 3) a[3] else "1"
pcol <- if (length(a) >= 4) a[4] else "emx_p"
mtag <- sub("_p$", "", pcol)
FDR  <- 0.05
Q_GRID  <- seq(0.50, 0.95, by = 0.05)          # ld_w quantile thresholds to sweep
C_GRID  <- seq(0.1, 1.0, by = 0.1)             # C-score gates

P   <- pool_group(sprintf("^adapt_bgs_chr[0-9]+_V%s_c%s_env%s\\.rds$", V, cc, env))
map <- flag_true_positive_QTNs(as.data.table(P$map))
th  <- score_thresholds(P$decay$decay_sum)
ldw <- P$ldw[map$marker]; pv <- map[[pcol]]
totQTN <- map[true_pos_QTN == TRUE, .N]
cat(sprintf("V%s_c%s_env%s [%s]: %d SNPs, %d chr, true_pos_QTN=%d\n",
    V, cc, env, mtag, nrow(map), uniqueN(map$Chr), totQTN))

## per-q called set (candidate & FDR-significant within candidates); no null
called_at_q <- function(q) {
  thr <- stats::quantile(ldw, q, na.rm = TRUE); cand <- which(ldw >= thr)
  qv <- stats::p.adjust(pv[cand], "BH"); map$marker[cand[qv < FDR]]
}
calls <- lapply(Q_GRID, called_at_q)
tab   <- table(unlist(calls))
Cs    <- data.table(marker = names(tab), C = as.numeric(tab) / length(Q_GRID))
map[, C := 0][Cs, on = "marker", C := i.C]
rmsc  <- rmsc_threshold(pv, ldw, grid = seq(0, 0.98, 0.01))
cat(sprintf("single-q* RMSC q*=%.2f | markers ever called across sweep: %d\n",
            rmsc$q_star, nrow(Cs)))

## (1) does C concentrate on QTN-linked SNPs?  bin C, show QTN linkage
map[, Cbin := cut(C, breaks = c(-.01, 0, .3, .6, .9, 1.01),
                  labels = c("0", "(0,.3]", "(.3,.6]", "(.6,.9]", "(.9,1]"))]
cat("\n=== C-score vs QTN linkage ===\n")
print(map[, .(n = .N, med_maxLD_QTN = round(median(max_LD_with_QTN), 3),
              frac_r2gt.5 = round(mean(max_LD_with_QTN > 0.5), 3),
              frac_true_pos = round(mean(true_pos_QTN), 4)), by = Cbin][order(Cbin)])

## (2) score: single-SNP, single-q* ld_w (RMSC), and C-gated at each tau_C
mk_single <- map$marker[p.adjust(pv, "BH") < FDR]
mk_ldw    <- { thr <- stats::quantile(ldw, rmsc$q_star); ci <- which(ldw >= thr)
              qv <- p.adjust(pv[ci], "BH"); map$marker[ci[qv < FDR]] }
gated <- lapply(C_GRID, function(ct) map[C >= ct, marker])
all_mk <- unique(c(mk_single, mk_ldw, unlist(gated)))
qtab   <- precompute_QTN_LD(P$GTs, map, all_mk, 2e6, cores = 4)
score  <- function(mk, tag) { if (!length(mk)) return(NULL)
  reg <- cluster_regions(mk, map, P$GTs, r2_link = 0.5, dcap = P_dcap)
  ev  <- evaluate_ORs_qtn(reg, map, qtab, th$r2min, th$dmax)
  data.table(set = tag, n_mk = length(mk), n_reg = length(reg), TP = ev$TP, FP = ev$FP,
             Recall = round(ev$Recall, 3), Precision = round(ev$Precision, 3)) }
P_dcap <- d_from_rho(median(P$decay$decay_sum$a_pred), 0.95)
res <- rbindlist(c(
  list(score(mk_single, "single-SNP"), score(mk_ldw, sprintf("ld_w q*=%.2f", rmsc$q_star))),
  lapply(seq_along(C_GRID), function(i) { r <- score(gated[[i]], sprintf("C>=%.1f", C_GRID[i]))
                                          if (!is.null(r)) r[, tau_C := C_GRID[i]]; r })), fill = TRUE)
cat("\n=== scoring: single-SNP vs single-q* ld_w vs C-gated (dedup-neutral) ===\n"); print(res)

## figure: PR path over tau_C, with single-SNP and single-q* references
cg <- res[!is.na(tau_C)]
refs <- res[is.na(tau_C)]
p <- ggplot(cg, aes(Recall, Precision)) +
  geom_path(color = "grey60") + geom_point(aes(color = tau_C), size = 3) +
  geom_point(data = refs, aes(shape = set), size = 3, color = "black") +
  scale_color_viridis_c(name = "tau_C") + scale_shape_manual(values = c(15, 17), name = NULL) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = sprintf("V%s_c%s_env%s [%s]: C-gated (q*-swept) vs single-q* vs single-SNP", V, cc, env, mtag),
       subtitle = "path = C-score gate tau_C; black = references; up-and-right = better") +
  theme_bw(base_size = 11)
ggsave(file.path(dir_fig, sprintf("cscore_V%s_c%s_env%s_%s.png", V, cc, env, mtag)), p, width = 7, height = 6, dpi = 150)
saveRDS(list(res = res, C = map[, .(marker, Chr, Pos, C, max_LD_with_QTN, true_pos_QTN)],
             rmsc_qstar = rmsc$q_star, params = list(V = V, c = cc, env = env, method = pcol)),
        file.path(dir_data, sprintf("cscore_V%s_c%s_env%s_%s.rds", V, cc, env, mtag)))
cat("wrote cscore figure + rds\n")
