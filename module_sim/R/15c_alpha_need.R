## module_sim/R/15c_alpha_need.R
## IS THE alpha SWEEP NEEDED? rho and q* span their ENTIRE natural domain [0,1], so
## sweeping them is assumption-free; alpha has no natural range, so folding a grid
## {.001,.01,.05,.1} into the C-score makes tau_C conditional on that arbitrary choice.
## Test: compare PR-AUC of the alpha-SWEPT C-score vs an alpha=0.05-FIXED C-score
## (C = frac of (rho,q*) cells where candidate & FDR<0.05, i.e. drop the alpha axis).
## If they match, drop alpha -> the C-score sweeps ONLY the two [0,1]-spanning knobs.
## Replicate-averaged over env1-5 (module rule), both methods, l_min=2.
## Run from LDscnR-paper/:  Rscript module_sim/R/15c_alpha_need.R [V c]
## Output (git-ignored): data/alpha_need_V*.rds, figures/alpha_need_V*.png

source("module_sim/R/_config.R")
suppressMessages(library(ggplot2))
a  <- commandArgs(trailingOnly = TRUE)
V  <- if (length(a) >= 1) a[1] else "2"
cc <- if (length(a) >= 2) a[2] else "1"
ENV <- 1:5
QSTAR <- seq(0, 0.95, by = 0.05); TAUC <- seq(0.05, 1.0, by = 0.05)
LMIN <- 2L; RHO_LD <- 0.75; RHO_D <- 0.95; DCAP <- 5e5
ALPHA_FULL <- c(0.001, 0.01, 0.05, 0.1); ALPHA_05 <- 0.05

## C-score over a chosen alpha set: frac of (rho,q*,alpha) cells candidate & FDR<alpha
cscore <- function(pv, LDW, RHO, markers, alphas) {
  ncell <- length(RHO) * length(QSTAR) * length(alphas)
  cnt <- integer(length(markers))
  for (rc in RHO) { lw <- LDW[, rc]
    for (q in QSTAR) { thr <- stats::quantile(lw, q, na.rm = TRUE); cand <- which(lw >= thr)
      if (!length(cand)) next; qv <- stats::p.adjust(pv[cand], "BH")
      for (al in alphas) { h <- cand[qv < al]; if (length(h)) cnt[h] <- cnt[h] + 1L } } }
  setNames(cnt / ncell, markers) }

one_env <- function(env) {
  files <- list.files(SIM_DATA, full.names = TRUE,
                      pattern = sprintf("^adapt_bgs_chr[0-9]+_V%s_c%s_env%s\\.rds$", V, cc, env))
  reps <- as.integer(sub(".*_chr([0-9]+)_.*", "\\1", basename(files))); files <- files[order(reps)]
  maps <- gts <- ldws <- decs <- vector("list", length(files))
  for (i in seq_along(files)) { d <- readRDS(files[i]); m <- as.data.table(d$map)
    m[, `:=`(Chr = paste0("R", i, "_", Chr), marker = paste0("R", i, "_", marker))]
    G <- d$GTs; colnames(G) <- paste0("R", i, "_", colnames(G))
    lw <- d$ld_ws; rownames(lw) <- m$marker
    ds <- as.data.table(d$LD_decay$decay_sum); ds[, Chr := paste0("R", i, "_", Chr)]
    maps[[i]] <- m; gts[[i]] <- G; ldws[[i]] <- lw; decs[[i]] <- ds }
  map <- flag_true_positive_QTNs(rbindlist(maps, fill = TRUE))
  GTs <- do.call(cbind, gts); LDW <- do.call(rbind, ldws)[map$marker, ]
  decay_sum <- rbindlist(decs, fill = TRUE); RHO <- colnames(LDW)
  th <- score_thresholds(decay_sum, rho_r2 = 0.75, rho_d = 0.95, dmax_cap = DCAP)
  clust <- function(mk) cluster_regions(mk, map, GTs, decay_sum = decay_sum,
                                        rho_ld = RHO_LD, rho_d = RHO_D, dcap_max = DCAP)
  out <- rbindlist(lapply(c("emx_p", "lfmm_p"), function(pcol) {
    pv <- map[[pcol]]
    Cf <- cscore(pv, LDW, RHO, map$marker, ALPHA_FULL)
    C5 <- cscore(pv, LDW, RHO, map$marker, ALPHA_05)
    qtab <- precompute_QTN_LD(GTs, map, unique(c(names(Cf)[Cf > 0], names(C5)[C5 > 0])), 2e6, cores = 2)
    prauc <- function(Cv) { pc <- rbindlist(lapply(TAUC, function(t) { mk <- names(Cv)[Cv >= t]
        reg <- if (length(mk)) clust(mk) else list()
        ev <- evaluate_ORs_qtn(reg[lengths(reg) >= LMIN], map, qtab, th$r2min, th$dmax)
        data.table(knob = t, recall = ev$Recall, precision = ev$Precision) }))
      pr_auc(pc$recall, pc$precision) }
    data.table(method = sub("_p$", "", pcol), env = env,
               PR_AUC_full = prauc(Cf), PR_AUC_a05 = prauc(C5)) }))
  cat(sprintf("env%s done\n", env)); out }

R <- rbindlist(lapply(ENV, one_env))
se <- function(x){ x <- x[is.finite(x)]; if (length(x) < 2) NA_real_ else sd(x)/sqrt(length(x)) }
agg <- R[, .(full = mean(PR_AUC_full), full_se = se(PR_AUC_full),
             a05 = mean(PR_AUC_a05),  a05_se = se(PR_AUC_a05),
             delta = mean(PR_AUC_full - PR_AUC_a05), delta_se = se(PR_AUC_full - PR_AUC_a05)), by = method]
cat("\n=== PR-AUC: alpha-swept vs alpha=0.05-fixed C-score (mean over 5 env, l_min=2) ===\n")
print(agg[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)])
cat("\n(delta = full - a05; if |delta| within delta_se, the alpha sweep is NOT needed)\n")

saveRDS(list(per_env = R, agg = agg), file.path(dir_data, sprintf("alpha_need_V%s_c%s.rds", V, cc)))
m <- melt(R, id.vars = c("method", "env"), measure.vars = c("PR_AUC_full", "PR_AUC_a05"),
          variable.name = "variant", value.name = "PR_AUC")
m[, variant := factor(variant, labels = c("alpha-swept", "alpha=0.05"))]
p <- ggplot(m, aes(variant, PR_AUC, group = env)) +
  geom_line(color = "grey70") + geom_point(aes(color = variant), size = 2) +
  facet_wrap(~method) + scale_color_manual(values = c("alpha-swept" = "#D62828", "alpha=0.05" = "#457B9D"), guide = "none") +
  labs(title = sprintf("V%s_c%s: does the alpha sweep help? PR-AUC per env (l_min=2)", V, cc),
       subtitle = "grey line joins the same env; overlapping clouds => alpha sweep redundant", x = NULL, y = "PR-AUC") +
  theme_bw(base_size = 11)
ggsave(file.path(dir_fig, sprintf("alpha_need_V%s_c%s.png", V, cc)), p, width = 8, height = 5, dpi = 150)
cat("wrote alpha_need figure\n")
