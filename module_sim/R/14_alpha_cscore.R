## module_sim/R/14_alpha_cscore.R
## Test the prediction: the C-score's benefit grows as the FDR threshold alpha
## gets more stringent (each single (rho,q*) run is then underpowered, so
## consistency across runs recovers signal no single run holds), and this is
## larger for EMMAX (calibrated/conservative) than LFMM (baseline-inflated =
## effectively lenient). alpha is a swept detection-stringency parameter (folded
## into C); l_min is a PRINCIPLED POST-C filter (drop size-1 regions).
##
## At each alpha, for each method: C_alpha = fraction of (rho,q*) cells where a
## SNP is candidate & FDR<alpha; gate on tau_C, cluster, keep regions >= LMIN loci,
## score (dedup-neutral); take best PR over tau_C. Compare to single-SNP (FDR<alpha,
## same LMIN filter). Prediction: PR_C - PR_single grows as alpha shrinks; EMMAX>LFMM.
## Run from LDscnR-paper/:  Rscript module_sim/R/14_alpha_cscore.R [V c env]
## Output (git-ignored): data/alpha_cscore_V*.rds, figures/alpha_cscore_V*.png

source("module_sim/R/_config.R")
suppressMessages(library(ggplot2))
a   <- commandArgs(trailingOnly = TRUE)
V   <- if (length(a) >= 1) a[1] else "2"
cc  <- if (length(a) >= 2) a[2] else "1"
env <- if (length(a) >= 3) a[3] else "1"
ALPHA  <- c(0.001, 0.005, 0.01, 0.05, 0.1)
QSTAR  <- seq(0, 0.95, by = 0.05)
C_GRID <- seq(0.1, 1.0, by = 0.1)
LMIN   <- 2L; R2LINK <- 0.4

## pool GTs + map + full ld_ws
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
th  <- score_thresholds(rbindlist(decs, fill = TRUE))
dcap <- d_from_rho(median(rbindlist(decs, fill = TRUE)$a_pred), 0.95)
RHO <- colnames(LDW); ncell <- length(RHO) * length(QSTAR)
cat(sprintf("V%s_c%s_env%s: %d SNPs, true_pos_QTN=%d | alpha=%s\n",
    V, cc, env, nrow(map), map[true_pos_QTN==TRUE,.N], paste(ALPHA, collapse=",")))

## C at fixed alpha across (rho, q*)
Cscore_alpha <- function(pv, alpha) {
  calls <- unlist(lapply(RHO, function(rc) { lw <- LDW[, rc]
    unlist(lapply(QSTAR, function(q) { thr <- stats::quantile(lw, q, na.rm = TRUE); cand <- which(lw >= thr)
      qv <- stats::p.adjust(pv[cand], "BH"); map$marker[cand[qv < alpha]] })) }))
  tb <- table(calls); C <- setNames(as.numeric(tb) / ncell, names(tb)); C
}
## score a marker set with the LMIN post-filter (drop size-1) using shared qtab
score_set <- function(mk, qtab) { if (!length(mk)) return(list(PR = 0, Precision = NA, Recall = 0, TP = 0, FP = 0))
  reg <- cluster_regions(mk, map, GTs, r2_link = R2LINK, dcap = dcap)
  reg <- reg[lengths(reg) >= LMIN]                              # principled post-C l_min filter
  ev  <- evaluate_ORs_qtn(reg, map, qtab, th$r2min, th$dmax)
  list(PR = ifelse(is.na(ev$PR), 0, ev$PR), Precision = ev$Precision, Recall = ev$Recall, TP = ev$TP, FP = ev$FP) }

res <- rbindlist(lapply(c("emx_p", "lfmm_p"), function(pcol) {
  pv <- map[[pcol]]; meth <- if (pcol == "emx_p") "EMMAX" else "LFMM"
  rbindlist(lapply(ALPHA, function(al) {
    C <- Cscore_alpha(pv, al)
    mk_single <- map$marker[p.adjust(pv, "BH") < al]
    gated <- lapply(C_GRID, function(ct) names(C)[C >= ct])
    allmk <- unique(c(mk_single, unlist(gated)))
    qtab  <- precompute_QTN_LD(GTs, map, allmk, 2e6, cores = 4)
    s_single <- score_set(mk_single, qtab)
    s_C <- lapply(gated, score_set, qtab = qtab)              # best over tau_C
    bi <- which.max(vapply(s_C, `[[`, numeric(1), "PR"))
    data.table(method = meth, alpha = al,
               PR_single = round(s_single$PR, 3), PR_C = round(s_C[[bi]]$PR, 3),
               dPR = round(s_C[[bi]]$PR - s_single$PR, 3), best_tauC = C_GRID[bi],
               Rec_single = round(s_single$Recall, 3), Rec_C = round(s_C[[bi]]$Recall, 3),
               Prec_single = round(s_single$Precision, 3), Prec_C = round(s_C[[bi]]$Precision, 3))
  }))
}))
cat("\n=== C-gated (best tau_C, >=2 post-filter) vs single-SNP, by alpha ===\n"); print(res)
saveRDS(res, file.path(dir_data, sprintf("alpha_cscore_V%s_c%s_env%s.rds", V, cc, env)))

## figure: improvement dPR vs alpha (stringent -> lenient), per method
p <- ggplot(res, aes(alpha, dPR, color = method)) +
  geom_line(linewidth = 0.7) + geom_point(size = 2.5) +
  scale_x_log10(breaks = ALPHA) + scale_color_manual(values = c(EMMAX = "#E1AF00", LFMM = "#3B9AB2")) +
  labs(title = sprintf("V%s_c%s_env%s: C-score benefit (PR_C - PR_single) vs alpha", V, cc, env),
       subtitle = "prediction: benefit grows as alpha tightens (weaker per-run signal); EMMAX > LFMM",
       x = "FDR threshold alpha (log; left = more stringent)", y = "PR improvement (C-gated - single-SNP)") +
  theme_bw(base_size = 11)
ggsave(file.path(dir_fig, sprintf("alpha_cscore_V%s_c%s_env%s.png", V, cc, env)), p, width = 8, height = 5, dpi = 150)
cat("wrote alpha_cscore figure\n")
