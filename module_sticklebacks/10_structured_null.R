## module_sticklebacks/10_structured_null.R
## STRUCTURED null tau_C calibration on the EMPIRICAL 3sp data -- the per-dataset
## calibration the sim (module_sim/18/18b) certified. Confounding here is GENETIC
## (marine/freshwater ecotype is the dominant structure axis, GRM corr ~0.95), so the
## structure basis is the GRM K itself (the sim used a spatial coordinate kernel):
##   surrogate y ~ MVN(0, K) then y_perp = resid(lm(y ~ eco))  (Gram-Schmidt).
## Same genetic autocorrelation as the real structure -> reproduces the structure-FP
## RATE; orthogonal to ecotype -> removes the true adaptation signal. Empirical FDR on
## the C-score then gives a truth-free tau_C FOR THIS DATASET.
## EMMAX only: fast EMMAX (eigendecompose K + rotate GTs once) makes B surrogates cheap;
## LFMM has no fast recompute (non-separable) -> its C landscape is reported in 09/11,
## its genome-wide spread being the known LD-correlated inflation.
## Run from LDscnR-paper/:  Rscript module_sticklebacks/10_structured_null.R [B]
## Output (git-ignored): module_sticklebacks/structured_null_3sp.rds + .png

suppressMessages({ library(LDscnR); library(data.table); library(ggplot2) })
source("module_sim/R/_config.R")
mod <- "/Users/petrikem/gitlab/LDscnR-paper/module_sticklebacks"
a <- commandArgs(trailingOnly = TRUE); B <- if (length(a) >= 1) as.integer(a[1]) else 100L
QSTAR <- seq(0, 0.95, by = 0.05); ALPHA_C <- 0.05; TAU <- seq(0.02, 1, by = 0.02)  # fixed alpha (null prices in anticonservativeness)

sr  <- readRDS(file.path(mod, "snp_stats_aligned.rds")); setDT(sr)
LDW <- readRDS("/Users/petrikem/gitlab/LDscnR-paper/3sp_data/ld_ws_3sp_MAF01.rds")[sr$marker, ]
RHO <- colnames(LDW); ncell <- length(RHO) * length(QSTAR) * length(ALPHA_C)
e <- new.env(); load("/Users/petrikem/gitlab/LDscnR-paper/LFMM_3sp/data/3sp_data.RData", envir = e)
GTs <- e$GTs_3sp; colnames(GTs) <- e$map_3sp$marker; GTs <- GTs[, sr$marker]
eco <- as.integer(e$pheno_3sp$ecotype == "Marine"); n <- length(eco)
K <- readRDS("/Users/petrikem/gitlab/LDscnR-paper/3sp_data/grm_null.rds")$K
stopifnot(nrow(K) == n, ncol(GTs) == nrow(sr))
pre <- fast_emmax_setup(GTs, K); rm(GTs); gc()          # rotate GTs ONCE; drop the big matrix

## genetic-structure surrogate: same covariance (K), orthogonal to observed ecotype
eK <- eigen(K, symmetric = TRUE); Lv <- pmax(eK$values, 0); Vk <- eK$vectors
gen_surrogate <- function() { y <- as.numeric(Vk %*% (sqrt(Lv) * rnorm(n))); as.numeric(resid(lm(y ~ eco))) }

Cscore <- function(pv) { cnt <- integer(nrow(sr))          # fast per-marker counter
  for (rc in RHO) { lw <- LDW[, rc]
    for (q in QSTAR) { thr <- stats::quantile(lw, q, na.rm = TRUE); cand <- which(lw >= thr)
      if (!length(cand)) next; qv <- stats::p.adjust(pv[cand], "BH")
      for (al in ALPHA_C) { h <- cand[qv < al]; if (length(h)) cnt[h] <- cnt[h] + 1L } } }
  cnt / ncell }
cat(sprintf("3sp EMMAX: %d SNPs, B=%d genetic-structure surrogates (K %dx%d)\n", nrow(sr), B, n, n))

C_obs <- Cscore(fast_emmax_p(pre, eco))                  # observed C from the SAME fast EMMAX
n_obs <- vapply(TAU, function(t) sum(C_obs >= t), numeric(1))
set.seed(1); null_cnt <- matrix(0, B, length(TAU))
for (b in seq_len(B)) { Cb <- Cscore(fast_emmax_p(pre, gen_surrogate()))
  null_cnt[b, ] <- vapply(TAU, function(t) sum(Cb >= t), numeric(1))
  if (b %% 10 == 0) cat("surrogate", b, "/", B, "\n") }
n_null <- colMeans(null_cnt); FDR <- pmin(1, n_null / pmax(n_obs, 1))
tau_at <- function(q){ ok <- which(FDR <= q & n_obs > 0); if (length(ok)) TAU[min(ok)] else NA_real_ }
res <- data.table(tau = TAU, n_obs = n_obs, n_null = round(n_null, 2), FDR = round(FDR, 3))
cat("\n=== 3sp EMMAX C-score empirical FDR (genetic-structure null) ===\n"); print(res[tau %in% seq(0.05, 0.5, 0.05)])
cat(sprintf("\n3sp EMMAX structured-null tau_C: FDR<=0.10 -> %.2f | FDR<=0.05 -> %.2f  (PER-SNP C, alpha=0.05, NO l_min filter)\n",
            tau_at(0.10), tau_at(0.05)))
## which SNPs/chr survive at the FDR<=0.05 tau_C
t05 <- tau_at(0.05); if (!is.na(t05)) { keep <- sr[C_obs >= t05]
  cat(sprintf("\nEMMAX survivors at tau_C=%.2f (FDR<=0.05): %d SNPs\n", t05, nrow(keep)))
  print(keep[, .N, by = Chr][order(-N)]) }
saveRDS(list(res = res, C_obs = C_obs, tau05 = tau_at(0.05), tau10 = tau_at(0.10)),
        file.path(mod, "structured_null_3sp.rds"))
p <- ggplot(res, aes(tau, FDR)) + geom_line(color = "#D62828", linewidth = 0.8) +
  geom_hline(yintercept = c(0.05, 0.10), linetype = 2, color = "grey40") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = sprintf("3sp EMMAX: C-score FDR vs tau_C (genetic-structure null, B=%d)", B),
       subtitle = "y ~ MVN(0,K) orthogonal to ecotype; dashed = FDR 0.05/0.10", x = "tau_C", y = "empirical FDR") +
  theme_bw(base_size = 11)
ggsave(file.path(mod, "structured_null_3sp.png"), p, width = 8, height = 5, dpi = 150)
cat("wrote structured_null_3sp figure\n")
