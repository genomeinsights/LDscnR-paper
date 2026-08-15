## module_sim/R/17_perm_null.R
## Permutation-null calibration of tau_C for EMMAX -- the truth-free rule for
## choosing the C-score cutoff. Individuals are shared across the pooled files
## (same landscape grid points, identical env), so we permute env ONCE per
## permutation and recompute per-file EMMAX (each file's own GCTA GRM, all
## markers; reproduces precomputed emx_p at r~0.96), pool -> null p -> null
## C-score (over rho x q* x alpha). Empirical FDR on C:
##   FDR(tau) = E[# null SNPs with C>=tau] / (# observed SNPs with C>=tau).
## Then check the FDR-controlled tau lands near the empirical PR-optimum (~0.35-0.5).
## Run from LDscnR-paper/:  Rscript module_sim/R/17_perm_null.R [V c env] [B]
## Output (git-ignored): data/perm_null_V*.rds, figures/perm_null_V*.png

source("module_sim/R/_config.R")
suppressMessages(library(ggplot2))
a   <- commandArgs(trailingOnly = TRUE)
V   <- if (length(a) >= 1) a[1] else "2"
cc  <- if (length(a) >= 2) a[2] else "1"
env <- if (length(a) >= 3) a[3] else "1"
B   <- if (length(a) >= 4) as.integer(a[4]) else 50L
QSTAR <- seq(0, 0.95, by = 0.05); ALPHA_C <- c(0.001, 0.01, 0.05, 0.1)
TAU <- seq(0.02, 1, by = 0.02)

gcta_grm <- function(X) { p <- colMeans(X) / 2; k <- p > 0 & p < 1; X <- X[, k, drop = FALSE]; p <- p[k]
  Z <- sweep(sweep(X, 2, 2*p, "-"), 2, sqrt(2*p*(1-p)), "/"); tcrossprod(Z) / ncol(Z) }

## load per file: GTs + GRM (once); pool map + ld_ws; env is shared
files <- list.files(SIM_DATA, full.names = TRUE,
                    pattern = sprintf("^adapt_bgs_chr[0-9]+_V%s_c%s_env%s\\.rds$", V, cc, env))
reps <- as.integer(sub(".*_chr([0-9]+)_.*", "\\1", basename(files))); files <- files[order(reps)]
Y <- NULL; GTl <- Kl <- vector("list", length(files)); maps <- ldws <- vector("list", length(files)); mk_i <- vector("list", length(files))
for (i in seq_along(files)) {
  d <- readRDS(files[i]); m <- as.data.table(d$map); m[, marker := paste0("R", i, "_", marker)][, Chr := paste0("R", i, "_", Chr)]
  if (is.null(Y)) Y <- d$env$env
  GTl[[i]] <- d$GTs; Kl[[i]] <- gcta_grm(d$GTs); lw <- d$ld_ws; rownames(lw) <- m$marker
  maps[[i]] <- m; ldws[[i]] <- lw; mk_i[[i]] <- m$marker
}
map <- flag_true_positive_QTNs(rbindlist(maps, fill = TRUE)); LDW <- do.call(rbind, ldws)[map$marker, ]
RHO <- colnames(LDW); ncell <- length(RHO) * length(QSTAR) * length(ALPHA_C)
cat(sprintf("V%s_c%s_env%s: %d SNPs, %d files, B=%d perms | C-cells=%d\n", V, cc, env, nrow(map), length(files), B, ncell))

## pooled EMMAX p for a phenotype Y (per-file EMMAX with per-file GRM)
emmax_pooled <- function(Yv) { pv <- numeric(0)
  for (i in seq_along(GTl)) { f <- emmax(Y = Yv, X = GTl[[i]], K = Kl[[i]]); p <- f$pval
    p <- setNames(p, mk_i[[i]]); pv <- c(pv, p) }; pv[map$marker] }
## C-score over (rho, q*, alpha) from a pooled p vector -> named C per SNP
Cscore <- function(pv) { calls <- unlist(lapply(RHO, function(rc) { lw <- LDW[, rc]
    unlist(lapply(QSTAR, function(q) { thr <- stats::quantile(lw, q, na.rm = TRUE); cand <- which(lw >= thr)
      qv <- stats::p.adjust(pv[cand], "BH"); unlist(lapply(ALPHA_C, function(al) map$marker[cand[qv < al]])) })) }))
  tb <- table(calls); C <- setNames(numeric(nrow(map)), map$marker); C[names(tb)] <- as.numeric(tb) / ncell; C }

## observed
p_obs <- emmax_pooled(Y); C_obs <- Cscore(p_obs)
cat(sprintf("observed EMMAX vs precomputed emx_p: cor(-log10p)=%.3f\n",
    cor(-log10(pmax(p_obs,1e-300)), -log10(pmax(map$emx_p,1e-300)))))
n_obs <- vapply(TAU, function(t) sum(C_obs >= t), numeric(1))

## null: permute env, recompute, count null SNPs >= tau (accumulate, don't store per-SNP)
set.seed(1); null_cnt <- matrix(0, nrow = B, ncol = length(TAU))
for (b in seq_len(B)) { Cb <- Cscore(emmax_pooled(sample(Y)))
  null_cnt[b, ] <- vapply(TAU, function(t) sum(Cb >= t), numeric(1))
  if (b %% 10 == 0) cat("perm", b, "/", B, "\n") }
n_null <- colMeans(null_cnt)
FDR <- pmin(1, n_null / pmax(n_obs, 1))

tau_at <- function(q) { ok <- which(FDR <= q & n_obs > 0); if (length(ok)) TAU[min(ok)] else NA_real_ }
res <- data.table(tau = TAU, n_obs = n_obs, n_null = round(n_null, 1), FDR = round(FDR, 3))
cat("\n=== C-score empirical FDR (permutation null) ===\n")
print(res[tau %in% seq(0.1, 1, 0.1)])
cat(sprintf("\nFDR-controlled tau_C: FDR<=0.10 -> %.2f | FDR<=0.05 -> %.2f | FDR<=0.01 -> %.2f\n",
    tau_at(0.10), tau_at(0.05), tau_at(0.01)))
cat("(empirical PR-optimum was tau_C ~ 0.35-0.5)\n")
saveRDS(list(res = res, C_obs = C_obs, TAU = TAU, tau05 = tau_at(0.05), tau10 = tau_at(0.10)),
        file.path(dir_data, sprintf("perm_null_V%s_c%s_env%s.rds", V, cc, env)))

p <- ggplot(res, aes(tau, FDR)) + geom_line(color = "#D62828", linewidth = 0.8) +
  geom_hline(yintercept = c(0.05, 0.10), linetype = 2, color = "grey40") +
  geom_vline(xintercept = c(0.35, 0.5), linetype = 3, color = "#457B9D") +
  annotate("text", x = 0.42, y = 0.9, label = "empirical\nPR-optimum", size = 3, color = "#457B9D") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = sprintf("V%s_c%s_env%s: C-score empirical FDR vs tau_C (permutation null, B=%d)", V, cc, env, B),
       subtitle = "dashed = FDR 0.05/0.10; dotted blue = empirical PR-optimum tau_C", x = "tau_C", y = "empirical FDR") +
  theme_bw(base_size = 11)
ggsave(file.path(dir_fig, sprintf("perm_null_V%s_c%s_env%s.png", V, cc, env)), p, width = 8, height = 5, dpi = 150)
cat("wrote perm_null figure\n")
