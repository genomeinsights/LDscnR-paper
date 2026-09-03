## =============================================================================
## null_structure_alignment.R -- do the surrogate nulls actually exercise the
## failure mode that population structure creates?
##
## The 3sp panel session finds that on their data a GRADIENT surrogate rejects in
## 92.5% of draws with a mean of 26 rejections, while my env_orth gives 0.338.
## Their question: are my surrogates orthogonal to the causal environment but
## also only weakly aligned with the DOMINANT genetic structure -- in which case
## my null is milder than it looks and the protective effect of filtering was
## measured in an easy regime.
##
## The surrogate environment vectors are not stored, only the p-values they
## produced. But alignment is recoverable from those: if a surrogate loads on the
## leading genetic axis, markers with large loadings on that axis get small
## p-values. So for each surrogate, correlate -log10(p) against |marker loading|
## on each leading PC. A null that does not exercise structure shows ~zero
## correlation; one that does shows a clear positive one.
##
## Env: INPUTS, ENVS, NPC
## =============================================================================
suppressMessages({library(data.table)})
IN   <- Sys.getenv("INPUTS", "/Volumes/Nemo/Nemo_sim/analysis_inputs")
ENVS <- as.integer(strsplit(Sys.getenv("ENVS", "1,2,3"), ",")[[1]])
NPC  <- as.integer(Sys.getenv("NPC", "3"))
## The 3sp panel session showed that a top-3 correlation is BLIND to diffuse
## structure: on their data gradient and nested surrogates are indistinguishable
## at PC1-3 (0.276 vs 0.290) and separate completely at PC1-10 (0.833 vs 0.456),
## and it is the PC1-10 value that tracks the rejection rate. So the leading
## SUBSPACE is measured here, not just the leading axes.
SUBS <- as.integer(strsplit(Sys.getenv("SUBSPACES", "3,10,20,50"), ",")[[1]])

out <- list(); subout <- list()
for (ev in ENVS) {
  pf <- file.path(IN, sprintf("panel_V2_c1_env%d.rds", ev)); if (!file.exists(pf)) next
  pan <- readRDS(pf)
  G <- pan$GTs                                  # individuals x markers
  keep <- which(is.finite(colSums(G)))
  Gc <- scale(G[, keep], center = TRUE, scale = FALSE)
  Gc[!is.finite(Gc)] <- 0
  ## PCs of individuals, then marker loadings on them
  K  <- tcrossprod(Gc) / ncol(Gc)
  ei <- eigen(K, symmetric = TRUE)
  pcs  <- ei$vectors[, seq_len(NPC), drop = FALSE]
  pmax_ <- max(SUBS)
  pcsA <- ei$vectors[, seq_len(pmax_), drop = FALSE]
  loadA <- abs(crossprod(Gc, pcsA))            # markers x max(SUBS)
  varex <- ei$values[seq_len(NPC)] / sum(ei$values[ei$values > 0])
  load <- abs(crossprod(Gc, pcs))               # markers x NPC
  envo <- pan$env_obs
  r_env <- suppressWarnings(cor(envo, pcs))     # observed env vs leading axes
  rm(pan, G, Gc, K, ei); gc(verbose = FALSE)

  for (bs in c("genetic","env_orth")) {
    qf <- file.path(IN, sprintf("pvals_V2_c1_env%d_emmax_%s_B100.rds", ev, bs))
    if (!file.exists(qf)) next
    z <- readRDS(qf); P <- z$p_perm[keep, , drop = FALSE]
    L <- -log10(pmax(P, .Machine$double.xmin))
    obs <- -log10(pmax(z$p_obs[keep], .Machine$double.xmin))
    rr <- sapply(seq_len(NPC), function(j)
            suppressWarnings(apply(L, 2, function(col) cor(col, load[, j], use = "complete.obs"))))
    r_obs <- sapply(seq_len(NPC), function(j) suppressWarnings(cor(obs, load[, j], use = "complete.obs")))
    ## R2 of the association signal on the leading SUBSPACE of marker loadings
    r2sub <- function(y) {
      y <- as.numeric(y); ok <- is.finite(y)
      vapply(SUBS, function(np) {
        X <- loadA[ok, seq_len(np), drop = FALSE]
        fit <- .lm.fit(cbind(1, X), y[ok])
        1 - sum(fit$residuals^2)/sum((y[ok] - mean(y[ok]))^2)
      }, numeric(1))
    }
    sub_surr <- apply(L, 2, r2sub)                   # SUBS x B
    sub_obs  <- r2sub(obs)
    sub_dt <- data.table(env = ev, basis = bs, n_pc = SUBS,
                         r2_surr_med = round(apply(matrix(sub_surr, nrow = length(SUBS)), 1, median), 4),
                         r2_surr_q90 = round(apply(matrix(sub_surr, nrow = length(SUBS)), 1, quantile, .9), 4),
                         r2_obs = round(sub_obs, 4))
    subout[[length(subout)+1]] <- sub_dt
    out[[length(out)+1]] <- data.table(env = ev, basis = bs,
      pc = seq_len(NPC), var_explained = round(varex, 3),
      env_obs_vs_pc = round(as.numeric(r_env), 3),
      surr_med = round(apply(rr, 2, median), 4),
      surr_q90 = round(apply(rr, 2, quantile, .9), 4),
      obs_align = round(r_obs, 4))
    rm(z, P, L); gc(verbose = FALSE)
  }
  cat(sprintf("  env%d done\n", ev))
}
res <- rbindlist(out)
OUT <- Sys.getenv("OUT", "module_sim_LDscnR/results/filter_then_test")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
fwrite(res, file.path(OUT, "null_structure_alignment.csv"))
sub <- rbindlist(subout)
fwrite(sub, file.path(OUT, "null_structure_subspace.csv"))
cat("\n=== R2 of the association signal on the leading GRM SUBSPACE ===\n")
print(dcast(sub, basis + n_pc ~ ., value.var = c("r2_surr_med","r2_surr_q90","r2_obs"),
            fun.aggregate = median)[order(basis, n_pc)])
cat(sprintf("  written: %s\n", file.path(OUT, "null_structure_alignment.csv")))
cat("\n=== alignment of surrogate p-values with leading genetic axes ===\n")
print(res)
cat("\n  surr_med: median across 100 surrogates of cor(-log10 p, |PC loading|)\n")
cat("  obs_align: the same for the OBSERVED data, as a reference point\n")
cat("  A null that exercises structure shows surr_med clearly above 0.\n")
