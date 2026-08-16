## module_sim/R/18_structured_null.R
## STRUCTURED null for tau_C calibration: instead of a naive phenotype permutation
## (which destroys the env<->structure confounding -> too permissive, see 17), draw
## null phenotypes that keep the observed env's SPATIAL AUTOCORRELATION but are
## ORTHOGONAL to it. Same autocorrelation -> reproduces the structure-FP RATE;
## orthogonal -> no true-QTN signal. Empirical FDR on C then gives a truth-free
## tau_C that (unlike 17) prices in structure confounding.
##
## Surrogate generator (EMPIRICALLY APPLICABLE -- needs only observed env + coords):
##   K_s = Gaussian kernel over the sampling coordinates; y ~ MVN(0, K_s);
##   y_perp = resid(lm(y ~ env)) (Gram-Schmidt orthogonalisation).
## For empirical data swap the coordinate kernel for the relevant structure basis
## (sampling coords if spatial confounding; genetic PCs / a low-rank structure
## kernel if the confounding is genetic) -- the rest is unchanged.
## Fast EMMAX (eigendecompose K once) is the future optimisation for large B / genome-wide.
## Run from LDscnR-paper/:  Rscript module_sim/R/18_structured_null.R [V c env] [B]
## Output (git-ignored): data/structured_null_V*.rds, figures/structured_null_V*.png

source("module_sim/R/_config.R")
suppressMessages(library(ggplot2))
a   <- commandArgs(trailingOnly = TRUE)
V   <- if (length(a) >= 1) a[1] else "2"
cc  <- if (length(a) >= 2) a[2] else "1"
env <- if (length(a) >= 3) a[3] else "1"
B   <- if (length(a) >= 4) as.integer(a[4]) else 50L
QSTAR <- seq(0, 0.95, by = 0.05); ALPHA_C <- c(0.001, 0.01, 0.05, 0.1); TAU <- seq(0.02, 1, by = 0.02)  # alpha SWEPT on sim (reviewer evidence)
gcta_grm <- function(X){ p<-colMeans(X)/2; k<-p>0&p<1; X<-X[,k,drop=F]; p<-p[k]
  Z<-sweep(sweep(X,2,2*p,"-"),2,sqrt(2*p*(1-p)),"/"); tcrossprod(Z)/ncol(Z) }

files <- list.files(SIM_DATA, full.names = TRUE,
                    pattern = sprintf("^adapt_bgs_chr[0-9]+_V%s_c%s_env%s\\.rds$", V, cc, env))
reps <- as.integer(sub(".*_chr([0-9]+)_.*", "\\1", basename(files))); files <- files[order(reps)]
Yobs <- COORDS <- NULL; PRE <- maps <- ldws <- mk_i <- vector("list", length(files))
for (i in seq_along(files)) { d <- readRDS(files[i]); m <- as.data.table(d$map)
  m[, marker := paste0("R", i, "_", marker)][, Chr := paste0("R", i, "_", Chr)]
  if (is.null(Yobs)) { Yobs <- d$env$env; COORDS <- cbind(d$env$x, d$env$y) }
  PRE[[i]] <- fast_emmax_setup(d$GTs, gcta_grm(d$GTs))     # eigendecompose + rotate GTs ONCE
  lw <- d$ld_ws; rownames(lw) <- m$marker
  maps[[i]] <- m; ldws[[i]] <- lw; mk_i[[i]] <- m$marker }
map <- flag_true_positive_QTNs(rbindlist(maps, fill = TRUE)); LDW <- do.call(rbind, ldws)[map$marker, ]
RHO <- colnames(LDW); ncell <- length(RHO) * length(QSTAR) * length(ALPHA_C); n <- length(Yobs)

## spatial surrogate generator: same autocorrelation, orthogonal to observed env
Dm <- as.matrix(dist(COORDS)); l <- median(Dm[lower.tri(Dm)])
Ks <- exp(-0.5 * (Dm / l)^2); eK <- eigen(Ks, symmetric = TRUE); Lv <- pmax(eK$values, 0); Vk <- eK$vectors
gen_surrogate <- function() { y <- Vk %*% (sqrt(Lv) * rnorm(n)); as.numeric(resid(lm(y ~ Yobs))) }

emmax_pooled <- function(Yv) { pv <- numeric(0)                        # fast EMMAX per file
  for (i in seq_along(PRE)) pv <- c(pv, setNames(fast_emmax_p(PRE[[i]], Yv), mk_i[[i]]))
  pv[map$marker] }
Cscore <- function(pv) { cnt <- integer(nrow(map))          # per-marker counter (fast; no unlist+table)
  for (rc in RHO) { lw <- LDW[, rc]
    for (q in QSTAR) { thr <- stats::quantile(lw, q, na.rm = TRUE); cand <- which(lw >= thr)
      if (!length(cand)) next; qv <- stats::p.adjust(pv[cand], "BH")
      for (al in ALPHA_C) { h <- cand[qv < al]; if (length(h)) cnt[h] <- cnt[h] + 1L } } }
  setNames(cnt / ncell, map$marker) }
cat(sprintf("V%s_c%s_env%s: %d SNPs, B=%d structured surrogates (spatial l=%.2f)\n", V, cc, env, nrow(map), B, l))

C_obs <- Cscore(emmax_pooled(Yobs)); n_obs <- vapply(TAU, function(t) sum(C_obs >= t), numeric(1))
set.seed(1); null_cnt <- matrix(0, B, length(TAU))
for (b in seq_len(B)) { Cb <- Cscore(emmax_pooled(gen_surrogate()))
  null_cnt[b, ] <- vapply(TAU, function(t) sum(Cb >= t), numeric(1))
  if (b %% 10 == 0) cat("surrogate", b, "/", B, "\n") }
n_null <- colMeans(null_cnt); FDR <- pmin(1, n_null / pmax(n_obs, 1))
tau_at <- function(q){ ok <- which(FDR <= q & n_obs > 0); if(length(ok)) TAU[min(ok)] else NA_real_ }
res <- data.table(tau = TAU, n_obs = n_obs, n_null = round(n_null,1), FDR = round(FDR,3))
cat("\n=== C-score empirical FDR (STRUCTURED spatial null) ===\n"); print(res[tau %in% seq(0.1,1,0.1)])
cat(sprintf("\nStructured-null tau_C: FDR<=0.10 -> %.2f | FDR<=0.05 -> %.2f  (naive-null was 0.06; PR-optimum ~0.35-0.5)\n",
    tau_at(0.10), tau_at(0.05)))
saveRDS(list(res=res, tau05=tau_at(0.05), tau10=tau_at(0.10)), file.path(dir_data, sprintf("structured_null_V%s_c%s_env%s.rds", V, cc, env)))
p <- ggplot(res, aes(tau, FDR)) + geom_line(color="#D62828", linewidth=0.8) +
  geom_hline(yintercept=c(0.05,0.10), linetype=2, color="grey40") +
  geom_vline(xintercept=c(0.35,0.5), linetype=3, color="#457B9D") +
  annotate("text", x=0.42, y=0.9, label="empirical\nPR-optimum", size=3, color="#457B9D") +
  coord_cartesian(ylim=c(0,1)) +
  labs(title=sprintf("V%s_c%s_env%s: C-score FDR vs tau_C (STRUCTURED spatial null, B=%d)", V, cc, env, B),
       subtitle="orthogonal-but-same-autocorrelation surrogates; dashed=FDR 0.05/0.10; dotted=PR-optimum", x="tau_C", y="empirical FDR") +
  theme_bw(base_size=11)
ggsave(file.path(dir_fig, sprintf("structured_null_V%s_c%s_env%s.png", V, cc, env)), p, width=8, height=5, dpi=150)
cat("wrote structured_null figure\n")
