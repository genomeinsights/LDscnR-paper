## =============================================================================
## module_sim_LDscnR / env_noise_null.R
##
## PK: add a residual component to the environment and test whether the
## validity/power dilemma relaxes.
##
## THE DILEMMA (mixed_basis_null.R): a surrogate has power against the real
## confounding only if it RESEMBLES the observed covariate, but resembling it
## makes the surrogate a replicate rather than a null. On bgs5 the environment is
## a deterministic spatial function -- REML puts vg/(vg+ve) = 1.00 -- so spatial
## surrogates are alternative spatial gradients and rediscover the same clusters.
##
## THE HYPOTHESIS THIS TESTS. A real environment carries measurement error and
## idiosyncratic variation. If the observed covariate is signal PLUS noise, then
## which units it discovers depends partly on its own noise realisation, and a
## surrogate with a different realisation should stop rediscovering them --
## locking falls -- while still carrying the structure that gives it power.
##
## Noise is added at the POPULATION level, because the environment is
## population-constant here and the real analogue is a per-population measurement
## error rather than per-individual. h2 is the share of the noisy covariate's
## variance that is the original signal, so h2 = 1 is the bgs5 environment as
## simulated.
##
## Reported per (h2, lambda): the observed neutral-chromosome tail, which is the
## TRUTH and moves with h2; the surrogate tail as a ratio to it (calibration);
## and the overlap between surrogate and observed discoveries against chance
## (validity). A usable null needs ratio ~ 1 AND overlap ~ 1x, and no setting has
## achieved both so far.
##
## RESULT: THE TEST IS INCONCLUSIVE AND THE DESIGN HAS A FLAW. Both are recorded
## because the flaw is easy to repeat.
##
##   h2    observed neutral tail   discoveries   locking measurable?
##   1.0          0.0576               16              yes
##   0.9          0.0183                0              no
##   0.7          0.0131                0              no
##   0.5          0.0052                0              no
##
## TEN PERCENT NOISE REMOVES EVERY DISCOVERY, so there is no observed set to
## measure locking against and the hypothesis cannot be tested this way. It is
## not merely that power falls: the observed neutral tail falls with it, 0.0576
## -> 0.0183 -> 0.0052, converging to and then below the uniform 0.0100. SIGNAL
## AND CONFOUNDING VANISH TOGETHER. The miscalibration that motivates the whole
## exercise exists because the simulated environment is noiseless and perfectly
## spatially structured, and it does not survive a covariate that is not.
##
## THE FLAW: vg/(vg+ve) stays at 1.00 at every h2. Noise was added at the
## population level, and with 80 populations against a 160-individual GRM of rank
## up to 159, a population-level vector lies largely in the GRM's span -- so REML
## absorbs the "residual" as genetic and no residual component is created. A real
## test needs individual-level noise, or more populations than the GRM can span.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/env_noise_null.R
## Env: CELL, TAG, ENV, NSIM, H2, LAMBDAS, SEED
## =============================================================================
suppressMessages({library(data.table); library(LDscnR)})
REMLE <- utils::getFromNamespace("emma.REMLE", "LDscnR")
EIGR  <- utils::getFromNamespace("emma.eigen.R.wo.Z", "LDscnR")
SIM   <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
CELL  <- Sys.getenv("CELL", "V0.5_c2"); TAG <- Sys.getenv("TAG", "nobgs")
ENV   <- as.integer(Sys.getenv("ENV", "1")); NSIM <- as.integer(Sys.getenv("NSIM", "60"))
H2    <- as.numeric(strsplit(Sys.getenv("H2", "1,0.9,0.7,0.5,0.3"), ",")[[1]])
LAMS  <- as.numeric(strsplit(Sys.getenv("LAMBDAS", "1,0.25,0.05"), ",")[[1]])
SEED  <- as.integer(Sys.getenv("SEED", "31"))

P <- list(); y0 <- NULL; GRM <- NULL; XY <- NULL; POP <- NULL; CT <- list()
for (i in 1:10) {
  x <- readRDS(sprintf("%s/adapt_%s_chr%d_%s_env%d.rds", SIM, TAG, i, CELL, ENV))
  m <- as.data.table(x$map)
  pr <- ld_prune_and_eMLG(GTs = x$GTs, stage1 = x$complexity_reduction$stage1,
        LD_decay = x$LD_decay, ld_w_col = "ld_w_095", ld_w_threshold = 0.025,
        score_threshold = 0.80, min_r2_rho = 0.5, distance_threshold = 1e5,
        compute_unflagged_eMLG = TRUE, cores = 1)
  E <- pr$eMLG; g <- as.data.table(pr$groups); gk <- g[g$group_id %in% colnames(E)]
  CT[[length(CT)+1]] <- merge(data.table(CL = paste0(i,"_",gk$group_id), Chr = as.character(gk$Chr)),
                              unique(m[, .(Chr = as.character(Chr), chr_type)]), by = "Chr")
  colnames(E) <- paste0(i, "_", colnames(E))
  if (is.null(y0)) { y0 <- as.numeric(x$env$env); GRM <- x$GRM
                     XY <- as.matrix(x$env[, .(x, y)]); POP <- x$env$pop }
  P[[length(P)+1]] <- list(prep = emmax_setup(E, x$GRM), cl = colnames(E))
}
CLall <- unlist(lapply(P, `[[`, "cl"))
ntrl  <- CLall %in% rbindlist(CT)[chr_type == "ntrl", CL]
scan1 <- function(v) unlist(lapply(P, function(z) as.numeric(emmax_fast(z$prep, v))))
tf    <- function(p) 100 * mean(p < 1e-4, na.rm = TRUE)
n <- length(y0); Xo <- matrix(1, n, 1); nrm <- function(M) M / mean(diag(M))
Kg <- nrm(GRM); Dm <- as.matrix(stats::dist(XY)); l <- stats::median(Dm[lower.tri(Dm)])
Ks <- nrm(exp(-0.5 * (Dm / l)^2))
pl <- sort(unique(POP)); idx <- lapply(pl, function(q) which(POP == q))

cat(sprintf("%s %s env%d | NSIM %d | uniform tail 0.0100\n\n", CELL, TAG, ENV, NSIM))
for (h2 in H2) {
  ## population-level noise, scaled so h2 is the signal's share of total variance
  set.seed(SEED)
  ep <- stats::rnorm(length(pl)); e <- numeric(n)
  for (a in seq_along(pl)) e[idx[[a]]] <- ep[a]
  e  <- e / stats::sd(e) * stats::sd(y0) * sqrt((1 - h2) / h2)
  yy <- y0 + e
  po <- scan1(yy); obs <- CLall[p.adjust(po, "BH") < 0.05]; ot <- tf(po[ntrl])
  ch <- length(obs) / length(CLall)
  reml_l <- LAMS[which.max(sapply(LAMS, function(L) {
    K <- L*Kg + (1-L)*Ks; REMLE(yy, Xo, K, eig.R = EIGR(K, Xo))$REML }))]
  cat(sprintf("h2 = %.1f | observed neutral tail %.4f | %d discoveries | REML lambda %.2f\n",
              h2, ot, length(obs), reml_l))
  for (L in LAMS) {
    K <- L*Kg + (1-L)*Ks; re <- REMLE(yy, Xo, K, eig.R = EIGR(K, Xo))
    Sig <- re$vg * K + re$ve * diag(n)
    ee <- eigen(Sig, symmetric = TRUE)
    A <- ee$vectors %*% diag(sqrt(pmax(ee$values, 0)), n)
    set.seed(SEED + 1); NE <- A %*% matrix(stats::rnorm(n * NSIM), n, NSIM)
    st <- numeric(NSIM); hits <- character(0)
    for (b in seq_len(NSIM)) {
      pn <- scan1(NE[, b]); st[b] <- tf(pn[ntrl])
      hits <- c(hits, CLall[p.adjust(pn, "BH") < 0.05])
    }
    ov <- if (length(hits)) mean(hits %in% obs) / ch else NA_real_
    cat(sprintf("   lambda %.2f  vg_share %.2f  surrogate tail %.4f (%.2fx truth)  locking %.0fx\n",
                L, re$vg/(re$vg + re$ve), mean(st), mean(st)/max(ot, 1e-9), ov))
  }
  cat("\n")
}
