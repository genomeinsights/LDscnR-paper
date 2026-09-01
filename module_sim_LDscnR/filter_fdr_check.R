## =============================================================================
## filter_fdr_check.R -- does BH INSIDE an ld_w-selected set still control FDR?
##
## filter_then_test.R measured POWER. That is necessary but not sufficient: a
## filter that concentrates false positives buys discoveries you cannot trust.
## ld_w is phenotype-blind, but it is NOT independent of the test statistic --
## both respond to MAF, and ld_w selects high-LD regions where tests are
## correlated. Either could distort BH inside the selected set.
##
## NULL: surrogate ENVIRONMENTS, not a naive phenotype shuffle. The surrogates
## preserve spatial and genetic structure (basis = env_orth / genetic / latent),
## which is essential here -- a naive shuffle destroys exactly the structure ld_w
## keys on and would make filtered BH look well behaved for the wrong reason.
##
## The filter needs no recomputation per draw: ld_w is computed from genotypes
## alone, so the SAME selection applies to the observed data and to all B
## surrogates. (2c's "recompute the filter per draw" caution binds for a
## phenotype-dependent filter; it does not bind for this one.)
##
## Under a surrogate there is no true association, so EVERY rejection is false.
## Two things are reported per selection:
##   reject_rate -- fraction of surrogates yielding >= 1 rejection. BH controls
##                  FDR, and with no true signal FDR = P(any rejection), so this
##                  should sit at or below ALPHA.
##   n_rej       -- mean rejections per surrogate, filtered vs genome-wide.
## An ld_w-filtered set that rejects MORE often than genome-wide is buying its
## power by inflating false positives, and the power result does not stand.
##
## Env: INPUTS, OUT, KS, ALPHA, ENGINE, BASIS, ENVS
## =============================================================================
suppressMessages({library(data.table)})

IN    <- Sys.getenv("INPUTS", "/Volumes/Nemo/Nemo_sim/analysis_inputs")
OUT   <- Sys.getenv("OUT", "module_sim_LDscnR/results/filter_then_test")
KS    <- as.integer(strsplit(Sys.getenv("KS", "500,1000,5000,20000,50000"), ",")[[1]])
ALPHA <- as.numeric(Sys.getenv("ALPHA", "0.05"))
ENG   <- Sys.getenv("ENGINE", "emmax")
BASES <- strsplit(Sys.getenv("BASIS", "genetic,env_orth"), ",")[[1]]
ENVS  <- as.integer(strsplit(Sys.getenv("ENVS", "1,2,3,4,5,6,7,8,9,10"), ",")[[1]])
LDCOL <- Sys.getenv("LDCOL", "rho_0.95")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

n_rej <- function(p, idx) {                      # rejections by BH within idx only
  q <- p.adjust(p[idx], "BH"); sum(q < ALPHA, na.rm = TRUE)
}

res <- list()
for (ev in ENVS) {
  pf <- file.path(IN, sprintf("panel_V2_c1_env%d.rds", ev))
  if (!file.exists(pf)) next
  pan <- readRDS(pf)
  ldw <- pan$ld_ws[, LDCOL]
  ## MAF from the panel genotypes, for the same-size control
  af  <- colMeans(pan$GTs, na.rm = TRUE) / 2
  maf <- pmin(af, 1 - af)
  rm(pan); gc(verbose = FALSE)
  set.seed(4242 + ev)
  sels <- list(ld_w = order(-ldw), MAF = order(-maf), random = sample.int(length(ldw)))

  for (bs in BASES) {
    qf <- file.path(IN, sprintf("pvals_V2_c1_env%d_%s_%s_B100.rds", ev, ENG, bs))
    if (!file.exists(qf)) next
    z  <- readRDS(qf)
    P  <- z$p_perm; B <- ncol(P); nm <- nrow(P)
    ## genomic inflation. These inputs were generated with gc_mode = "none", so
    ## any anticonservatism may be uncorrected inflation rather than a filter effect.
    lam_o <- z$lambda_obs; lam_s <- mean(z$lambda_surr, na.rm = TRUE); gcm <- z$gc_mode
    stopifnot(nm == length(ldw))

    ## genome-wide baseline, per surrogate
    gw <- vapply(seq_len(B), function(b) n_rej(P[, b], seq_len(nm)), integer(1))
    res[[length(res)+1]] <- data.table(env = ev, basis = bs, engine = ENG, method = "genome_wide",
        k = nm, B = B, reject_rate = mean(gw > 0), mean_rej = mean(gw), max_rej = max(gw),
        obs_rej = n_rej(z$p_obs, seq_len(nm)), lambda_obs = lam_o, lambda_surr = lam_s, gc_mode = gcm)

    for (kk in KS) {
      if (kk >= nm) next
      for (mth in names(sels)) {
        idx <- head(sels[[mth]], kk)
        v <- vapply(seq_len(B), function(b) n_rej(P[, b], idx), integer(1))
        res[[length(res)+1]] <- data.table(env = ev, basis = bs, engine = ENG, method = mth,
            k = kk, B = B, reject_rate = mean(v > 0), mean_rej = mean(v), max_rej = max(v),
            obs_rej = n_rej(z$p_obs, idx), lambda_obs = lam_o, lambda_surr = lam_s, gc_mode = gcm)
      }
    }
    cat(sprintf("  env%-3d %-9s done (gw reject rate %.2f, lambda_obs %.3f, lambda_surr %.3f, gc %s)\n",
                ev, bs, mean(gw > 0), lam_o, lam_s, gcm))
    rm(z, P); gc(verbose = FALSE)
  }
}
out <- rbindlist(res)
stopifnot(nrow(out) > 0)
fwrite(out, file.path(OUT, sprintf("filter_fdr_check_%s.csv", ENG)))

cat("\n=== FALSE-POSITIVE BEHAVIOUR UNDER SURROGATE ENVIRONMENTS ===\n")
cat(sprintf("  alpha = %.2f. With no true signal, reject_rate should be <= alpha.\n", ALPHA))
for (bs in unique(out$basis)) {
  cat(sprintf("\n  --- basis = %s ---\n", bs))
  s <- out[basis == bs, .(n_env = .N, reject_rate = mean(reject_rate),
                          mean_rej = mean(mean_rej), obs_rej = mean(obs_rej)),
           by = .(method, k)][order(k, method)]
  print(s)
  gwr <- out[basis == bs & method == "genome_wide"]
  cat("\n  paired vs genome-wide, same environment (positive = filter rejects MORE often):\n")
  for (kk in sort(unique(out[method != "genome_wide"]$k))) {
    for (mth in c("ld_w","MAF","random")) {
      a <- out[basis == bs & method == mth & k == kk]
      if (!nrow(a)) next
      m <- merge(a[, .(env, r = reject_rate)], gwr[, .(env, rg = reject_rate)], by = "env")
      d <- m$r - m$rg; nz <- d[d != 0]
      cat(sprintf("    k=%-6d %-7s  mean rate %.3f vs %.3f   diff %+.3f%s\n", kk, mth,
          mean(m$r), mean(m$rg), mean(d),
          if (length(nz)) sprintf("  (higher in %d/%d env)", sum(nz > 0), length(nz)) else "  (all tied)"))
    }
  }
}
cat(sprintf("\n  written: %s\n", file.path(OUT, sprintf("filter_fdr_check_%s.csv", ENG))))
