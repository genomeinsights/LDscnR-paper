## module_sim/R/21_estimate.R
## AD-HOC estimation from the cache (20_build_cache.R) -- NO genotypes, NO re-pooling,
## NO re-clustering. Everything the C-score chain explores (any rho/q*/alpha grid,
## tau_C sweep, l_min, single-SNP alpha, PR-AUC) is a lookup on the cached
## {map, LDW, th, qtab, edges}. Demonstrated here by REPRODUCING the alpha-need PR-AUC
## (15c) in seconds instead of minutes. Source this file to get the helpers, or run it.
## Run from LDscnR-paper/:  Rscript module_sim/R/21_estimate.R [V c]

source("module_sim/R/_config.R")

## ---- ad-hoc estimators (cache in, numbers out) ------------------------------
load_cache <- function(V, cc, env)
  readRDS(file.path(dir_data, sprintf("cache_V%s_c%s_env%s.rds", V, cc, env)))

## PR-AUC of the C-score under ANY alpha set, sweeping tau_C, at one l_min -- ad-hoc.
prauc_cscore <- function(cache, pcol, alpha_c, tauc = seq(0.05, 1, 0.05), lmin = 2L) {
  Cv <- cscore_count(cache$map[[pcol]], cache$LDW, cache$grid$RHO, cache$grid$QSTAR, alpha_c)
  pc <- rbindlist(lapply(tauc, function(t) {
    mk <- names(Cv)[Cv >= t]; reg <- if (length(mk)) cluster_from_cache(mk, cache$edges) else list()
    ev <- evaluate_ORs_qtn(reg[lengths(reg) >= lmin], cache$map, cache$qtab, cache$th$r2min, cache$th$dmax)
    data.table(knob = t, recall = ev$Recall, precision = ev$Precision) }))
  pr_auc(pc$recall, pc$precision)
}
## structured-null observed side is also ad-hoc (C-score only); the surrogate loop
## still needs cached fast-EMMAX 'pre' (not stored here -- see module_sim/18).

if (sys.nframe() == 0L) {   # run as a script: reproduce the alpha-need result from cache
  a  <- commandArgs(trailingOnly = TRUE)
  V  <- if (length(a) >= 1) a[1] else "2"
  cc <- if (length(a) >= 2) a[2] else "1"; ENV <- as.character(1:5)
  ALPHA_FULL <- c(0.001, 0.01, 0.05, 0.1); ALPHA_05 <- 0.05
  t0 <- proc.time()[3]
  R <- rbindlist(lapply(ENV, function(env) { ca <- load_cache(V, cc, env)
    rbindlist(lapply(c("emx_p", "lfmm_p"), function(p) data.table(method = sub("_p$", "", p), env = env,
      full = prauc_cscore(ca, p, ALPHA_FULL), a05 = prauc_cscore(ca, p, ALPHA_05)))) }))
  se <- function(x){ x <- x[is.finite(x)]; if (length(x) < 2) NA_real_ else sd(x)/sqrt(length(x)) }
  agg <- R[, .(full = mean(full), full_se = se(full), a05 = mean(a05), a05_se = se(a05),
               delta = mean(full - a05), delta_se = se(full - a05)), by = method]
  cat(sprintf("=== alpha-need FROM CACHE (mean over 5 env, l_min=2), %.1fs ===\n", proc.time()[3] - t0))
  print(agg[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)])
  cat("\n(compare to 15c which recomputes from genotypes in minutes)\n")
}
