## =====================================================================
## module_C2 / R/07_sim_validation.R    [Question 7]
##
## Do TRUE simulated regions get higher operating-grid stability than FALSE ones?
##
## Feasibility (checked before writing anything heavy):
##   * the 10-chromosome Nemo bundles for V2_c1, env1-5 are already local under
##     module_sim_LDscnR/data/cache/ (GTs, map, LD_decay) -- nothing to re-simulate;
##   * truth is already carried in the map (`true_QTN`, `focal_QTN`,
##     `bp_to_focal_QTN`, `max_LD_with_QTN`), so region truth is assignable with the
##     package's own flag_true_qtns/qtn_ld_table/evaluate_ors -- the upstream
##     pipeline is NOT touched;
##   * the pooled spatial structured null (C_obs + B=100 C_surr) is cached in
##     module_sim_LDscnR/results/nullE_V2_c1_env{1..5}.rds -- no new nulls needed.
## So only the clustering + location-matched p-values are recomputed here, over a
## 2.6k-marker universe (cheap).
##
## Design. Per env: fix anchors at the sim operating point, score D / Q / S_G over
## the sim's own (tau_C, l_min) grid, label each anchor TRUE if it tags a true QTN
## (the package's LD-aware rule), and compare distributions + PR-AUC of the
## stability ranking against maxC and best q_R baselines. REPLICATE-AVERAGED over
## env1-5 (env1 alone is known to fluke on these sims).
##
## Rscript module_C2/R/07_sim_validation.R
## =====================================================================
source("module_C2/R/00_helpers.R")
suppressMessages(library(LDscnR))

SIM_CACHE <- "module_sim_LDscnR/data/cache"
SIM_RES   <- "module_sim_LDscnR/results"
V <- "2"; CC <- "1"; TAG <- "nobgs"; ENVS <- 1:5
## sim-canonical clustering + grid (module_sim_LDscnR/run_sim_LDscnR.R PAR)
S_RHO_LD <- 0.75; S_RHO_D <- 0.95; S_DCAP <- 5e5
S_TAUS <- seq(0.05, 1.00, by = 0.05); S_LMINS <- c(1L, 2L, 4L, 8L)
S_NCELL <- length(S_TAUS) * length(S_LMINS)
OP_TAU <- 0.05; OP_LMIN <- 2L                 # permissive anchor point: keeps false
                                              # regions in the set, which is the point
CACHE <- file.path(C2$CACHE, "sim_grid.rds")

pool_cell <- function(env) {
  files <- list.files(SIM_CACHE, full.names = TRUE,
    pattern = sprintf("^adapt_%s_chr[0-9]+_V%s_c%s_env%s\\.rds$", TAG, V, CC, env))
  files <- files[order(as.integer(sub(".*_chr([0-9]+)_.*", "\\1", basename(files))))]
  maps <- gts <- decs <- vector("list", length(files))
  for (i in seq_along(files)) {
    d <- readRDS(files[i]); m <- as.data.table(d$map)
    m[, `:=`(Chr = paste0("R", i, "_", Chr), marker = paste0("R", i, "_", marker))]
    G <- d$GTs; colnames(G) <- m$marker
    ds <- as.data.table(d$LD_decay$decay_sum); ds[, Chr := paste0("R", i, "_", Chr)]
    maps[[i]] <- m; gts[[i]] <- G; decs[[i]] <- ds
  }
  map <- flag_true_qtns(rbindlist(maps, fill = TRUE))
  list(map = map, GTs = do.call(cbind, gts)[, map$marker],
       decay_sum = rbindlist(decs, fill = TRUE), nfile = length(files))
}

if (file.exists(CACHE)) {
  SIM <- readRDS(CACHE); c2_msg("[0] loaded cached sim grid\n")
} else {
 SIM <- list()
 for (env in ENVS) {
  nf <- file.path(SIM_RES, sprintf("nullE_V%s_c%s_env%s.rds", V, CC, env))
  if (!file.exists(nf)) { c2_msg("[!] missing %s -- skipping env%d\n", nf, env); next }
  P  <- pool_cell(env); nb <- readRDS(nf); B <- length(nb$C_surr)
  th <- score_thresholds(P$decay_sum, rho_r2 = S_RHO_LD, rho_d = S_RHO_D, dmax_cap = S_DCAP)
  edges <- ld_edges(nb$universe, P$GTs, P$map[, .(marker, Chr, Pos)], P$decay_sum,
                    rho_ld = S_RHO_LD, dcap = S_DCAP)
  ## LD-aware truth table: which markers tag a true QTN
  qtab <- qtn_ld_table(P$GTs, P$map, nb$universe, 2e6, cores = 4L)
  Dl <- list(edges = edges, mpos = setNames(P$map$Pos, P$map$marker),
             mchr = setNames(as.integer(factor(P$map$Chr)), P$map$marker))
  ## ---- anchors at the sim operating point --------------------------
  ar <- ld_regions(names(nb$C_obs)[nb$C_obs >= OP_TAU], edges)
  ar <- ar[lengths(ar) >= OP_LMIN]
  ## truth: a region is TRUE if evaluate_ors counts it as tagging a true QTN
  ## the package's own per-region verdict (dedup-neutral): a region is TRUE when
  ## it is assigned a focal QTN that is a true positive. `.diagnose_ors` is
  ## internal, so fall back to the aggregate check if it moves.
  dg <- LDscnR:::.diagnose_ors(ar, P$map, qtab, th$r2min, th$dmax)
  istrue <- dg[order(CL_id)]$is_TP
  ## regions dropped by dedup (a duplicate call on an already-claimed true QTN)
  ## are neither TP nor FP in the package's scoring -- exclude them here too.
  keepr  <- !(dg[order(CL_id)]$dropped_by_dedup %in% TRUE &
              dg[order(CL_id)]$candidate_qtn_is_true_positive %in% TRUE)
  c2_msg("[%d] env%d: %d anchors (tau=%.2f l_min=%d) ; %d TRUE / %d FALSE ; B=%d\n",
         env, env, length(ar), OP_TAU, OP_LMIN, sum(istrue), sum(!istrue), B)
  ## ---- grid pass ----------------------------------------------------
  rows <- list()
  for (tau in S_TAUS) {
    oc <- c2_cluster(nb$C_obs, tau, Dl, keep_markers = TRUE)
    S  <- rbindlist(lapply(seq_len(B), function(b) {
            s <- c2_cluster(nb$C_surr[[b]], tau, Dl); if (nrow(s)) s[, b := b] else NULL }),
          fill = TRUE)
    if (!nrow(S)) S <- data.table(size=integer(), chr=integer(), lo=numeric(), hi=numeric(),
                                  score=numeric(), maxC=numeric(), b=integer())
    nO <- nrow(oc$tab)
    ms <- matrix(0, length(ar), max(nO, 1L))
    if (nO) for (j in seq_len(nO)) for (i in seq_along(ar))
      ms[i, j] <- c2_match_score(ar[[i]], oc$mk[[j]], "recover")
    for (lm in S_LMINS) {
      keep <- which(oc$tab$size >= lm)
      if (!length(keep)) { rows[[length(rows)+1L]] <- data.table(
        anchor = seq_along(ar), tau = tau, lmin = lm, detected = FALSE, sig = FALSE); next }
      O <- oc$tab[keep]; pq <- c2_emp_pq(O, S[size >= lm], B)
      sub <- ms[, keep, drop = FALSE]
      best <- vapply(seq_along(ar), function(i) { w <- which(sub[i, ] >= 0.5)
        if (!length(w)) NA_integer_ else w[which.max(sub[i, w])] }, integer(1))
      rows[[length(rows)+1L]] <- data.table(
        anchor = seq_along(ar), tau = tau, lmin = lm, detected = !is.na(best),
        sig = !is.na(best) & !is.na(pq$q_R[best]) & pq$q_R[best] < C2$FDR)
    }
  }
  G <- rbindlist(rows)
  sc <- G[, .(D = sum(detected)/S_NCELL, Q = sum(sig)/S_NCELL), by = anchor]
  sc[, `:=`(env = env, is_true = istrue[anchor], keep = keepr[anchor],
            Q_over_D = ifelse(D > 0, Q/D, 0), size = lengths(ar)[anchor],
            maxC = vapply(ar, function(r) max(nb$C_obs[r]), 0)[anchor])]
  sc <- sc[keep == TRUE]                      # dedup-neutral, as in evaluate_ors()
  SIM[[as.character(env)]] <- sc
 }
 saveRDS(SIM, CACHE)
}

A <- rbindlist(SIM)
if (!nrow(A)) stop("no sim envs available")

## ---- distributions true vs false ------------------------------------
ds <- A[, .(n = .N, D = mean(D), Q = mean(Q), QD = mean(Q_over_D)), by = .(env, is_true)]
c2_msg("\n[1] mean stability by truth, per env:\n"); print(ds)
avg <- ds[, .(n = sum(n), D = mean(D), Q = mean(Q), QD = mean(QD)), by = is_true]
c2_msg("\n[1] replicate-averaged over %d envs:\n", uniqueN(ds$env)); print(avg)
w <- A[, tryCatch(stats::wilcox.test(Q ~ is_true)$p.value, error = function(e) NA_real_), by = env]
c2_msg("[1] per-env Wilcoxon Q ~ truth p: %s\n", paste(signif(w$V1, 3), collapse = " "))

## ---- PR-AUC of each ranking -----------------------------------------
pr_auc_of <- function(score, truth) {
  o <- order(-score); t <- truth[o]
  tp <- cumsum(t); prec <- tp / seq_along(t); rec <- tp / sum(t)
  if (!sum(t)) return(NA_real_)
  sum(diff(c(0, rec)) * prec)
}
prt <- A[, .(PR_stability_Q = pr_auc_of(Q, is_true),
             PR_detection_D = pr_auc_of(D, is_true),
             PR_QoverD      = pr_auc_of(Q_over_D, is_true),
             PR_maxC        = pr_auc_of(maxC, is_true),
             PR_size        = pr_auc_of(size, is_true),
             base_rate      = mean(is_true)), by = env]
c2_msg("\n[2] PR-AUC by ranking (per env):\n"); print(prt)
mn <- prt[, lapply(.SD, mean, na.rm = TRUE), .SDcols = patterns("PR_|base")]
se <- prt[, lapply(.SD, function(x) stats::sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x)))),
          .SDcols = patterns("PR_|base")]
c2_msg("\n[2] replicate-averaged PR-AUC (mean +- SE over %d envs):\n", nrow(prt))
print(rbind(mean = round(unlist(mn), 3), SE = round(unlist(se), 3)))
fwrite(prt, file.path(C2$RES, "sim_validation_prauc.csv"))
fwrite(A,   file.path(C2$RES, "sim_anchor_scores.csv"))
c2_msg("\n[2] wrote results/sim_validation_prauc.csv + sim_anchor_scores.csv\n")

## ---- figure ----------------------------------------------------------
suppressMessages(library(ggplot2))
A[, truth := ifelse(is_true, "true", "false")]
f6 <- ggplot(A, aes(truth, Q, fill = truth)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6) +
  facet_wrap(~ paste0("env", env), nrow = 1) +
  scale_fill_manual(values = c(`false` = "grey75", `true` = "#F21A00"), guide = "none") +
  labs(x = NULL, y = "significance stability  Q",
       title = "Simulation truth: operating-grid stability of true vs false anchor regions",
       subtitle = sprintf("Nemo V2_c1, %d envs, spatial structured null B=100, anchors at tau_C=%.2f l_min=%d",
                          uniqueN(A$env), OP_TAU, OP_LMIN)) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 7))
ggsave(file.path(C2$FIG, "fig6_sim_truth.png"), f6, width = 8, height = 3.5, dpi = 170)
c2_msg("[3] wrote figures/fig6_sim_truth.png\n")
