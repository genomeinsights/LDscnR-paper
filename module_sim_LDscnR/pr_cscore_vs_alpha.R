## =====================================================================
## module_sim_LDscnR / pr_cscore_vs_alpha.R
##
## The paper's central comparison: does the consistency C-score beat plain
## single-SNP association at finding causal loci?
##
## Both arms are treated IDENTICALLY except for the marker-selection rule:
##
##   C-score     markers with C >= tau        (tau swept)
##   alpha       markers with BH q < alpha    (alpha swept)
##
## Everything downstream is shared -- the same LD-edge graph, the same clustering,
## the same l_min filter, the same evaluate_ors() truth matching with thresholds
## from the panel's own decay fit. So any difference is the selection rule, not
## the machinery around it.
##
## THE EDGE GRAPH MUST SPAN BOTH ARMS. The saved fits' edges cover only the
## C > 0 universe; BH can select markers outside it, and those would silently fail
## to cluster. The graph is therefore rebuilt over the union of both arms'
## candidate markers -- cheap now that ld_regions() uses igraph.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/pr_cscore_vs_alpha.R <scan_dir> [outdir]
## Env: PANEL_DIR, ENGINE (default emmax), LMINS (default 1,2,3,5),
##      MAX_TAU (default 40), P_COL (default emx_p)
## =====================================================================
suppressMessages({ library(data.table); library(LDscnR); library(ggplot2) })
`%||%` <- function(a, b) if (is.null(a)) b else a

a <- commandArgs(trailingOnly = TRUE)
if (!length(a)) stop("usage: pr_cscore_vs_alpha.R <scan_dir> [outdir]")
SCAN_DIR <- a[1]; OUT <- if (length(a) >= 2) a[2] else SCAN_DIR
PANEL_DIR <- Sys.getenv("PANEL_DIR", "/Volumes/Nemo/Nemo_sim/analysis_inputs")
ENGINE <- Sys.getenv("ENGINE", "emmax")
LMINS <- as.integer(strsplit(Sys.getenv("LMINS", "1,2,3,5"), ",")[[1]])
MAX_TAU <- as.integer(Sys.getenv("MAX_TAU", "40"))
## 1e5, not 5e5: the stage-2 partition in the bundles moved to
## distance_threshold = 1e5 on 2026-08-29, and the scoring geometry has to match
## it or regions are formed on one distance scale and the partition on another.
## It is also load-bearing rather than nominal: d(rho=0.95) is 636-845 kb on these
## cells, so a 5e5 cap BINDS and is what actually sets the window -- region
## formation is cap-governed, not decay-governed, at either value.
RHO_LD <- 0.75; RHO_D <- 0.95; DCAP <- 1e5
## alpha grid spanning the usable range on a log scale -- 7 points as in
## run_sim_LDscnR.R is too sparse to integrate a PR curve over
ALPHAS <- sort(unique(c(10^seq(-6, log10(0.5), length.out = 34), 0.05)))
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

scans <- list.files(SCAN_DIR, pattern = sprintf("^scan_.*_%s_.*[.]rds$", ENGINE), full.names = TRUE)
cellof <- function(f) sub(sprintf("^scan_(V[0-9.]+_c[0-9.]+_env[0-9]+)_%s_.*$", ENGINE), "\\1", basename(f))
cells <- unique(vapply(scans, cellof, character(1)))
cat(sprintf("[1] engine %s | %d cells | l_min %s | %d alpha, <=%d tau\n",
            ENGINE, length(cells), paste(LMINS, collapse=","), length(ALPHAS), MAX_TAU)); flush.console()

res <- list(); k <- 0L
for (cell in cells) {
  pf <- file.path(PANEL_DIR, sprintf("panel_%s.rds", cell))
  if (!file.exists(pf)) next
  panel <- readRDS(pf); map <- as.data.table(panel$map)
  if (!"true_pos_QTN" %in% names(map)) map <- flag_true_qtns(map)
  n_true <- sum(map$true_pos_QTN %in% TRUE); if (!n_true) next
  th <- score_thresholds(as.data.table(panel$decay_sum), rho_r2 = RHO_LD, rho_d = RHO_D, dmax_cap = DCAP)

  sc <- scans[vapply(scans, cellof, character(1)) == cell]
  f <- sc[grepl("env_orth", sc)][1]; if (is.na(f)) f <- sc[1]
  x <- readRDS(f); C <- x$null$C_obs

  ## the two arms' candidate pools
  p_obs <- map[[Sys.getenv("P_COL", "emx_p")]]
  q <- stats::p.adjust(p_obs, "BH")
  c_pool <- names(C)[which(C > 0)]
  a_pool <- map$marker[which(q < max(ALPHAS))]
  uni <- unique(c(c_pool, a_pool))
  if (!length(uni)) { cat(sprintf("  [skip] %s: no candidates in either arm\n", cell)); next }

  t0 <- Sys.time()
  edges <- ld_edges(uni, panel$GTs, map[, .(marker, Chr, Pos)], panel$decay_sum,
                    rho_ld = RHO_LD, dcap = DCAP)
  qtab <- qtn_ld_table(panel$GTs, map, uni, 2e6, cores = 1)
  cat(sprintf("  %s: %d QTN | C-pool %d, alpha-pool %d, union %d | setup %.1f s\n",
              cell, n_true, length(c_pool), length(a_pool), length(uni),
              as.numeric(Sys.time()-t0, units="secs"))); flush.console()

  tv <- sort(unique(C[C > 0]))
  if (length(tv) > MAX_TAU) tv <- as.numeric(stats::quantile(tv, seq(0, 1, length.out = MAX_TAU), type = 1))

  score_at <- function(mk, arm, knob) {
    if (!length(mk)) return(NULL)
    r_all <- ld_regions(mk, edges)
    rbindlist(lapply(LMINS, function(L) {
      r <- r_all[lengths(r_all) >= L]
      ev <- if (length(r)) evaluate_ors(r, map, qtab, th$r2min, th$dmax)
            else list(TP=0L, FP=0L, Precision=NA_real_, Recall=0)
      data.table(cell, arm, knob, l_min = L, n_regions = length(r),
                 TP = ev$TP, FP = ev$FP, precision = ev$Precision, recall = ev$Recall,
                 n_true = n_true)
    }))
  }
  for (tau in unique(tv)) { k <- k+1L; res[[k]] <- score_at(names(C)[which(C >= tau)], "C-score", tau) }
  for (al in ALPHAS)      { k <- k+1L; res[[k]] <- score_at(map$marker[which(q < al)], "alpha", al) }
}
pr <- rbindlist(res)
fwrite(pr, file.path(OUT, sprintf("pr_cscore_vs_alpha_%s.csv", ENGINE)))

auc <- pr[!is.na(precision), .(PR_AUC = pr_auc(recall, precision)), by = c("cell","arm","l_min")]
w <- dcast(auc, cell + l_min ~ arm, value.var = "PR_AUC")
setnames(w, c("C-score","alpha"), c("C_score","alpha"), skip_absent = TRUE)
cat(sprintf("\n=== PR-AUC, paired within genome (%d genomes) ===\n", uniqueN(auc$cell)))
agg <- w[, .(genomes = .N,
             C_score = sprintf("%.3f+-%.3f", mean(C_score, na.rm=TRUE), sd(C_score, na.rm=TRUE)/sqrt(.N)),
             alpha   = sprintf("%.3f+-%.3f", mean(alpha,   na.rm=TRUE), sd(alpha,   na.rm=TRUE)/sqrt(.N)),
             mean_diff = round(mean(C_score - alpha, na.rm=TRUE), 3),
             SE_diff = round(sd(C_score - alpha, na.rm=TRUE)/sqrt(.N), 3),
             C_wins = sum(C_score > alpha, na.rm=TRUE),
             p_paired = round(tryCatch(stats::t.test(w[l_min==.BY$l_min, C_score - alpha])$p.value,
                                       error=function(e) NA_real_), 4)), by = l_min]
setorderv(agg, "l_min"); print(agg)
fwrite(w, file.path(OUT, sprintf("pr_auc_paired_%s.csv", ENGINE)))
cat("\n[2] wrote pr_cscore_vs_alpha and pr_auc_paired csv\n")
