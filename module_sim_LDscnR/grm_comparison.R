## =====================================================================
## module_sim_LDscnR / grm_comparison.R
##
## How the EMMAX GRM (its marker set) affects the C-score benchmark on the Nemo
## sims -- the sim side of the GRM question we settled on 3sp. Compares, per (V,c)
## cell and replicate-averaged over env1-5:
##   * "chain"  -- the complexity-reduction pruned set (ld_prune_and_eMLG, ~45%),
##   * "ld_w<b" -- LD-INDEPENDENT markers below background LD (the 3sp-style rule).
## Reports, per GRM: gif (calibration), pooled maxC, and PR-AUC vs truth on the
## ADAPTIVE tau grid (the sweep = the distinct observed C values, so it covers
## exactly where regions appear -- a fixed [0.05,1] grid is meaningless when a
## cell's C-score maxes out low). Also draws the ld_w_095 Manhattan that explains
## the mechanism.
##
## FINDINGS (V2_c1, env1-5; see README):
##   * gif: chain ~1.02-1.08 (calibrated); ld_w<b ~0.91-0.96 (mildly DEFLATED).
##   * PR-AUC: chain beats ld_w<b by ~0.06-0.10 at l_min=1-2, but the gap CLOSES
##     by l_min>=4-8 (l_min filter compensates for the deflation).
##   * ld_w<b over-corrects on the sims because ~96% of markers (the QTN included,
##     since sim QTN are LOW-LD) fall below b -> the GRM keeps them. But they are
##     low-LD, so the deflation is gentle, not catastrophic. This is the OPPOSITE
##     of 3sp, where the signal is HIGH-LD and ld_w<b correctly excludes it.
##   * Transferable lesson: on the sims the complexity chain (gif~=1) is the better
##     GRM; gif ~= 1 is a good check against GLOBAL deflation (the sim failure mode)
##     but does NOT catch REGIONAL signal absorption (the 3sp failure mode, where
##     the chain has gif 1.061 yet collapses C). See module_sticklebacks_LDscnR.
##
## Run from the LDscnR-paper root (reads regen bundles from $SIM_DATA):
##   Rscript module_sim_LDscnR/grm_comparison.R [V c]
## =====================================================================

suppressMessages({ library(data.table); library(ggplot2); library(LDscnR) })

## ---- 0. config -------------------------------------------------------
a  <- commandArgs(trailingOnly = TRUE)
V  <- if (length(a) >= 1) a[1] else "2"
CC <- if (length(a) >= 2) a[2] else "1"
SIM_DATA <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_nobgs")
TAG    <- "nobgs"

## ENVS is derived from what is on disk, not hardcoded, so adding replicates does
## not mean editing this file. Only COMPLETE cells are used: pooling needs all
## CHR_N chromosomes of an env, and a partial cell would silently pool a smaller
## genome. Override with SIM_ENVS=1,2,3.
CHR_N <- 10L
discover_envs <- function(dir, tag, V, CC, chr_n = CHR_N) {
  fs <- list.files(dir, pattern = sprintf("^adapt_%s_chr[0-9]+_V%s_c%s_env[0-9]+[.]rds$", tag, V, CC))
  if (!length(fs)) stop("no bundles matching V", V, "_c", CC, " (", tag, ") in ", dir)
  e <- as.integer(sub(".*_env([0-9]+)[.]rds$", "\\1", fs)); tab <- table(e)
  full <- as.integer(names(tab)[tab == chr_n]); short <- setdiff(as.integer(names(tab)), full)
  if (length(short)) message(sprintf("  [envs] incomplete cells skipped: %s",
    paste(sprintf("env%d (%d/%d chr)", short, as.integer(tab[as.character(short)]), chr_n), collapse = ", ")))
  sort(full)
}
ENVS <- { .e <- Sys.getenv("SIM_ENVS", "")
          if (nzchar(.e)) as.integer(strsplit(.e, ",")[[1]]) else discover_envs(SIM_DATA, TAG, V, CC) }
if (!length(ENVS)) stop("no complete env cells in ", SIM_DATA)
message(sprintf("  [envs] using %d env cell(s): %s", length(ENVS), paste(ENVS, collapse = ",")))
OUTFIG <- "module_sim_LDscnR/figures"; OUTRES <- "module_sim_LDscnR/results"
CORES  <- as.integer(Sys.getenv("SIM_CORES", "1"))          # env-level parallelism over ENVS
QTAB_C <- if (CORES > 1L) 1L else 4L
if (CORES > 1L) { data.table::setDTthreads(1L); Sys.setenv(OMP_NUM_THREADS = "1") }  # avoid fork thread oversubscription                        # avoid nested threads under mclapply
PAR <- list(qstar = seq(0, 0.95, by = 0.05), alpha = c(0.001, 0.01, 0.05, 0.1),
            lmin = c(1L, 2L, 4L, 8L), rho_ld = 0.75, rho_d = 0.95, dcap = 5e5, max_tau = 50L)

gcta_grm <- function(X) { p <- colMeans(X) / 2; k <- p > 0 & p < 1; X <- X[, k, drop = FALSE]; p <- p[k]
  Z <- sweep(sweep(X, 2, 2 * p, "-"), 2, sqrt(2 * p * (1 - p)), "/"); tcrossprod(Z) / ncol(Z) }
## per-file EMMAX from a GRM marker set (+ genomic control if gif > 1.1)
emx_p <- function(G, mk, y) { K <- gcta_grm(G[, mk, drop = FALSE]); pv <- emmax_fast(emmax_setup(G, K), y)
  n <- nrow(G); Fv <- stats::qf(pv, 1, n - 2, lower.tail = FALSE)
  gif <- stats::median(Fv) / stats::qf(0.5, 1, n - 2, lower.tail = FALSE)
  if (gif > 1.1) { Fv <- Fv / gif; pv <- stats::pf(Fv, 1, n - 2, lower.tail = FALSE) }
  list(p = pv, gif = gif, n = length(mk)) }
## The regen bundles renamed this field: older ones populate `pruned_markers`
## and leave `grm_markers` empty, current ones do the reverse. Accept either, and
## fail loudly rather than returning an empty set -- an empty marker set would
## build a GRM from nothing and report it as a result.
grm_set <- function(d) {
  mk <- if (length(d$grm_markers)) d$grm_markers else d$pruned_markers
  if (!length(mk)) stop("bundle carries neither grm_markers nor pruned_markers")
  mk
}
## the GRM marker sets under comparison (extend here to add variants)
grm_markers <- list(
  chain  = function(d, G, lw, b) grm_set(d),                     # complexity-reduction pruned set
  ldw_b  = function(d, G, lw, b) colnames(G)[which(lw < b)])      # ld_w_095 < background LD

## ---- 1. per-env: pool, per-GRM EMMAX, C-score, PR-AUC ----------------
per_env <- function(env) {
  files <- list.files(SIM_DATA, full.names = TRUE,
                      pattern = sprintf("^adapt_%s_chr[0-9]+_V%s_c%s_env%s\\.rds$", TAG, V, CC, env))
  files <- files[order(as.integer(sub(".*_chr([0-9]+)_.*", "\\1", basename(files))))]
  maps <- gts <- ldws <- decs <- vector("list", length(files)); gifs <- vector("list", length(grm_markers))
  for (i in seq_along(files)) {
    d <- readRDS(files[i]); m <- as.data.table(d$map); G <- d$GTs; y <- d$env$env
    lwc <- if ("rho_0.95" %in% colnames(d$ld_ws)) "rho_0.95" else "0.95"; lw <- d$ld_ws[, lwc]
    b <- stats::median(as.data.table(d$LD_decay$decay_sum)$b)
    for (g in names(grm_markers)) { r <- emx_p(G, grm_markers[[g]](d, G, lw, b), y)
      m[, (paste0("p_", g)) := r$p]; gifs[[g]] <- c(gifs[[g]], r$gif) }
    m[, `:=`(Chr = paste0("R", i, "_", Chr), marker = paste0("R", i, "_", marker))]
    colnames(G) <- paste0("R", i, "_", colnames(G))
    lwm <- d$ld_ws; rownames(lwm) <- m$marker; ds <- as.data.table(d$LD_decay$decay_sum); ds[, Chr := paste0("R", i, "_", Chr)]
    maps[[i]] <- m; gts[[i]] <- G; ldws[[i]] <- lwm; decs[[i]] <- ds
  }
  map <- flag_true_qtns(rbindlist(maps, fill = TRUE)); GTs <- do.call(cbind, gts)[, map$marker]
  LDW <- do.call(rbind, ldws)[map$marker, ]; decay_sum <- rbindlist(decs, fill = TRUE)
  th  <- score_thresholds(decay_sum, rho_r2 = PAR$rho_ld, rho_d = PAR$rho_d, dmax_cap = PAR$dcap)
  Cof <- function(g) ld_cscore(map[[paste0("p_", g)]], LDW, alpha = PAR$alpha, rho = colnames(LDW), qstar = PAR$qstar)
  Cs  <- lapply(stats::setNames(names(grm_markers), names(grm_markers)), Cof)
  ## edges + truth-match computed ONCE over the C>0 universe (no per-cell LD work)
  uni   <- unique(unlist(lapply(Cs, function(C) names(C)[C > 0])))
  edges <- ld_edges(uni, GTs, map[, .(marker, Chr, Pos)], decay_sum, rho_ld = PAR$rho_ld, dcap = PAR$dcap)
  qtab  <- qtn_ld_table(GTs, map, uni, 2e6, cores = QTAB_C)
  ## ADAPTIVE tau grid = the distinct observed C values (both GRMs), thinned to max_tau
  tv <- sort(unique(unlist(lapply(Cs, function(C) C[C > 0]))))
  if (length(tv) > PAR$max_tau) tv <- as.numeric(stats::quantile(tv, seq(0, 1, length.out = PAR$max_tau), type = 1))
  TAUC <- unique(tv)
  prauc <- function(C) { sc <- function(mk) { reg <- if (length(mk)) ld_regions(mk, edges) else list()
      rbindlist(lapply(PAR$lmin, function(lm) { ev <- evaluate_ors(reg[lengths(reg) >= lm], map, qtab, th$r2min, th$dmax)
        data.table(l_min = lm, recall = ev$Recall, precision = ev$Precision) })) }
    rbindlist(lapply(TAUC, function(t) sc(names(C)[C >= t])))[, .(PR_AUC = pr_auc(recall, precision)), by = l_min] }
  cat(sprintf("env%d: %s ; |tau|=%d\n", env,
      paste(sprintf("%s gif=%.3f maxC=%.3f", names(Cs), sapply(names(Cs), function(g) mean(gifs[[g]])),
                    sapply(Cs, max)), collapse = " | "), length(TAUC))); utils::flush.console()
  rbindlist(lapply(names(Cs), function(g) prauc(Cs[[g]])[, `:=`(method = g, env = env,
    gif = mean(gifs[[g]]))]))
}

cat(sprintf("V%s_c%s (%s): GRM comparison over env %s\n", V, CC, TAG, paste(ENVS, collapse = ",")))
auc <- rbindlist(if (CORES > 1L) parallel::mclapply(ENVS, per_env, mc.cores = CORES) else lapply(ENVS, per_env))
summ <- auc[, .(gif = round(mean(gif), 3), PR_AUC = round(mean(PR_AUC, na.rm = TRUE), 3),
                SE = round(stats::sd(PR_AUC, na.rm = TRUE) / sqrt(sum(!is.na(PR_AUC))), 3)),
            by = .(method, l_min)][order(method, l_min)]
fwrite(summ, file.path(OUTRES, "grm_comparison_prauc.csv"))
cat("\n=== PR-AUC (mean +/- SE, adaptive tau) ===\n"); print(summ)

## ---- 2. ld_w_095 Manhattan (mechanism figure) ------------------------
## why the GRMs behave differently: sim QTN sit BELOW b (kept by ld_w<b), the
## high-ld_w blocks are neutral structure (R5 inversion, R9). Env 1.
mmaps <- list(); bv <- c()
for (i in seq_along(list.files(SIM_DATA, pattern = sprintf("^adapt_%s_chr[0-9]+_V%s_c%s_env1\\.rds$", TAG, V, CC)))) {
  d <- readRDS(sprintf("%s/adapt_%s_chr%d_V%s_c%s_env1.rds", SIM_DATA, TAG, i, V, CC))
  m <- as.data.table(d$map); m[, `:=`(Chr = paste0("R", i), marker = paste0("R", i, "_", marker))]
  mmaps[[i]] <- m[, .(marker, Chr, Pos, ld_w_095, type)]; bv <- c(bv, stats::median(as.data.table(d$LD_decay$decay_sum)$b))
}
mm <- rbindlist(mmaps); b <- mean(bv); mm[, Chr := factor(Chr, levels = paste0("R", seq_along(mmaps)))]
g <- ld_manhattan(mm[, .(marker, Chr, Pos)], stats::setNames(mm$ld_w_095, mm$marker),
                  value_label = "ld_w (rho=0.95)", qtn = mm[type == "QTN", marker], hline = b,
                  title = sprintf("V%s_c%s env1: ld_w_095 (b=%.3f dashed; + = QTN) -- ld_w<b keeps %.0f%%",
                                  V, CC, b, 100 * mean(mm$ld_w_095 < b, na.rm = TRUE)), point_size = 0.7)
ggsave(file.path(OUTFIG, "grm_ldw095_manhattan.png"), g, width = 13, height = 4.2, dpi = 135)
cat(sprintf("\nwrote %s and %s\n", file.path(OUTRES, "grm_comparison_prauc.csv"),
            file.path(OUTFIG, "grm_ldw095_manhattan.png")))
