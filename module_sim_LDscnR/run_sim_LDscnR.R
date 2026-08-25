## =====================================================================
## module_sim_LDscnR / run_sim_LDscnR.R
##
## Benchmark the LD-aware C-score outlier-region method on the Nemo simulations,
## using only exported LDscnR functions -- the sim counterpart of
## module_sticklebacks_LDscnR/run_3sp_LDscnR.R. Unlike 3sp, the sims carry TRUTH
## (true_pos_QTN), so the analysis is a TP/FP benchmark, not region discovery.
##
## Story:
##   * the C-score (integrated over rho x q* x alpha) beats single-SNP association
##     on PR-AUC, and the >=2-SNP region filter is the key precision lever;
##   * the structure-aware null gives a data-driven (tau_C, l_min) whose PR sits on
##     the sweep's high-PR ridge -- i.e. calibration lands in the right place;
##   * the tau_C x l_min truth-PR ridge coincides with the null-quiet zone.
##
## ALWAYS replicate-average env1-5 (mean +- SE) -- env1 alone repeatedly flukes.
##
## Produces (into figures/ and results/):
##   Fig 1  fig1_pr_auc.png            PR curves: C-score vs single-SNP (per l_min), env-averaged
##   Fig 2  fig2_heatmaps.png          mean truth-PR + mean structured-null regions over tau_C x l_min
##   Fig 3  fig3_manhattan.png         C-score Manhattan with true-QTN crosses at the operating point
##   results/pr_auc_sim.csv, heatmap_sim.rds, operating_point_sim.csv
##
## RUN from the LDscnR-paper root:  Rscript module_sim_LDscnR/run_sim_LDscnR.R [V c]
## Reads the regen bundles (module_sim_LDscnR/regen_sim_data.R) from $SIM_DATA.
## =====================================================================

suppressMessages({
  library(data.table); library(ggplot2); library(patchwork); library(LDscnR)
})

## ---- 0. configuration -------------------------------------------------
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
ENVS <- as.character(ENVS)   # used as %s in filename patterns
OUTFIG <- "module_sim_LDscnR/figures"; OUTRES <- "module_sim_LDscnR/results"
CORES  <- as.integer(Sys.getenv("SIM_CORES", "1"))
QTAB_C <- if (CORES > 1L) 1L else 4L
if (CORES > 1L) { data.table::setDTthreads(1L); Sys.setenv(OMP_NUM_THREADS = "1") }  # avoid fork thread oversubscription

PAR <- list(
  alpha   = c(0.001, 0.01, 0.05, 0.1),          # alpha axis of the C-score (swept -> integrated)
  qstar   = seq(0, 0.95, by = 0.05),
  rho_ld  = 0.75, rho_d = 0.95, dcap = 5e5,      # clustering (sim-calibrated)
  B       = 100L, seed = 1L, basis = "spatial",  # structure-null: spatial autocorrelation (Nemo)
  fdr     = 0.05, lmin_q = 0.99, lmin_tau = 0.05,
  tau_grid  = seq(0.05, 1.0, by = 0.05),
  lmin_grid = c(1L, 2L, 4L, 8L),
  alpha_single = c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1)  # single-SNP BH cutoffs
)
gcta_grm <- function(X) { p <- colMeans(X)/2; k <- p>0 & p<1; X <- X[,k,drop=FALSE]; p <- p[k]
  Z <- sweep(sweep(X,2,2*p,"-"),2,sqrt(2*p*(1-p)),"/"); tcrossprod(Z)/ncol(Z) }

## ---- 1. pool one (V, c, env) cell ------------------------------------
## Pools 10 chromosome files (R1_.. prefixes) and keeps each file's SAVED GRM +
## marker set so the structured null runs per-file EMMAX on the identical kinship.
pool_cell <- function(env) {
  files <- list.files(SIM_DATA, full.names = TRUE,
                      pattern = sprintf("^adapt_%s_chr[0-9]+_V%s_c%s_env%s\\.rds$", TAG, V, CC, env))
  files <- files[order(as.integer(sub(".*_chr([0-9]+)_.*", "\\1", basename(files))))]
  maps <- gts <- ldws <- decs <- prep <- mk_i <- vector("list", length(files)); Yobs <- coords <- NULL
  for (i in seq_along(files)) {
    d <- readRDS(files[i]); m <- as.data.table(d$map)
    m[, `:=`(Chr = paste0("R", i, "_", Chr), marker = paste0("R", i, "_", marker))]
    G <- d$GTs; colnames(G) <- m$marker
    lw <- d$ld_ws; rownames(lw) <- m$marker
    ds <- as.data.table(d$LD_decay$decay_sum); ds[, Chr := paste0("R", i, "_", Chr)]
    prep[[i]] <- emmax_setup(G, d$GRM); mk_i[[i]] <- m$marker      # per-file engine on the SAVED GRM
    maps[[i]] <- m; gts[[i]] <- G; ldws[[i]] <- lw; decs[[i]] <- ds
    if (is.null(Yobs)) { Yobs <- d$env$env; coords <- cbind(d$env$x, d$env$y) }
  }
  map <- flag_true_qtns(rbindlist(maps, fill = TRUE))
  list(map = map, GTs = do.call(cbind, gts)[, map$marker],
       LDW = do.call(rbind, ldws)[map$marker, ], decay_sum = rbindlist(decs, fill = TRUE),
       prep = prep, mk_i = mk_i, Yobs = Yobs, coords = coords)
}

## per-file EMMAX for a phenotype, p-values concatenated in pooled marker order
emmax_pooled <- function(P, Yv) unlist(lapply(seq_along(P$prep), function(i)
  stats::setNames(emmax_fast(P$prep[[i]], Yv), P$mk_i[[i]])))[P$map$marker]

cscore <- function(pv, LDW) ld_cscore(pv, LDW, alpha = PAR$alpha, rho = colnames(LDW), qstar = PAR$qstar)

## ---- 2. structured null (pooled, spatial surrogate) ------------------
## One surrogate phenotype drawn from the shared spatial autocorrelation, run
## through every chromosome's own EMMAX, pooled -> one surrogate C. Reproduces the
## structure-driven false-positive rate; l_min-free (calibrate later). Cached.
null_bundle <- function(P, env) {
  cache <- file.path(OUTRES, sprintf("null_V%s_c%s_env%s.rds", V, CC, env))
  if (file.exists(cache)) return(readRDS(cache))
  n <- length(P$Yobs)
  Dm <- as.matrix(stats::dist(P$coords)); l <- stats::median(Dm[lower.tri(Dm)])
  eK <- eigen(exp(-0.5 * (Dm / l)^2), symmetric = TRUE); Lv <- pmax(eK$values, 0); Vk <- eK$vectors
  surr <- function() as.numeric(stats::resid(stats::lm((Vk %*% (sqrt(Lv) * stats::rnorm(n))) ~ P$Yobs)))
  sparseC <- function(pv) { C <- cscore(pv, P$LDW); C[C > 0] }
  C_obs <- cscore(emmax_pooled(P, P$Yobs), P$LDW)
  set.seed(PAR$seed); C_surr <- vector("list", PAR$B)
  for (b in seq_len(PAR$B)) C_surr[[b]] <- sparseC(emmax_pooled(P, surr()))
  bundle <- structure(list(C_obs = C_obs, C_surr = C_surr,
                           universe = unique(c(names(C_obs)[C_obs > 0], unlist(lapply(C_surr, names)))),
                           basis = PAR$basis, B = PAR$B), class = "ld_null")
  saveRDS(bundle, cache); bundle
}

## ---- 3. per-env compute: C-score, null, clustering, truth eval --------
per_env <- function(env) {
  P <- pool_cell(env)
  th   <- score_thresholds(P$decay_sum, rho_r2 = PAR$rho_ld, rho_d = PAR$rho_d, dmax_cap = PAR$dcap)
  nb   <- null_bundle(P, env)
  C    <- nb$C_obs                                     # observed C (coherent with the null engine)
  single <- lapply(PAR$alpha_single, function(al) P$map$marker[p.adjust(P$map$emx_p, "BH") < al])
  ## edges + truth-match over the lit-up universe ONLY (observed C>0, surrogate C>0,
  ## and the single-SNP sets) -- never all markers (a full r^2 graph over ~300k
  ## markers is ~8 GB/chromosome and needless: clustering only ever gates subsets).
  uni  <- unique(c(nb$universe, unlist(single)))
  edges <- ld_edges(uni, P$GTs, P$map[, .(marker, Chr, Pos)], P$decay_sum,
                    rho_ld = PAR$rho_ld, dcap = PAR$dcap)
  qtab <- qtn_ld_table(P$GTs, P$map, uni, 2e6, cores = QTAB_C)
  score <- function(mk) { reg <- if (length(mk)) ld_regions(mk, edges) else list()
    rbindlist(lapply(PAR$lmin_grid, function(lm) {
      ev <- evaluate_ors(reg[lengths(reg) >= lm], P$map, qtab, th$r2min, th$dmax)
      data.table(l_min = lm, recall = ev$Recall, precision = ev$Precision) })) }
  ## PR sweeps: C-score (tau_C) vs single-SNP (alpha)
  pr_C <- rbindlist(lapply(PAR$tau_grid, function(t) score(names(C)[C >= t])[, `:=`(curve = "C-score", knob = t)]))
  pr_S <- rbindlist(lapply(seq_along(PAR$alpha_single), function(i)
    score(single[[i]])[, `:=`(curve = "single-SNP", knob = PAR$alpha_single[i])]))
  ## structured-null operating point
  op_lmin <- calibrate_lmin(nb, edges, tau = PAR$lmin_tau, q = PAR$lmin_q)
  op_tau  <- calibrate_tauc(nb, edges, l_min = op_lmin, fdr = PAR$fdr, tau_grid = PAR$tau_grid)
  ## tau_C x l_min surfaces: truth-PR (observed) + mean null region count
  surf <- rbindlist(lapply(PAR$lmin_grid, function(lm) {
    f <- null_fdr(nb, edges, PAR$tau_grid, lm)
    pr <- vapply(PAR$tau_grid, function(t) { ev <- evaluate_ors(
      { r <- ld_regions(names(C)[C >= t], edges); r[lengths(r) >= lm] }, P$map, qtab, th$r2min, th$dmax); ev$PR }, 0)
    data.table(tau = PAR$tau_grid, lmin = lm, truthPR = pr, nullreg = f$n_null, fdr = f$fdr) }))
  list(env = env, pr = rbind(pr_C, pr_S), op = data.table(env = env, tau_C = op_tau, l_min = op_lmin),
       surf = surf, C = C, map = P$map, edges = edges, qtab = qtab, th = th)
}

cat(sprintf("V%s_c%s (%s): computing env %s ...\n", V, CC, TAG, paste(ENVS, collapse = ",")))
res <- if (CORES > 1L) parallel::mclapply(ENVS, per_env, mc.cores = CORES) else lapply(ENVS, per_env)

## ---- 4. Fig 1: env-averaged PR-AUC benchmark -------------------------
prA <- rbindlist(lapply(res, function(r) r$pr[, .(PR_AUC = pr_auc(recall, precision)), by = .(curve, l_min)][, env := r$env]))
prS <- prA[, .(PR_AUC = mean(PR_AUC), SE = sd(PR_AUC)/sqrt(.N)), by = .(curve, l_min)]
fwrite(prS, file.path(OUTRES, "pr_auc_sim.csv"))
g1 <- ggplot(prS, aes(factor(l_min), PR_AUC, fill = curve)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = PR_AUC - SE, ymax = PR_AUC + SE), position = position_dodge(0.8), width = 0.2) +
  labs(x = expression(l[min]), y = "PR-AUC (mean +/- SE, env1-5)",
       title = sprintf("V%s_c%s: C-score vs single-SNP", V, CC)) +
  theme_bw(base_size = 11) + theme(strip.background = element_blank())
ggsave(file.path(OUTFIG, "fig1_pr_auc.png"), g1, width = 7, height = 5, dpi = 130)

## ---- 5. Fig 2: env-averaged tau_C x l_min heatmaps -------------------
surf <- rbindlist(lapply(res, `[[`, "surf"))[, .(truthPR = mean(truthPR), nullreg = mean(nullreg)),
                                             by = .(tau, lmin)]
saveRDS(surf, file.path(OUTRES, "heatmap_sim.rds"))
zis <- grDevices::hcl.colors(100, "Zissou 1")
hm <- function(fill, title, ...) ggplot(surf, aes(factor(tau), factor(lmin), fill = .data[[fill]])) +
  geom_tile() + scale_fill_gradientn(colors = zis, ...) +
  labs(x = expression(tau[C]), y = expression(l[min]), title = title) +
  theme_minimal(base_size = 8) + theme(axis.text.x = element_text(angle = 90, vjust = .5))
g2 <- hm("truthPR", "mean truth PR (env1-5)") / hm("nullreg", "mean # structured-null regions (log1p)", trans = "log1p")
ggsave(file.path(OUTFIG, "fig2_heatmaps.png"), g2, width = 9, height = 6, dpi = 130)

## ---- 6. Fig 3: TP/FP Manhattan at the operating point (env1) ---------
r1 <- res[[1]]; op <- r1$op
regs <- if (!is.na(op$tau_C)) { rr <- ld_regions(names(r1$C)[r1$C >= op$tau_C], r1$edges); rr[lengths(rr) >= op$l_min] } else list()
qtn <- r1$map[true_pos_QTN == TRUE, marker]
g3 <- ld_manhattan(r1$map[, .(marker, Chr, Pos)], r1$C, value_label = "C-score", regions = regs, qtn = qtn,
                   title = sprintf("V%s_c%s env1: C-score at operating point (tau_C=%s, l_min=%s); + = true QTN",
                                   V, CC, format(op$tau_C), format(op$l_min)))
ggsave(file.path(OUTFIG, "fig3_manhattan.png"), g3, width = 12, height = 4, dpi = 130)

fwrite(rbindlist(lapply(res, `[[`, "op")), file.path(OUTRES, "operating_point_sim.csv"))
cat("wrote 3 figures + tables\n")
