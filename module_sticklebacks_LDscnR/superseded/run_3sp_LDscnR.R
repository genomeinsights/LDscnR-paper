## =====================================================================
## module_sticklebacks_LDscnR / run_3sp_LDscnR.R
##
## Clean re-implementation of the 3-spine stickleback (3sp) LD-aware
## outlier-region analysis using ONLY exported LDscnR functions.
##
## Story (saturated genome, no neutral floor):
##   * single-SNP EMMAX association is dead on 3sp (nothing FDR-significant),
##   * yet the LD-aware C-score resolves the known marine<->freshwater
##     architecture (Chr1 inversion, Chr4 Eda, Chr7, Chr20, Chr17, ...),
##   * the structure-aware null certifies these regions (FDR ~ 0 everywhere:
##     there is no null-danger corner -- the challenge is ranking many real
##     regions, not signal-vs-structure),
##   * so regions are ranked threshold-free by cross-parameter stability.
##
## Produces (into figures/ and results/):
##   Fig 1  fig1_manhattan_C_vs_q.png   C-score vs -log10(q), regions highlighted
##   Fig 2  fig2_heatmaps_tau_lmin.png  observed / null / FDR over tau_C x l_min
##   Fig 3  fig3_region_ranking.png     regions recoloured by stability & persist_tau
##   results/region_stability_3sp.csv   the ranked region table
##   results/operating_point_3sp.csv    the zero-null max-cluster region set
##
## RUN from the LDscnR-paper root:  Rscript module_sticklebacks_LDscnR/run_3sp_LDscnR.R
## =====================================================================

suppressMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(LDscnR)                      # devtools::load_all() during development
})

## ---- 0. configuration -------------------------------------------------
## single self-contained bundle from regen_3sp_data.R (GTs, map, eco, ld_ws,
## LD_decay, complexity_reduction, GRM, emx, lfmm). Regenerate it first.
BUNDLE <- "module_sticklebacks_LDscnR/data/3sp_LDscnR_data.rds"
OUTFIG <- "module_sticklebacks_LDscnR/figures"
OUTRES <- "module_sticklebacks_LDscnR/results"

PAR <- list(
  alpha   = 0.05,                       # fixed within-candidate FDR (C-score)
  qstar   = seq(0, 0.95, by = 0.05),    # C-score stringency grid
  rho_ld  = 0.75,                       # LD-decay-relative r^2 link (matches the sims)
  dcap    = 1e5,                        # hard 500 kb gap-split
  B       = 100L,                       # structure-null surrogates
  seed    = 1L,
  fdr      = 0.05,                       # target region-level FDR for calibrate_tauc
  lmin_q   = 0.99, lmin_tau = 0.05,      # calibrate_lmin quantile + reference tau
  tau_grid  = seq(0.02, 0.50, by = 0.02),        # tau_C axis (heatmaps + stability)
  lmin_grid = c(1L, 2L, 3L, 5L, 10L, 15L, 20L)   # l_min axis
)

## ---- 1. load the regenerated bundle ----------------------------------
d    <- readRDS(BUNDLE)
GTs  <- d$GTs                                                # individuals x SNPs (MAF>0.1)
map  <- as.data.table(d$map)[, .(marker, Chr, Pos)]
chr_lev <- paste0("Chr", 1:21); map[, Chr := factor(Chr, levels = chr_lev)]
LDW  <- d$ld_ws[map$marker, ]                               # SNP x rho ld_w matrix
decs <- as.data.table(d$LD_decay$decay_sum)                 # per-chr LD decay (b, c, a_pred)
K    <- d$GRM                                               # shared kinship (sim-consistent pruning)
eco  <- d$eco                                               # phenotype (Marine = 1)
cat(sprintf("[1] %d individuals x %d markers ; Marine = %d ; GRM from %d LD-independent markers\n",
            nrow(GTs), ncol(GTs), sum(eco), length(d$grm_markers)))

## ---- 2. EMMAX scan, C-score, structure-aware null --------------------
## Observed and null share ONE engine (fast EMMAX with the saved GRM); the bundle
## carries the observed C-score (C_obs) and the B sparse surrogate C-scores.
prep <- emmax_setup(GTs, K)
pval <- emmax_fast(prep, eco); names(pval) <- map$marker      # per-SNP association p
q    <- p.adjust(pval, "BH"); nlq <- setNames(-log10(q), map$marker)
cat(sprintf("[2] single-SNP EMMAX: min q = %.3f ; SNPs with q<=0.05 = %d\n", min(q), sum(q <= 0.05)))

## the B=100 null is the expensive step -- cache it (with C + edges) so the
## downstream calibration/figures can be re-run cheaply. Delete the cache to force
## a recompute (e.g. after changing the bundle, GRM, or null parameters).
cache <- file.path(OUTRES, "null_cache_3sp.rds")
if (file.exists(cache)) {
  cc <- readRDS(cache); null <- cc$null; edges <- cc$edges
  cat("[2] loaded cached structured null (delete ", cache, " to recompute)\n", sep = "")
} else {
  null  <- structured_null(eco, GTs, K, LDW, basis = "genetic", B = PAR$B,
                           alpha = PAR$alpha, qstar = PAR$qstar, prep = prep, seed = PAR$seed)
  edges <- ld_edges(null$universe, GTs, map, decs, rho_ld = PAR$rho_ld, dcap = PAR$dcap)
  saveRDS(list(null = null, edges = edges), cache)
}
C <- null$C_obs                                              # full per-SNP C-score
cat(sprintf("[2] C-score: %d SNPs with C>0, %d with C>=0.3\n", sum(C > 0), sum(C >= 0.3)))

## ---- 3. null calibration + operating point ---------------------------
## grid of observed / null region counts + region-level FDR over tau_C x l_min
grid <- rbindlist(lapply(PAR$lmin_grid, function(L) {
  f <- null_fdr(null, edges, PAR$tau_grid, L)
  data.table(tau = f$tau_C, lmin = L, obs = f$n_obs, nullreg = f$n_null, fdr = f$fdr)
}))
fwrite(grid, file.path(OUTRES, "fdr_grid_3sp.csv"))

## operating point: the package's data-driven calibration -- l_min from the null's
## max-region-size (calibrate_lmin) and the smallest tau_C whose region-level FDR is
## at/below the target (calibrate_tauc). Robust where no cell is exactly null-silent.
op_lmin <- calibrate_lmin(null, edges, tau = PAR$lmin_tau, q = PAR$lmin_q)
op_tau  <- calibrate_tauc(null, edges, l_min = op_lmin, fdr = PAR$fdr, tau_grid = PAR$tau_grid)
if (is.na(op_tau)) {                                # fallback: loosest cell meeting FDR at any l_min
  ok <- grid[fdr <= PAR$fdr & obs > 0][order(-obs)]
  if (nrow(ok)) { op_tau <- ok$tau[1]; op_lmin <- ok$lmin[1] }
}
cat(sprintf("[3] operating point: tau_C = %s, l_min = %s (FDR target %.2f)\n",
            format(op_tau), format(op_lmin), PAR$fdr))

## the region set at that operating point
regs <- if (!is.na(op_tau)) { r <- ld_regions(names(C)[C >= op_tau], edges); r[lengths(r) >= op_lmin] } else list()
POS <- setNames(map$Pos, map$marker); CH <- setNames(as.character(map$Chr), map$marker)
op_tab <- if (length(regs)) rbindlist(lapply(seq_along(regs), function(i) { m <- regs[[i]]
  data.table(Chr = CH[m[1]], start = min(POS[m]), end = max(POS[m]),
             span_kb = round((max(POS[m]) - min(POS[m])) / 1e3, 1),
             n_snp = length(m), maxC = round(max(C[m]), 3),
             min_q = signif(min(q[m]), 2)) }))[order(Chr, start)] else data.table()
fwrite(op_tab, file.path(OUTRES, "operating_point_3sp.csv"))
cat(sprintf("[3] %d regions at the operating point\n", length(regs)))

## ---- 4. threshold-free region ranking (cross-parameter stability) ----
if (length(regs)) {
  stab <- ld_region_stability(C, edges = edges, regions = regs,
                              tau_grid = PAR$tau_grid, l_min_grid = PAR$lmin_grid,
                              base_tau = op_tau, base_lmin = op_lmin, map = map)
  fwrite(as.data.table(stab), file.path(OUTRES, "region_stability_3sp.csv"))
  print(stab)
} else stab <- data.table()

## ---- 5. figures ------------------------------------------------------
## Fig 1 : C-score vs -log10(q), operating-point regions highlighted (same colours)
g1a <- ld_manhattan(map, C,   value_label = "C-score", regions = regs,
                    title = "3sp EMMAX C-score (LD-aware regions highlighted)", point_size = 1.1)
g1b <- ld_manhattan(map, nlq, value_label = expression(-log[10](q)), regions = regs,
                    title = "3sp EMMAX single-SNP association (same regions; nothing survives q<=0.05)",
                    point_size = 1.1)
ggsave(file.path(OUTFIG, "fig1_manhattan_C_vs_q.png"), g1a / g1b, width = 12, height = 7, dpi = 120)

## Fig 2 : observed / null / FDR heatmaps over tau_C x l_min
zis <- grDevices::hcl.colors(100, "Zissou 1")
hm <- function(fill, title, ...) ggplot(grid, aes(factor(tau), factor(lmin), fill = .data[[fill]])) +
  geom_tile() + scale_fill_gradientn(colors = zis, ...) +
  labs(x = expression(tau[C]), y = expression(l[min]), title = title) +
  theme_minimal(base_size = 8) + theme(axis.text.x = element_text(angle = 90, vjust = .5))
g2 <- hm("obs", "3sp EMMAX: # observed outlier regions (log1p)", trans = "log1p") /
      hm("nullreg", "mean # regions under structured NULL (log1p)", trans = "log1p") /
      hm("fdr", "region-level FDR (null / observed)", limits = c(0, 1))
ggsave(file.path(OUTFIG, "fig2_heatmaps_tau_lmin.png"), g2, width = 9, height = 9, dpi = 115)

## Fig 3 : regions recoloured by the two threshold-free ranks
##   persist_tau  = continuous LD-region consistency
##   stability    = fraction of tau_C x l_min cells the region survives
col_persist <- col_stab <- setNames(rep(NA_real_, nrow(map)), map$marker)
for (i in seq_along(regs)) {
  col_persist[regs[[i]]] <- stab$persist_tau[stab$region == i]
  col_stab[regs[[i]]]    <- stab$stability[stab$region == i]
}
g3a <- ld_manhattan(map, C, value_label = "C-score", colour = col_persist,
                    colour_label = "persist_tau", title = "coloured by persist_tau (LD-region consistency)",
                    point_size = 1.2)
g3b <- ld_manhattan(map, C, value_label = "C-score", colour = col_stab,
                    colour_label = "stability", title = "coloured by grid-survival stability (rank only)",
                    point_size = 1.2)
ggsave(file.path(OUTFIG, "fig3_region_ranking.png"), g3a / g3b, width = 12, height = 7, dpi = 120)

cat("\n[5] wrote 3 figures to ", OUTFIG, " and tables to ", OUTRES, "\n", sep = "")
