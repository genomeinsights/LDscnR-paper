## =====================================================================
## module_sticklebacks_LDscnR / null_fdr_landscape.R
##
## The tau_C x l_min operating landscape (observed # regions / mean # regions under
## the null / region-level FDR = null/observed) computed for ALL FOUR structure
## nulls, not just the genetic MVN. The observed panel is null-independent (C_obs is
## the same EMMAX C-score); the null-region-count and FDR panels are null-specific.
## Shows how the FDR-safe operating region differs by null: the genetic MVN is silent
## (FDR ~ 0 everywhere), the permutation nulls fire in the small-cluster corner, and
## the spatial null is broadly inflated.
##
## Run from the LDscnR-paper root:
##   Rscript module_sticklebacks_LDscnR/null_fdr_landscape.R
## Writes figures/null_fdr_landscape_4nulls.png + results/null_fdr_landscape.csv
## =====================================================================
suppressMessages({ library(data.table); library(LDscnR); library(ggplot2); library(patchwork) })

BND    <- "module_sticklebacks_LDscnR/data/3sp_LDscnR_data.rds"
RES    <- "module_sticklebacks_LDscnR/results"; OUTFIG <- "module_sticklebacks_LDscnR/figures"
NULLS  <- c("genetic (MVN)"      = "null_uncapped_3sp",  "global perm."        = "null_popperm_3sp",
            "regional perm."     = "null_regionperm_3sp","spatial (MVN kernel)"= "null_spatial_3sp")
TAUS   <- seq(0.02, 0.5, by = 0.02); LMINS <- c(1,2,3,5,10,15,20)
RHO_LD <- 0.60; DCAP <- 1e5; OP_TAU <- 0.05; OP_LMIN <- 3L
ZISSOU <- c("#3B9AB2","#78B7C5","#EBCC2A","#E1AF00","#F21A00")

d  <- readRDS(BND); map <- as.data.table(d$map); decay <- as.data.table(d$LD_decay$decay_sum)
nl <- lapply(NULLS, function(f) readRDS(file.path(RES, paste0(f, ".rds"))))
C_obs <- nl[[1]]$C_obs; BCAP <- 100L      # cap surrogates used for the null-mean (stable, faster)
cat(sprintf("[1] grid %d tau x %d l_min x %d nulls ; per-null edge graphs\n",
            length(TAUS), length(LMINS), length(nl))); flush.console()

## region counts by l_min for a C-vector at threshold tau, given a null's edge graph
counts <- function(C, tau, edges) { mk <- names(C)[C >= tau]
  if (!length(mk)) return(rep(0L, length(LMINS)))
  s <- lengths(ld_regions(mk, edges)); vapply(LMINS, function(l) sum(s >= l), integer(1)) }

## Build edges ON EACH NULL'S OWN UNIVERSE (contains the observed C>0 markers, so the
## observed clustering -- an induced subgraph -- is identical across nulls, but the graph
## is far smaller than the union). Compute observed + mean-null region counts per null.
grid <- rbindlist(lapply(names(nl), function(nm) { n <- nl[[nm]]
  e_nm <- ld_edges(n$universe, d$GTs, map[, .(marker, Chr, Pos)], decay, rho_ld = RHO_LD, dcap = DCAP)
  surr <- n$C_surr[seq_len(min(BCAP, length(n$C_surr)))]
  cat(sprintf("   null: %-20s universe=%5d ; %d surrogates\n", nm, length(n$universe), length(surr))); flush.console()
  rbindlist(lapply(TAUS, function(t) {
    obs_c <- counts(C_obs, t, e_nm)
    M <- vapply(surr, function(C) counts(C, t, e_nm), integer(length(LMINS)))
    data.table(null = nm, tau = t, lmin = LMINS, obs = obs_c, nullreg = rowMeans(M)) })) }))
grid[, fdr := pmin(nullreg / pmax(obs, 1), 1)]
grid[, null := factor(null, levels = names(nl))]
obs <- unique(grid[null == names(nl)[1], .(tau, lmin, obs)])   # shared observed panel
fwrite(grid, file.path(RES, "null_fdr_landscape.csv"))

## ---- plot -------------------------------------------------------------------------
hm <- function(DT, fill, title, trans = "identity", limits = NULL) {
  ggplot(DT, aes(factor(tau), factor(lmin, levels = LMINS), fill = .data[[fill]])) +
    geom_tile() +
    geom_point(data = data.table(tau = OP_TAU, lmin = OP_LMIN), inherit.aes = FALSE,
               aes(factor(tau), factor(lmin, levels = LMINS)), shape = 4, size = 2, stroke = 1.1) +
    scale_fill_gradientn(colours = ZISSOU, name = NULL, trans = trans, limits = limits) +
    labs(x = expression(tau[C]), y = expression(l[min]), title = title) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, size = 4.5),
          axis.text.y = element_text(size = 6), plot.title = element_text(size = 8, face = "bold"),
          legend.key.width = unit(0.4, "lines"))
}
g_obs  <- hm(obs, "obs", "Observed # regions (shared across nulls)", trans = "log1p")
g_null <- lapply(names(nl), function(nm) hm(grid[null == nm], "nullreg",
            sprintf("mean # null regions: %s", nm), trans = "log1p"))
g_fdr  <- lapply(names(nl), function(nm) hm(grid[null == nm], "fdr",
            sprintf("region-level FDR: %s", nm), limits = c(0, 1)))

design <- g_obs /
  (g_null[[1]] | g_null[[2]] | g_null[[3]] | g_null[[4]]) /
  (g_fdr[[1]]  | g_fdr[[2]]  | g_fdr[[3]]  | g_fdr[[4]]) +
  plot_layout(heights = c(0.9, 1, 1)) +
  plot_annotation(
    title = "3sp EMMAX: null-region count and region-level FDR over the tau_C x l_min grid, for four structure nulls",
    subtitle = "x = operating point (tau=0.05, l_min=3). Top panel (observed) is null-independent; the null-count and FDR rows are null-specific.")
ggsave(file.path(OUTFIG, "null_fdr_landscape_4nulls.png"), design, width = 15, height = 10, dpi = 160)
cat("[2] wrote figures/null_fdr_landscape_4nulls.png\n")
