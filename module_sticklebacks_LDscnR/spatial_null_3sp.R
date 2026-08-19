## =====================================================================
## module_sticklebacks_LDscnR / spatial_null_3sp.R
##
## Build the SPATIAL-autocorrelation structure null for 3sp -- the empirical
## counterpart of the sim spatial null (module_sim_LDscnR). The surrogate
## phenotype is a Gaussian process over the population GPS coordinates,
## residualised on the observed ecotype, so it carries the SAME spatial
## autocorrelation as the environment but is orthogonal to it. Scanned by
## emmax_fast on the SAME GRM as the observed ecotype (kinship still corrected).
##
## This complements the other two 3sp nulls:
##   * basis="genetic"  (null_uncapped_3sp.rds) -- MVN from the GRM; self-absorbed
##     by the kinship correction, so it produces almost no null peaks (a floor).
##   * pop-permutation  (null_popperm_3sp.rds)  -- shuffles ecotype among pops;
##     preserves population structure, ignores geography.
##   * basis="spatial"  (this)                  -- preserves spatial autocorrelation
##     of the environment; asks whether observed peaks exceed spatially-structured
##     but ecotype-orthogonal variation.
##
## 3sp sampling is continental (lat 47-75, lon -21..37), so spatial structure is
## real. Coordinates are per individual (GPS_N/E), aligned to the GTs rows.
##
## Run from the LDscnR-paper root (heavy, ~35-40 min at B=200):
##   Rscript module_sticklebacks_LDscnR/spatial_null_3sp.R [B]
## =====================================================================

suppressMessages({ library(data.table); library(LDscnR) })

a    <- commandArgs(trailingOnly = TRUE)
B    <- if (length(a) >= 1) as.integer(a[1]) else 200L
BUNDLE <- "module_sticklebacks_LDscnR/data/3sp_LDscnR_data.rds"
RAW    <- "~/gitlab/LD-scaling-genome-scans/empirical_data/3sp/3sp_data.RData"
OUTF   <- "module_sticklebacks_LDscnR/results/null_spatial_3sp.rds"

d  <- readRDS(BUNDLE); map <- as.data.table(d$map)
e  <- new.env(); load(path.expand(RAW), envir = e); ph <- as.data.table(e$pheno_3sp)
## align raw phenotype rows to the GTs rows used in the bundle
stopifnot(nrow(ph) == nrow(d$GTs), all(d$eco == as.integer(ph$ecotype == "Marine")))
coords <- as.matrix(ph[, .(GPS_N_updated, GPS_E_updated)])
cat(sprintf("[1] %d individuals ; %d pops ; coords lat %.1f..%.1f lon %.1f..%.1f\n",
            nrow(coords), length(unique(ph$pop_ID)),
            min(coords[,1]), max(coords[,1]), min(coords[,2]), max(coords[,2])))

QSTAR <- seq(0, 0.95, by = 0.05); ALPHA <- 0.05
prep  <- emmax_setup(d$GTs, d$GRM)
cat(sprintf("[2] building spatial null (GP over GPS coords, resid on ecotype), B=%d ...\n", B)); flush.console()
null <- structured_null(d$eco, d$GTs, d$GRM, d$ld_ws, basis = "spatial", coords = coords,
                        B = B, alpha = ALPHA, qstar = QSTAR, prep = prep, seed = 1L)
saveRDS(null, OUTF)
cat(sprintf("[3] C_obs: %d C>0 ; universe=%d ; saved %s\n",
            sum(null$C_obs > 0), length(null$universe), OUTF))

## quick calibrated look (same knobs as the pop-perm report)
edges <- ld_edges(null$universe, d$GTs, map[, .(marker, Chr, Pos)],
                  as.data.table(d$LD_decay$decay_sum), rho_ld = 0.60, dcap = 5e5)
op_l  <- calibrate_lmin(null, edges, tau = 0.05, q = 0.99)
tau_e <- calibrate_tauc(null, edges, l_min = op_l, fdr = 0.05, tau_grid = seq(0.05, 1, 0.05))
cat(sprintf("\n=== SPATIAL null: calibrated tau_C=%.3f, l_min=%d ===\n", tau_e, op_l))
cat("null_fdr grid (region-level FDR at l_min=3):\n")
print(as.data.table(null_fdr(null, edges, c(0.05, 0.1, 0.2), 3L)))
