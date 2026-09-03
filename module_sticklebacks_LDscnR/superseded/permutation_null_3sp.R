## =====================================================================
## module_sticklebacks_LDscnR / permutation_null_3sp.R
##
## Ecotype-PERMUTATION structure null for 3sp, stratified or global. Shuffles the
## marine/freshwater label among populations (preserving population structure and
## breaking the ecotype<->genotype link), then scans each surrogate with
## emmax_fast on the SAME GRM as the observed ecotype and reduces to a C-score.
##
##   strata = "none"          global among-population permutation (all 35 pops
##                            in one pool) -- reproduces null_popperm_3sp.rds.
##   strata = "pop_locality"  REGIONAL permutation: shuffle ecotype only AMONG
##                            pops WITHIN each locality (Baltic, NorthSea,
##                            Norwegian, WB_Sea -- each has both ecotypes), so
##                            each region's ecotype composition is fixed and only
##                            the within-region marine/freshwater contrast is
##                            tested. Strictly more conservative than global:
##                            a parallel locus (same freshwater allele across
##                            regions) is NOT rebuilt by within-region shuffling,
##                            whereas a locus riding regional/clinal structure is.
##
## Run from the LDscnR-paper root (heavy, ~35-40 min at B=200):
##   Rscript module_sticklebacks_LDscnR/permutation_null_3sp.R [strata] [B]
##   defaults: strata=pop_locality  B=200
## =====================================================================

suppressMessages({ library(data.table); library(LDscnR) })

a      <- commandArgs(trailingOnly = TRUE)
STRATA <- if (length(a) >= 1) a[1] else "pop_locality"
B      <- if (length(a) >= 2) as.integer(a[2]) else 200L
BUNDLE <- "module_sticklebacks_LDscnR/data/3sp_LDscnR_data.rds"
RAW    <- "~/gitlab/LD-scaling-genome-scans/empirical_data/3sp/3sp_data.RData"
TAG    <- if (STRATA == "none") "popperm" else "regionperm"
BASIS  <- if (STRATA == "none") "pop_perm" else "region_perm"
OUTF   <- sprintf("module_sticklebacks_LDscnR/results/null_%s_3sp.rds", TAG)

d  <- readRDS(BUNDLE); map <- as.data.table(d$map)
e  <- new.env(); load(path.expand(RAW), envir = e); ph <- as.data.table(e$pheno_3sp)
stopifnot(nrow(ph) == nrow(d$GTs), all(d$eco == as.integer(ph$ecotype == "Marine")))

## population-level table: one ecotype + one stratum per pop
ph[, stratum := if (STRATA == "none") "ALL" else get(STRATA)]
pop_tab <- unique(ph[, .(pop_ID, ecotype, stratum)])
if (anyDuplicated(pop_tab$pop_ID)) stop("a pop maps to >1 ecotype/stratum")
cat(sprintf("[1] strata=%s ; %d pops in %d stratum/strata:\n", STRATA, nrow(pop_tab), uniqueN(pop_tab$stratum)))
print(dcast(pop_tab, stratum ~ ecotype, fun.aggregate = length, value.var = "pop_ID"))

## permute ecotype among pops WITHIN each stratum, map back to individuals
permute_within <- function() {
  pt <- copy(pop_tab)
  pt[, eperm := sample(ecotype), by = stratum]      # shuffle labels inside each stratum
  as.integer(pt$eperm[match(ph$pop_ID, pt$pop_ID)] == "Marine")
}

QSTAR <- seq(0, 0.95, by = 0.05); ALPHA <- 0.05
prep  <- emmax_setup(d$GTs, d$GRM)
Cof   <- function(y) { C <- ld_cscore(emmax_fast(prep, y), d$ld_ws, alpha = ALPHA, qstar = QSTAR); names(C) <- map$marker; C }
C_obs <- Cof(d$eco)
cat(sprintf("[2] C_obs: %d C>0 ; building %s null B=%d ...\n", sum(C_obs > 0), BASIS, B)); flush.console()
set.seed(1); C_surr <- vector("list", B)
for (b in 1:B) { Cs <- Cof(permute_within()); C_surr[[b]] <- Cs[Cs > 0]
  if (b %% 40 == 0) { cat("  perm", b, "/", B, "\n"); flush.console() } }

null <- structure(list(C_obs = C_obs, C_surr = C_surr,
  universe = unique(c(names(C_obs)[C_obs > 0], unlist(lapply(C_surr, names)))),
  basis = BASIS, B = B, strata = STRATA), class = "ld_null")
saveRDS(null, OUTF)
cat(sprintf("[3] saved %s ; universe=%d\n", OUTF, length(null$universe)))

## quick calibrated look
edges <- ld_edges(null$universe, d$GTs, map[, .(marker, Chr, Pos)],
                  as.data.table(d$LD_decay$decay_sum), rho_ld = 0.60, dcap = 1e5)
op_l  <- calibrate_lmin(null, edges, tau = 0.05, q = 0.99)
tau_e <- calibrate_tauc(null, edges, l_min = op_l, fdr = 0.05, tau_grid = seq(0.05, 1, 0.05))
cat(sprintf("\n=== %s null: calibrated tau_C=%.3f, l_min=%d ===\n", BASIS, tau_e, op_l))
