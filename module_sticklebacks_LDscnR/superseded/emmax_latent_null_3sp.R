## =====================================================================
## module_sticklebacks_LDscnR / emmax_latent_null_3sp.R
##
## The EMMAX side of the symmetric 2x2 null design. EMMAX already has its home-field
## K-MVN null (null_uncapped_3sp.rds); this adds the CROSS null -- phenotypes drawn
## from LFMM's home-field structure model (the top-K genotype-PC latent subspace) and
## scanned by EMMAX on the same GRM. It is the mirror of running K-MVN through LFMM:
## together the four cells ({EMMAX,LFMM} x {K-MVN,latent}) show how each engine
## behaves against its OWN structure model versus the OTHER's, so the comparison is
## not stacked in favour of either.
##
## Construction is identical to structured_null(basis="genetic") except the surrogate
## covariance is U diag(d^2) U^T (top-K genotype PCs) instead of the full-rank K. The
## SAME basis file feeds the LFMM runner, so both engines see the identical latent
## space.  y_surr = U %*% (d * rnorm(K)), residualised on the observed ecotype.
##
## Run from the LDscnR-paper root (fast -- EMMAX, not LFMM):
##   Rscript module_sticklebacks_LDscnR/emmax_latent_null_3sp.R [B]     # default B=200
## Writes results/null_latent_3sp.rds (an ld_null bundle for region_empirical_pvals.R).
## =====================================================================
suppressMessages({ library(data.table); library(LDscnR) })

## args: [B] [basis_file]   (basis_file default = raw all-SNP latent basis)
a    <- commandArgs(trailingOnly = TRUE)
B    <- if (length(a) >= 1) as.integer(a[1]) else 200L
BND  <- "module_sticklebacks_LDscnR/data/3sp_LDscnR_data.rds"
LAT  <- if (length(a) >= 2) a[2] else "module_sticklebacks_LDscnR/data/3sp_latent_basis.rds"
btag <- sub("\\.rds$", "", sub(".*3sp_latent_basis", "", basename(LAT)))   # "", "_thin250", ...
OUTF <- sprintf("module_sticklebacks_LDscnR/results/null_latent%s_3sp.rds", btag)
QSTAR <- seq(0, 0.95, by = 0.05); ALPHA <- 0.05; SEED <- 1L

d  <- readRDS(BND); map <- as.data.table(d$map)
pc <- readRDS(LAT); U <- pc$U; dd <- pc$d; K <- pc$K
y  <- d$eco; n <- length(y)
stopifnot(nrow(U) == n)

prep <- emmax_setup(d$GTs, d$GRM)                        # same kinship correction as every EMMAX null
gen  <- function() as.numeric(stats::resid(stats::lm(as.numeric(U %*% (dd * stats::rnorm(K))) ~ y)))
sparseC <- function(pv) { C <- ld_cscore(pv, d$ld_ws, ALPHA, colnames(d$ld_ws), QSTAR); C[C > 0] }

cat(sprintf("[1] EMMAX latent null: n=%d, top-%d PC basis (var expl %s), B=%d\n",
            n, K, paste(sprintf("%.1f%%", 100 * pc$var_explained), collapse=" "), B)); flush.console()
C_obs <- ld_cscore(emmax_fast(prep, y), d$ld_ws, ALPHA, colnames(d$ld_ws), QSTAR)
set.seed(SEED)
C_surr <- vector("list", B)
for (b in seq_len(B)) C_surr[[b]] <- sparseC(emmax_fast(prep, gen()))
universe <- unique(c(names(C_obs)[C_obs > 0], unlist(lapply(C_surr, names))))
null <- structure(list(C_obs = C_obs, C_surr = C_surr, universe = universe,
                       basis = "latent", B = B,
                       params = list(alpha = ALPHA, rho = colnames(d$ld_ws), qstar = QSTAR)),
                  class = "ld_null")
saveRDS(null, OUTF)
cat(sprintf("[2] C_obs: %d C>0 ; universe=%d ; median null C>0/surrogate = %.0f (vs observed %d) ; saved %s\n",
            sum(C_obs > 0), length(universe),
            stats::median(vapply(C_surr, length, integer(1))), sum(C_obs > 0), OUTF))
cat("   score it: Rscript module_sticklebacks_LDscnR/region_empirical_pvals.R 0.05 3 0.60", OUTF, "\n")
