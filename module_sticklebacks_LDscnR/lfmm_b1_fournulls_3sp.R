## =====================================================================
## module_sticklebacks_LDscnR / lfmm_b1_fournulls_3sp.R
##
## Quick DIRECTIONAL test of LFMM inflation on 3sp: run LFMM once on the observed
## ecotype and once on ONE surrogate from each of the four null constructions
## (genetic MVN, global permutation, regional permutation, spatial MVN kernel), and
## compare the surrogate "background" (number of C>0 markers and >=3-SNP regions) to
## the observed. For EMMAX the genetic surrogate was silent (~0 peaks); if LFMM's
## surrogates instead look like the observed (high background), the observed regions
## will not clear the null and B=1 already tells us the direction.
##
## 5 whole-genome LFMM scans (~5 GB, several min each), run in parallel over cores.
## Run from the LDscnR-paper root:
##   LFMM_CORES=5 Rscript module_sticklebacks_LDscnR/lfmm_b1_fournulls_3sp.R
## Writes results/lfmm_b1_fournulls_Cs.rds + prints the background table.
## =====================================================================
suppressMessages({ library(data.table); library(LEA); library(LDscnR); library(parallel) })
CORES <- as.integer(Sys.getenv("LFMM_CORES", "5"))
if (CORES > 1L) { data.table::setDTthreads(1L); Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1") }
BND <- "module_sticklebacks_LDscnR/data/3sp_LDscnR_data.rds"
RAW <- "~/gitlab/LD-scaling-genome-scans/empirical_data/3sp/3sp_data.RData"
K_LFMM <- 5L; QSTAR <- seq(0, 0.95, by = 0.05); ALPHA <- 0.05; TAU <- 0.05; LMIN <- 3L; RHO_LD <- 0.60; DCAP <- 1e5

d  <- readRDS(BND); map <- as.data.table(d$map); eco <- d$eco; n <- length(eco)
e  <- new.env(); load(path.expand(RAW), envir = e); ph <- as.data.table(e$pheno_3sp)
stopifnot(nrow(ph) == n, all(eco == as.integer(ph$ecotype == "Marine")))
pop <- ph$pop_ID
pe  <- unique(ph[, .(pop = pop_ID, ecotype)])
ple <- unique(ph[, .(pop = pop_ID, ecotype, loc = pop_locality)])
coords <- as.matrix(ph[, .(GPS_N_updated, GPS_E_updated)])

## ---- four null phenotype generators (same as the EMMAX four-null analysis) ----------
gen_mvn <- function(K) { eK <- eigen(K, symmetric = TRUE); Lv <- pmax(eK$values, 0); Vk <- eK$vectors
  as.numeric(stats::resid(stats::lm(as.numeric(Vk %*% (sqrt(Lv) * stats::rnorm(n))) ~ eco))) }
spatial_kernel <- function() { Dm <- as.matrix(stats::dist(coords)); l <- stats::median(Dm[lower.tri(Dm)]); exp(-0.5 * (Dm / l)^2) }
perm_global <- function() { pt <- copy(pe);  pt[, ep := sample(ecotype)];          as.integer(pt$ep[match(pop, pt$pop)] == "Marine") }
perm_region <- function() { pt <- copy(ple); pt[, ep := sample(ecotype), by = loc]; as.integer(pt$ep[match(pop, pt$pop)] == "Marine") }

set.seed(1)
phenos <- list(observed      = eco,
               genetic       = gen_mvn(d$GRM),
               global_perm   = perm_global(),
               regional_perm = perm_region(),
               spatial       = gen_mvn(spatial_kernel()))

lfmm_scan <- function(y) {
  tmp <- tempfile("lfmm3sp_"); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  gf <- file.path(tmp, "g.lfmm"); ef <- file.path(tmp, "e.env")
  write.lfmm(d$GTs, gf); write.env(y, ef)
  proj <- lfmm2(gf, ef, K = K_LFMM)
  suppressWarnings(lfmm2.test(proj, gf, ef, genomic.control = TRUE, full = TRUE))$pvalues
}
Cof <- function(pv) { C <- ld_cscore(pv, d$ld_ws, alpha = ALPHA, qstar = QSTAR); names(C) <- map$marker; C }

cat(sprintf("[1] 5 whole-genome LFMM scans (observed + 4 null types) on %d cores...\n", CORES)); flush.console()
t0 <- Sys.time()
Cs <- mclapply(phenos, function(y) Cof(lfmm_scan(y)), mc.cores = CORES, mc.preschedule = FALSE)
names(Cs) <- names(phenos)
saveRDS(Cs, "module_sticklebacks_LDscnR/results/lfmm_b1_fournulls_Cs.rds")
cat(sprintf("[2] scans done in %.1f min\n", as.numeric(Sys.time() - t0, units = "mins"))); flush.console()

## ---- background table: C>0 markers and >=3-SNP regions per type ---------------------
uni   <- unique(unlist(lapply(Cs, function(C) names(C)[C > 0])))
edges <- ld_edges(uni, d$GTs, map[, .(marker, Chr, Pos)], as.data.table(d$LD_decay$decay_sum), rho_ld = RHO_LD, dcap = DCAP)
nreg  <- function(C) { mk <- names(C)[C >= TAU]; if (!length(mk)) return(0L); sum(lengths(ld_regions(mk, edges)) >= LMIN) }
tab   <- data.table(type = names(Cs),
                    n_Cgt0 = vapply(Cs, function(C) sum(C > 0), integer(1)),
                    n_regions = vapply(Cs, nreg, integer(1)))
obs_r <- tab[type == "observed", n_regions]; obs_c <- tab[type == "observed", n_Cgt0]
tab[, frac_of_obs_regions := round(n_regions / max(obs_r, 1), 2)]
cat("\n=== LFMM background per null type (tau=0.05, l_min=3) ===\n"); print(tab)
cat(sprintf("\nINTERPRETATION: observed = %d regions / %d C>0. If a surrogate's counts approach these,\n", obs_r, obs_c))
cat("the null reproduces the observed 'signal' -> LFMM inflated, regions will not clear it.\n")
cat("If surrogates are near zero (as EMMAX's genetic null was), the observed LFMM regions are real.\n")
