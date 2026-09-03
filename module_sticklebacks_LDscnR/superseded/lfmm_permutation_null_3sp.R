## =====================================================================
## module_sticklebacks_LDscnR / lfmm_permutation_null_3sp.R
##
## LFMM permutation structure-null for 3sp -- the direct test of whether the LFMM
## scan is inflated on this real genome. Builds the SAME among-population ecotype
## permutation null used for EMMAX (permutation_null_3sp.R), but with a full LFMM
## (K=5, genomic control) re-run per surrogate instead of the fast emmax engine.
## If LFMM's observed regions are real they will clear this null (as EMMAX's did);
## if LFMM is inflated by residual structure, the permuted scans will reproduce its
## peaks and the observed regions will NOT clear it.
##
## EXPENSIVE: each surrogate is a whole-genome LFMM fit+test on 117 x ~790k SNPs.
## Parallelised over B via mclapply. Observed LFMM is RE-RUN here (not the bundle's
## reused values) so observed and null share one engine.
##
## Run from the LDscnR-paper root:
##   LFMM_CORES=8 Rscript module_sticklebacks_LDscnR/lfmm_permutation_null_3sp.R [B]
##   defaults: B=100, cores from LFMM_CORES (else 4)
## Writes results/null_lfmm_perm_3sp.rds (an ld_null bundle for region_empirical_pvals.R).
## =====================================================================
suppressMessages({ library(data.table); library(LEA); library(LDscnR); library(parallel) })

a     <- commandArgs(trailingOnly = TRUE); B <- if (length(a) >= 1) as.integer(a[1]) else 100L
CORES <- as.integer(Sys.getenv("LFMM_CORES", "4"))
if (CORES > 1L) { data.table::setDTthreads(1L); Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1") }
BND   <- "module_sticklebacks_LDscnR/data/3sp_LDscnR_data.rds"
RAW   <- "~/gitlab/LD-scaling-genome-scans/empirical_data/3sp/3sp_data.RData"
OUTF  <- "module_sticklebacks_LDscnR/results/null_lfmm_perm_3sp.rds"
K_LFMM <- 5L; QSTAR <- seq(0, 0.95, by = 0.05); ALPHA <- 0.05

d  <- readRDS(BND); map <- as.data.table(d$map)
e  <- new.env(); load(path.expand(RAW), envir = e); ph <- as.data.table(e$pheno_3sp)
stopifnot(nrow(ph) == nrow(d$GTs), all(d$eco == as.integer(ph$ecotype == "Marine")))
pop <- ph$pop_ID; pop_ecotype <- unique(ph[, .(pop = pop_ID, ecotype)])
if (anyDuplicated(pop_ecotype$pop)) stop("one ecotype per pop expected")
permute_among_pops <- function() { pt <- copy(pop_ecotype); pt[, ep := sample(ecotype)]
  as.integer(pt$ep[match(pop, pt$pop)] == "Marine") }

## one whole-genome LFMM scan (K=5, GC) for a phenotype vector -> per-SNP p, ordered as map
lfmm_scan <- function(y) {
  tmp <- tempfile("lfmm3sp_"); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  gf <- file.path(tmp, "g.lfmm"); ef <- file.path(tmp, "e.env")
  write.lfmm(d$GTs, gf); write.env(y, ef)
  proj <- lfmm2(gf, ef, K = K_LFMM)
  pv <- suppressWarnings(lfmm2.test(proj, gf, ef, genomic.control = TRUE, full = TRUE))
  pv$pvalues
}
Cof <- function(pv) { C <- ld_cscore(pv, d$ld_ws, alpha = ALPHA, qstar = QSTAR); names(C) <- map$marker; C }

cat(sprintf("[1] 3sp LFMM permutation null: %d indiv x %d SNPs, %d pops, B=%d, cores=%d\n",
            nrow(d$GTs), ncol(d$GTs), length(unique(pop)), B, CORES)); utils::flush.console()
t0 <- Sys.time()
C_obs <- Cof(lfmm_scan(d$eco))
cat(sprintf("[2] observed LFMM: %d C>0 ; one scan = %.1f min\n",
            sum(C_obs > 0), as.numeric(Sys.time() - t0, units = "mins"))); utils::flush.console()

set.seed(1); perms <- lapply(seq_len(B), function(b) permute_among_pops())
C_surr <- mclapply(seq_len(B), function(b) { Cs <- Cof(lfmm_scan(perms[[b]])); Cs[Cs > 0] },
                   mc.cores = CORES, mc.preschedule = FALSE)
bad <- !vapply(C_surr, is.numeric, logical(1)); if (any(bad)) { cat(sprintf(" dropped %d failed surrogate(s)\n", sum(bad))); C_surr <- C_surr[!bad] }
null <- structure(list(C_obs = C_obs, C_surr = C_surr,
  universe = unique(c(names(C_obs)[C_obs > 0], unlist(lapply(C_surr, names)))),
  basis = "lfmm_perm", B = length(C_surr)), class = "ld_null")
saveRDS(null, OUTF)
cat(sprintf("[3] done in %.1f min ; %d surrogates ; universe=%d ; saved %s\n",
            as.numeric(Sys.time() - t0, units = "mins"), length(C_surr), length(null$universe), OUTF))
cat("   score it with: Rscript module_sticklebacks_LDscnR/region_empirical_pvals.R 0.05 3 0.60", OUTF, "\n")
