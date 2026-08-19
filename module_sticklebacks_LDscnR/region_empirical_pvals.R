## =====================================================================
## module_sticklebacks_LDscnR / region_empirical_pvals.R
##
## Region-level empirical p-values for 3sp from the AMONG-POPULATION ecotype
## permutation null (null_popperm_3sp.rds).
##
## Rationale. A global count FDR from this null is too conservative: ecotype is
## nested in population, so a locus that drives marine<->freshwater adaptation in
## only PART of the geographic range is charged to "structure". But such a peak
## needs THOSE particular populations to be freshwater -- so permuting ecotype
## among populations rarely reproduces it AT THAT LOCATION. That gives a
## location-resolved test: for each observed region, how often does a permutation
## produce an overlapping peak at least as strong?
##
##   emp_p(region) = (1 + #perms with an overlapping null region of score >= obs) / (1 + B)
##   overlap_freq  = fraction of perms producing ANY region overlapping the location
##
## Region score = summed C over the region's markers (mass = size x strength).
## Peaks that the structure-preserving permutation seldom rebuilds at their locus
## get a low emp_p and are the robust, possibly limited-range, adaptation signals.
##
## Run from the LDscnR-paper root (fast; the heavy null is already computed):
##   Rscript module_sticklebacks_LDscnR/region_empirical_pvals.R [tau_C] [l_min] [rho_ld] [null.rds]
##   defaults: tau_C=0.05  l_min=3  rho_ld=0.60  null=null_popperm_3sp.rds
##   MVN/genetic null: ... 0.05 3 0.60 module_sticklebacks_LDscnR/results/null_uncapped_3sp.rds
## Outputs are tagged by the null's basis, so the two nulls do not overwrite.
## =====================================================================

suppressMessages({ library(data.table); library(LDscnR) })

a      <- commandArgs(trailingOnly = TRUE)
TAU    <- if (length(a) >= 1) as.numeric(a[1]) else 0.05
LMIN   <- if (length(a) >= 2) as.integer(a[2]) else 3L
RHO_LD <- if (length(a) >= 3) as.numeric(a[3]) else 0.60
NULLF  <- if (length(a) >= 4) a[4] else "module_sticklebacks_LDscnR/results/null_popperm_3sp.rds"
DCAP   <- 5e5
BUNDLE <- "module_sticklebacks_LDscnR/data/3sp_LDscnR_data.rds"
OUTRES <- "module_sticklebacks_LDscnR/results"; OUTFIG <- "module_sticklebacks_LDscnR/figures"
for (p in c(OUTRES, OUTFIG)) if (!dir.exists(p)) dir.create(p, recursive = TRUE)

## ---- 1. data + null --------------------------------------------------
d    <- readRDS(BUNDLE); map <- as.data.table(d$map)
null <- readRDS(NULLF)
BASIS <- if (!is.null(null$basis)) null$basis else "null"
B    <- length(null$C_surr)
cat(sprintf("[1] %s null B=%d ; C_obs: %d C>0 ; tau=%.3f l_min=%d rho_ld=%.2f ; universe=%d\n",
            BASIS, B, sum(null$C_obs > 0), TAU, LMIN, RHO_LD, length(null$universe)))
edges <- ld_edges(null$universe, d$GTs, map[, .(marker, Chr, Pos)],
                  as.data.table(d$LD_decay$decay_sum), rho_ld = RHO_LD, dcap = DCAP)
mpos <- stats::setNames(map$Pos, map$marker)
mchr <- stats::setNames(as.character(map$Chr), map$marker)

## ---- 2. cluster a C-vector into scored genomic regions ---------------
regstats <- function(Cvec) {
  mk <- names(Cvec)[Cvec >= TAU]
  if (!length(mk)) return(data.table(Chr = character(), lo = numeric(), hi = numeric(),
                                     score = numeric(), size = integer(), maxC = numeric()))
  regs <- ld_regions(mk, edges); regs <- regs[lengths(regs) >= LMIN]
  if (!length(regs)) return(data.table(Chr = character(), lo = numeric(), hi = numeric(),
                                        score = numeric(), size = integer(), maxC = numeric()))
  rbindlist(lapply(regs, function(r) data.table(
    Chr = unname(mchr[r[1]]), lo = min(mpos[r]), hi = max(mpos[r]),
    score = sum(Cvec[r]), size = length(r), maxC = max(Cvec[r]))))
}

obs <- regstats(null$C_obs)
if (!nrow(obs)) stop("no observed regions at this (tau, l_min)")
obs[, id := .I]
cat(sprintf("[2] observed regions: %d  (chromosomes: %s)\n",
            nrow(obs), paste(sort(unique(obs$Chr)), collapse = ", ")))

## null regions per permutation
null_regs <- lapply(null$C_surr, regstats)

## ---- 3. location-matched empirical p per observed region -------------
## for each perm: best score among null regions overlapping the observed locus
## NB: args are o_* (not lo/hi) so they do not collide with nr's columns inside
## the data.table `[`, where bare lo/hi bind to the columns of nr.
overlap_best <- function(o_chr, o_lo, o_hi, nr) {
  if (!nrow(nr)) return(0)
  h <- nr[Chr == o_chr & lo <= o_hi & hi >= o_lo]   # same chromosome, position ranges intersect
  if (!nrow(h)) 0 else max(h$score)
}
emp <- rbindlist(lapply(seq_len(nrow(obs)), function(i) {
  best <- vapply(null_regs, function(nr) overlap_best(obs$Chr[i], obs$lo[i], obs$hi[i], nr), numeric(1))
  data.table(id = obs$id[i],
             overlap_freq = mean(best > 0),
             null_max_overlap = max(best),
             emp_p = (1 + sum(best >= obs$score[i])) / (1 + B))
}))
out <- merge(obs, emp, by = "id")[order(emp_p, -score)]
out[, emp_q := stats::p.adjust(emp_p, "fdr")]
out[, `:=`(lo_Mb = round(lo / 1e6, 3), hi_Mb = round(hi / 1e6, 3))]

## ---- 4. report + save ------------------------------------------------
cat("\n=== observed regions, ranked by location-matched empirical p ===\n")
print(out[, .(id, Chr, lo_Mb, hi_Mb, size, maxC = round(maxC, 3), score = round(score, 2),
              overlap_freq = round(overlap_freq, 3), emp_p = round(emp_p, 4), emp_q = round(emp_q, 4))])
cat(sprintf("\nregions with emp_p < 0.05 : %d / %d ; emp_q < 0.05 : %d\n",
            sum(out$emp_p < 0.05), nrow(out), sum(out$emp_q < 0.05)))
cat(sprintf("median overlap_freq (null rebuilds the locus): %.2f\n", stats::median(out$overlap_freq)))
fn <- sprintf("region_emp_pvals_%s_tau%.2f_lmin%d_rho%.2f.csv", BASIS, TAU, LMIN, RHO_LD)
fwrite(out, file.path(OUTRES, fn))
cat(sprintf("[4] wrote %s/%s\n", OUTRES, fn))
