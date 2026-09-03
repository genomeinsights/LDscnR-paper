## =============================================================================
## module_sim_LDscnR / neutral_chr_control.R
##
## PK: "one control left that is for free -- the neutral chromosomes." Every
## bundle holds one QTN chromosome and one `ntrl` chromosome with ZERO QTN, and
## the neutral half carries 53.4% of markers. So a discovery there is an
## unambiguous false positive: no dmax, no r2 threshold, no tagging convention.
##
## IT VALIDATES THE CONVENTION AND THEN GOES PAST IT.
##   neutral-chromosome FDP   0.68 - 0.84   (no convention)
##   tagging-based FDP        0.57 - 0.76   (dmax + r2 convention)
## They agree, so the tagging convention was not doing hidden work. Note the
## neutral estimate is a LOWER bound: a unit on a QTN chromosome can be
## significant through long-range LD with the QTN, which inflates the
## false-positive rate there relative to a neutral chromosome.
##
## THE HEADLINE, AND IT IS NOT COMFORTABLE. Cluster-level BH at a nominal 5% is
## running at roughly 70% FDP. Single-SNP is worse (~81%, from precision 0.187),
## so the 2x precision advantage is real -- but NEITHER analysis is anywhere near
## its nominal level on these data, and that has to be said plainly rather than
## reported as a ratio.
##
## WHERE THE FAILURE IS -- NOT WHERE IT LOOKED. Three diagnostics, in order:
##
##   MARKER-level p on neutral chromosomes: gif 0.91-0.99, 0.78-1.22% below 0.01,
##   0.07-0.25% below 0.001. The association test is FINE; genomic inflation is
##   not the cause.
##
##   CLUSTER-level p on neutral chromosomes: the BODY is calibrated (consensus
##   gif 0.998) but the TAIL is 1.8-3x too heavy -- 0.185% below 0.001 against an
##   expected 0.100% for consensus, 0.300% for best, 0.178% for Simes. BH lives
##   entirely in the tail, so a 2x tail excess is enough to wreck it.
##
##   best-SNP has gif 1.334, confirming that the minimum of correlated member
##   p-values is anti-conservative -- which is why BH could never use it and why
##   a permutation null rescues it: the null carries the same excess.
##
## AND THE TAIL EXCESS EXPLAINS THE "ONE CELL WINS" PHENOMENON, which was
## previously unexplained and which the sensitivity grid showed was NOT an ld_w
## artefact:
##
##   cell        tail (% p < 0.001)   precision
##   V0.5_c1          0.096            0.58      <- calibrated, expected 0.100
##   V2_c1            0.119            0.52
##   V1_c1.5          0.187            0.21
##   V0.5_c2          0.338            0.13      <- 3.4x too heavy
##
## Precision tracks null calibration, not detection power. The cells where
## clustering appears to fail are the cells where the cluster-level null is
## anti-conservative in the tail. THAT MAKES THE STRUCTURE-AWARE NULL A FIX
## RATHER THAN A FEATURE, because it calibrates the tail empirically instead of
## assuming uniformity.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/neutral_chr_control.R
## Env: DIR (structure_null results), PATTERN (which basis), OUT
## =============================================================================
suppressMessages(library(data.table))
DIR <- Sys.getenv("DIR", "module_sim_LDscnR/results/structure_null")
PAT <- Sys.getenv("PATTERN", "spatial")
CT  <- as.data.table(readRDS(file.path(DIR, "cl_chrtype.rds")))
PCOL <- c(cons = "p_cons", best = "p_best", simes = "p_simes")

one <- function(f) {
  z <- readRDS(f); U <- z$units; key <- z$summary[1]
  ct <- CT[cell == key$cell & tag == key$tag & env == key$env, .(CL, chr_type)]
  U  <- merge(U, ct, by = "CL", all.x = TRUE)
  fp <- list(); gi <- list()
  for (rt in names(PCOL)) {
    pv <- U[[PCOL[rt]]]
    ## calibration on the neutral half, where the null is exact by construction
    pn <- pv[U$chr_type == "ntrl" & is.finite(pv)]
    if (length(pn)) gi[[length(gi)+1]] <- data.table(route = rt, n = length(pn),
      gif = stats::median(stats::qchisq(pn, 1, lower.tail = FALSE)) / stats::qchisq(0.5, 1),
      pct_lt_01 = 100 * mean(pn < 0.01), pct_lt_001 = 100 * mean(pn < 0.001))
    for (fl in c(1, 2, 5, 8)) {
      idx <- which(U$n_loci >= fl & is.finite(pv)); if (!length(idx)) next
      sig <- idx[p.adjust(pv[idx], "BH") < 0.05]
      fp[[length(fp)+1]] <- data.table(route = rt, floor = fl, R = length(sig),
        R_ntrl = sum(U$chr_type[sig] == "ntrl", na.rm = TRUE),
        f_ntrl = mean(U$chr_type[idx] == "ntrl", na.rm = TRUE))
    }
  }
  list(fp = rbindlist(fp)[, `:=`(cell = key$cell, tag = key$tag, env = key$env)],
       gi = rbindlist(gi)[, `:=`(cell = key$cell, tag = key$tag, env = key$env)])
}
Z  <- lapply(list.files(DIR, pattern = PAT, full.names = TRUE), one)
FP <- rbindlist(lapply(Z, `[[`, "fp")); GI <- rbindlist(lapply(Z, `[[`, "gi"))
fwrite(FP, file.path(DIR, "neutral_chr_control.csv"))
fwrite(GI, file.path(DIR, "neutral_chr_calibration.csv"))

A <- FP[, .(R = mean(R), R_ntrl = mean(R_ntrl), f = mean(f_ntrl)), by = .(route, floor)]
A[, FDP_neutral := (R_ntrl / f) / pmax(R, 1e-9)]
cat("== neutral-chromosome FDP (no convention; a lower bound)\n")
print(A[order(route, floor), .(route, floor, R = round(R, 1), R_ntrl = round(R_ntrl, 2),
                               FDP_neutral = round(FDP_neutral, 3))])
cat("\n== cluster-level p on neutral chromosomes (expected: gif 1.00, 1.00%, 0.100%)\n")
print(GI[, .(gif = round(stats::median(gif), 3), pct_lt_0.01 = round(mean(pct_lt_01), 2),
             pct_lt_0.001 = round(mean(pct_lt_001), 3)), by = route])
cat("\n== tail excess by simulation cell, consensus route\n")
print(GI[route == "cons", .(pct_lt_0.001 = round(mean(pct_lt_001), 3)), by = cell][order(pct_lt_0.001)])
