## Merge three cross-based (pedigree) threespine stickleback recombination maps
## onto one set of bins in gasAcu1 coordinates.
##   Glazer et al. 2015 G3   : FTC and BEPA F2 crosses (bins already tile gasAcu1)
##   Roesti et al. 2013 MolEc: Appendix S4, 1872 markers, one F2 cross
## Run from the directory holding glazer2015_rate_bins.tsv and roesti_AppendixS4.txt.
suppressMessages(library(data.table))

## ---- Roesti: local rate per marker, then intervals in gasAcu1 coordinates ----
## Precondition on the Roesti input. Checked, not merely documented: Wiley is behind
## Cloudflare, which answers a scripted curl with HTTP 200 and a ~6 KB HTML interstitial
## rather than an error, so a failed fetch looks like a truncated file rather than an
## access block. See README.md, "Re-fetching the source files".
ROESTI_SRC  <- "source/roesti_AppendixS4.txt"
ROESTI_SIZE <- 120011
ROESTI_MD5  <- "255a8f02adc2bed9993c5284edfe3e0c"
if (!file.exists(ROESTI_SRC))
  stop(sprintf("MISSING: %s\n  source/ is gitignored, so a fresh clone will not have it.\n  README.md gives the DOI and says why a scripted fetch will not work.", ROESTI_SRC))
.sz <- file.size(ROESTI_SRC); .md5 <- unname(tools::md5sum(ROESTI_SRC))
if (!identical(.md5, ROESTI_MD5))
  stop(sprintf(paste0("WRONG CONTENT: %s\n  expected %d bytes, md5 %s\n  got      %d bytes, md5 %s\n  %s"),
       ROESTI_SRC, ROESTI_SIZE, ROESTI_MD5, .sz, .md5,
       if (.sz < 50000) "That size is in the range of a Cloudflare interstitial -- you have an HTML block page, not the file. Fetch it through a browser."
       else "Size plausible but content differs: revised upstream, or edited locally. Do not build on it silently."))
cat(sprintf("source file verified: %s (%d bytes)\n", ROESTI_SRC, .sz))

ro <- fread(ROESTI_SRC)
setnames(ro, c("chromosome_BroadS1", "position_BroadS1",
               "chromosome_reassembled", "position_reassembled"),
             c("old_chr", "old_pos", "new_chr", "new_pos"))

## slope of cM against physical position, over the two flanking markers, taken in
## the REASSEMBLED order (the order Roesti corrected; BroadS1 has misassemblies)
setorder(ro, new_chr, new_pos)
ro[, rate_ROESTI := {
  n <- .N
  lo <- pmax(1L, seq_len(n) - 1L); hi <- pmin(n, seq_len(n) + 1L)
  d_cm <- abs(cM[hi] - cM[lo]); d_bp <- abs(new_pos[hi] - new_pos[lo])
  fifelse(d_bp > 0, d_cm / (d_bp / 1e6), NA_real_)
}, by = new_chr]

len <- ro[, .(cM = max(cM) - min(cM), Mb = (max(new_pos) - min(new_pos)) / 1e6), by = new_chr]
cat(sprintf("Roesti: %d markers, %d linkage groups, %.1f cM total, %.2f cM/Mb genome-wide\n",
            nrow(ro), nrow(len), sum(len$cM), sum(len$cM) / sum(len$Mb)))
cat("  (published: 1872 markers, 3.11 cM/Mb)\n")

## carry each marker's rate back to its gasAcu1 position; markers the old assembly
## left unplaced (old_chr == "Un") cannot be looked up by gasAcu1 coordinate
ro_a <- ro[old_chr != "Un" & !is.na(rate_ROESTI)]
ro_a[, old_chr := as.integer(old_chr)]
setorder(ro_a, old_chr, old_pos)
## tile each gasAcu1 chromosome: interval boundaries at midpoints between markers
ro_iv <- ro_a[, {
  n <- .N
  mid <- (old_pos[-n] + old_pos[-1]) / 2
  .(start = as.integer(c(1, mid + 1)), end = as.integer(c(mid, old_pos[n])),
    rate_ROESTI = rate_ROESTI)
}, by = old_chr][end > start]
cat(sprintf("Roesti: %d markers anchored in gasAcu1 -> %d intervals\n", nrow(ro_a), nrow(ro_iv)))

## ---- Glazer bins as the common frame ----
gl <- fread("glazer2015_rate_bins.tsv")
gl <- gl[old_chr != "Un"][, old_chr := as.integer(old_chr)]
setorder(gl, old_chr, old_start)
setnames(gl, c("rate_FTC_cMperMb", "rate_BEPA_cMperMb"), c("rate_FTC", "rate_BEPA"))

## overlap-weighted mean of the Roesti intervals within each Glazer bin
setkey(ro_iv, old_chr, start, end)
ov <- foverlaps(gl[, .(old_chr, start = old_start, end = old_end, marker)],
                ro_iv, type = "any", nomatch = NULL)
ov[, w := pmin(end, i.end) - pmax(start, i.start) + 1]
ro_bin <- ov[, .(rate_ROESTI = sum(rate_ROESTI * w) / sum(w)), by = marker]
gl <- merge(gl, ro_bin, by = "marker", all.x = TRUE)
setorder(gl, old_chr, old_start)

## ---- consensus: normalise each cross to its own genome-wide mean, then average ----
crosses <- c("rate_FTC", "rate_BEPA", "rate_ROESTI")
gwm <- sapply(crosses, function(k) mean(gl[[k]], na.rm = TRUE))
cat("\nbin-mean rate per cross (cM/Mb): ",
    paste(sprintf("%s=%.2f", sub("rate_", "", crosses), gwm), collapse = "  "), "\n")
for (k in crosses) gl[[sub("rate_", "rel_", k)]] <- round(gl[[k]] / gwm[[k]], 4)
rel <- as.matrix(gl[, .(rel_FTC, rel_BEPA, rel_ROESTI)])
set(gl, j = "n_crosses",     value = rowSums(!is.na(rel)))
set(gl, j = "rel_consensus", value = round(rowMeans(rel, na.rm = TRUE), 4))
## rel_sd: plain spread of the relative rates. Kept as a description, but it is NOT a
## quality filter -- it is an SD of untransformed rates, so it scales with the mean
## (cor with rel_consensus about +0.64) and thresholding it selects LOW-RECOMBINATION
## bins rather than well-measured ones. Filtering on it truncates the top of the rate
## range and will attenuate any correlation against recombination by range restriction.
set(gl, j = "rel_sd", value = round(apply(rel, 1, sd, na.rm = TRUE), 4))

## rank_sd: disagreement measured on each cross's own genome-wide ORDERING of bins.
## That is the scale this map is meant to be read on (the relative landscape), and it
## is free of the magnitude confound because every cross contributes a uniform [0,1]
## percentile whatever its absolute cM/Mb calibration.
pctile <- apply(rel, 2, function(x) {
  p <- rep(NA_real_, length(x)); ok <- !is.na(x)
  p[ok] <- (rank(x[ok]) - 0.5) / sum(ok); p
})
set(gl, j = "rank_sd", value = round(apply(pctile, 1, sd, na.rm = TRUE), 4))

## rank_sd removes the magnitude scaling but is still not orthogonal to rate: percentiles
## are bounded, so bins at the very top or bottom of every cross's ordering have no room
## to disagree. Thresholding it keeps the extremes and drops the middle.
## `disagree` fixes that by ranking rank_sd WITHIN rate decile. It is uniform on [0,1] by
## construction and therefore carries no information about the rate itself -- 0 = the
## crosses agree about this bin relative to other bins of similar rate, 1 = they do not.
ok <- gl$n_crosses >= 2 & is.finite(gl$rank_sd) & is.finite(gl$rel_consensus)
set(gl, j = "disagree", value = NA_real_)
gl[ok, dec_ := cut(rel_consensus, quantile(rel_consensus, 0:10 / 10),
                   include.lowest = TRUE)]
gl[ok, disagree := round((rank(rank_sd) - 0.5) / .N, 4), by = dec_]
gl[, dec_ := NULL]
set(gl, j = "concordant", value = !is.na(gl$disagree) & gl$disagree <= 0.75)
## put the consensus on an absolute scale using Roesti's published genome-wide mean,
## the most conservative of the three (F2 map length inflates with marker number)
ROESTI_PUBLISHED <- 3.11
gl[, rate_consensus_cMperMb := round(rel_consensus * ROESTI_PUBLISHED, 4)]

cat("\npairwise correlation of local rate across bins scored by both crosses:\n")
for (p in list(c("rate_FTC","rate_BEPA"), c("rate_FTC","rate_ROESTI"), c("rate_BEPA","rate_ROESTI"))) {
  ok <- complete.cases(gl[[p[1]]], gl[[p[2]]])
  cat(sprintf("  %-12s vs %-12s r = %.3f  (n = %d)\n", sub("rate_","",p[1]), sub("rate_","",p[2]),
              cor(gl[[p[1]]][ok], gl[[p[2]]][ok]), sum(ok)))
}
cat("\nbins by number of crosses contributing:\n"); print(gl[, .N, by = n_crosses][order(-n_crosses)])
cat(sprintf("bins flagged concordant (>=2 crosses, disagree <= 0.75): %d of %d\n",
            sum(gl$concordant, na.rm = TRUE), nrow(gl)))
## the flag must not be a proxy for rate: check the confound both ways
cat(sprintf("confound check   cor(rel_sd , rel_consensus) = %+.3f   <- why rel_sd is not a filter\n",
            cor(gl$rel_sd, gl$rel_consensus, use = "complete.obs")))
cat(sprintf("                 cor(rank_sd, rel_consensus) = %+.3f   <- bounded, still not orthogonal\n",
            cor(gl$rank_sd, gl$rel_consensus, use = "complete.obs")))
cat(sprintf("                 cor(disagree, rel_consensus) = %+.3f   <- the criterion actually used\n",
            cor(gl$disagree, gl$rel_consensus, use = "complete.obs")))
dec <- gl[!is.na(rate_consensus_cMperMb)][, dec := cut(rate_consensus_cMperMb,
            quantile(rate_consensus_cMperMb, 0:10/10), include.lowest = TRUE, labels = FALSE)]
cat("       % of bins retained by `concordant`, by rate decile (1 = lowest rate):\n         ",
    paste(sprintf("%d:%.0f%%", 1:10, 100 * dec[, mean(concordant), by = dec][order(dec)]$V1),
          collapse = "  "), "\n")

out <- gl[, .(marker, old_chr, old_start, old_end, new_chr, new_start, new_end, bin_length,
              cM_FTC, cM_BEPA,
              rate_FTC, rate_BEPA, rate_ROESTI = round(rate_ROESTI, 4),
              rel_FTC, rel_BEPA, rel_ROESTI,
              n_crosses, rel_consensus, rel_sd, rank_sd, disagree, concordant,
              rate_consensus_cMperMb)]
fwrite(out, "stickleback_recomb_3crosses_gasAcu1.tsv", sep = "\t")
cat(sprintf("\nwrote stickleback_recomb_3crosses_gasAcu1.tsv: %d bins\n", nrow(out)))

## ---- coverage against the 3sp SNP set ----
sp <- "/Users/petrikem/gitlab/LD-scaling-genome-scans/empirical_data/3sp/SNP_res_3sp.rds"
if (file.exists(sp)) {
  sr <- as.data.table(readRDS(sp)); sr[, chr_i := as.integer(sub("Chr", "", Chr))]
  b <- copy(out); setkey(b, old_chr, old_start, old_end)
  h <- foverlaps(sr[, .(old_chr = chr_i, old_start = Pos, old_end = Pos)], b, type = "within")
  cat(sprintf("3sp SNPs %d; in a bin %.1f%%; with a consensus rate %.1f%%; all three crosses %.1f%%\n",
              nrow(sr), 100*mean(!is.na(h$marker)),
              100*mean(!is.na(h$rate_consensus_cMperMb)), 100*mean(h$n_crosses %in% 3)))
}
