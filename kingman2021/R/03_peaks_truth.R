## kingman2021/R/03_peaks_truth.R
## Turn the published EcoPeak / TempoPeak BEDs into the near-truth region table used to
## score the C-score outlier-region pipeline, and tag every SNP of a cohort with which
## peak (if any) it falls in.
##
## The pairing that matters: EcoPeaks were called on exactly the cohort we analyse, so
##   cohort c155_global -> gasAcu1-4.c155.*  peaks
##   cohort c150_pacNW  -> gasAcu1-4.c150.*  peaks
## Using the other cohort's peaks as truth would be a mismatch.
##
## Run from LDscnR-paper/ AFTER 02_build_rds.R:
##   Rscript kingman2021/R/03_peaks_truth.R c155_global
## Output (git-ignored): kingman2021/data/kingman2021_<cohort>_truth.rds

suppressMessages({ library(data.table) })

args   <- commandArgs(trailingOnly = TRUE)
COHORT <- if (length(args) >= 1) args[1] else "c155_global"
TAG    <- if (grepl("c150", COHORT)) "c150" else "c155"

DATA <- "/Users/petrikem/gitlab/LD-scaling-genome-scans/empirical_data/kingman2021"
MOD  <- "/Users/petrikem/gitlab/LDscnR-paper/kingman2021/data"
TR   <- file.path(DATA, "tracks")

read_peaks <- function(f, set) {
  b <- fread(file.path(TR, f), header = FALSE, sep = "\t")
  setnames(b, 1:3, c("Chr", "start", "end"))
  b[, `:=`(set = set,
           p_snp = if (ncol(b) >= 4) as.numeric(sub(".*=", "", V4)) else NA_real_,
           p_win = if (ncol(b) >= 5) as.numeric(sub(".*=", "", V5)) else NA_real_)]
  b[, .(Chr, start, end, set, p_snp, p_win)]
}

peaks <- rbindlist(list(
  read_peaks(sprintf("gasAcu1-4.%s.specific.50kb.final.peaks.bed",  TAG), "eco_specific"),
  read_peaks(sprintf("gasAcu1-4.%s.sensitive.50kb.final.peaks.bed", TAG), "eco_sensitive"),
  read_peaks("CH_SC_LB.specific.final.bed",  "tempo_specific"),
  read_peaks("CH_SC_LB.sensitive.final.bed", "tempo_sensitive")
))
peaks[, width := end - start]

cat(sprintf("cohort %s -> %s peak sets\n", COHORT, TAG))
print(peaks[, .(n = .N, Mb = round(sum(width) / 1e6, 2),
                median_kb = round(median(width) / 1e3, 1)), by = set])

## ---- tag the cohort's SNPs --------------------------------------------------------
rds <- file.path(MOD, sprintf("kingman2021_%s.rds", COHORT))
if (file.exists(rds)) {
  map <- as.data.table(readRDS(rds)$map)
  ## setkey below re-sorts by Chr as a *string*; keep the original row order so that
  ## snp_tags stays row-aligned with the map inside kingman2021_<cohort>.rds
  map[, row_id := .I]
  map[, `:=`(start = Pos - 1L, end = Pos)]
  setkey(map, Chr, start, end)
  for (s in unique(peaks$set)) {
    pk <- peaks[set == s]; setkey(pk, Chr, start, end)
    ov <- foverlaps(map, pk, type = "any", mult = "first", nomatch = NA)
    set(map, j = s, value = !is.na(ov$set))
  }
  cat("\nSNPs inside each peak set (of ", nrow(map), " total):\n", sep = "")
  for (s in unique(peaks$set))
    cat(sprintf("  %-16s %8d SNPs (%.2f%%)\n", s, sum(map[[s]]), 100 * mean(map[[s]])))
  setorder(map, row_id)                       # back to the map's own row order
  map[, c("start", "end", "row_id") := NULL]
} else {
  cat(sprintf("\nNOTE: %s not built yet -- saving peak table only.\n", basename(rds)))
  map <- NULL
}

out <- file.path(MOD, sprintf("kingman2021_%s_truth.rds", COHORT))
saveRDS(list(peaks = peaks, snp_tags = map, cohort = COHORT, peak_tag = TAG), out)
cat(sprintf("\nwrote %s\n", out))
