## =============================================================================
## module_3sp/R/16_ldw_peaks.R
##
## Do outlier regions coincide with local peaks in ld_w? (Methods claim:
## "Where such a region is associated with a selection signal, the outlier
## region is expected to coincide with a local peak in ld_w.")
##
## Interval sets: the independently identified marine-freshwater regions
## (EcoPeaks; independent of LDscnR) and the reported consensus and Simes
## regions. Three tests:
##   1. Rotation: mean ld_w of the markers inside each interval, averaged over
##      intervals, against 10,000 within-chromosome random placements of every
##      interval with its span preserved.
##   2. Local peak: mean ld_w inside each interval versus its flanks
##      (FLANK bp either side, interval excluded); sign test.
##   3. Floor control (reported regions only): significant Stage-1 units are
##      by construction at least SIZE_FLOOR markers, and large units sit in
##      high-LD stretches. So significant units are compared with the
##      non-significant units that also passed the floor, adjusting for size.
## Output: out/16_ldw_peaks/ldw_peaks.tsv, ldw_units.tsv
## =============================================================================
suppressMessages({library(data.table)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_3sp"), "R", "00_config.R"))
STAGE <- "16_ldw_peaks"
say("=== %s ===\n\n", STAGE)
N_ROT <- 10000L; FLANK <- 250e3

b  <- readRDS(file.path(stage_dir("02_bundle"), "bundle.rds"))
sc <- readRDS(file.path(stage_dir("03_EMMAX"), "scan.rds"))
mk <- as.data.table(b$map)[!is.na(ld_w_095), .(Chr = as.character(Chr), Pos, ldw = ld_w_095)][order(Chr, Pos)]
chr_idx <- split(mk, by = "Chr", keep.by = FALSE)
for (ch in names(chr_idx)) chr_idx[[ch]][, cs := cumsum(ldw)]

## sum and count of ld_w over [from, to] for vectors of intervals on one chromosome
win <- function(d, from, to) {
  lo <- findInterval(from - 1, d$Pos); hi <- findInterval(to, d$Pos)
  cs <- c(0, d$cs)
  list(s = cs[hi + 1] - cs[lo + 1], n = hi - lo)
}

sets <- list(
  marine_freshwater = as.data.table(sc$ecopeaks)[, .(Chr = chr, from = start, to = end)],
  consensus = sc$consensus$test$regions[, .(Chr, from, to)],
  simes     = sc$simes$test$regions[, .(Chr, from, to)]
)

set.seed(1L)
res <- rbindlist(lapply(names(sets), function(nm) {
  iv <- copy(sets[[nm]])[Chr %chin% names(chr_idx)]
  iv[, span := to - from]
  ## observed per-interval mean ld_w
  iv[, c("s", "n") := {
    d <- chr_idx[[.BY$Chr]]; w <- win(d, from, to); list(w$s, w$n)
  }, by = Chr]
  iv <- iv[n > 0]
  iv[, m := s / n]
  obs <- mean(iv$m)
  ## rotation null: uniform start on the same chromosome, span preserved
  null <- numeric(N_ROT)
  per_chr <- split(iv, by = "Chr")
  acc <- matrix(NA_real_, N_ROT, nrow(iv)); col <- 0L
  for (ch in names(per_chr)) {
    d <- chr_idx[[ch]]; x <- per_chr[[ch]]; lo <- min(d$Pos); hi <- max(d$Pos)
    for (k in seq_len(nrow(x))) {
      col <- col + 1L
      st <- runif(N_ROT, lo, max(lo, hi - x$span[k]))
      w <- win(d, st, st + x$span[k])
      acc[, col] <- ifelse(w$n > 0, w$s / pmax(w$n, 1), NA_real_)
    }
  }
  null <- rowMeans(acc, na.rm = TRUE)
  ## local peak: inside versus flanks
  iv[, flank := {
    d <- chr_idx[[.BY$Chr]]
    l <- win(d, pmax(from - FLANK, 0), from - 1); r <- win(d, to + 1, to + FLANK)
    (l$s + r$s) / pmax(l$n + r$n, 1)
  }, by = Chr]
  k_up <- sum(iv$m > iv$flank); n_pk <- nrow(iv)
  data.table(set = nm, n_intervals = nrow(sets[[nm]]), n_with_markers = n_pk,
             mean_ldw = obs, null_mean = mean(null), fold = obs / mean(null),
             p_rotation = (1 + sum(null >= obs)) / (N_ROT + 1),
             genome_mean = mean(mk$ldw),
             median_inside_over_flank = median(iv$m / iv$flank),
             n_above_flank = k_up, p_sign = binom.test(k_up, n_pk, alternative = "greater")$p.value)
}))
print(res)

## ---- 3. floor control: significant vs eligible non-significant units --------
cl <- as.data.table(b$stage1$clusters)[n_snps >= SIZE_FLOOR]
cl[, `:=`(from = vapply(members, function(x) min(as.numeric(sub(".*:", "", x))), 0),
          to   = vapply(members, function(x) max(as.numeric(sub(".*:", "", x))), 0))]
ldw_of <- setNames(as.data.table(b$map)$ld_w_095, as.data.table(b$map)$marker)
cl[, ldw := vapply(members, function(x) mean(ldw_of[x], na.rm = TRUE), 0)]
units <- rbindlist(lapply(c("consensus", "simes"), function(st) {
  u <- copy(sc[[st]]$test$units)[, statistic := st]
  merge(u, cl[, .(Chr = as.character(Chr), from, to, n_markers = n_snps, ldw)],
        by = c("Chr", "from", "to", "n_markers"), all.x = TRUE)
}))
stopifnot(!anyNA(units$ldw))
ctrl <- units[, {
  f <- lm(log(ldw) ~ log(n_markers) + significant)
  co <- summary(f)$coefficients["significantTRUE", ]
  .(n_sig = sum(significant), n_nonsig = sum(!significant),
    mean_ldw_sig = mean(ldw[significant]), mean_ldw_nonsig = mean(ldw[!significant]),
    p_wilcox = wilcox.test(ldw[significant], ldw[!significant])$p.value,
    size_adj_ratio = exp(co[["Estimate"]]), p_size_adj = co[["Pr(>|t|)"]])
}, by = statistic]
print(ctrl)

out <- stage_dir(STAGE); dir.create(out, showWarnings = FALSE, recursive = TRUE)
fwrite(res, file.path(out, "ldw_peaks.tsv"), sep = "\t")
fwrite(ctrl, file.path(out, "ldw_units.tsv"), sep = "\t")
say("wrote %s\n", out)

## ---- 4. reverse direction: how many ld_w peaks lie in outlier regions? -------
## Peaks = runs of markers at or above a genome-wide ld_w percentile, split at
## gaps > GAP bp, with >= 3 markers. Overlap with the marine-freshwater regions
## (and with reported consensus regions) is compared with random placement of
## each peak on its own chromosome (200 replicates).
LEN <- as.data.table(sc$chrom_lengths)[, .(Chr = as.character(chr), len)]
overlapping <- function(p, iv) {
  iv <- copy(iv); setkey(iv, Chr, from, to)
  unique(foverlaps(p[, .(Chr, from, to, pid)], iv, nomatch = NULL)$pid)
}
peaks <- rbindlist(lapply(c(0.95, 0.99), function(q) rbindlist(lapply(c(1e4, 5e4), function(gap) {
  thr <- quantile(mk$ldw, q)
  hi <- mk[ldw >= thr][, grp := cumsum(c(TRUE, diff(Pos) > gap | Chr[-1] != Chr[-.N]))]
  pk <- hi[, .(from = min(Pos), to = max(Pos), n = .N), by = .(Chr, grp)][n >= 3][, pid := .I]
  set.seed(1L)
  rnd <- replicate(200, {
    r <- pk[LEN, on = "Chr", nomatch = NULL][, s := runif(.N, 0, len - (to - from))][, `:=`(from = s, to = s + (to - from))]
    length(overlapping(r, sets$marine_freshwater)) / nrow(r)
  })
  data.table(percentile = q, gap_bp = gap, n_peaks = nrow(pk), median_span_kb = median(pk$to - pk$from) / 1e3,
             pct_in_marine_freshwater = 100 * length(overlapping(pk, sets$marine_freshwater)) / nrow(pk),
             pct_expected_random = 100 * mean(rnd),
             pct_in_reported_consensus = 100 * length(overlapping(pk, sets$consensus)) / nrow(pk))
}))))
print(peaks)
fwrite(peaks, file.path(out, "ldw_peaks_in_regions.tsv"), sep = "\t")
say("wrote %s\n", file.path(out, "ldw_peaks_in_regions.tsv"))

## ---- 5. low-ld_w outlier regions --------------------------------------------
## For each interval, the percentile of its mean ld_w among 2,000 random
## placements of the same span on the same chromosome. "Low" = below the
## median of its own random placements, i.e. no elevation at all.
interval_pctl <- function(iv) {
  iv <- copy(iv)[Chr %chin% names(chr_idx)]
  iv[, c("m", "pctl") := {
    d <- chr_idx[[.BY$Chr]]; lo <- min(d$Pos); hi <- max(d$Pos)
    res <- mapply(function(f, t) {
      w <- win(d, f, t); if (w$n == 0) return(c(NA_real_, NA_real_))
      st <- runif(2000, lo, max(lo, hi - (t - f))); r <- win(d, st, st + (t - f))
      c(w$s / w$n, mean((r$s / pmax(r$n, 1))[r$n > 0] < w$s / w$n))
    }, from, to)
    list(res[1, ], res[2, ])
  }, by = Chr][!is.na(m)]
}
set.seed(2L)
mf <- interval_pctl(sets$marine_freshwater)
low <- rbindlist(lapply(c("consensus", "simes"), function(st) {
  rg <- interval_pctl(sets[[st]])
  e <- copy(sets$marine_freshwater); setkey(e, Chr, from, to); r2 <- copy(rg); setkey(r2, Chr, from, to)
  on <- unique(foverlaps(r2, e, nomatch = NULL)[, paste(Chr, i.from, i.to)])
  rg[, on_mf := paste(Chr, from, to) %chin% on]
  data.table(set = st, n = nrow(rg), n_low = sum(rg$pctl < 0.5), n_on_mf = sum(rg$on_mf),
             n_on_mf_low = sum(rg$on_mf & rg$pctl < 0.5))
}))
low <- rbind(data.table(set = "marine_freshwater", n = nrow(mf), n_low = sum(mf$pctl < 0.5),
                        n_on_mf = NA_integer_, n_on_mf_low = NA_integer_,
                        pct_below_p75 = 100 * mean(mf$pctl < 0.75), pct_above_p90 = 100 * mean(mf$pctl > 0.9)),
             low, fill = TRUE)
low[, pct_low := 100 * n_low / n]
print(low)
fwrite(low, file.path(out, "ldw_low_regions.tsv"), sep = "\t")
say("wrote %s\n", file.path(out, "ldw_low_regions.tsv"))

## ---- 6. would an ld_w-restricted analysis have missed reported regions? -----
## For each reported region, whether any marker of its significant Stage-1
## units reaches the genome-wide 75th/90th/95th ld_w percentile. Regions with
## none would be lost entirely if analyses were restricted to high-ld_w markers.
ldw_by_marker <- as.data.table(b$map)[!is.na(ld_w_095), .(marker, Pos, ldw = ld_w_095)]
cl_pos <- as.data.table(b$stage1$clusters)[n_snps >= SIZE_FLOOR]
cl_pos[, `:=`(from = vapply(members, function(x) min(as.numeric(sub(".*:", "", x))), 0),
              to   = vapply(members, function(x) max(as.numeric(sub(".*:", "", x))), 0))]
filt <- rbindlist(lapply(c("consensus", "simes"), function(st) {
  u <- merge(sc[[st]]$test$units[significant == TRUE],
             cl_pos[, .(Chr = as.character(Chr), from, to, n_markers = n_snps, members)],
             by = c("Chr", "from", "to", "n_markers"))
  stopifnot(nrow(u) == sum(sc[[st]]$test$units$significant))
  mm <- merge(u[, .(marker = unlist(members)), by = Chr], ldw_by_marker, by = "marker")
  rg <- copy(sc[[st]]$test$regions)[, rid := .I]; setkey(rg, Chr, from, to)
  x <- foverlaps(mm[, .(Chr, from = Pos, to = Pos, ldw)], rg[, .(Chr, from, to, rid)], nomatch = NULL)
  rbindlist(lapply(c(0.75, 0.90, 0.95), function(q) {
    thr <- quantile(ldw_by_marker$ldw, q)
    k <- x[, .(any_hi = any(ldw >= thr)), by = rid]
    data.table(set = st, percentile = q, ldw_threshold = thr, n_regions = nrow(rg),
               n_without_high_ldw_marker = sum(!k$any_hi))
  }))
}))
filt[, pct_without := 100 * n_without_high_ldw_marker / n_regions]
print(filt)
fwrite(filt, file.path(out, "ldw_filter_missed_regions.tsv"), sep = "\t")
say("wrote %s\n", file.path(out, "ldw_filter_missed_regions.tsv"))
