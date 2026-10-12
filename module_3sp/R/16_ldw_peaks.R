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
