## =====================================================================
## patch_bundle_va.R -- correct Va / sum_Va in EXISTING bundles, in place.
##
## get_va() used the allele effect a where the SUBSTITUTION effect is 2a, so
## stored Va is a quarter of its true value (see get_va() for the phenotype
## reconstruction that established the coding). p_Va is a ratio and was never
## affected, so flag_true_qtns() and every truth definition built on it are
## unchanged -- this touches the absolute Va and sum_Va only.
##
## Rebuilding the bundles instead was not an option: the decay fit subsamples
## pairs, so a rebuilt bundle would carry different ld_ws, GRM and p-values and
## would silently invalidate every panel/pvals file already delivered from it.
## Va is recomputed from the bundle's OWN GTs and allelic_values, so nothing
## else in the file changes.
##
## Refuses to touch a bundle whose p_Va would move -- if that ever happens the
## assumption above is wrong and the file should be looked at, not rewritten.
##
##   Rscript patch_bundle_va.R <bundle-dir> [--dry-run]
## =====================================================================
suppressMessages(library(data.table))
a   <- commandArgs(trailingOnly = TRUE)
DIR <- if (length(a) >= 1) a[1] else stop("usage: patch_bundle_va.R <bundle-dir> [--dry-run]")
DRY <- "--dry-run" %in% a

fs <- list.files(DIR, pattern = "[.]rds$", full.names = TRUE)
cat(sprintf("%s %d bundle(s) in %s\n", if (DRY) "DRY RUN over" else "patching", length(fs), DIR))
n_ok <- n_skip <- n_bad <- 0L

for (f in fs) {
  d <- readRDS(f); m <- as.data.table(d$map)
  qtn <- which(m$type == "QTN")
  if (!length(qtn) || all(is.na(m$Va))) { n_skip <- n_skip + 1L; next }

  ## unname(): colMeans() carries marker names, and comparing a named vector with
  ## an unnamed one makes all.equal() report a difference that is purely cosmetic
  p   <- unname(colMeans(d$GTs[, m$marker[qtn], drop = FALSE]) / 2)
  aa  <- m$allelic_values[qtn]
  new <- 2 * p * (1 - p) * (2 * aa)^2
  old <- m$Va[qtn]

  ## the factor must be exactly 4, and p_Va must not move
  ratio <- new / old
  p_old <- unname(old / sum(old)); p_new <- unname(new / sum(new))
  if (any(is.finite(ratio) & abs(ratio - 4) > 1e-6) ||
      !isTRUE(all.equal(p_old, p_new, tolerance = 1e-10))) {
    cat(sprintf("  !! %s: unexpected change, left untouched\n", basename(f))); n_bad <- n_bad + 1L; next
  }
  if (!DRY) {
    m[qtn, Va := unname(new)][qtn, sum_Va := sum(new)][qtn, p_Va := unname(new / sum(new))]
    d$map <- m
    saveRDS(d, f)
  }
  n_ok <- n_ok + 1L
}
cat(sprintf("%s: %d | no QTN to patch: %d | refused: %d\n",
            if (DRY) "would patch" else "patched", n_ok, n_skip, n_bad))
