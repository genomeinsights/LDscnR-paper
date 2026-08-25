## =====================================================================
## module_sim_LDscnR / check_inputs.R
##
## Validate a (panel, pvals) pair against the contract analyse_one_dataset.R
## expects. Run this BEFORE handing data over: every check here is something
## that would otherwise fail deep inside the analysis, or -- worse -- not fail at
## all and quietly produce wrong regions.
##
##   Rscript module_sim_LDscnR/check_inputs.R <panel.rds> <pvals.rds>
##
## Exit status 0 = pass, 1 = fail. Warnings do not fail the run.
## =====================================================================
suppressMessages({ library(data.table) })
`%||%` <- function(a, b) if (is.null(a)) b else a

a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 2) stop("usage: check_inputs.R <panel.rds> <pvals.rds>")
ERR <- 0L; WARN <- 0L
ok   <- function(m) cat(sprintf("  [ ok ] %s\n", m))
bad  <- function(m) { cat(sprintf("  [FAIL] %s\n", m)); ERR <<- ERR + 1L }
warn <- function(m) { cat(sprintf("  [warn] %s\n", m)); WARN <<- WARN + 1L }

panel <- readRDS(a[1]); P <- readRDS(a[2])

cat("== panel:", basename(a[1]), "==\n")
for (f in c("GTs", "map", "ld_ws", "decay_sum"))
  if (is.null(panel[[f]])) bad(sprintf("missing `%s`", f))
if (ERR) { cat("\nFAILED: panel is missing required fields\n"); quit(status = 1L) }

GTs <- panel$GTs; map <- as.data.table(panel$map); ld_ws <- panel$ld_ws; ds <- as.data.table(panel$decay_sum)
nm <- ncol(GTs); ni <- nrow(GTs)

if (!all(c("marker", "Chr", "Pos") %in% names(map))) bad("map needs marker, Chr, Pos") else
  ok(sprintf("map has marker/Chr/Pos (%d rows)", nrow(map)))
if (nrow(map) != nm) bad(sprintf("map has %d rows but GTs has %d columns", nrow(map), nm)) else
  ok(sprintf("GTs is %d individuals x %d markers", ni, nm))
if (is.null(colnames(GTs))) bad("GTs has no colnames (marker IDs)") else
  if (!identical(as.character(colnames(GTs)), as.character(map$marker)))
    bad("colnames(GTs) is not identical to map$marker, in order") else ok("colnames(GTs) == map$marker, in order")

if (is.null(rownames(ld_ws))) bad("ld_ws has no rownames") else
  if (!identical(as.character(rownames(ld_ws)), as.character(map$marker)))
    bad("rownames(ld_ws) is not identical to map$marker, in order") else ok("rownames(ld_ws) == map$marker, in order")
if (is.null(colnames(ld_ws))) bad("ld_ws has no colnames (the rho windows)") else
  ok(sprintf("ld_ws has %d rho window(s): %s", ncol(ld_ws),
             paste(utils::head(colnames(ld_ws), 4), collapse = ", ")))
if (anyNA(ld_ws)) warn(sprintf("ld_ws has %d NA values", sum(is.na(ld_ws))))

## decay_sum drives the per-chromosome r^2 link. A chromosome missing from it, or
## carrying an NA fit, does NOT error: ld_edges() silently falls back to r2 = 0.5,
## far stricter than a typical fitted value (~0.27), so that chromosome quietly
## clusters almost nothing. This is the check most worth having.
if (!all(c("Chr", "b") %in% names(ds))) bad("decay_sum needs at least Chr and b") else {
  miss <- setdiff(unique(as.character(map$Chr)), as.character(ds$Chr))
  if (length(miss)) bad(sprintf("decay_sum is missing %d chromosome(s): %s -- ld_edges() would silently use r2 = 0.5 for these",
                                length(miss), paste(utils::head(miss, 5), collapse = ", ")))
  else ok(sprintf("decay_sum covers all %d chromosomes in map", uniqueN(map$Chr)))
  cc <- if ("c" %in% names(ds)) ds$c else rep(1, nrow(ds))
  nabad <- sum(is.na(ds$b) | is.na(cc))
  if (nabad) bad(sprintf("%d chromosome(s) have an NA decay fit -- same silent r2 = 0.5 fallback", nabad))
  else ok("no NA decay fits")
}

has_truth <- "true_QTN" %in% names(map)
cat(sprintf("  [info] truth column: %s\n",
            if (has_truth) sprintf("present, %d true_QTN -- evaluation section will run", sum(map$true_QTN %in% TRUE))
            else "absent -- treated as real data, evaluation skipped"))

cat("\n== pvals:", basename(a[2]), "==\n")
for (f in c("p_obs", "p_perm")) if (is.null(P[[f]])) bad(sprintf("missing `%s`", f))
if (ERR) { cat("\nFAILED\n"); quit(status = 1L) }

p_obs <- P$p_obs; pp <- P$p_perm
if (length(p_obs) != nm) bad(sprintf("p_obs has %d values, expected %d", length(p_obs), nm)) else
  ok(sprintf("p_obs has %d values", length(p_obs)))
if (is.null(names(p_obs))) warn("p_obs is unnamed -- legal, but names are the only guard against misordering") else
  if (!identical(as.character(names(p_obs)), as.character(map$marker)))
    bad("names(p_obs) is not identical to map$marker, in order") else ok("names(p_obs) == map$marker, in order")
rng <- range(p_obs, na.rm = TRUE)
if (rng[1] < 0 || rng[2] > 1) {
  bad(sprintf("p_obs outside [0,1]: %.3g to %.3g -- are these p-values, not F or -log10p?", rng[1], rng[2]))
} else ok(sprintf("p_obs in [0,1] (%.3g to %.3g)", rng[1], rng[2]))
if (anyNA(p_obs)) warn(sprintf("p_obs has %d NA -- tolerated (they never become hits), but check they are expected", sum(is.na(p_obs))))

isM <- is.matrix(pp)
B <- if (isM) ncol(pp) else length(pp)
if (isM) {
  if (nrow(pp) != nm) bad(sprintf("p_perm matrix has %d rows, expected %d", nrow(pp), nm)) else
    ok(sprintf("p_perm is a %d x %d markers-by-B matrix", nrow(pp), B))
  if (!is.null(rownames(pp)) && !identical(as.character(rownames(pp)), as.character(map$marker)))
    bad("rownames(p_perm) is not identical to map$marker, in order")
} else if (is.list(pp)) {
  len <- vapply(pp, length, integer(1))
  if (any(len != nm)) bad(sprintf("%d surrogate(s) have the wrong length (first bad: %d values)",
                                  sum(len != nm), len[which(len != nm)[1]]))
  else ok(sprintf("p_perm is a list of %d surrogate vectors, all length %d", B, nm))
} else bad("p_perm must be a markers-by-B matrix or a list of vectors")

if (B < 20L) {
  bad(sprintf("B = %d is too small: the smallest attainable region p is 1/(1+B) = %.3f", B, 1/(1+B)))
} else if (B < 100L) {
  warn(sprintf("B = %d gives a p-floor of %.4f; B >= 100 recommended", B, 1/(1+B)))
} else ok(sprintf("B = %d, p-floor = %.4f", B, 1/(1+B)))

cat(sprintf("  [info] engine: %s | basis: %s\n", P$engine %||% "(unset)", P$basis %||% "(unset)"))
if (is.null(P$engine) || is.null(P$basis))
  warn("engine/basis labels are unset -- they are carried into every printed output; set them")

cat(sprintf("\n%s  (%d error(s), %d warning(s))\n",
            if (ERR) "FAILED" else "PASSED", ERR, WARN))
quit(status = if (ERR) 1L else 0L)
