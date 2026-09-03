## =====================================================================
## module_C2 / R/10_support.R    [Questions 3, 4A, 6]
##
## Detection support of reference-point regions across each null-admissible grid.
##
##   D_r = #{admissible cells recovering r} / |G_adm|
##
## The denominator is ALL admissible cells, including admissible cells in which no
## empirical region is detected. The conditional "usable cells" denominator is
## computed only as a comparator.
##
## Matching is by MARKER MEMBERSHIP: an anchor r is recovered in a cell when some
## called region R satisfies |A n R| / |A| >= thr (rec25 / rec50 / rec75), with a
## reciprocal variant |A n R| / |R| >= thr as a sensitivity check.
## Splits    : the anchor's markers land in several called regions; the best single
##             region must clear thr on its own, so a fully split anchor fails.
## Mergers   : the anchor is swallowed by a larger component; recovery = 1, so a
##             merge counts as detection (the reciprocal variant penalises it).
## Boundaries: shifts are absorbed -- only shared MEMBERSHIP counts, never spans.
##
## Rscript module_C2/R/10_support.R
## =====================================================================
source("module_C2/R/00_helpers.R")
core <- readRDS(file.path(C2$CACHE, "grid_core.rds")); B <- core$B
KEEP <- readRDS(file.path(C2$CACHE, "admissible_cells.rds"))
NL <- fread(file.path(C2$RES, "null_grid_diagnostics.csv"))
NL[, tau := C2$TAUS[match(round(tau, 6), round(C2$TAUS, 6))]]
FULL  <- CJ(tau = C2$TAUS, lmin = C2$LMINS)
NCELL <- nrow(FULL)

## ---- reference-point regions (coordinate display only) --------------
anc <- c2_anchors(core)                      # tau_C = 0.05, l_min = 3 -> 17 regions
c2_msg("[0] reference point (%.2f, %d) supplies %d displayed regions; sizes %s\n",
       C2$REF_TAU, C2$OP_LMIN, nrow(anc$tab), paste(sort(anc$tab$size), collapse = " "))
c2_msg("    NB the reference point is NOT a grid cell and does not define admissibility.\n")

## ---- detection + significance over an arbitrary cell set ------------
support_grid <- function(cells, thr, recip = FALSE) {
  A <- anc$mk; alab <- names(A); out <- list()
  for (tk in unique(as.character(cells$tau))) {
    ct <- core$by_tau[[tk]]; tau <- as.numeric(tk)
    lms <- cells[as.character(tau) == tk]$lmin
    nO <- nrow(ct$obs)
    rec <- rcp <- matrix(0, length(A), max(nO, 1L))
    if (nO) for (j in seq_len(nO)) { Rj <- ct$mk[[j]]
      for (i in seq_along(A)) { n <- length(intersect(A[[i]], Rj))
        if (n) { rec[i, j] <- n / length(A[[i]]); rcp[i, j] <- n / length(Rj) } } }
    for (lm in lms) {
      keep <- which(ct$obs$size >= lm)
      if (!length(keep)) { out[[length(out) + 1L]] <- data.table(
        tau = tau, lmin = lm, label = alab, detected = FALSE, sig = FALSE,
        n_match = 0L, n_tested = 0L); next }
      O <- ct$obs[keep]; pq <- c2_emp_pq(O, ct$surr[b <= B & size >= lm], B)
      sr <- rec[, keep, drop = FALSE]; sp <- rcp[, keep, drop = FALSE]
      hit <- sr >= thr & (if (recip) sp >= thr else TRUE)
      best <- vapply(seq_along(A), function(i) { w <- which(hit[i, ])
        if (!length(w)) NA_integer_ else w[which.max(sr[i, w])] }, integer(1))
      out[[length(out) + 1L]] <- data.table(
        tau = tau, lmin = lm, label = alab, detected = !is.na(best),
        sig = !is.na(best) & !is.na(pq$q_R[best]) & pq$q_R[best] < C2$FDR,
        n_match = rowSums(hit), n_tested = nrow(O))
    }
  }
  rbindlist(out)
}

## ---- compute once per (grid, thr) -----------------------------------
GRIDS <- list(full = FULL, P20 = KEEP$P20, P10 = KEEP$P10, P05 = KEEP$P05,
              U10 = KEEP$U10, C001 = KEEP$C001)
THRS  <- c(rec25 = 0.25, rec50 = 0.50, rec75 = 0.75)
CACHE <- file.path(C2$CACHE, "support.rds")
if (file.exists(CACHE)) { SUP <- readRDS(CACHE); c2_msg("[1] loaded cached support\n") } else {
  SUP <- list()
  for (gn in names(GRIDS)) for (tn in names(THRS)) {
    if (!nrow(GRIDS[[gn]])) next
    SUP[[paste(gn, tn, sep = ".")]] <- support_grid(GRIDS[[gn]], THRS[[tn]])
    c2_msg("[1] %s x %s done (%d cells)\n", gn, tn, nrow(GRIDS[[gn]]))
  }
  SUP[["full.rec50recip"]] <- support_grid(FULL, 0.50, recip = TRUE)
  saveRDS(SUP, CACHE)
}
sc <- function(key, n) SUP[[key]][, .(D = sum(detected) / n, Q = sum(sig) / n), by = label]

## ---- 2. does significance differ from detection? --------------------
c2_msg("\n[2] detection vs BH-significance, every grid x threshold:\n")
dq <- rbindlist(lapply(names(SUP), function(k) { g <- SUP[[k]]
  data.table(config = k, cells_detected = sum(g$detected), cells_significant = sum(g$sig),
             detected_not_sig = sum(g$detected & !g$sig)) }))
print(dq)
c2_msg("    -> detected-but-not-significant cells, total across ALL configurations: %d\n",
       sum(dq$detected_not_sig))
c2_msg("       Significance support is REDUNDANT with detection support; it is reported\n")
c2_msg("       below for completeness but is not independent evidence.\n")

## ---- 3. support tables ----------------------------------------------
res <- list()
for (gn in names(GRIDS)) { n <- nrow(GRIDS[[gn]]); if (!n) next
  for (tn in names(THRS)) { k <- paste(gn, tn, sep = ".")
    s <- sc(k, n); s[, `:=`(grid = gn, thr = tn, n_cells = n)]; res[[k]] <- s } }
SS <- rbindlist(res)
fwrite(SS, file.path(C2$RES, "support_by_grid_threshold.csv"))

W <- dcast(SS[thr == "rec50"], label ~ grid, value.var = "D")
setcolorder(W, c("label", "full", "P20", "P10", "P05", "U10", "C001"))
setorder(W, -full)
c2_msg("\n[3] detection support D_r at rec50, by admissible grid (denominator = ALL cells in that grid):\n")
print(W[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
c2_msg("\n    grid sizes: %s\n", paste(sprintf("%s=%d", names(GRIDS),
       vapply(GRIDS, nrow, 0L)), collapse = " "))

## zero-support regions
c2_msg("\n[3] regions with ZERO support, by grid (rec50):\n")
for (gn in setdiff(names(GRIDS), character(0))) { if (!nrow(GRIDS[[gn]])) next
  z <- W$label[W[[gn]] == 0]
  c2_msg("    %-5s : %2d of %d  %s\n", gn, length(z), nrow(W),
         if (length(z)) paste(substr(z, 1, 22), collapse = " | ") else "") }

## ---- 4. the old conditional denominator, as comparator --------------
sg <- readRDS(file.path(C2$CACHE, "grid_scored.rds"))
usable <- sg$resF[n_sig > 0, .N]
Wc <- sc("full.rec50", NCELL)[, .(label, D_fullgrid = D)]
Wc[, D_conditional := SUP[["full.rec50"]][, sum(sig), by = label]$V1 / usable]
c2_msg("\n[4] old conditional denominator |U| = %d vs full grid |G| = %d:\n", usable, NCELL)
c2_msg("    Spearman = %.3f (a pure rescaling by |G|/|U| = %.2f; ranking identical)\n",
       stats::cor(Wc$D_fullgrid, Wc$D_conditional, method = "spearman"), NCELL / usable)

## ---- 5. required comparisons: rank correlation + top overlap --------
vec <- function(gn, tn) { s <- sc(paste(gn, tn, sep = "."), nrow(GRIDS[[gn]]))
                          stats::setNames(s$D, s$label) }
base <- vec("full", "rec50")
cmp <- rbindlist(lapply(names(GRIDS), function(gn) { if (!nrow(GRIDS[[gn]])) return(NULL)
  rbindlist(lapply(names(THRS), function(tn) { v <- vec(gn, tn); a <- c2_agree(base, v, k = 5L)
    data.table(grid = gn, thr = tn, n_cells = nrow(GRIDS[[gn]]),
               spearman_vs_full_rec50 = round(a$spearman, 3),
               top5_overlap = round(a$top_k_overlap, 2),
               n_zero = sum(v == 0), mean_D = round(mean(v), 4)) })) }))
## reciprocal variant
vr <- SUP[["full.rec50recip"]][, .(D = sum(detected) / NCELL), by = label]
ar <- c2_agree(base, stats::setNames(vr$D, vr$label), k = 5L)
cmp <- rbind(cmp, data.table(grid = "full", thr = "rec50recip", n_cells = NCELL,
             spearman_vs_full_rec50 = round(ar$spearman, 3), top5_overlap = round(ar$top_k_overlap, 2),
             n_zero = sum(vr$D == 0), mean_D = round(mean(vr$D), 4)))
fwrite(cmp, file.path(C2$RES, "support_comparisons.csv"))
c2_msg("\n[5] rank agreement with full-grid rec50:\n"); print(cmp)

## ---- 6. reference-point sensitivity ---------------------------------
## Do conclusions change if a NEIGHBOURING reference point supplies the regions?
c2_msg("\n[6] reference-point sensitivity (regions redefined at neighbouring points,\n")
c2_msg("    then scored on the SAME P20 admissible grid, rec50):\n")
REFS <- list(c(0.05, 3), c(0.04, 3), c(0.06, 3), c(0.05, 2), c(0.05, 5))
rs <- list()
for (r in REFS) {
  a2 <- c2_anchors(core, tau = r[1], lmin = as.integer(r[2]))
  old <- anc; anc <<- a2
  s <- support_grid(KEEP$P20, 0.50)[, .(D = sum(detected) / nrow(KEEP$P20)), by = label]
  anc <<- old
  rs[[length(rs) + 1L]] <- data.table(ref = sprintf("(%.2f, %d)", r[1], r[2]),
    n_regions = nrow(a2$tab), mean_D = round(mean(s$D), 4),
    n_zero = sum(s$D == 0), max_D = round(max(s$D), 3),
    top_region = s$label[which.max(s$D)])
}
RS <- rbindlist(rs); fwrite(RS, file.path(C2$RES, "reference_point_sensitivity.csv"))
print(RS)

## ---- 7. validation: are lenient-only regions in high-p_null cells? --
## For each anchor, the mean p_null_any of the cells in which it is detected.
## This CORROBORATES the null rule; it was not used to build it.
g <- SUP[["full.rec50"]]
gm <- merge(g[detected == TRUE], NL[, .(tau, lmin, p_null_any)], by = c("tau", "lmin"))
vv <- gm[, .(cells_detected = .N, mean_p_null_of_detecting_cells = round(mean(p_null_any), 3),
             frac_in_P20 = round(mean(paste(tau, lmin) %in% KEEP$P20[, paste(tau, lmin)]), 3)), by = label]
setorder(vv, -mean_p_null_of_detecting_cells)
fwrite(vv, file.path(C2$RES, "lenient_cell_validation.csv"))
c2_msg("\n[7] validation -- where are the detecting cells? (rec50, full grid)\n")
print(vv)
zeroP20 <- W$label[W$P20 == 0]
c2_msg("\n    anchors with zero P20 support are detected in cells of mean p_null = %.3f;\n",
       vv[label %in% zeroP20, mean(mean_p_null_of_detecting_cells)])
c2_msg("    anchors retaining P20 support: mean p_null = %.3f\n",
       vv[!label %in% zeroP20, mean(mean_p_null_of_detecting_cells)])
c2_msg("\n[8] wrote support_by_grid_threshold.csv, support_comparisons.csv,\n")
c2_msg("    reference_point_sensitivity.csv, lenient_cell_validation.csv\n")
