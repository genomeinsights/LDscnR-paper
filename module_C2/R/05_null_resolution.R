## =====================================================================
## module_C2 / R/05_null_resolution.R    [Question 6]
##
## Monte Carlo resolution of the per-cell empirical p-values and of BH within a
## cell. Reports B, the floor 1/(B+1), ties at the floor, whether a single region
## at the floor can pass BH given the cell's hypothesis count, and the rank
## stability of the stability score under subsampled B.
##
## Rscript module_C2/R/05_null_resolution.R
## =====================================================================
source("module_C2/R/00_helpers.R")
core <- readRDS(file.path(C2$CACHE, "grid_core.rds")); B <- core$B
anc  <- c2_anchors(core)
sg   <- readRDS(file.path(C2$CACHE, "grid_scored.rds"))
NCELL <- length(C2$TAUS) * length(C2$LMINS)

L <- sg$Lfull
FLOOR <- 1 / (B + 1)
c2_msg("[1] B = %d ; p floor = 1/(B+1) = %.5f ; (prototype cap B=100 -> %.5f)\n",
       B, FLOOR, 1 / 101)

## ---- per-cell BH behaviour, and the role of ties at the floor -------
## A naive bound would say: since every p >= 1/(B+1), a cell testing m regions can
## reject at all only if 1/(B+1) <= alpha/m, i.e. m <= alpha(B+1). That bound is
## WRONG here, because it assumes distinct p-values. BH is a step-up procedure: if
## k regions are TIED at the floor p0, it rejects all of them as soon as
## p0 <= alpha * k / m. Ties at the floor rescue each other, so the binding
## quantity is the FRACTION of the cell's regions sitting at the floor, not m.
cells <- L[!is.na(emp_p), .(m = .N, n_floor = sum(emp_p <= FLOOR + 1e-12),
                            min_p = min(emp_p), n_sig = sum(q_R < C2$FDR)), by = .(tau, lmin)]
cells[, `:=`(frac_floor = n_floor / m,
             naive_bound_ok = m <= C2$FDR * (B + 1),          # the WRONG bound
             tie_bound_ok   = FLOOR <= C2$FDR * n_floor / m)] # the correct one
c2_msg("[2] naive (distinct-p) bound m <= alpha*(B+1) = %.1f would block %d of %d productive cells;\n",
       C2$FDR * (B + 1), cells[!(naive_bound_ok), .N], nrow(cells))
c2_msg("    but %d of those actually DO produce significant regions, because their p-values\n",
       cells[!(naive_bound_ok) & n_sig > 0, .N])
c2_msg("    are tied at the floor and BH rejects tied sets together.\n")
c2_msg("[2] tie-aware condition (p0 <= alpha*k/m) predicts the cell's outcome correctly in %d of %d cells\n",
       cells[tie_bound_ok == (n_sig > 0), .N], nrow(cells))
c2_msg("    -> the operative driver of per-cell significance is the FRACTION of regions at the\n")
c2_msg("       p-value floor (median %.2f, range %.2f-%.2f), i.e. how many regions the null\n",
       stats::median(cells$frac_floor), min(cells$frac_floor), max(cells$frac_floor))
c2_msg("       fails to rebuild AT ALL -- not the strength of any individual region.\n")
c2_msg("[2] a cell where EVERY region is at the floor rejects every one of them: %d such cells\n",
       cells[frac_floor == 1, .N])
print(cells[, .(cells = .N, m_median = as.integer(stats::median(m)), m_max = max(m),
                median_frac_floor = round(stats::median(frac_floor), 2),
                cells_with_sig = sum(n_sig > 0)), by = lmin][order(lmin)])
fwrite(cells, file.path(C2$RES, "null_resolution.csv"))

## ---- rank stability under subsampled B ------------------------------
## rescore the anchors using only the first b surrogates
set.seed(1)
BS <- c(25L, 50L, 100L, 150L, B)
ref <- NULL; rows <- list()
for (b in BS) {
  g <- c2_anchor_grid(core, anc, b, rule = "recover", thr = 0.5)
  s <- g[, .(S_G = sum(sig) / NCELL, D = sum(detected) / NCELL), by = label]
  v <- stats::setNames(s$S_G, s$label)
  if (is.null(ref)) ref <- v
  a <- c2_agree(ref, v, k = 5L)
  rows[[length(rows) + 1L]] <- data.table(B = b, mean_S_G = round(mean(v), 4),
                                          n_pos = sum(v > 0),
                                          spearman_vs_Bmin = round(a$spearman, 3),
                                          top5_vs_Bmin = round(a$top_k_overlap, 2))
  c2_msg("[4] B = %3d : mean S_G = %.4f ; anchors with S_G>0 = %2d\n", b, mean(v), sum(v > 0))
  if (b == B) full <- v
}
sub <- rbindlist(rows)
## agreement of each B against the FULL B (the quantity that matters)
for (i in seq_len(nrow(sub))) {
  g <- c2_anchor_grid(core, anc, sub$B[i], rule = "recover", thr = 0.5)
  v <- stats::setNames(g[, sum(sig) / NCELL, by = label]$V1, g[, sum(sig) / NCELL, by = label]$label)
  a <- c2_agree(full, v, k = 5L)
  sub[i, `:=`(spearman_vs_full = round(a$spearman, 3), top5_vs_full = round(a$top_k_overlap, 2))]
}
print(sub)
fwrite(sub, file.path(C2$RES, "null_resolution_subsampling.csv"))
c2_msg("\n[4] wrote results/null_resolution.csv + null_resolution_subsampling.csv\n")
