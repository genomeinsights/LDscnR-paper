## =====================================================================
## module_C2 / R/02_compare_denominators.R    [Question 2]
##
## Compare denominator definitions on the FIXED anchor loci (so the locus set is
## constant across variants and only the denominator changes):
##   A  S^U = #{usable cells where l significant} / |U|        (current proposal)
##   B  S^G = #{all cells where l significant}    / |G|        (leading candidate)
##   C  A = |U|/|G| reported alongside S^U                     (availability)
##   D  D_l (detection), Q_l (significance), Q_l/D_l           (decomposition)
## and the prototype's chromosome-level grid-derived score as a reference point.
##
## NB depends on cache/anchor_grid.rds -- run 03_anchor_loci.R first.
## Rscript module_C2/R/02_compare_denominators.R
## =====================================================================
source("module_C2/R/00_helpers.R")
AG <- file.path(C2$CACHE, "anchor_grid.rds")
if (!file.exists(AG)) stop("run module_C2/R/03_anchor_loci.R first (builds ", AG, ")")
ag <- readRDS(AG); B <- ag$B
g  <- ag$G$rec50                                   # recommended matching rule
NCELL <- length(C2$TAUS) * length(C2$LMINS)

## ---- usable cells: cells yielding >= 1 significant region genome-wide ----
sg <- readRDS(file.path(C2$CACHE, "grid_scored.rds"))
cellsF <- sg$resF                                  # full-B per-cell counts
usable <- cellsF[n_sig > 0, paste(tau, lmin)]
NUSE <- length(usable); AVAIL <- NUSE / NCELL
c2_msg("[1] |G| = %d ; |U| = %d ; availability A = %.3f (B = %d)\n", NCELL, NUSE, AVAIL, B)

g[, cell := paste(tau, lmin)][, in_U := cell %in% usable]

## ---- the variants ----------------------------------------------------
V <- g[, .(
  n_sig_all   = sum(sig),
  n_sig_U     = sum(sig & in_U),
  n_det_all   = sum(detected),
  S_G         = sum(sig) / NCELL,                      # B: full-grid significance
  S_U         = sum(sig & in_U) / max(NUSE, 1L),       # A: current conditional
  D           = sum(detected) / NCELL,                 # D: detection stability
  Q           = sum(sig) / NCELL                       # D: significance stability
), by = label]
V[, Q_over_D := ifelse(D > 0, Q / D, NA_real_)]
V[, A := AVAIL]
setorder(V, -S_G)
fwrite(V, file.path(C2$RES, "C2_variant_comparison.csv"))
c2_msg("\n[2] denominator variants on the %d fixed anchor loci:\n", nrow(V))
print(V[, .(label, S_U = round(S_U, 3), S_G = round(S_G, 3), D = round(D, 3),
            Q = round(Q, 3), Q_over_D = round(Q_over_D, 3))])

## ---- rank agreement --------------------------------------------------
vec <- function(x) stats::setNames(V[[x]], V$label)
c2_msg("\n[3] rank agreement between variants (Spearman; top-5 overlap):\n")
pairs <- list(c("S_U", "S_G"), c("S_G", "D"), c("S_G", "Q_over_D"), c("D", "Q_over_D"),
              c("S_U", "Q_over_D"))
agr <- rbindlist(lapply(pairs, function(p) {
  a <- c2_agree(vec(p[1]), vec(p[2]), k = 5L)
  data.table(a = p[1], b = p[2], spearman = round(a$spearman, 3),
             top5_overlap = round(a$top_k_overlap, 2))
}))
print(agr)

## ---- S^U vs S^G: where do they differ, and why ----------------------
V[, delta := S_U - S_G]
c2_msg("\n[4] S^U inflates every locus by the factor |G|/|U| = %.2f (identical rankings\n", NCELL / NUSE)
c2_msg("    when every significant cell is usable, which it is by definition:\n")
c2_msg("    S^U == S_G * |G|/|U| exactly for all loci: %s\n",
       all(abs(V$S_U - V$S_G * NCELL / NUSE) < 1e-12))
c2_msg("    => the conditional denominator is a pure RESCALING here: it changes the\n")
c2_msg("       reported magnitude (and so its interpretability) but not the ranking.\n")

## ---- the single-usable-cell pathology, made concrete ----------------
c2_msg("\n[5] pathology check: a locus significant in exactly 1 usable cell scores\n")
c2_msg("    S^U = 1/%d = %.4f now, but would score S^U = 1.0 if |U| = 1.\n", NUSE, 1 / NUSE)
c2_msg("    S^G for that same locus is 1/%d = %.4f regardless of how dead the grid is.\n",
       NCELL, 1 / NCELL)

## ---- reference: the prototype's grid-derived chromosome score -------
## Done for BOTH matching rules, because they disagree sharply: the permissive
## rule lets an anchor claim the whole merged chromosome-scale component at low
## tau_C, which silently reproduces the prototype's per-chromosome statistic.
rep_f <- file.path(C2$RES, "current_C2_reproduction.csv")
if (file.exists(rep_f)) {
  pr <- fread(rep_f)
  out <- list()
  for (nm in names(ag$G)) {
    gg <- ag$G[[nm]]
    ss <- gg[, .(S_G = sum(sig) / NCELL), by = label]
    ss[, chr := as.integer(sub(".*Chr([0-9]+):.*", "\\1", label))]
    m <- merge(ss, pr[, .(chr, proto_stability = stability)], by = "chr", all.x = TRUE)
    rho <- stats::cor(m$S_G, m$proto_stability, method = "spearman", use = "complete.obs")
    out[[nm]] <- data.table(rule = nm, spearman_vs_prototype = round(rho, 3),
                            mean_S_G = round(mean(ss$S_G), 4))
    c2_msg("[6] rule %-6s : Spearman(anchor S_G, prototype per-chromosome score) = %+.3f\n", nm, rho)
  }
  cmp2 <- rbindlist(out)
  fwrite(cmp2, file.path(C2$RES, "anchor_vs_prototype.csv"))
  c2_msg("    (%d distinct prototype values are shared among the %d anchors -- 7 of the 17\n",
         uniqueN(pr$stability[pr$chr %in% unique(as.integer(sub(".*Chr([0-9]+):.*", "\\1", V$label)))]), nrow(V))
  c2_msg("     anchors sit on a chromosome that carries another anchor, so the prototype\n")
  c2_msg("     cannot separate them at all.)\n")
}

## ---- how much does the matching rule move the anchor ranking? -------
c2_msg("\n[7] anchor S_G ranking by matching rule:\n")
RS <- lapply(ag$G, function(gg) { s <- gg[, .(S_G = sum(sig) / NCELL), by = label]
                                  stats::setNames(s$S_G, s$label) })
for (a in names(RS)) for (b in names(RS)) if (a < b) {
  x <- c2_agree(RS[[a]], RS[[b]], k = 5L)
  c2_msg("    %-6s vs %-6s : rho = %+.3f ; top5 overlap = %.2f\n", a, b, x$spearman, x$top_k_overlap)
}
c2_msg("\n[7] top-5 anchors under each rule:\n")
for (a in names(RS)) c2_msg("    %-6s : %s\n", a,
  paste(names(sort(RS[[a]], decreasing = TRUE))[1:5], collapse = " | "))
c2_msg("\n[7] wrote results/C2_variant_comparison.csv, anchor_vs_prototype.csv\n")
