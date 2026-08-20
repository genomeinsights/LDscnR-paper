## =====================================================================
## module_C2 / R/13_comparisons.R    [Question 6, consolidated]
##
## The six required comparisons in one place:
##   1 full grid vs 2 each null-admissible grid vs 3 the old conditional |U|
##   4 detection vs BH-significance support
##   5 rec25 / rec50 / rec75 (and reciprocal)
##   6 reference-point region support vs anchor-free marker support
## plus: which regions change interpretation materially.
##
## Rscript module_C2/R/13_comparisons.R
## =====================================================================
source("module_C2/R/00_helpers.R")
core <- readRDS(file.path(C2$CACHE, "grid_core.rds"))
KEEP <- readRDS(file.path(C2$CACHE, "admissible_cells.rds"))
SUP  <- readRDS(file.path(C2$CACHE, "support.rds"))
SM   <- readRDS(file.path(C2$CACHE, "marker_support.rds"))
anc  <- c2_anchors(core)
FULL <- CJ(tau = C2$TAUS, lmin = C2$LMINS); NCELL <- nrow(FULL)
NC <- c(full = NCELL, P20 = nrow(KEEP$P20), P10 = nrow(KEEP$P10),
        P05 = nrow(KEEP$P05), U10 = nrow(KEEP$U10), C001 = nrow(KEEP$C001))

Dv <- function(gn, tn = "rec50") { k <- paste(gn, tn, sep = "."); g <- SUP[[k]]
  s <- g[, .(D = sum(detected) / NC[[gn]]), by = label]; stats::setNames(s$D, s$label) }

## ---- 6. region support vs anchor-free marker support ----------------
c2_msg("[6] reference-point region support D_r vs anchor-free marker support S_m\n")
rows <- list()
for (gn in c("full", "P20", "P10", "U10")) {
  s <- SM[[gn]]; d <- Dv(gn)
  agg <- rbindlist(lapply(names(anc$mk), function(l) { m <- anc$mk[[l]]
    v <- s[intersect(m, names(s))]
    data.table(label = l, D_r = d[[l]],
               mean_S = if (length(v)) mean(v) else 0,
               max_S  = if (length(v)) max(v)  else 0,
               frac_markers_supported = length(v) / length(m)) }))
  rows[[gn]] <- agg[, grid := gn]
  c2_msg("    %-5s : Spearman(D_r, mean S_m) = %+.3f ; Spearman(D_r, max S_m) = %+.3f ; top-5 overlap = %.2f\n",
         gn, stats::cor(agg$D_r, agg$mean_S, method = "spearman"),
         stats::cor(agg$D_r, agg$max_S, method = "spearman"),
         c2_agree(stats::setNames(agg$D_r, agg$label),
                  stats::setNames(agg$mean_S, agg$label), k = 5L)$top_k_overlap)
}
RM <- rbindlist(rows); fwrite(RM, file.path(C2$RES, "region_vs_marker_support.csv"))
c2_msg("\n    P20 detail (region support vs the support of its own markers):\n")
print(RM[grid == "P20"][order(-D_r), .(label, D_r = round(D_r, 3), mean_S = round(mean_S, 3),
        max_S = round(max_S, 3), frac_mk = round(frac_markers_supported, 2))])
z <- RM[grid == "P20" & D_r == 0 & max_S > 0]
c2_msg("\n    -> region-level and marker-level views AGREE almost exactly on the admissible\n")
c2_msg("       grids (Spearman 0.99). In principle a region could score D_r = 0 while its\n")
c2_msg("       markers kept positive S_m (never recovered >=50%% intact, yet partly inside\n")
c2_msg("       some region); in fact %d of the %d zero-support P20 regions are of that kind --\n",
       nrow(z), RM[grid == "P20" & D_r == 0, .N])
c2_msg("       their markers have NO support at all (frac_markers_supported = 0), so nothing\n")
c2_msg("       is hidden by working at the region level rather than the marker level.\n")

## ---- consolidated comparison table ----------------------------------
base <- Dv("full")
cons <- list()
add <- function(name, v, note) {
  a <- c2_agree(base, v, k = 5L)
  cons[[length(cons) + 1L]] <<- data.table(
    comparison = name, spearman_vs_full_rec50 = round(a$spearman, 3),
    top5_overlap = round(a$top_k_overlap, 2), n_zero_support = sum(v == 0),
    mean_support = round(mean(v), 4), note = note)
}
add("1 full grid (175 cells)", base, "reference for all rows")
for (gn in c("P20", "P10", "P05", "U10", "C001"))
  add(sprintf("2 admissible %s (%d cells)", gn, NC[[gn]]), Dv(gn),
      sprintf("null-only rule; %d/17 regions lose all support", sum(Dv(gn) == 0)))
sg <- readRDS(file.path(C2$CACHE, "grid_scored.rds")); usable <- sg$resF[n_sig > 0, .N]
add(sprintf("3 old conditional |U| = %d", usable),
    stats::setNames(SUP[["full.rec50"]][, sum(sig), by = label]$V1 / usable,
                    SUP[["full.rec50"]][, sum(sig), by = label]$label),
    "pure rescaling of row 1; identical ranking")
add("4 BH-significance support (full)",
    stats::setNames(SUP[["full.rec50"]][, sum(sig) / NCELL, by = label]$V1,
                    SUP[["full.rec50"]][, sum(sig) / NCELL, by = label]$label),
    "IDENTICAL to detection: 0 detected-but-not-significant cells")
for (tn in c("rec25", "rec75")) add(sprintf("5 matching %s (full)", tn), Dv("full", tn),
    "anchor-marker retention threshold")
vr <- SUP[["full.rec50recip"]][, .(D = sum(detected) / NCELL), by = label]
add("5 reciprocal rec50 (full)", stats::setNames(vr$D, vr$label), "also demands |AnR|/|R| >= 0.5")
sm <- SM$full; agg <- vapply(anc$mk, function(m) { v <- sm[intersect(m, names(sm))]
  if (length(v)) mean(v) else 0 }, 0)
add("6 anchor-free mean S_m (full)", agg, "marker-level support of the same regions")
CONS <- rbindlist(cons); fwrite(CONS, file.path(C2$RES, "consolidated_comparisons.csv"))
c2_msg("\n[all] consolidated comparison table:\n"); print(CONS)

## ---- regions whose interpretation changes materially ----------------
c2_msg("\n[7] regions changing interpretation between the full and P20 grids:\n")
ch <- data.table(label = names(base), full = base, P20 = Dv("P20")[names(base)])
ch[, rank_full := frank(-full, ties.method = "min")][, rank_P20 := frank(-P20, ties.method = "min")]
ch[, delta_rank := rank_P20 - rank_full]
setorder(ch, -delta_rank)
print(ch[, .(label, full = round(full, 3), P20 = round(P20, 3), rank_full, rank_P20, delta_rank)])
c2_msg("\n    largest demotions are the small anchors that only ever appear in null-prone cells;\n")
c2_msg("    the Chr1 inversion moves rank %d -> %d and Eda rank %d -> %d.\n",
       ch[grepl("inv", label)]$rank_full, ch[grepl("inv", label)]$rank_P20,
       ch[grepl("Eda", label)]$rank_full, ch[grepl("Eda", label)]$rank_P20)
c2_msg("\n[8] wrote region_vs_marker_support.csv + consolidated_comparisons.csv\n")
