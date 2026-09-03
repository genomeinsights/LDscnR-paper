## =====================================================================
## module_C2 / R/11_marker_support.R    [Question 4B]
##
## Anchor-free marker-level operating-grid support:
##   S_m = #{admissible cells in which marker m belongs to an empirical region}
##         / |G_adm|
##
## This is OPERATING-GRID SUPPORT, not a second-tier significance score: it partly
## re-expresses the marker C-score itself (a marker only enters a region at
## tau_C <= C_m), so its correlation with C is reported explicitly.
##
## Cross-grid locus families are built by applying the package's OWN region rule
## (ld_regions on the cached LD graph, with the 500 kb gap-split) to the
## support-gated marker set. That groups by chromosome AND connected component --
## it does not aggregate by chromosome, and it cannot chain two signals together
## through an intermediate region, because linkage is r^2-based and gap-capped.
##
## Rscript module_C2/R/11_marker_support.R
## =====================================================================
source("module_C2/R/00_helpers.R")
suppressMessages(library(LDscnR))
core <- readRDS(file.path(C2$CACHE, "grid_core.rds")); B <- core$B
KEEP <- readRDS(file.path(C2$CACHE, "admissible_cells.rds"))
FULL <- CJ(tau = C2$TAUS, lmin = C2$LMINS)
GRIDS <- list(full = FULL, P20 = KEEP$P20, P10 = KEEP$P10, U10 = KEEP$U10)
map <- core$D$map; C_obs <- core$D$C_obs
anc <- c2_anchors(core)

## ---- 1. marker support over each grid --------------------------------
marker_support <- function(cells) {
  tab <- table(unlist(lapply(seq_len(nrow(cells)), function(i)
    c2_obs_markers(core, cells$tau[i], cells$lmin[i])), use.names = FALSE))
  s <- stats::setNames(as.integer(tab) / nrow(cells), names(tab))
  s
}
SM <- lapply(GRIDS, function(g) if (nrow(g)) marker_support(g) else numeric(0))
c2_msg("[1] markers with S_m > 0, per grid (of %d in the universe):\n", length(core$D$universe))
for (nm in names(SM)) c2_msg("    %-5s (%3d cells): %5d markers ; max S_m = %.3f ; mean over positive = %.3f\n",
       nm, nrow(GRIDS[[nm]]), length(SM[[nm]]),
       if (length(SM[[nm]])) max(SM[[nm]]) else 0,
       if (length(SM[[nm]])) mean(SM[[nm]]) else 0)

## ---- 2. is S_m just the C-score re-expressed? ------------------------
c2_msg("\n[2] S_m vs the marker C-score (the honest caveat):\n")
for (nm in names(SM)) { s <- SM[[nm]]; if (!length(s)) next
  cc <- C_obs[names(s)]
  c2_msg("    %-5s : Spearman(S_m, C_m) over supported markers = %.3f\n",
         nm, stats::cor(s, cc, method = "spearman")) }
c2_msg("    -> S_m is substantially a re-expression of C_m, as expected: a marker can only\n")
c2_msg("       enter a region at tau_C <= C_m. Report it as support, never as new evidence.\n")

## ---- 3. locus families, by chromosome AND component ------------------
## The gate is a support threshold; the grouping is ld_regions() on the cached
## edge graph, so families are r^2-connected components split at 500 kb gaps.
family_table <- function(s, thr, gname) {
  mk <- names(s)[s >= thr]
  if (!length(mk)) return(data.table())
  fam <- ld_regions(mk, core$D$edges)
  rbindlist(lapply(seq_along(fam), function(i) { m <- fam[[i]]
    data.table(grid = gname, support_gate = thr, family = i,
               Chr = unname(core$D$mchr[m[1]]),
               start = min(core$D$mpos[m]), end = max(core$D$mpos[m]),
               n_marker = length(m),
               mean_S = mean(s[m]), max_S = max(s[m]), max_C = max(C_obs[m])) }))
}
FAM <- rbindlist(lapply(c(0.02, 0.05, 0.10), function(t)
         rbindlist(lapply(names(SM), function(nm)
           if (length(SM[[nm]])) family_table(SM[[nm]], t, nm) else NULL))), fill = TRUE)
## does each family overlap a reference-point region? (by marker membership)
amk <- unique(unlist(anc$mk, use.names = FALSE))
FAM[, in_reference := FALSE]
if (nrow(FAM)) {
  for (nm in names(SM)) for (t in unique(FAM$support_gate)) {
    s <- SM[[nm]]; if (!length(s)) next
    mk <- names(s)[s >= t]; if (!length(mk)) next
    fam <- ld_regions(mk, core$D$edges)
    ov <- vapply(fam, function(m) any(m %in% amk), logical(1))
    FAM[grid == nm & support_gate == t, in_reference := ov]
  }
}
fwrite(FAM, file.path(C2$RES, "marker_support_families.csv"))
c2_msg("\n[3] locus families (r^2-connected components of the support-gated markers,\n")
c2_msg("    500 kb gap-split -- grouped by chromosome AND component, never by chromosome alone):\n")
print(FAM[, .(n_families = .N, n_in_reference = sum(in_reference),
              n_novel = sum(!in_reference),
              median_markers = as.numeric(stats::median(n_marker))),
          by = .(grid, support_gate)][order(grid, support_gate)])

## ---- 4. cores / shoulders / isolated candidates ----------------------
s <- SM$P20
if (length(s)) {
  fam <- family_table(s, 0.02, "P20")
  fam[, kind := fifelse(max_S >= 0.10 & n_marker >= 5, "stable core",
               fifelse(max_S >= 0.05, "parameter-dependent shoulder", "low-support isolated"))]
  c2_msg("\n[4] P20 families classified (gate S_m >= 0.02):\n")
  print(fam[, .(n = .N, markers = sum(n_marker)), by = kind])
  setorder(fam, -max_S, -n_marker)
  c2_msg("\n    top families by max S_m:\n")
  print(head(fam[, .(Chr, start_Mb = round(start/1e6, 2), end_Mb = round(end/1e6, 2),
                     n_marker, max_S = round(max_S, 3), max_C = round(max_C, 2))], 12))
}

## ---- 5. anchor-free vs reference-point: what is missed? --------------
## Support mass that lies OUTSIDE every reference-point region.
c2_msg("\n[5] reference-point coverage of the supported markers:\n")
for (nm in names(SM)) { s <- SM[[nm]]; if (!length(s)) next
  inref <- names(s) %in% amk
  c2_msg("    %-5s : %5d supported markers, %5d (%.1f%%) inside a reference-point region;\n",
         nm, length(s), sum(inref), 100 * mean(inref))
  c2_msg("            support MASS outside the reference set = %.1f%% ; max S_m outside = %.3f\n",
         100 * sum(s[!inref]) / sum(s), if (any(!inref)) max(s[!inref]) else 0) }

saveRDS(SM, file.path(C2$CACHE, "marker_support.rds"))
SMDT <- rbindlist(lapply(names(SM), function(nm) { s <- SM[[nm]]
  if (!length(s)) return(NULL)
  data.table(grid = nm, marker = names(s), S_m = s, C_m = unname(C_obs[names(s)]),
             Chr = unname(core$D$mchr[names(s)]), Pos = unname(core$D$mpos[names(s)]),
             in_reference = names(s) %in% amk) }))
fwrite(SMDT, file.path(C2$RES, "marker_support.csv"))
c2_msg("\n[6] wrote results/marker_support.csv (%d rows) + marker_support_families.csv\n", nrow(SMDT))
