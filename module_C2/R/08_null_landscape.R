## =====================================================================
## module_C2 / R/08_null_landscape.R    [Question 1]
##
## Characterise NULL behaviour over the complete 25 x 7 = 175 cell grid, using
## all B = 200 population-permutation surrogates and the same LD graph / region
## rules as the empirical data.
##
## The grid is tau_C = seq(0.02, 0.50, 0.02) x l_min = {1,2,3,5,10,15,20}.
## tau_C = 0.05 is a REFERENCE coordinate only (appended to the cache for anchor
## construction); it is reported separately and never enters a grid denominator.
##
## Empirical counts are computed for DIAGNOSIS and are never used to define
## admissibility -- that is done from null quantities alone in 09_.
##
## Rscript module_C2/R/08_null_landscape.R
## =====================================================================
source("module_C2/R/00_helpers.R")
core <- readRDS(file.path(C2$CACHE, "grid_core.rds")); B <- core$B
NCELL <- length(C2$TAUS) * length(C2$LMINS)
UNI <- length(core$D$universe)

## ---- guard: the grid is exactly what we say it is -------------------
## compare on the cache KEYS: seq() values do not round-trip through
## as.numeric(as.character()), so a numeric %in% test spuriously fails.
cached <- names(core$by_tau)
gridk  <- as.character(C2$TAUS)
stopifnot(all(gridk %in% cached))
extra  <- as.numeric(setdiff(cached, gridk))
c2_msg("[0] grid = %d tau x %d l_min = %d cells ; B = %d ; universe = %d markers\n",
       length(C2$TAUS), length(C2$LMINS), NCELL, B, UNI)
c2_msg("[0] cached-but-OFF-GRID tau values (reference only, excluded from all grid\n")
c2_msg("    denominators): %s\n", paste(sprintf("%.2f", extra), collapse = ", "))

## ---- 1. per-cell null diagnostics -----------------------------------
NL <- rbindlist(lapply(C2$TAUS, function(t)
        rbindlist(lapply(C2$LMINS, function(l) c2_null_cell(core, t, l, B)))))
stopifnot(nrow(NL) == NCELL)

## empirical burden, for diagnosis only
OB <- rbindlist(lapply(C2$TAUS, function(t) rbindlist(lapply(C2$LMINS, function(l) {
  ct <- core$by_tau[[as.character(t)]]; keep <- which(ct$obs$size >= l)
  data.table(tau = t, lmin = l, n_obs_regions = length(keep),
             n_obs_markers = if (length(keep)) sum(ct$obs$size[keep]) else 0L) }))))
NL <- merge(NL, OB, by = c("tau", "lmin"))
NL[, obs_coverage := n_obs_markers / UNI]
setorder(NL, tau, lmin)
fwrite(NL, file.path(C2$RES, "null_grid_diagnostics.csv"))

## the reference point, reported apart
REF <- rbindlist(lapply(C2$LMINS, function(l) c2_null_cell(core, C2$REF_TAU, l, B)))
REFo <- rbindlist(lapply(C2$LMINS, function(l) { ct <- core$by_tau[[as.character(C2$REF_TAU)]]
  keep <- which(ct$obs$size >= l)
  data.table(tau = C2$REF_TAU, lmin = l, n_obs_regions = length(keep)) }))
REF <- merge(REF, REFo, by = c("tau", "lmin"))
fwrite(REF, file.path(C2$RES, "null_reference_point.csv"))

## ---- 2. is there a distinct lenient corner / transition? ------------
c2_msg("\n[1] p_null_any over the grid: min %.3f, median %.3f, max %.3f\n",
       min(NL$p_null_any), stats::median(NL$p_null_any), max(NL$p_null_any))
c2_msg("    cells with p_null_any == 0 : %d of %d (%.0f%%)\n",
       NL[p_null_any == 0, .N], NCELL, 100 * NL[p_null_any == 0, .N] / NCELL)
c2_msg("    cells with p_null_any == 1 : %d\n", NL[p_null_any == 1, .N])
c2_msg("\n[1] p_null_any by l_min (across tau_C):\n")
print(NL[, .(min = round(min(p_null_any), 3), median = round(stats::median(p_null_any), 3),
             max = round(max(p_null_any), 3),
             n_zero = sum(p_null_any == 0)), by = lmin][order(lmin)])
c2_msg("\n[1] p_null_any at the most lenient l_min = 1, across tau_C:\n")
print(NL[lmin == 1, .(tau, p_null_any = round(p_null_any, 3),
                      mean_null_regions = round(mean_null_regions, 2),
                      mean_null_coverage = round(mean_null_coverage, 4))])

## where does p_null_any cross the candidate tolerances, per l_min?
c2_msg("\n[2] transition: largest tau_C at which p_null_any still exceeds each tolerance\n")
c2_msg("    (i.e. the lenient region is tau_C <= this value); NA = never exceeds it\n")
tr <- rbindlist(lapply(c(0.05, 0.10, 0.20), function(eps)
  NL[, .(eps = eps, last_tau_above = { w <- tau[p_null_any > eps]
         if (length(w)) max(w) else NA_real_ },
         n_cells_above = sum(p_null_any > eps)), by = lmin]))
print(dcast(tr, lmin ~ eps, value.var = c("last_tau_above", "n_cells_above")))

## is the exclusion contiguous in tau at each l_min?
c2_msg("\n[3] contiguity of {p_null_any > eps} in tau_C, per l_min:\n")
for (eps in c(0.05, 0.10, 0.20)) {
  bad <- NL[p_null_any > eps]
  ok <- all(vapply(split(bad, bad$lmin), function(d) {
    i <- match(sort(d$tau), C2$TAUS); all(diff(i) == 1L) && min(i) == 1L }, logical(1)))
  c2_msg("    eps = %.2f : excluded set is a contiguous low-tau prefix at every l_min: %s\n", eps, ok)
}

## ---- 3. null vs empirical burden ------------------------------------
c2_msg("\n[4] Spearman(p_null_any, n_obs_regions) = %.3f  [DIAGNOSTIC ONLY]\n",
       stats::cor(NL$p_null_any, NL$n_obs_regions, method = "spearman"))
c2_msg("    Spearman(mean_null_coverage, obs_coverage) = %.3f\n",
       stats::cor(NL$mean_null_coverage, NL$obs_coverage, method = "spearman"))
c2_msg("    -> null and empirical burden co-vary, but admissibility below uses the\n")
c2_msg("       null columns ONLY; the empirical columns are never consulted.\n")

## ---- 4. figures -----------------------------------------------------
suppressMessages(library(ggplot2))
th <- theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 7),
        plot.subtitle = element_text(size = 7, colour = "grey30"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5))
hm <- function(col, title, sub, lab, digits = 2) {
  ggplot(NL, aes(factor(tau), factor(lmin, levels = C2$LMINS), fill = .data[[col]])) +
    geom_tile(colour = "white", linewidth = 0.2) +
    geom_text(aes(label = ifelse(.data[[col]] > 0, formatC(.data[[col]], format = "f", digits = digits), "")),
              size = 1.7, colour = "grey15") +
    scale_fill_gradientn(colours = C2$ZISSOU, name = lab) +
    labs(x = expression(tau[C]), y = expression(l[min]), title = title, subtitle = sub) + th
}
g1 <- hm("p_null_any", "Null-region production: Pr(a surrogate yields >= 1 region)",
         sprintf("population-permutation null, B = %d surrogates; this is the ONLY quantity used to define admissibility", B),
         expression(p[null]), 2)
g2 <- hm("mean_null_regions", "Mean number of null regions per surrogate",
         "same grid, same LD graph and region rules as the observed data", "mean", 1)
g3 <- hm("mean_null_coverage", "Mean fraction of the marker universe covered by null regions",
         sprintf("universe = %d markers", UNI), "coverage", 3)
g4 <- hm("n_obs_regions", "Observed region count (DIAGNOSTIC -- not used to define admissibility)",
         "shown on the identical grid for comparison with the null panels above", "regions", 0)
ggsave(file.path(C2$FIG, "fig7_null_p_any.png"),     g1, width = 11, height = 3.6, dpi = 170)
ggsave(file.path(C2$FIG, "fig8_null_mean_regions.png"), g2, width = 11, height = 3.6, dpi = 170)
ggsave(file.path(C2$FIG, "fig9_null_coverage.png"),  g3, width = 11, height = 3.6, dpi = 170)
ggsave(file.path(C2$FIG, "fig10_obs_regions.png"),   g4, width = 11, height = 3.6, dpi = 170)

## null vs empirical burden
g5 <- ggplot(NL, aes(p_null_any, n_obs_regions + 1, colour = factor(lmin, levels = C2$LMINS))) +
  geom_point(size = 1.3) + scale_y_log10() +
  scale_colour_viridis_d(name = expression(l[min]), option = "C", end = 0.9) +
  geom_vline(xintercept = c(0.05, 0.10, 0.20), linetype = c(3, 2, 4), colour = "grey40") +
  labs(x = expression(Pr[0](N[regions] > 0)), y = "observed regions + 1 (log scale)",
       title = "Null-region probability vs empirical region burden, per grid cell",
       subtitle = "dashed/dotted lines = candidate tolerances eps = 0.05, 0.10, 0.20. The association is corroborating evidence only: admissibility is defined from the x axis alone.") +
  th + theme(axis.text.x = element_text(angle = 0))
ggsave(file.path(C2$FIG, "fig11_null_vs_obs.png"), g5, width = 7.5, height = 4.5, dpi = 170)

c2_msg("\n[5] wrote results/null_grid_diagnostics.csv, null_reference_point.csv\n")
c2_msg("[5] wrote figures fig7-fig11\n")
