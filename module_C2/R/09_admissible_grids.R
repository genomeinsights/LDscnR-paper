## =====================================================================
## module_C2 / R/09_admissible_grids.R    [Question 2]
##
## Candidate NULL-ONLY admissibility rules over the 175-cell grid:
##   P(eps)   : p_null_any <= eps                      (point estimate)
##   U(eps)   : upper Clopper-Pearson bound <= eps     (honest about B = 200)
##   C(c)     : mean null marker coverage <= c
##   M(m)     : mean null region count <= m
## No empirical quantity is consulted anywhere in this script.
##
## Rscript module_C2/R/09_admissible_grids.R
## =====================================================================
source("module_C2/R/00_helpers.R")
NL <- fread(file.path(C2$RES, "null_grid_diagnostics.csv"))
NCELL <- length(C2$TAUS) * length(C2$LMINS)
stopifnot(nrow(NL) == NCELL)
## snap tau back onto the exact seq() doubles: a CSV round-trip yields 0.06, while
## seq(0.02, 0.50, 0.02) holds 0.06000000000000000005, so match()/== silently fail.
NL[, tau := C2$TAUS[match(round(tau, 6), round(C2$TAUS, 6))]]
stopifnot(!anyNA(NL$tau))

RULES <- list(
  P05 = quote(p_null_any  <= 0.05), P10 = quote(p_null_any <= 0.10), P20 = quote(p_null_any <= 0.20),
  U05 = quote(p_null_hi   <= 0.05), U10 = quote(p_null_hi  <= 0.10), U20 = quote(p_null_hi  <= 0.20),
  C001 = quote(mean_null_coverage <= 0.001), C0005 = quote(mean_null_coverage <= 0.0005),
  M1   = quote(mean_null_regions  <= 1.0),   M05   = quote(mean_null_regions  <= 0.5))

## ---- shape / connectivity of a retained set --------------------------
shape <- function(keep) {
  d <- NL[keep]
  if (!nrow(d)) return(list(n = 0L, n_comp = 0L, lmins = "(none)",
                            contig_suffix = NA, rect = NA))
  ## 4-neighbour connectivity on the (tau index, l_min index) lattice
  ti <- match(round(d$tau, 6), round(C2$TAUS, 6)); li <- match(d$lmin, C2$LMINS)
  stopifnot(!anyNA(ti), !anyNA(li))
  idx <- paste(ti, li); comp <- seq_len(nrow(d))
  for (i in seq_len(nrow(d))) for (j in seq_len(nrow(d))) if (i < j) {
    if (abs(ti[i] - ti[j]) + abs(li[i] - li[j]) == 1L) {
      a <- comp[i]; b <- comp[j]; if (a != b) comp[comp == b] <- a } }
  ## is it a contiguous HIGH-tau suffix within each retained l_min?
  contig <- all(vapply(split(ti, li), function(v) { v <- sort(v)
    all(diff(v) == 1L) && max(v) == length(C2$TAUS) }, logical(1)))
  list(n = nrow(d), n_comp = uniqueN(comp),
       lmins = paste(sort(unique(d$lmin)), collapse = ","),
       contig_suffix = contig,
       rect = nrow(d) == uniqueN(d$tau) * uniqueN(d$lmin))
}

out <- list(); KEEP <- list()
for (nm in names(RULES)) {
  keep <- NL[, eval(RULES[[nm]])]
  KEEP[[nm]] <- NL[keep, .(tau, lmin)]
  s <- shape(keep)
  out[[nm]] <- data.table(rule = nm, expr = deparse(RULES[[nm]]),
                          n_adm = s$n, frac_adm = round(s$n / NCELL, 3),
                          n_components = s$n_comp, lmins_retained = s$lmins,
                          contiguous_high_tau = s$contig_suffix, rectangular = s$rect,
                          max_p_null_retained = if (s$n) round(max(NL[keep]$p_null_any), 3) else NA_real_)
}
AD <- rbindlist(out)
fwrite(AD, file.path(C2$RES, "admissible_grid_rules.csv"))
saveRDS(KEEP, file.path(C2$CACHE, "admissible_cells.rds"))
c2_msg("[1] candidate null-only admissibility rules (%d-cell grid):\n", NCELL)
print(AD)

## ---- is exclusion concentrated in the lenient corner? ---------------
c2_msg("\n[2] which l_min values survive each tolerance (tau_C plays almost no role):\n")
for (nm in c("P05", "P10", "P20")) {
  d <- KEEP[[nm]]
  if (!nrow(d)) { c2_msg("    %-4s : EMPTY grid\n", nm); next }
  z <- d[, .(n_tau = .N, tau_min = min(tau), tau_max = max(tau)), by = lmin][order(lmin)]
  c2_msg("    %-4s : %s\n", nm, paste(sprintf("l_min=%d (%d tau, %.2f-%.2f)",
         z$lmin, z$n_tau, z$tau_min, z$tau_max), collapse = " ; "))
}
c2_msg("\n[2] Spearman(p_null_any, l_min) = %.3f  vs  Spearman(p_null_any, tau_C) = %.3f\n",
       stats::cor(NL$p_null_any, NL$lmin, method = "spearman"),
       stats::cor(NL$p_null_any, NL$tau, method = "spearman"))
c2_msg("    -> null-region production is governed by l_min, not by tau_C. There is no\n")
c2_msg("       'lenient tau_C corner': p_null_any at l_min=1 falls only from %.2f to %.2f\n",
       NL[lmin == 1 & tau == min(C2$TAUS)]$p_null_any, NL[lmin == 1 & tau == max(C2$TAUS)]$p_null_any)
c2_msg("       across the entire tau_C range, and never reaches any candidate tolerance.\n")

## ---- is there a natural transition to cut at? -----------------------
ps <- sort(NL$p_null_any)
gaps <- diff(ps)
c2_msg("\n[3] natural-break check on the sorted p_null_any values:\n")
c2_msg("    range %.3f-%.3f ; largest gap = %.3f between %.3f and %.3f\n",
       min(ps), max(ps), max(gaps), ps[which.max(gaps)], ps[which.max(gaps) + 1])
c2_msg("    second largest gap = %.3f ; median gap = %.4f\n",
       sort(gaps, decreasing = TRUE)[2], stats::median(gaps))
c2_msg("    -> %s\n", if (max(gaps) > 5 * stats::quantile(gaps, 0.9))
       "a distinct break exists" else
       "NO distinct break: p_null_any is a smooth continuum, so ANY tolerance is arbitrary")

## ---- Monte-Carlo uncertainty at the boundary ------------------------
c2_msg("\n[4] finite-B uncertainty (B = 200; Clopper-Pearson):\n")
for (eps in c(0.05, 0.10, 0.20)) {
  pt <- NL[p_null_any <= eps, .N]; up <- NL[p_null_hi <= eps, .N]
  amb <- NL[p_null_any <= eps & p_null_hi > eps, .N]
  c2_msg("    eps = %.2f : point rule keeps %3d cells, upper-bound rule keeps %3d;\n", eps, pt, up)
  c2_msg("               %3d cells (%.0f%% of those kept) are statistically indistinguishable\n",
         amb, if (pt) 100 * amb / pt else 0)
  c2_msg("               from failing the rule at B = %d.\n", 200)
}
c2_msg("    typical CI width at the boundary: %.3f\n",
       NL[abs(p_null_any - 0.10) < 0.03, mean(p_null_hi - p_null_lo)])
c2_msg("    -> at B = 200 a cell with k = %d/200 has CI [%.3f, %.3f]; distinguishing\n", 20,
       stats::binom.test(20, 200)$conf.int[1], stats::binom.test(20, 200)$conf.int[2])
c2_msg("       eps = 0.05 from eps = 0.10 is at the edge of this resolution.\n")

## ---- figure: the admissible grids ------------------------------------
suppressMessages(library(ggplot2))
PL <- rbindlist(lapply(c("P05", "P10", "P20", "U10"), function(nm) {
  d <- copy(NL); d[, rule := nm]
  d[, admissible := paste(tau, lmin) %in% KEEP[[nm]][, paste(tau, lmin)]]; d }))
g <- ggplot(PL, aes(factor(tau), factor(lmin, levels = C2$LMINS), fill = p_null_any)) +
  geom_tile(aes(colour = admissible), linewidth = 0.4) +
  facet_wrap(~ rule, ncol = 1) +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "white"),
                      name = "admissible", labels = c("excluded", "retained")) +
  scale_fill_gradientn(colours = C2$ZISSOU, name = expression(p[null])) +
  labs(x = expression(tau[C]), y = expression(l[min]),
       title = "Null-only admissible grids (black outline = retained)",
       subtitle = "Retention is governed almost entirely by l_min; tau_C barely matters. P = point estimate, U = upper Clopper-Pearson bound.") +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 8),
        plot.subtitle = element_text(size = 7, colour = "grey30"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5))
ggsave(file.path(C2$FIG, "fig12_admissible_grids.png"), g, width = 10, height = 8, dpi = 170)
c2_msg("\n[5] wrote results/admissible_grid_rules.csv + figures/fig12_admissible_grids.png\n")
