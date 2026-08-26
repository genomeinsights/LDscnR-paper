## =====================================================================
## module_sim_LDscnR / pr_curves.R
##
## Precision-recall curves over tau_C, one per l_min, scored against truth with
## evaluate_ors(). This is the benchmark framing rather than the discovery one:
## tau_C is SWEPT as the operating knob, so no significance filter is applied and
## NO SURROGATES ARE NEEDED. C_obs is a full-length vector over every marker on
## every chromosome, and the saved edge graph spans every marker with C > 0, so
## any tau can be clustered from the saved fits at negligible cost.
##
## C_obs does not depend on the surrogate basis (verified: identical between the
## env_orth and genetic fits), so there is ONE curve set per cell per ENGINE.
##
## The tau grid is ADAPTIVE -- the distinct observed C values, thinned -- because
## a fixed [0.05, 1] grid is meaningless when a cell's C-score maxes out low, and
## the distinct values are exactly where the region set changes.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/pr_curves.R <scan_dir> [outdir]
## Env: PANEL_DIR, ENGINE (emmax|lfmm, default emmax), MAX_TAU (default 60),
##      LMINS (default 1,2,3,5,10,20)
## =====================================================================
suppressMessages({ library(data.table); library(LDscnR); library(ggplot2) })
`%||%` <- function(a, b) if (is.null(a)) b else a

a <- commandArgs(trailingOnly = TRUE)
if (!length(a)) stop("usage: pr_curves.R <scan_dir> [outdir]")
SCAN_DIR <- a[1]; OUT <- if (length(a) >= 2) a[2] else SCAN_DIR
PANEL_DIR <- Sys.getenv("PANEL_DIR", "/Volumes/Nemo/Nemo_sim/analysis_inputs")
ENGINE <- Sys.getenv("ENGINE", "emmax")
MAX_TAU <- as.integer(Sys.getenv("MAX_TAU", "60"))
LMINS <- as.integer(strsplit(Sys.getenv("LMINS", "1,2,3,5,10,20"), ",")[[1]])
RHO_R2 <- 0.75; RHO_D <- 0.95; DCAP <- 5e5
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

## one fit per cell is enough -- prefer env_orth, fall back to whatever exists
scans <- list.files(SCAN_DIR, pattern = sprintf("^scan_.*_%s_.*[.]rds$", ENGINE), full.names = TRUE)
cellof <- function(f) sub(sprintf("^scan_(V[0-9.]+_c[0-9.]+_env[0-9]+)_%s_.*$", ENGINE), "\\1", basename(f))
cells <- unique(vapply(scans, cellof, character(1)))
cat(sprintf("[1] engine %s | %d cell(s) | l_min: %s\n", ENGINE, length(cells), paste(LMINS, collapse = ","))); flush.console()

res <- list(); k <- 0L
for (cell in cells) {
  pf <- file.path(PANEL_DIR, sprintf("panel_%s.rds", cell))
  if (!file.exists(pf)) next
  panel <- readRDS(pf); map <- as.data.table(panel$map)
  if (!"true_pos_QTN" %in% names(map)) map <- flag_true_qtns(map)
  n_true <- sum(map$true_pos_QTN %in% TRUE)
  th <- score_thresholds(as.data.table(panel$decay_sum), rho_r2 = RHO_R2, rho_d = RHO_D, dmax_cap = DCAP)

  sc <- scans[vapply(scans, cellof, character(1)) == cell]
  f <- sc[grepl("env_orth", sc)][1]; if (is.na(f)) f <- sc[1]
  x <- readRDS(f); C <- x$null$C_obs
  pos <- names(C)[which(C > 0)]
  if (!length(pos) || !n_true) { cat(sprintf("  [skip] %s: %d markers with C>0, %d QTN\n", cell, length(pos), n_true)); next }
  qtab <- qtn_ld_table(panel$GTs, map, pos, 2e6, cores = 1)

  ## adaptive grid: the distinct positive C values, thinned
  tv <- sort(unique(C[C > 0]))
  if (length(tv) > MAX_TAU) tv <- as.numeric(stats::quantile(tv, seq(0, 1, length.out = MAX_TAU), type = 1))
  tv <- unique(tv)

  for (tau in tv) {
    mk <- names(C)[which(C >= tau)]
    if (!length(mk)) next
    r_all <- ld_regions(mk, x$edges)
    for (L in LMINS) {
      r <- r_all[lengths(r_all) >= L]
      ev <- if (length(r)) evaluate_ors(r, map, qtab, th$r2min, th$dmax)
            else list(TP = 0L, FP = 0L, Precision = NA_real_, Recall = 0)
      k <- k + 1L
      res[[k]] <- data.table(cell, engine = ENGINE, tau = tau, l_min = L,
                             n_regions = length(r), TP = ev$TP, FP = ev$FP,
                             precision = ev$Precision, recall = ev$Recall, n_true = n_true)
    }
  }
  cat(sprintf("  %s: %d QTN, %d tau values\n", cell, n_true, length(tv))); flush.console()
}
pr <- rbindlist(res)
fwrite(pr, file.path(OUT, sprintf("pr_curves_%s.csv", ENGINE)))

## PR-AUC per cell per l_min, then replicate-averaged
auc <- pr[!is.na(precision), .(PR_AUC = pr_auc(recall, precision)), by = c("cell", "l_min")]
agg <- auc[, .(cells = .N, PR_AUC = round(mean(PR_AUC, na.rm = TRUE), 3),
               SE = round(stats::sd(PR_AUC, na.rm = TRUE) / sqrt(.N), 3)), by = l_min]
setorderv(agg, "l_min")
cat(sprintf("\n=== PR-AUC by l_min, replicate-averaged over %d cells (%s) ===\n", uniqueN(auc$cell), ENGINE))
print(agg)

## env-averaged curve: interpolate precision onto a common recall grid
rg <- seq(0, 1, by = 0.02)
curve <- pr[!is.na(precision), {
  o <- order(recall); rr <- recall[o]; pp <- precision[o]
  .(recall = rg, precision = stats::approx(rr, pp, xout = rg, rule = 2, ties = "ordered")$y)
}, by = c("cell", "l_min")][, .(precision = mean(precision, na.rm = TRUE),
                                 se = stats::sd(precision, na.rm = TRUE) / sqrt(.N)), by = c("l_min", "recall")]
fwrite(curve, file.path(OUT, sprintf("pr_curve_mean_%s.csv", ENGINE)))

g <- ggplot(curve, aes(recall, precision, colour = factor(l_min))) +
  geom_ribbon(aes(ymin = precision - se, ymax = precision + se, fill = factor(l_min)),
              alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.8) +
  scale_colour_viridis_d(name = expression(l[min]), end = 0.9) +
  scale_fill_viridis_d(guide = "none", end = 0.9) +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
  labs(x = "Recall", y = "Precision",
       title = sprintf("Precision-recall over tau_C by region-size filter (%s)", ENGINE),
       subtitle = sprintf("evaluate_ors against true_pos_QTN; mean +- SE over %d simulated genomes", uniqueN(auc$cell))) +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(), strip.background = element_blank())
ggsave(file.path(OUT, sprintf("pr_curves_%s.png", ENGINE)), g, width = 7, height = 5, dpi = 170)
cat(sprintf("\n[2] wrote pr_curves_%s.csv, pr_curve_mean_%s.csv, pr_curves_%s.png\n", ENGINE, ENGINE, ENGINE))
