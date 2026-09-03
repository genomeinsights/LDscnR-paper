## =====================================================================
## module_sim_LDscnR / null_operating_grid.R
##
## The tau_C x l_min OPERATING grid, averaged over every environment analysed
## with permutations. Distinct from the C-score's rho x q* INTEGRATION grid
## (see rho_q_integration_grid.R): that one is the 400 analyses C averages over,
## this one is what happens downstream once C exists and you must choose where
## to cut.
##
## Three panels, and the third is the one that decides anything:
##   regions surviving q_R < 0.05   the null side
##   regions containing a QTN       the truth side
##   precision                      the trade between them
##
## READ l_min, NOT tau. Across a six-fold change in tau at fixed l_min the
## region count moves little and the QTN count barely at all; l_min does the
## work. That matches the earlier finding that this grid is really an l_min rule.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/null_operating_grid.R
## Env: SIM_INPUTS, NULL_DIR, OUT
## =====================================================================
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(LDscnR)})
A   <- Sys.getenv("SIM_INPUTS", "/Volumes/Nemo/Nemo_sim/analysis_inputs")
D   <- Sys.getenv("NULL_DIR",   "module_sim_LDscnR/results/nulls_V2_c1")
OUT <- Sys.getenv("OUT",        "module_sim_LDscnR/figures")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

fs <- list.files(D, pattern = "^scan_.*[.]rds$", full.names = TRUE)
if (!length(fs)) stop("no scans in ", D)
cat(sprintf("  %d scans\n", length(fs)))

one <- function(f) {
  x <- readRDS(f)
  e <- as.integer(sub(".*_env([0-9]+)_.*", "\\1", basename(f)))
  b <- sub(".*_emmax_(.+)_B[0-9]+\\.rds$", "\\1", basename(f))
  panel <- readRDS(file.path(A, sprintf("panel_V2_c1_env%d.rds", e)))
  map <- flag_true_qtns(as.data.table(panel$map))
  qtn <- map[true_pos_QTN %in% TRUE, .(Chr = as.character(Chr), Pos)]
  if (!nrow(qtn)) return(NULL)
  ls_ <- as.data.table(x$c2$landscape); C <- x$null$C_obs; ed <- x$edges
  tr <- rbindlist(lapply(seq_len(nrow(ls_)), function(i) {
    tau <- ls_$tau[i]; L <- ls_$l_min[i]
    mk <- names(C)[which(C >= tau)]
    if (!length(mk)) return(data.table(tau = tau, l_min = L, n = 0L, hit = 0L))
    ra <- ld_regions(mk, ed); ra <- ra[lengths(ra) >= L]
    if (!length(ra)) return(data.table(tau = tau, l_min = L, n = 0L, hit = 0L))
    co <- rbindlist(lapply(ra, function(m) { mm <- map[marker %in% m]
      data.table(chr = as.character(mm$Chr[1]), lo = min(mm$Pos), hi = max(mm$Pos)) }))
    data.table(tau = tau, l_min = L, n = nrow(co),
               hit = sum(vapply(seq_len(nrow(co)), function(k)
                 any(qtn$Chr == co$chr[k] & qtn$Pos >= co$lo[k] & qtn$Pos <= co$hi[k]),
                 logical(1)))) }))
  merge(cbind(env = e, basis = b, ls_), tr, by = c("tau", "l_min"))
}
g <- rbindlist(lapply(fs, one), fill = TRUE)
fwrite(g, file.path(D, "operating_grid_all_env.csv"))
m <- g[, .(n_env = uniqueN(env), n_sig = mean(n_sig), hit = mean(hit),
           precision = mean(fifelse(n > 0, hit/n, NA_real_), na.rm = TRUE)),
       by = .(basis, tau, l_min)]
cat(sprintf("  %d environments x %d grid cells x %d bases\n",
            uniqueN(g$env), uniqueN(g[, .(tau, l_min)]), uniqueN(g$basis)))

hm <- function(v, lab, pal, dec) ggplot(m, aes(factor(tau), factor(l_min), fill = .data[[v]])) +
  geom_tile(colour = "white", linewidth = .4) +
  geom_text(aes(label = formatC(.data[[v]], format = "f", digits = dec)),
            size = 2.6, colour = "grey15") +
  facet_wrap(~ basis) +
  scale_fill_distiller(palette = pal, direction = 1, na.value = "grey92", name = NULL) +
  labs(x = expression(tau[C]), y = expression(l[min]), title = lab) +
  theme_bw(base_size = 9) +
  theme(strip.background = element_blank(), panel.grid = element_blank())

p <- (hm("n_sig", "regions surviving q_R < 0.05  (null side)", "Blues", 1) /
      hm("hit",   "regions containing a detectable QTN  (truth side)", "Greens", 2) /
      hm("precision", "precision = QTN-containing / total", "Purples", 2)) +
  plot_annotation(
    title = sprintf("tau_C x l_min operating grid, mean over %d environments", uniqueN(g$env)),
    subtitle = "V2_c1, EMMAX. Two null bases: genetic = EMMAX's own kinship model, env_orth = method-agnostic.",
    theme = theme(plot.subtitle = element_text(size = 8, colour = "grey35")))
f <- file.path(OUT, "operating_grid_tau_lmin_allenv.png")
ggsave(f, p, width = 11, height = 11, dpi = 170)
cat("  wrote", f, "\n\n")
cat("=== at tau = 0.05, averaged over environments ===\n")
print(m[tau == 0.05, .(basis, l_min, regions = round(n_sig,1), QTN = round(hit,2),
                       precision = round(precision,3))][order(basis, l_min)])
