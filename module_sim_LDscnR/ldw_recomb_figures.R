## =====================================================================
## module_sim_LDscnR / ldw_recomb_figures.R
##
## What the local-LD statistic ld_w actually tracks, on one pooled genome.
## Three figures, all on the replicates concatenated into a single 20-chromosome
## pseudo-genome (10 bundle files x 2 chromosomes each, prefixed R1_..R10_):
##
##   1. ldw_manhattan_*      ld_w along the genome, coloured by r2 with the FOCAL QTN
##   2. ldw_recomb500kb_*    ld_w with a recombination track over it
##   3. ldw_vs_recomb*       ld_w against recombination rate, and the same
##                           coloured by r2 with the focal QTN
##
## FOCAL, not nearest. map$max_LD_with_QTN is max r2 over the QTN on that
## chromosome and map$focal_QTN is the QTN achieving it (see focal_QTN() in the
## parse pipeline: max.col over the r2 matrix). On chromosomes with >=2 QTN the
## focal QTN is NOT the nearest one for 56% of markers, so the two labels are
## not interchangeable.
##
## Two presentation choices, both forced by the data rather than taste:
##   - the recombination track is a BINNED MEDIAN. rec_rate is extremely skewed
##     (median ~0.9, max ~16700), so rescaling raw per-marker values puts the
##     whole track on the floor with a couple of spikes. 100 kb bins were also
##     unusable -- ~4300 noisy bins that buried the points -- hence 500 kb.
##   - the scatter uses log1p(rec_rate) (defined at the exact zeros the map
##     contains) and log10(ld_w). Spearman is reported alongside because it is
##     invariant to both, so the relationship cannot be a transform artefact.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/ldw_recomb_figures.R
## Env: SIM_DATA, CELL (default V0.5_c1_env1), BIN (default 5e5), OUT
## =====================================================================
suppressMessages({ library(data.table); library(ggplot2); library(patchwork); library(LDscnR) })
SIM  <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
CELL <- Sys.getenv("CELL", "V0.5_c1_env1")
BIN  <- as.numeric(Sys.getenv("BIN", "5e5"))
OUT  <- Sys.getenv("OUT", "module_sim_LDscnR/figures")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

pool_map <- function(tag) {
  ff <- list.files(SIM, full.names = TRUE,
    pattern = sprintf("^adapt_%s_chr[0-9]+_%s[.]rds$", tag, gsub("\\.", "[.]", CELL)))
  if (!length(ff)) return(NULL)
  ff <- ff[order(as.integer(sub(".*_chr([0-9]+)_.*", "\\1", basename(ff))))]
  m <- rbindlist(lapply(seq_along(ff), function(i) {
    x <- readRDS(ff[i]); mm <- as.data.table(x$map)
    mm[, `:=`(Chr = paste0("R", i, "_", Chr), marker = paste0("R", i, "_", marker))]
    mm[, .(marker, Chr, Pos, ld_w_095, max_LD_with_QTN, rec_rate, true_QTN)]
  }), fill = TRUE)
  m <- m[is.finite(ld_w_095)]
  ## R1..R10 numerically; the default factor order puts R10 before R1
  lv <- unique(m$Chr)
  lv <- lv[order(as.integer(sub("^R([0-9]+)_.*", "\\1", lv)), sub("^R[0-9]+_", "", lv))]
  m[, `:=`(Chr = factor(Chr, levels = lv), x = Pos / 1e6)][]
}

for (tag in c("bgs", "nobgs")) {
  m <- pool_map(tag); if (is.null(m)) { cat("  [skip]", tag, "\n"); next }
  ymax <- max(m$ld_w_095, na.rm = TRUE)
  qt   <- m[true_QTN %in% TRUE]
  lab_focal <- expression(r^2~"with focal QTN")

  ## 1. ld_w along the genome, coloured by LD with the focal QTN
  p1 <- ld_manhattan(map = m, value = "ld_w_095",
         value_label = expression(ld[w]~"("*rho*"=0.95)"),
         colour = "max_LD_with_QTN", colour_label = lab_focal, qtn = qt$marker,
         title = sprintf("%s / %s -- local LD, coloured by LD with the focal (highest-r2) causal variant (%d markers, %d QTN)",
                         CELL, tag, nrow(m), nrow(qt))) +
    theme(strip.background = element_blank())
  ggsave(file.path(OUT, sprintf("ldw_manhattan_%s_%s.png", CELL, tag)), p1,
         width = 15, height = 5, dpi = 150)

  ## 2. the same, with a recombination track
  trk  <- m[is.finite(rec_rate), .(rec = median(rec_rate, na.rm = TRUE)),
            by = .(Chr, bin = floor(Pos / BIN))]
  rmax <- max(trk$rec, na.rm = TRUE)
  trk[, `:=`(x = (bin + 0.5) * BIN / 1e6, y = rec / rmax * ymax)]
  p2 <- ggplot(m, aes(x, ld_w_095)) +
    geom_point(colour = "grey68", size = 0.5, alpha = 0.9) +
    geom_line(data = trk, aes(x, y), colour = "black", linewidth = 0.4) +
    geom_point(data = qt, aes(x, ld_w_095), shape = 3, colour = "black",
               size = 2.4, stroke = 0.7) +
    facet_wrap(~ Chr, nrow = 1, scales = "free_x") +
    scale_y_continuous(name = expression(ld[w]~"("*rho*"=0.95)"),
      sec.axis = sec_axis(~ . / ymax * rmax, name = "rec. rate (binned median)")) +
    labs(x = "Position (Mbp)",
         title = sprintf("%s / %s -- local LD and the recombination environment (%g kb bins)",
                         CELL, tag, BIN / 1e3)) +
    theme_bw(base_size = 10) +
    theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 6))
  ggsave(file.path(OUT, sprintf("ldw_recomb%gkb_%s_%s.png", BIN/1e3, CELL, tag)), p2,
         width = 16, height = 5, dpi = 150)

  ## 3. ld_w against recombination: per marker, and at window scale
  d <- m[is.finite(rec_rate) & ld_w_095 > 0]
  rs <- suppressWarnings(cor(d$rec_rate, d$ld_w_095, method = "spearman"))
  pa <- ggplot(d, aes(log1p(rec_rate), log10(ld_w_095))) +
    geom_bin2d(bins = 70) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), colour = "black",
                linewidth = 0.6, se = FALSE) +
    scale_fill_viridis_c(trans = "log10", name = "markers") +
    labs(x = "log1p(recombination rate)", y = expression(log[10]~ld[w]),
         title = sprintf("per marker (n = %s)", format(nrow(d), big.mark = ",")),
         subtitle = sprintf("Spearman rho = %.3f", rs)) +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())
  b <- m[is.finite(rec_rate), .(rec = median(rec_rate, na.rm = TRUE),
                                ldw = median(ld_w_095, na.rm = TRUE), n = .N),
         by = .(Chr, bin = floor(Pos / BIN))][n >= 20]
  rb <- suppressWarnings(cor(b$rec, b$ldw, method = "spearman"))
  pb <- ggplot(b, aes(log1p(rec), log10(ldw))) +
    geom_point(colour = "grey45", size = 0.9, alpha = 0.7) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), colour = "black",
                linewidth = 0.6, se = TRUE) +
    labs(x = "log1p(median recombination rate)", y = expression(log[10]~median~ld[w]),
         title = sprintf("%g kb windows (n = %d)", BIN/1e3, nrow(b)),
         subtitle = sprintf("Spearman rho = %.3f", rb)) +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())
  ggsave(file.path(OUT, sprintf("ldw_vs_recomb_%s_%s.png", CELL, tag)),
         pa + pb + plot_annotation(title = sprintf("%s / %s -- local LD against recombination rate", CELL, tag),
                                   theme = theme(plot.title = element_text(size = 12))),
         width = 11, height = 4.6, dpi = 150)

  ## 4. the same plane, coloured by LD with the focal QTN
  dd <- d[is.finite(max_LD_with_QTN)]
  setorder(dd, max_LD_with_QTN)      # high r2 drawn last
  dd[, `:=`(X = log1p(rec_rate), Y = log10(ld_w_095))]
  qa <- ggplot(dd, aes(X, Y, colour = max_LD_with_QTN)) +
    geom_point(size = 0.45, alpha = 0.75) +
    geom_point(data = dd[true_QTN %in% TRUE], shape = 3, colour = "black",
               size = 2.6, stroke = 0.8) +
    scale_colour_viridis_c(name = lab_focal, limits = c(0, 1)) +
    labs(x = "log1p(recombination rate)", y = expression(log[10]~ld[w]),
         title = sprintf("per marker (n = %s)", format(nrow(dd), big.mark = ",")),
         subtitle = "crosses = the causal variants; focal = the QTN a marker has its HIGHEST r2 with") +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())
  ## mean per 2-D cell: unlike the panel above, this does not depend on draw order
  qb <- ggplot(dd, aes(X, Y, z = max_LD_with_QTN)) +
    stat_summary_2d(bins = 55, fun = mean) +
    scale_fill_viridis_c(name = expression(mean~r^2~"(focal)"), limits = c(0, 1)) +
    labs(x = "log1p(recombination rate)", y = expression(log[10]~ld[w]),
         title = "mean LD with focal QTN per cell",
         subtitle = "same axes; colour is an average, not an overplot") +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())
  ggsave(file.path(OUT, sprintf("ldw_vs_recomb_qtnLD_%s_%s.png", CELL, tag)),
         qa + qb + plot_annotation(
           title = sprintf("%s / %s -- ld_w vs recombination, coloured by LD with the focal QTN", CELL, tag),
           theme = theme(plot.title = element_text(size = 12))),
         width = 12, height = 4.8, dpi = 150)

  hi <- dd[max_LD_with_QTN > 0.8]
  cat(sprintf("  %-5s | %s markers | per-marker rho %+.3f | %gkb-window rho %+.3f | r2>0.8 median log1p(rec) %.2f vs %.2f genome-wide\n",
              tag, format(nrow(d), big.mark = ","), rs, BIN/1e3, rb,
              median(hi$X), median(dd$X)))
}
