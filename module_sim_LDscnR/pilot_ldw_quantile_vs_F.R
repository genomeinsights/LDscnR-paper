## =====================================================================
## module_sim_LDscnR / pilot_ldw_quantile_vs_F.R
##
## Does association strength track local LD, and do the two carry DIFFERENT
## information? Plots the F statistic against each marker's ld_w quantile,
## coloured by r2 with the focal QTN, for the map-rebuild pilots and the
## matching bgs5 cell.
##
## The quantile, not ld_w itself, is on x: ld_w's scale differs about 2x
## between bgs5 and the pilots (median 0.0085 against 0.017-0.022), so raw
## values are not comparable across panels and the quantile is.
##
## FOCAL QTN = the QTN a marker has its HIGHEST r2 with, which is not the
## nearest one for ~56% of markers. See focal_QTN() in the parse pipeline.
##
## What to read, and what not to:
##   - the cloud is broadly FLAT, i.e. ld_w and F are close to independent.
##     That is the point -- it is why ld_w adds information to a p-value
##     rather than duplicating it.
##   - high-r2 markers concentrate at the right edge in every panel. That is
##     the C-score's premise, and it holds on the rebuilt maps.
##   - bgs5's high-r2 markers extend much further LEFT than the pilots'. Its
##     33-44% zero-recombination loci blur the quantile; the pilots carry
##     1.5-2%.
##   - n = 2 burn-in replicates, ONE environment (env2), ONE cell. env2 is an
##     extreme draw for ld_w PEAKS (see pilot_decay/README_ldw_landscapes.md);
##     the bulk distribution shown here is the stable part, the extremes are not.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/pilot_ldw_quantile_vs_F.R
## Env: SIM_ROOT, REF_CELL (default V1_c1.5_env2), OUT
## =====================================================================
suppressMessages({ library(data.table); library(ggplot2) })
ROOT <- Sys.getenv("SIM_ROOT", "/Volumes/Nemo/Nemo_sim")
REF  <- Sys.getenv("REF_CELL", "V1_c1.5_env2")
OUT  <- Sys.getenv("OUT", "module_sim_LDscnR/figures")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

SETS <- c(bgs5 = "regen_sim_data_bgs5", pilot_26 = "regen_sim_data_pilot_26",
          pilotC_26 = "regen_sim_data_pilotC_26", pilot_53 = "regen_sim_data_pilot_53",
          pilotC_53 = "regen_sim_data_pilotC_53")

load_set <- function(dir, nm) {
  ff <- list.files(file.path(ROOT, dir), pattern = "[.]rds$", full.names = TRUE)
  if (nm == "bgs5") ff <- grep(sprintf("_%s[.]rds$", gsub("\\.", "[.]", REF)), ff, value = TRUE)
  if (!length(ff)) return(NULL)
  rbindlist(lapply(ff, function(f) {
    m <- as.data.table(readRDS(f)$map)
    m[, .(dataset = nm, arm = if (grepl("_nobgs_", basename(f))) "nobgs" else "bgs",
          ld_w = ld_w_095, maxLD = max_LD_with_QTN,
          emmax = emx_F, lfmm = lfmm_F, true_QTN)]
  }), fill = TRUE)
}

d <- rbindlist(lapply(names(SETS), function(nm) load_set(SETS[[nm]], nm)), fill = TRUE)
d <- d[is.finite(ld_w)]
d[, ldw_q := (frank(ld_w, ties.method = "average") - 0.5) / .N, by = .(dataset, arm)]

l <- melt(d, id.vars = c("dataset", "arm", "ld_w", "ldw_q", "maxLD", "true_QTN"),
          measure.vars = c("emmax", "lfmm"), variable.name = "engine", value.name = "F")
l <- l[is.finite(F)]
l[, dataset := factor(dataset, levels = names(SETS))]
l[, row := factor(paste(engine, arm, sep = " / "),
                  levels = c("emmax / nobgs", "emmax / bgs", "lfmm / nobgs", "lfmm / bgs"))]
setorder(l, maxLD)   # high r2 drawn last, so it is not buried

p <- ggplot(l, aes(ldw_q, F, colour = maxLD)) +
  geom_point(size = .32, alpha = .65) +
  geom_point(data = l[true_QTN %in% TRUE], shape = 3, colour = "black", size = 2, stroke = .7) +
  facet_grid(row ~ dataset) +
  scale_colour_viridis_c(name = expression(r^2~"with focal QTN"), limits = c(0, 1)) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 3, 10, 30)) +
  scale_x_continuous(breaks = c(0, .5, 1), labels = c("0", "0.5", "1")) +
  labs(x = expression("quantile of "*ld[w]*" (within dataset x arm)"),
       y = "F statistic (log1p scale)",
       title = "Association strength against local-LD quantile, rebuilt maps vs bgs5",
       subtitle = paste(sprintf("%s, two burn-in replicates pooled per panel. Crosses = causal variants.", REF),
                        "Colour is r2 with the FOCAL (highest-r2) QTN, not the nearest.", sep = "\n")) +
  theme_bw(base_size = 10) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "top", plot.subtitle = element_text(size = 8, colour = "grey35"))

f <- file.path(OUT, "pilot_ldw_quantile_vs_F.png")
ggsave(f, p, width = 13, height = 9, dpi = 150)
cat(sprintf("  wrote %s | %d markers x %d datasets x %d arms x %d engines\n",
            f, nrow(l) / 2, uniqueN(l$dataset), uniqueN(l$arm), uniqueN(l$engine)))
