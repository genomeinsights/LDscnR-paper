## =====================================================================
## module_sim_LDscnR / region_size_distribution.R
##
## How big are the called clusters, for BOTH selection rules, with NO l_min
## filter? The l_min = 3 operating point removes three quarters of them
## (1738 -> 458) and the entire left mode with it, so the unfiltered view is the
## one that shows what the rules actually produce.
##
## THE TWO RULES SHARE A NOISE FLOOR AND DIFFER IN SIGNAL. QTN-free clusters are
## near-identical between them -- 62% singletons, median 1, max ~107-143 -- so
## the C-score is not producing a different KIND of junk, just more of the same
## (1628 against 766). The difference is in the QTN-containing clusters: 110
## against 90, median size 35 against 12, 40% over 50 markers against 22%.
##
## That is why l_min works better for the C-score than for BH: the signal and
## noise distributions are further apart under C, so a size filter discriminates
## better, which is also why C's advantage GROWS as l_min tightens.
##
## And it reframes tau_C. On its own the threshold admits a great deal of noise
## -- the modal C-score cluster is a single isolated marker and essentially none
## of those contain a QTN. The rule that works is tau_C AND an l_min, and the
## second does more of the work than the first.
##
## Two panels: SNPs per cluster on a log axis, split by whether the cluster
## contains a detectable QTN; and clusters per dataset, which is the other thing
## "number of regions" can mean.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/region_size_distribution.R
## Env: NULL_DIR, OUT
## =====================================================================
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
D   <- Sys.getenv("NULL_DIR", "module_sim_LDscnR/results/nulls_V2_c1")
OUT <- Sys.getenv("OUT", "module_sim_LDscnR/figures")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
## ALL clusters at tau_C = 0.05, with NO l_min filter -- the regions_*.csv files
## are written at the l_min = 3 operating point, which removes the entire left
## mode. Rebuilt from the scans so singletons and pairs are included: they are
## most of the distribution and leaving them out hides the shape.
A   <- Sys.getenv("SIM_INPUTS", "/Volumes/Nemo/Nemo_sim/analysis_inputs")
TAU <- as.numeric(Sys.getenv("TAU", "0.05"))
suppressMessages(library(LDscnR))
fs <- list.files(D, pattern = "^scan_.*[.]rds$", full.names = TRUE)
d <- rbindlist(lapply(fs, function(f) {
  x <- readRDS(f)
  e <- as.integer(sub(".*_env([0-9]+)_.*", "\\1", basename(f)))
  b <- sub(".*_emmax_(.+)_B[0-9]+\\.rds$", "\\1", basename(f))
  panel <- readRDS(file.path(A, sprintf("panel_V2_c1_env%d.rds", e)))
  map <- flag_true_qtns(as.data.table(panel$map))
  qtn <- map[true_pos_QTN %in% TRUE, .(Chr = as.character(Chr), Pos)]
  if (!nrow(qtn)) return(NULL)
  ## BOTH rules, clustered identically and with no l_min filter, so the two
  ## distributions are directly comparable: only the selection differs.
  C  <- x$null$C_obs
  pv <- readRDS(file.path(A, sprintf("pvals_V2_c1_env%d_emmax_%s_B100.rds", e, b)))
  q  <- p.adjust(pv$p_obs, "BH")
  sets <- list(`C-score (tau 0.05)` = names(C)[which(C >= TAU)],
               `BH alpha 0.05`      = names(q)[which(q < 0.05)])
  rbindlist(lapply(names(sets), function(rule) {
    mk <- sets[[rule]]; if (!length(mk)) return(NULL)
    ra <- ld_regions(mk, x$edges)                    # NO l_min filter
    rbindlist(lapply(ra, function(m) { mm <- map[marker %in% m]
      ch <- as.character(mm$Chr[1]); lo <- min(mm$Pos); hi <- max(mm$Pos)
      data.table(env = e, basis = b, rule = rule, size = length(m),
                 has_qtn = any(qtn$Chr == ch & qtn$Pos >= lo & qtn$Pos <= hi)) })) })) }),
  fill = TRUE)
d <- d[!is.na(size)]
d[, qtn := fifelse(has_qtn %in% TRUE, "contains a QTN", "no QTN")]

cat(sprintf("  %d regions over %d datasets\n", nrow(d), uniqueN(d[, .(env, basis)])))
cat(sprintf("  size: min %.0f, median %.0f, mean %.1f, q90 %.0f, q99 %.0f, max %.0f\n",
            min(d$size), median(d$size), mean(d$size),
            quantile(d$size, .90), quantile(d$size, .99), max(d$size)))
cat(sprintf("  the largest 1%% of regions hold %.0f%% of all clustered markers\n",
            100 * sum(d[size >= quantile(size, .99)]$size) / sum(d$size)))

p1 <- ggplot(d, aes(size, fill = qtn)) +
  facet_wrap(~ rule, ncol = 1, scales = "free_y") +
  geom_histogram(bins = 40, colour = "white", linewidth = .2, position = "stack") +
  scale_x_log10(breaks = c(1, 3, 10, 30, 100, 300, 1000)) +
  scale_fill_manual(values = c("contains a QTN" = "#2E7156", "no QTN" = "#B4BCC0"), name = NULL) +
  labs(x = "SNPs in region (log scale)", y = "regions",
       title = "All clusters, no l_min filter -- both selection rules",
       subtitle = "Same clustering, same truth rule, same datasets. Only the marker selection differs.") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        plot.subtitle = element_text(size = 8, colour = "grey35"))

n <- d[, .(regions = .N), by = .(env, basis, rule)]
p2 <- ggplot(n, aes(regions)) +
  facet_wrap(~ rule, ncol = 1, scales = "free") +
  geom_histogram(bins = 18, fill = "#1F6F8B", colour = "white", linewidth = .2) +
  labs(x = "regions called per dataset", y = "datasets",
       title = "Regions per dataset",
       subtitle = sprintf("%d datasets, median %.0f regions, range %.0f-%.0f",
                          nrow(n), median(n$regions), min(n$regions), max(n$regions))) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 8, colour = "grey35"))

ggsave(file.path(OUT, "cluster_size_distribution_all.png"),
       (p1 | p2) + plot_layout(widths = c(1.6, 1)) +
       plot_annotation(title = "V2_c1, EMMAX -- what each rule selects before any l_min filter",
         theme = theme(plot.title = element_text(size = 12))),
       width = 13, height = 8, dpi = 170)
cat("  wrote region_size_distribution.png\n")
print(d[, .(clusters = .N, singletons = paste0(round(100*mean(size==1)),"%"),
            median_size = median(size), max_size = max(size),
            pct_ge_50 = round(100*mean(size >= 50))), by = .(rule, qtn)][order(rule, qtn)])
