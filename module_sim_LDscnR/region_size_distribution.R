## =====================================================================
## module_sim_LDscnR / region_size_distribution.R
##
## How big are the called regions? The distribution is heavily right-skewed --
## most regions are a handful of markers and a few are enormous -- which is the
## reason s_R ranks as well as it does (see region_ranking_tests.R: s_R IS size)
## and the reason a mean region size is a poor summary.
##
## Two panels: SNPs per region on a log axis, split by whether the region
## contains a detectable QTN; and regions per dataset, which is the other thing
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
fs <- list.files(D, pattern = "^regions_.*[.]csv$", full.names = TRUE)
d <- rbindlist(lapply(fs, function(f) { r <- fread(f)
  r[, `:=`(env = as.integer(sub(".*_env([0-9]+)_.*", "\\1", basename(f))),
           basis = sub(".*_emmax_(.+)_B[0-9]+\\.csv$", "\\1", basename(f)))] }), fill = TRUE)
d <- d[!is.na(size)]
d[, qtn := fifelse(has_qtn %in% TRUE, "contains a QTN", "no QTN")]

cat(sprintf("  %d regions over %d datasets\n", nrow(d), uniqueN(d[, .(env, basis)])))
cat(sprintf("  size: min %.0f, median %.0f, mean %.1f, q90 %.0f, q99 %.0f, max %.0f\n",
            min(d$size), median(d$size), mean(d$size),
            quantile(d$size, .90), quantile(d$size, .99), max(d$size)))
cat(sprintf("  the largest 1%% of regions hold %.0f%% of all clustered markers\n",
            100 * sum(d[size >= quantile(size, .99)]$size) / sum(d$size)))

p1 <- ggplot(d, aes(size, fill = qtn)) +
  geom_histogram(bins = 40, colour = "white", linewidth = .2, position = "stack") +
  scale_x_log10(breaks = c(3, 10, 30, 100, 300, 1000)) +
  scale_fill_manual(values = c("contains a QTN" = "#2E7156", "no QTN" = "#B4BCC0"), name = NULL) +
  labs(x = "SNPs in region (log scale)", y = "regions",
       title = "Region size is heavily right-skewed",
       subtitle = sprintf("median %.0f SNPs, max %.0f. Log axis -- on a linear one the tail is invisible.",
                          median(d$size), max(d$size))) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        plot.subtitle = element_text(size = 8, colour = "grey35"))

n <- d[, .(regions = .N), by = .(env, basis)]
p2 <- ggplot(n, aes(regions)) +
  geom_histogram(bins = 18, fill = "#1F6F8B", colour = "white", linewidth = .2) +
  labs(x = "regions called per dataset", y = "datasets",
       title = "Regions per dataset",
       subtitle = sprintf("%d datasets, median %.0f regions, range %.0f-%.0f",
                          nrow(n), median(n$regions), min(n$regions), max(n$regions))) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 8, colour = "grey35"))

ggsave(file.path(OUT, "region_size_distribution.png"),
       (p1 | p2) + plot_layout(widths = c(1.6, 1)) +
       plot_annotation(title = "V2_c1, EMMAX, tau_C 0.05 / l_min 3 -- what the called regions look like",
         theme = theme(plot.title = element_text(size = 12))),
       width = 13, height = 5, dpi = 170)
cat("  wrote region_size_distribution.png\n")
print(d[, .(regions = .N, median_size = median(size), max_size = max(size),
            pct_ge_50 = round(100*mean(size >= 50)), with_QTN = sum(has_qtn %in% TRUE)), by = qtn])
