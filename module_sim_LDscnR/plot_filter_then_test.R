## Figure: filter-then-test vs the alpha-based genome-wide scan.
## The "genome-wide" line IS the conventional analysis (BH at 0.05 over all
## units), so the comparison to alpha is the dashed reference in every panel.
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
R   <- "module_sim_LDscnR/results/filter_then_test"
OUT <- "module_sim_LDscnR/figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

d <- fread(file.path(R, "filter_then_test_clusters_emx.csv"))
## facet label ordered NUMERICALLY, not as a string (10, 50, 100)
wlev <- paste0(sort(unique(d$window_kb)), " kb window")
d[, wlab := factor(paste0(window_kb, " kb window"), levels = wlev)]
key <- c("cell","tag","env","window_kb")
gw  <- d[method == "genome_wide", c(key,"n_tp","n_sig"), with = FALSE]
setnames(gw, c("n_tp","n_sig"), c("gw_tp","gw_sig"))
f   <- merge(d[method != "genome_wide"], gw, by = key)

lab <- c(ld_w = "ld_w", size = "cluster size", MAF = "MAF", random = "random")
f[, method := factor(lab[method], levels = lab)]
pal <- c("ld_w" = "#1F6F8B", "cluster size" = "#2E7156", "MAF" = "#C1622F", "random" = "grey60")

## panel A -- true discoveries against selection size
a <- f[, .(tp = median(n_tp), lo = quantile(n_tp,.25), hi = quantile(n_tp,.75)),
       by = .(method, k, window_kb, wlab)]
gwl <- gw[, .(tp = median(gw_tp)), by = window_kb]
gwl[, wlab := factor(paste0(window_kb, " kb window"), levels = wlev)]
pA <- ggplot(a, aes(k, tp, colour = method, fill = method)) +
  geom_hline(data = gwl, aes(yintercept = tp), linetype = "22", colour = "grey25") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .13, colour = NA) +
  geom_line(linewidth = .7) + geom_point(size = 1.4) +
  scale_x_log10(breaks = sort(unique(a$k)), labels = function(z) format(z, big.mark = ",")) +
  scale_colour_manual(values = pal) + scale_fill_manual(values = pal) +
  facet_wrap(~ wlab, nrow = 1) +
  labs(x = "units selected (k), log scale", y = "true discoveries",
       title = "A  Filtering beats the genome-wide scan, but only if the filter is LD-structure-based",
       subtitle = "median over 80 panels (IQR shaded). Dashed line = genome-wide BH at 0.05, i.e. the conventional alpha analysis.") +
  theme_bw(base_size = 9) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom", legend.title = element_blank(),
        plot.subtitle = element_text(colour = "grey30", size = 7.5))

## panel B -- paired win rate against the same alpha baseline
b <- f[, {x <- n_tp - gw_tp; nz <- x[x != 0]
          .(win = if (length(nz)) 100*sum(nz > 0)/length(nz) else NA_real_, n = length(nz))},
       by = .(method, k, window_kb, wlab)]
pB <- ggplot(b, aes(k, win, colour = method)) +
  geom_hline(yintercept = 50, linetype = "22", colour = "grey25") +
  geom_line(linewidth = .7) + geom_point(size = 1.4) +
  scale_x_log10(breaks = sort(unique(b$k)), labels = function(z) format(z, big.mark = ",")) +
  scale_colour_manual(values = pal) + ylim(0, 100) +
  facet_wrap(~ wlab, nrow = 1) +
  labs(x = "units selected (k), log scale", y = "% of panels beating genome-wide",
       title = "B  Paired against the alpha analysis, panel by panel",
       subtitle = "50% (dashed) = no better than the conventional scan. Paired within panel; never a difference of medians.") +
  theme_bw(base_size = 9) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none", plot.subtitle = element_text(colour = "grey30", size = 7.5))

## panel C -- the control that matters: ld_w against cluster size
cc <- merge(f[method == "ld_w", c(key,"k","n_tp","wlab"), with=FALSE],
            f[method == "cluster size", c(key,"k","n_tp"), with=FALSE],
            by = c(key,"k"), suffixes = c("_ldw","_size"))
cs <- cc[, {x <- n_tp_ldw - n_tp_size; nz <- x[x != 0]
            .(win = if (length(nz)) 100*sum(nz>0)/length(nz) else NA_real_,
              p = if (length(nz)) binom.test(sum(nz>0), length(nz))$p.value else NA_real_)},
         by = .(k, window_kb, wlab)]
pC <- ggplot(cs, aes(k, win)) +
  geom_hline(yintercept = 50, linetype = "22", colour = "grey25") +
  geom_line(colour = "#1F6F8B", linewidth = .7) +
  geom_point(aes(shape = p < 0.05), colour = "#1F6F8B", size = 1.8) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), labels = c("n.s.","p < 0.05")) +
  scale_x_log10(breaks = sort(unique(cs$k)), labels = function(z) format(z, big.mark = ",")) +
  ylim(0, 100) + facet_wrap(~ wlab, nrow = 1) +
  labs(x = "units selected (k), log scale", y = "% of panels where ld_w > size",
       title = "C  The control that changes the story: ld_w against cluster size, same k",
       subtitle = "Mostly indistinguishable. Once ld_w has built the clusters, it adds little to choosing among them.") +
  theme_bw(base_size = 9) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom", legend.title = element_blank(),
        plot.subtitle = element_text(colour = "grey30", size = 7.5))

ggsave(file.path(OUT, "filter_then_test_summary.png"), pA / pB / pC,
       width = 10, height = 11, dpi = 190)
cat(sprintf("  written: %s\n", file.path(OUT, "filter_then_test_summary.png")))
