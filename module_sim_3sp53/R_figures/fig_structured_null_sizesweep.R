## module_sim_3sp53/R_figures/fig_structured_null_sizesweep.R
##
## Visualize R/12_structured_null.R's minimum-stage-1-unit-size sweep for
## the group null. PK: "can we say that larger stage-2 clusters are more
## likely to be true positives? We should sweep for minimum cluster size."
##
## [!] Shows BOTH the naive pooled mean (dashed, mean_surrogate/max(n_obs,1)
## averaged over ALL combos) and the corrected mean (solid, restricted to
## combos with n_obs>0 at that floor) side by side -- PK's 2026-09-09
## question ("are you accounting for the fact that the number of false
## positives also drop by randomly removing clusters from the outlier
## list?") showed the naive version is an artefact: the eligible pool of
## large Stage-1 units collapses with floor, so more and more combos hit
## n_obs=0 and the ratio stops being an FDR estimate (36.5% of combos at
## floor=2, 95.2% at floor=50). The gap between the two lines below IS
## that artefact, drawn explicitly rather than left implicit.
suppressMessages({library(data.table); library(ggplot2)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53"), "R", "00_config.R"))
say("=== fig_structured_null_sizesweep ===\n\n")

SN_PATH <- file.path(PATHS$module, "results", "structured_null_summary.rds")
if (!file.exists(SN_PATH)) stop("R/13_structnull_pool.R has not produced: ", SN_PATH)
d <- readRDS(SN_PATH)
sw <- d$size_sweep[scheme == "group"]

naive <- sw[, .(realised_fdr = mean(realised_fdr, na.rm = TRUE), variant = "naive (pooled, all combos)"),
            by = size_floor]
corrected <- sw[, .(realised_fdr = mean(realised_fdr[n_obs > 0], na.rm = TRUE),
                     n = sum(n_obs > 0), frac_pos = mean(n_obs > 0),
                     variant = "corrected (n_obs>0 only)"), by = size_floor]
say("[1] naive vs corrected at each floor:\n")
print(merge(naive[, .(size_floor, naive_fdr = realised_fdr)],
            corrected[, .(size_floor, corrected_fdr = realised_fdr, frac_pos)], by = "size_floor"))

both <- rbindlist(list(naive, corrected), fill = TRUE)
both[, variant := factor(variant, levels = c("naive (pooled, all combos)", "corrected (n_obs>0 only)"))]

lab <- corrected[, .(size_floor, y = realised_fdr, lab = sprintf("n=%d\n(%.0f%%)", n, 100 * frac_pos))]

p <- ggplot(both, aes(size_floor, realised_fdr, linetype = variant)) +
  geom_line(colour = "#1565C0", linewidth = 0.7) +
  geom_point(colour = "#1565C0", size = 2) +
  geom_text(data = lab, aes(size_floor, y, label = lab), inherit.aes = FALSE,
            vjust = -0.6, size = 2.8, colour = "grey30") +
  scale_x_continuous(breaks = c(2, 3, 5, 10, 20, 50)) +
  scale_linetype_manual(values = c("naive (pooled, all combos)" = "dashed",
                                    "corrected (n_obs>0 only)" = "solid"), name = NULL) +
  expand_limits(y = 0.9) +
  labs(x = "minimum stage-1 unit size (markers)", y = "group-null realised FDR",
       title = "Minimum-cluster-size sweep: naive vs. n_obs-corrected",
       subtitle = "dashed = original pooled mean (mean_surrogate/max(n_obs,1) over ALL 2800 combos) -- collapses toward 0\nsimply because most combos stop having any real discovery at large floors, not because big clusters are cleaner.\nsolid = restricted to combos with a genuine discovery at that floor; labels show n (of 2800) and the surviving fraction.") +
  theme_bw(11) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "top")

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "fig_structured_null_sizesweep.pdf")
ggsave(OUT, p, width = 8, height = 6, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 8, height = 6, dpi = 200)
say("\n[2] wrote %s (+ .png)\n", OUT)
