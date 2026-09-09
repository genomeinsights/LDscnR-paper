## module_sim_3sp53/R_figures/fig_structured_null.R
##
## Visualize R/12_structured_null.R's three calibration-null schemes
## (group/mvn/spatial) at the base size floor, by cell x tag. PK: "run the
## structured nulls... pull the same numbers as module_3sp/9sp."
##
## [!] Reports realised_fdr CONDITIONAL ON n_obs>0 (see the 2026-09-09
## correction in README.md and the caveat comment on R/12's
## .fdr_by_floor()) -- max(n_obs,1) makes the ratio meaningless once a
## combo has zero genuine observed discoveries, which is common enough
## (37% of combos even at the base floor) that the unconditional pooled
## mean is a materially different, misleading number. Point size encodes
## the number of genuine-discovery combos behind each mean (out of 200
## per cell x tag x scheme) so a thin/underpowered estimate is visible,
## not hidden.
##
## [!] SE errorbars added 2026-09-09 (PK: "add SE's if possible") --
## computed from a SINGLE-LEVEL aggregation straight to (tag,cell,scheme)
## pooling across BOTH arms and all (rep,env) combos at once (400 rows
## before the n_obs>0 filter), rather than the previous two-stage "mean
## of two arm-level means" -- simpler, and a well-defined n/SE needs one
## aggregation level, not two nested means.
suppressMessages({library(data.table); library(ggplot2)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53"), "R", "00_config.R"))
say("=== fig_structured_null ===\n\n")

.se <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) return(NA_real_); stats::sd(x) / sqrt(length(x)) }

SN_PATH <- file.path(PATHS$module, "results", "structured_null_summary.rds")
if (!file.exists(SN_PATH)) stop("R/13_structnull_pool.R has not produced: ", SN_PATH)
d <- readRDS(SN_PATH)
sw <- d$size_sweep[size_floor == 2]  ## base floor, matches the headline `summary` numbers

CELL_LEVELS <- c("V0.5_c1", "V1_c1", "V2_c1", "V0.5_c1.5", "V1_c1.5", "V2_c1.5", "V0.5_c2")
SCHEME_LEVELS <- c("mvn", "group", "spatial")
SCHEME_LABELS <- c(mvn = "mvn (negative control)", group = "group (5-population)", spatial = "spatial (individual)")
SCHEME_COLOURS <- c(mvn = "#546E7A", group = "#1565C0", spatial = "#C62828")

agg <- sw[, .(realised_fdr = mean(realised_fdr[n_obs > 0], na.rm = TRUE), se = .se(realised_fdr[n_obs > 0]),
              n_pos = sum(n_obs > 0), n_total = .N),
          by = .(tag, cell, scheme)]  ## pools both arms + all rep x env combos directly
agg[, cell := factor(cell, levels = CELL_LEVELS)]
agg[, scheme := factor(scheme, levels = SCHEME_LEVELS)]
agg[, tag := factor(tag, levels = c("nobgs", "bgs"))]

say("[1] %d (tag,cell,scheme) points, n_pos range %d-%d of %d combos each\n",
    nrow(agg), min(agg$n_pos), max(agg$n_pos), max(agg$n_total))

p <- ggplot(agg, aes(cell, realised_fdr, colour = scheme)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_errorbar(aes(ymin = realised_fdr - se, ymax = realised_fdr + se),
               width = 0.2, position = position_dodge(width = 0.5), alpha = 0.6) +
  geom_point(aes(size = n_pos), position = position_dodge(width = 0.5), alpha = 0.85) +
  facet_grid(scheme ~ tag, scales = "free_y",
             labeller = labeller(scheme = SCHEME_LABELS)) +
  scale_colour_manual(values = SCHEME_COLOURS, guide = "none") +
  scale_size_continuous(name = "genuine-discovery\ncombos (of 400)", range = c(1, 5)) +
  labs(x = NULL, y = "realised FDR (mean +/- SE, surrogate sig. units / observed, n_obs>0 only)",
       title = "Structured-null calibration check, full 1400-combo grid",
       subtitle = "group is conservative as predicted by env ICC; mvn negative control holds; spatial reveals a dispersal-dependent\nisolation-by-distance confound. Point size = how many of 400 (both arms x 200 rep-env) combos had a discovery.") +
  theme_bw(11) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "right")

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "fig_structured_null.pdf")
ggsave(OUT, p, width = 8, height = 8, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 8, height = 8, dpi = 200)
say("\n[2] wrote %s (+ .png)\n", OUT)
