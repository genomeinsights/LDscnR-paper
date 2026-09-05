## module_sim/R_figures/fig_fbeta.R
##
## F_beta = (1+b^2)*Precision*Recall / (b^2*Precision + Recall) across a
## range of beta, one line per method (arm). PK: "add it as well" -- beta
## here is the F-beta weighting parameter (beta<1 weights precision more,
## beta>1 weights recall more), NOT the Type II error rate fig_pooled_pr.R's
## "beta" row means -- a genuinely different quantity sharing the name; that
## row is left as-is (1-Recall), this is a new, separate figure.
##
## RESTRICTED TO HIGH-DISPERSAL CELLS (PK, reading fig_pooled_pr.R's
## Precision*Recall row: "for everything but high dispersal the results are
## useless regardless, so no point in comparing there"). Medians taken
## across the two tags x three high-disp cells (V0.5_c1, V1_c1, V2_c1) --
## six (tag,cell) panels per arm per beta, mirroring the reference figure's
## own "medians over N panels per arm" convention.
suppressMessages({library(data.table); library(ggplot2)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim"), "R", "00_config.R"))
say("=== fig_fbeta ===\n\n")

POOL_PATH <- file.path(PATHS$out, "05_pool", "pooled_pr.rds")
if (!file.exists(POOL_PATH)) stop("R/05_pool.R has not produced: ", POOL_PATH)
pooled <- readRDS(POOL_PATH)$pooled

HIGH_DISP_CELLS <- c("V0.5_c1", "V1_c1", "V2_c1")
pr <- pooled[cell %in% HIGH_DISP_CELLS, .(tag, cell, arm, Precision, Recall)]
say("[0] %d (tag,cell,arm) panels, cells restricted to high dispersal: %s\n",
    nrow(pr), paste(HIGH_DISP_CELLS, collapse = ", "))

## [!] FIXED 2026-09-06 (external audit, item 8): was "EMMAX representative" --
## the arm uses the consensus dosage across the cluster's markers, not a
## single representative SNP.
##
## [!] ADDED 2026-09-06 (PK: "single SNP included/excluded ... different
## colors"): emmax_snp/lfmm_snp score every significant marker as its own
## region (singletons INCLUDED); emmax_snp_clustered/lfmm_snp_clustered only
## count a marker inside a real >=2-marker Stage-1 unit (singletons EXCLUDED,
## the pre-fix behaviour kept as a comparator -- R/04_score.R). Colours pair
## by hue: saturated = included, pastel = excluded.
ARM_LEVELS <- c("emmax_consensus", "emmax_simes", "lfmm_simes",
                "emmax_snp", "emmax_snp_clustered", "lfmm_snp", "lfmm_snp_clustered")
ARM_LABELS <- c(emmax_consensus = "EMMAX consensus", emmax_simes = "EMMAX Simes", lfmm_simes = "LFMM Simes",
                emmax_snp = "EMMAX single-SNP (incl. singletons)",
                emmax_snp_clustered = "EMMAX single-SNP (excl. singletons)",
                lfmm_snp = "LFMM single-SNP (incl. singletons)",
                lfmm_snp_clustered = "LFMM single-SNP (excl. singletons)")
ARM_COLOURS <- c(emmax_consensus = "#1565C0", emmax_simes = "#26A69A", lfmm_simes = "#7B1FA2",
                 emmax_snp = "#F9A825", emmax_snp_clustered = "#FFCC80",
                 lfmm_snp = "#C0392B", lfmm_snp_clustered = "#EF9A9A")

## same beta grid as the reference figure (log-spaced, 0.25 to 4)
BETA_GRID <- c(0.25, 0.35, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 4)
f_beta <- function(P, R, b) (1 + b^2) * P * R / (b^2 * P + R)

fb <- pr[, {
  vals <- sapply(BETA_GRID, function(b) f_beta(Precision, Recall, b))
  list(beta = BETA_GRID, F_beta = vals)
}, by = .(tag, cell, arm)]

say("[1] median F_beta over %d (tag,cell) panels per (arm,beta)\n", uniqueN(pr, by = c("tag", "cell")))
med <- fb[, .(F_beta = median(F_beta)), by = .(arm, beta)]
med[, arm := factor(arm, levels = ARM_LEVELS, labels = ARM_LABELS[ARM_LEVELS])]
names(ARM_COLOURS) <- ARM_LABELS[names(ARM_COLOURS)]

## [!] FIXED 2026-09-06: subtitle used to hardcode "6 (tag x cell) panels"
## (2 tags x 3 high-dispersal cells) -- wrong once a cell is missing (bgs5
## currently has only 2 of the 3 target high-dispersal cells). Now built
## from the actual count.
n_panels <- uniqueN(pr, by = c("tag", "cell"))
p <- ggplot(med, aes(beta, F_beta, colour = arm)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  scale_x_log10(breaks = BETA_GRID) +
  scale_colour_manual(values = ARM_COLOURS, name = "method") +
  labs(x = expression(beta ~ "  (< 1 weights precision, > 1 weights recall; log scale)"),
      y = expression("median " * F[beta]),
      title = "Which analysis is better depends entirely on beta",
      subtitle = bquote(F[beta] == frac((1+beta^2)*PR, beta^2*P + R) ~ "  -- medians over" ~ .(n_panels) ~ "(tag x cell) panels per method, high-dispersal cells only")) +
  theme_bw(11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "fig_fbeta.pdf")
ggsave(OUT, p, width = 9, height = 5.5, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 9, height = 5.5, dpi = 200)
say("\n[2] wrote %s (+ .png)\n", OUT)
