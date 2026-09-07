## =============================================================================
## module_9sp/R_figures/figure_permutation_comparison.R
##
## MAIN FIGURE (PK, 2026-09-07): the 3sp-vs-9sp contrast as a compact figure
## rather than only tables -- observed discovery count against the permutation
## surrogate distribution, one panel per (species, statistic). Communicates
## the same numbers as Supplementary Table 6 (module_9sp/R/03_EMMAX.R,
## module_3sp/R/03_EMMAX.R), but the point of a permutation test -- is the
## observed count somewhere a well-populated null would rarely land -- is
## immediate from a histogram in a way a p-value column is not.
##
## FREE X-AXIS PER PANEL, deliberately (matching figure_manhattan.R's free
## y-axis convention): the four surrogate distributions span wildly different
## ranges (3sp consensus max 141, 9sp consensus max 357), and forcing a shared
## axis would flatten 3sp's story (observed sitting in an almost-empty tail)
## to illegibility to accommodate 9sp's much wider spread.
##
## LINEAR, NOT LOG, y-axis: these count distributions are heavily right-skewed
## (3sp consensus: mean 3.476, but a max of 141), so a log axis would need
## zero-count-bin handling that obscures the exact visual the figure is for --
## whether the observed line falls in a dense or a near-empty part of the
## histogram is legible directly on a linear axis, and is the whole story.
## =============================================================================
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_9sp"), "R", "00_config.R"))
STAGE <- "figure_permutation_comparison"
say("=== %s ===\n\n", STAGE)

sc3 <- readRDS(path.expand("~/gitlab/LDscnR-paper/module_3sp/out/03_EMMAX/scan.rds"))
sc9 <- readRDS(file.path(PATHS$out, "03_EMMAX", "scan.rds"))

## one row per (species, statistic) -- keeps the four null objects and their labels together
## rather than four parallel variables that have to be kept in sync by hand.
PANELS <- list(
  list(species = "Three-spined stickleback", stat = "Consensus", null = sc3$consensus$null, colour = "#2C7FB8"),
  list(species = "Three-spined stickleback", stat = "Simes",     null = sc3$simes$null,     colour = "#2C7FB8"),
  list(species = "Nine-spined stickleback",  stat = "Consensus", null = sc9$consensus$null, colour = "#D95F02"),
  list(species = "Nine-spined stickleback",  stat = "Simes",     null = sc9$simes$null,     colour = "#D95F02")
)

mk_panel <- function(p) {
  surr <- p$null$surrogates; obs <- p$null$observed; pval <- p$null$p
  B <- length(surr)
  null_mean <- mean(surr)
  D <- data.table(x = surr)
  ## binwidth scaled to each panel's own range so the histogram shape is legible whether
  ## the range is ~0-140 (3sp) or ~0-360 (9sp) -- a fixed binwidth would over- or
  ## under-resolve one pair or the other.
  bw <- max(1, round(max(surr, obs) / 40))
  ## headroom past the last bin's OWN right edge (not just past max(surr,obs)), so a bin
  ## that straddles max(surr,obs) is never truncated by scale_x_continuous()'s limit.
  x_max <- (ceiling(max(surr, obs) / bw) + 1) * bw
  ## single text box (not two floating labels at each line) so the null-mean label never
  ## collides with the observed label even when the two vlines sit close together in x
  ## (9sp: observed 55 vs null mean 67.36 are only ~12 apart on a 0-360 axis) -- colour-
  ## and linetype-coded lines let the reader connect each number in the box to its own line.
  ##
  ## FIXED at the top-right corner, not anchored to x = obs: all four panels are
  ## right-skewed (tall bars near 0, thinning out towards x_max), so the top-right corner
  ## is reliably empty in every panel -- unlike anchoring on obs, which put the box on top
  ## of the histogram's own bulk once obs and the null mean sit close together within it
  ## (9sp panels: the observed/null-mean pair sits inside the bulk of the distribution,
  ## not off in an empty tail the way 3sp's does).
  lbl <- sprintf("observed = %d\nnull mean = %.2f\np = %s", obs, null_mean,
                 if (pval < 0.001) "<0.001" else sprintf("%.3f", pval))
  ggplot(D, aes(x)) +
    geom_histogram(binwidth = bw, boundary = 0, fill = "grey75", colour = "white", linewidth = 0.15) +
    geom_vline(xintercept = null_mean, colour = "grey40", linetype = "dashed", linewidth = 0.7) +
    geom_vline(xintercept = obs, colour = p$colour, linewidth = 1) +
    annotate("text", x = x_max, y = Inf, label = lbl,
             colour = p$colour, hjust = 1.05, vjust = 1.3, size = 3.2, fontface = "bold") +
    scale_x_continuous(limits = c(0, x_max), expand = c(0.01, 0)) +
    ## PLAIN ASCII hyphen with spaces, not "--" (rendered as two literal hyphens -- ggplot's
    ## plain-text rendering, unlike LaTeX, never turns "--" into a dash) and not a Unicode
    ## en-dash either (the cairo_pdf device's default font substituted it with "..." on this
    ## machine -- a font/glyph fallback issue, not worth chasing when a hyphen is unambiguous).
    labs(x = "BH-significant Stage-1 units", y = "Number of permutation draws",
        title = p$species,
        subtitle = sprintf("%s (%s permutations)", p$stat, format(B, big.mark = ","))) +
    theme_bw(11) + theme(panel.grid.minor = element_blank(),
                          plot.title = element_text(size = 11, face = "bold"),
                          plot.subtitle = element_text(size = 9.5, colour = "grey30"))
}

panels <- lapply(PANELS, mk_panel)
## tag_levels = "a": panel labels a-d in reading order (top-left, top-right, bottom-left,
## bottom-right), matching PANELS' own order -- the standard multi-panel convention, not
## embedded in each panel's own title string.
FIG <- ((panels[[1]] | panels[[2]]) / (panels[[3]] | panels[[4]])) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 13, face = "bold"))

OUT_DIR <- file.path(PATHS$out, STAGE); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_PDF <- file.path(PATHS$figures, "permutation_comparison_3sp_9sp.pdf")
ggsave(OUT_PDF, FIG, width = 9, height = 6.5, device = cairo_pdf)
ggsave(file.path(PATHS$figures, "permutation_comparison_3sp_9sp.png"), FIG, width = 9, height = 6.5, dpi = 220)
say("[1] wrote %s\n", OUT_PDF)
write_receipt(STAGE,
  inputs = c(path.expand("~/gitlab/LDscnR-paper/module_3sp/out/03_EMMAX/_receipt.rds"),
             file.path(PATHS$out, "03_EMMAX", "_receipt.rds")),
  params = list(), outputs = OUT_PDF)
say("    receipt: %s\n", receipt_path(STAGE))
