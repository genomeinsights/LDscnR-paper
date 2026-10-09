## final_analysis/R_figures/figureS_qtn_detection.R
##
## Supplementary figure for the review request (Frederic, 2026-09-25): detection
## of every analysed QTN (MAF > 0.10) in the 600 primary runs by its share of
## additive variance (the detectability rule), its absolute additive variance and
## its MAF, plus the distribution of variance shares. Reads
## results/simulation_qtn_detection_c1.tsv (R/28). Bins: equal-count bins for
## absolute Va; fixed bins for the share and MAF; Wilson 95% intervals.
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
MODULE_ROOT <- path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis")
d <- fread(file.path(MODULE_ROOT, "results", "simulation_qtn_detection_c1.tsv"))

METH <- c(emmax_snp_region = "Marker-wise EMMAX", emmax_simes_region = "EMMAX Simes",
          emmax_consensus_region = "EMMAX consensus", lfmm_snp_region = "Marker-wise LFMM",
          lfmm_simes_region = "LFMM Simes")
COL <- c("Marker-wise EMMAX" = "#7a7a7a", "EMMAX Simes" = "#2a78d6", "EMMAX consensus" = "#eb6834",
         "Marker-wise LFMM" = "#7a7a7a", "LFMM Simes" = "#2a78d6")
LTY <- c("Marker-wise EMMAX" = "solid", "EMMAX Simes" = "solid", "EMMAX consensus" = "solid",
         "Marker-wise LFMM" = "22", "LFMM Simes" = "22")
wilson <- function(k, n, z = 1.96) {
  p <- k / n; den <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / den; half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  list(lo = ctr - half, hi = ctr + half)
}
long <- melt(d, id.vars = c("qtn_marker", "MAF", "Va", "p_Va", "tag", "cell", "rep", "env"),
             measure.vars = names(METH), variable.name = "method", value.name = "detected")
long[, method := factor(METH[as.character(method)], levels = METH)]
summ <- function(bin_var, x_var) {
  s <- long[, .(x = median(get(x_var)), k = sum(detected), n = .N), by = c(bin_var, "method")]
  s[, `:=`(p = k / n, lo = wilson(k, n)$lo, hi = wilson(k, n)$hi)][]
}
d_va_q <- unique(quantile(d$Va, seq(0, 1, 1 / 8)))
long[, va_bin := cut(Va, d_va_q, include.lowest = TRUE)]
long[, share_bin := cut(p_Va, c(0, .01, .02, .05, .1, .2, .5, 1), include.lowest = TRUE)]
long[, maf_bin := cut(MAF, c(.1, .2, .3, .4, .5), include.lowest = TRUE)]

th <- theme_bw(10) + theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
                           legend.position = "bottom", legend.title = element_blank())
pd <- position_dodge(width = 0.06)
curve <- function(s, xlab, logx = TRUE) {
  p <- ggplot(s, aes(x, p, colour = method, linetype = method)) +
    geom_line(linewidth = 0.6, position = pd) + geom_point(size = 1.4, position = pd) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.3, position = pd, linetype = "solid") +
    scale_colour_manual(values = COL) + scale_linetype_manual(values = LTY) +
    scale_y_continuous(limits = c(0, 1)) + labs(x = xlab, y = "Proportion of QTN detected") + th
  if (logx) p + scale_x_log10() else p
}
pA <- ggplot(d, aes(p_Va)) + geom_histogram(bins = 30, fill = "grey55", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0.05, linetype = "dashed") + scale_x_log10() +
  labs(x = "Share of the chromosome's additive variance", y = "Number of QTN") + th
pB <- curve(summ("share_bin", "p_Va"), "Share of additive variance (bin median)") +
  geom_vline(xintercept = 0.05, linetype = "dashed")
pC <- curve(summ("va_bin", "Va"), "Additive variance of the QTN (bin median)")
pD <- curve(summ("maf_bin", "MAF"), "Minor-allele frequency (bin median)", logx = FALSE)
fig <- (pA + pB) / (pC + pD) + plot_annotation(tag_levels = "A") + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

FIG_DIR <- file.path(MODULE_ROOT, "figures"); dir.create(FIG_DIR, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "figureS_qtn_detection.pdf")
ggsave(OUT, fig, width = 9, height = 7.5, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), fig, width = 9, height = 7.5, dpi = 200)
fwrite(rbindlist(list(summ("share_bin", "p_Va")[, axis := "share"], summ("va_bin", "Va")[, axis := "Va"],
                      summ("maf_bin", "MAF")[, axis := "MAF"]), use.names = TRUE, fill = TRUE),
       file.path(MODULE_ROOT, "results", "simulation_qtn_detection_binned_c1.tsv"), sep = "\t")
cat("wrote", OUT, "\n")
