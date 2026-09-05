## module_sim/R_figures/fig_bgs_recomb.R
##
## Does the BGS effect on Fst depend on local recombination rate? PK:
## "measure the difference between low and high recombination regions to
## see how bgs affects patterns of genetic differentiation as a function of
## recombination rate." Two panels from R/09_bgs_recomb_pool.R:
##   A. raw mean Fst by recombination-rate bin, one line per tag -- the
##      classic linked-selection signature would show bgs pulling further
##      above nobgs as recombination falls.
##   B. the paired BGS effect itself (log2(bgs/nobgs), matched at
##      cell,rep,env,bin before pooling) by bin, with its 0-line -- whether
##      that gap is real and where it sits, directly.
## Pooled across all 7 cells (mechanism question, not a per-scenario one --
## see 09_bgs_recomb_pool.R's header).
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim"), "R", "00_config.R"))
say("=== fig_bgs_recomb ===\n\n")

POOL_PATH <- file.path(PATHS$out, "09_bgs_recomb_pool", "bgs_recomb_summary.rds")
if (!file.exists(POOL_PATH)) stop("R/09_bgs_recomb_pool.R has not produced: ", POOL_PATH)
x <- readRDS(POOL_PATH)
fst_by_tag_bin <- x$fst_by_tag_bin
bgs_effect_by_bin <- x$bgs_effect_by_bin

BIN_LEVELS <- c("zero", "low", "medium", "high")
fst_by_tag_bin[, recomb_bin := factor(recomb_bin, levels = BIN_LEVELS)]
fst_by_tag_bin[, tag := factor(tag, levels = c("nobgs", "bgs"))]
bgs_effect_by_bin[, recomb_bin := factor(recomb_bin, levels = BIN_LEVELS)]

say("[0] %d (tag,bin) Fst points ; %d bin-level BGS effects\n", nrow(fst_by_tag_bin), nrow(bgs_effect_by_bin))

pA <- ggplot(fst_by_tag_bin, aes(recomb_bin, Fst, colour = tag, group = tag)) +
  geom_line(alpha = 0.5, linewidth = 0.5) +
  geom_pointrange(aes(ymin = Fst - Fst_SE, ymax = Fst + Fst_SE), size = 0.4) +
  scale_colour_manual(values = c(nobgs = "#1565C0", bgs = "#C0392B")) +
  labs(x = NULL, y = "Fst", title = "A. Fst by recombination-rate bin") +
  theme_bw(11) + theme(panel.grid.minor = element_blank(), legend.position = "top")

pB <- ggplot(bgs_effect_by_bin, aes(recomb_bin, log2_Fst)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_pointrange(aes(ymin = log2_Fst - log2_Fst_SE, ymax = log2_Fst + log2_Fst_SE),
                  colour = "#7B1FA2", size = 0.4) +
  labs(x = "local recombination-rate bin", y = expression(log[2] * (bgs / nobgs)),
      title = "B. paired BGS effect on Fst") +
  theme_bw(11) + theme(panel.grid.minor = element_blank())

p <- pA / pB +
  plot_annotation(title = "Does BGS's effect on differentiation depend on recombination rate?",
                  subtitle = "pooled across all 7 cells; recombination-rate bins cut from each rep's own map, shared identically between bgs and nobgs")

FIG_DIR <- file.path(PATHS$module, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_DIR, "fig_bgs_recomb.pdf")
ggsave(OUT, p, width = 7.5, height = 8, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 7.5, height = 8, dpi = 200)
say("\n[1] wrote %s (+ .png)\n", OUT)
