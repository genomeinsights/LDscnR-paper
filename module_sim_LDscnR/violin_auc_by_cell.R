## =====================================================================
## module_sim_LDscnR / violin_auc_by_cell.R
##
## PR-AUC across the bgs5 (V, c) design: C-score against BH alpha, as violins.
##
## c2 is EXCLUDED by default. It is the one cell whose background LD is far off
## the rest (b = 0.42 against 0.03-0.12) and which also carries ~3x the QTN per
## chromosome, so it varies two things at once and is not comparable with the
## others as a design point. Set KEEP_C2=1 to include it.
##
## The V x c grid is deliberately left as a grid even though it is not full:
## bgs5 has V0.5_c1, V1_c1.5, V2_c1 (+ the excluded V0.5_c2), so empty panels
## are real gaps in the design rather than missing data, and hiding them would
## make the design look complete when it is not.
##
## V is selection VARIANCE, so larger V = weaker selection = lower power.
## The x axis is therefore ordered by decreasing power left to right.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/violin_auc_by_cell.R [outdir]
## Env: RESULTS (dir holding the per-cell result dirs), KEEP_C2
## =====================================================================
suppressMessages({ library(data.table); library(ggplot2) })
a <- commandArgs(trailingOnly = TRUE)
OUT <- if (length(a)) a[1] else "module_sim_LDscnR/figures"
RES <- Sys.getenv("RESULTS", "module_sim_LDscnR/results")
KEEP_C2 <- nzchar(Sys.getenv("KEEP_C2", ""))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

fs <- list.files(RES, pattern = "^bgs_vs_nobgs_prauc[.]csv$",
                 recursive = TRUE, full.names = TRUE)
if (!length(fs)) stop("no bgs_vs_nobgs_prauc.csv under ", RES)
d <- unique(rbindlist(lapply(fs, fread), fill = TRUE))
cat(sprintf("  read %d file(s), %d rows, %d cells\n", length(fs), nrow(d), uniqueN(d$cell)))

## cell -> (V, c, env); cell strings look like V0.5_c1_env7
d[, V   := as.numeric(sub("^V([0-9.]+)_c.*",  "\\1", cell))]
d[, cc  := as.numeric(sub("^V[0-9.]+_c([0-9.]+)_env.*", "\\1", cell))]
d[, env := as.integer(sub(".*_env",  "", cell))]
if (!KEEP_C2) { n0 <- nrow(d); d <- d[cc != 2]
  cat(sprintf("  dropped c2: %d -> %d rows\n", n0, nrow(d))) }
if (!nrow(d)) stop("nothing left after filtering")

long <- melt(d, id.vars = c("cell","V","cc","env","tag","engine","l_min"),
             measure.vars = c("PR_AUC_C","PR_AUC_alpha"),
             variable.name = "method", value.name = "PR_AUC")
long[, method := factor(fifelse(method == "PR_AUC_C", "C-score", "BH alpha"),
                        levels = c("BH alpha", "C-score"))]
long <- long[is.finite(PR_AUC)]

## method is the COLOUR and tag the x grouping, not the other way round: the
## comparison that matters is C against alpha, and a fill contrast within a
## group is easier to read than one across adjacent groups.
p <- ggplot(long, aes(tag, PR_AUC, fill = method)) +
  geom_violin(position = position_dodge(0.8), width = 0.85,
              alpha = 0.55, colour = NA, scale = "width", trim = TRUE) +
  geom_boxplot(position = position_dodge(0.8), width = 0.13,
               outlier.shape = NA, alpha = 0.9, linewidth = 0.3) +
  stat_summary(aes(group = method), fun = mean, geom = "point",
               position = position_dodge(0.8), shape = 23, size = 1.8,
               fill = "white", stroke = 0.4) +
  facet_grid(cc ~ V, labeller = labeller(
    V  = function(x) sprintf("V = %s  (%s selection)", x,
           c("0.5" = "strong", "1" = "medium", "2" = "weak")[as.character(x)]),
    cc = function(x) sprintf("c = %s", x))) +
  scale_fill_manual(values = c("BH alpha" = "#B4B8BC", "C-score" = "#2C7FB8"), name = NULL) +
  labs(x = NULL, y = "PR-AUC",
       title = "PR-AUC by selection regime: C-score against BH alpha",
       subtitle = paste("bgs5, 10 environments per cell, both engines and l_min 1/3 pooled.",
                        "Diamond = mean. V is selection VARIANCE, so larger V = weaker selection.",
                        "Empty panels are gaps in the design, not missing runs.", sep = "\n")) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "top", plot.subtitle = element_text(size = 8, colour = "grey30"))

f <- file.path(OUT, "violin_auc_by_cell.png")
ggsave(f, p, width = 10, height = 6.5, dpi = 150)
cat("  wrote", f, "\n\n")

cat("=== mean PR-AUC by cell x method (replicate-averaged, +/- SE) ===\n")
se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
print(long[, .(n = .N, mean = round(mean(PR_AUC), 3), se = round(se(PR_AUC), 3)),
           by = .(V, cc, method)][order(V, cc, method)])
cat("\n=== C - alpha, paired within row (the quantity of interest) ===\n")
print(d[, .(n = .N, diff = round(mean(PR_AUC_C - PR_AUC_alpha, na.rm = TRUE), 3),
            se = round(se(PR_AUC_C - PR_AUC_alpha), 3),
            C_wins = sum(PR_AUC_C > PR_AUC_alpha, na.rm = TRUE)),
        by = .(V, cc)][order(V, cc)])
