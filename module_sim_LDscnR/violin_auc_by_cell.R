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
## facet_wrap, not facet_grid: the V x c design is not full (bgs5 has V0.5_c1,
## V1_c1.5, V2_c1 + the excluded V0.5_c2), and a grid spends half the figure on
## empty panels. Wrap shows only the cells that exist and absorbs new ones as
## they are simulated. The panel label carries V, c and what V means, so the
## gaps are still legible from the labels rather than from blank space.
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
sel <- c("0.5" = "strong", "1" = "medium", "2" = "weak")
long[, panel := sprintf("V = %s, c = %s\n(%s selection)", V, cc,
                        ifelse(is.na(sel[as.character(V)]), "?", sel[as.character(V)]))]
long[, panel := factor(panel, levels = unique(panel[order(V, cc)]))]

## Colour is engine x method, with the two crossed deliberately in the palette:
## HUE carries alpha vs C (grey vs blue) because that is the comparison being
## read, and SHADE carries emmax vs lfmm as the secondary split. x stays
## bgs/nobgs so every contrast of interest is within a panel.
long[, grp := interaction(engine, method, sep = " / ", lex.order = TRUE)]
PAL <- c("emmax / BH alpha" = "#C6C9CC", "lfmm / BH alpha" = "#7F868C",
         "emmax / C-score"  = "#7FBCDD", "lfmm / C-score"  = "#1F5F8B")
p <- ggplot(long, aes(tag, PR_AUC, fill = grp)) +
  geom_violin(position = position_dodge(0.8), width = 0.85,
              alpha = 0.55, colour = NA, scale = "width", trim = TRUE) +
  geom_boxplot(position = position_dodge(0.8), width = 0.13,
               outlier.shape = NA, alpha = 0.9, linewidth = 0.3) +
  stat_summary(aes(group = grp), fun = mean, geom = "point",
               position = position_dodge(0.8), shape = 23, size = 1.8,
               fill = "white", stroke = 0.4) +
  facet_wrap(~ panel, nrow = 1) +
  scale_fill_manual(values = PAL, name = NULL) +
  labs(x = NULL, y = "PR-AUC",
       title = "PR-AUC by selection regime: C-score against BH alpha",
       subtitle = paste("bgs5, 10 environments per cell, both engines and l_min 1/3 pooled.",
                        "Diamond = mean. V is selection VARIANCE, so larger V = weaker selection.",
                        "Hue = BH alpha vs C-score; shade = emmax vs lfmm.", sep = "\n")) +
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
