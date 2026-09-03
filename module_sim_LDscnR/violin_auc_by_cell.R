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
d_all <- copy(d)   # every cell, for the numeric tests below
if (!KEEP_C2) { n0 <- nrow(d); d <- d[cc != 2]
  cat(sprintf("  dropped c2 FROM THE FIGURE: %d -> %d rows\n", n0, nrow(d))) }
if (!nrow(d)) stop("nothing left after filtering")

## Prefer the COMMON-SUPPORT AUCs when the CSV has them. pr_auc() integrates to
## whatever recall a sweep reached, and the two arms do not reach the same
## place: measured here, BH alpha attains higher max recall than C in nearly
## every row (C reaches further in only ~2-5% of them), despite C sweeping a
## 1.6-2.7x larger candidate set. So the raw comparison is integrated over a
## wider domain for alpha and is biased AGAINST C. The *_cs columns truncate
## both arms at the recall both attain, which is the like-for-like number.
## Falls back to the raw columns for older CSVs that lack them.
CS <- all(c("PR_AUC_C_cs","PR_AUC_alpha_cs") %in% names(d)) &&
      any(is.finite(d$PR_AUC_C_cs))
mv <- if (CS) c("PR_AUC_C_cs","PR_AUC_alpha_cs") else c("PR_AUC_C","PR_AUC_alpha")
cat(sprintf("  metric: %s\n", if (CS) "common-support AUC (both arms truncated at shared max recall)"
                               else "RAW AUC (no _cs columns present)"))
long <- melt(d, id.vars = c("cell","V","cc","env","tag","engine","l_min"),
             measure.vars = mv, variable.name = "method", value.name = "PR_AUC")
long[, method := factor(fifelse(grepl("_C", method), "C-score", "BH alpha"),
                        levels = c("BH alpha", "C-score"))]
long <- long[is.finite(PR_AUC)]
## Panels are ordered by MEASURED difficulty (mean alpha PR-AUC), not by V.
## V is selection variance, so larger V is weaker selection -- but power here is
## governed by c, not V: at V = 0.5 the alpha AUC runs 0.539 (c=1) to 0.031
## (c=2), a 17-fold range at identical V, and V = 2 is EASIER than V = 1.
## Ordering by V would therefore read as a difficulty gradient that does not
## exist. The label carries the measured value so the axis states itself.
diff_by <- long[method == "BH alpha", .(a = mean(PR_AUC, na.rm = TRUE)), by = .(V, cc)]
long <- merge(long, diff_by, by = c("V", "cc"), all.x = TRUE)
long[, panel := sprintf("V = %s, c = %s\nalpha AUC = %.2f", V, cc, a)]
long[, panel := factor(panel, levels = unique(panel[order(-a)]))]

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
  labs(x = NULL, y = if (CS) "PR-AUC (common support)" else "PR-AUC",
       title = "PR-AUC by selection regime: C-score against BH alpha",
       subtitle = paste(sprintf("bgs5, 10 environments per cell, both engines and l_min 1/3 pooled. %s.",
                           if (CS) "Both arms integrated to the recall both attain" else "Raw AUCs"),
                        "Panels ordered by measured difficulty (mean BH-alpha PR-AUC), easiest first.",
                        "V is selection variance; power is governed by c, not V, so V is not the difficulty axis.",
                        "Hue = BH alpha vs C-score; shade = emmax vs lfmm.", sep = "\n")) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "top", plot.subtitle = element_text(size = 8, colour = "grey30"))

f <- file.path(OUT, "violin_auc_by_cell.png")
ggsave(f, p, width = 10, height = 6.5, dpi = 150)
cat("  wrote", f, "\n\n")

se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
cat("=== mean PR-AUC by cell x method ===\n")
print(long[, .(n = .N, alpha_AUC = round(a[1], 3), mean = round(mean(PR_AUC), 3)),
           by = .(V, cc, method)][order(-alpha_AUC, method)])

## SEs are computed over GENOMES, not rows. A genome is (cell, env, tag); the
## 2 engines x 2 l_min rows within it are four measurements OF that genome, not
## four independent draws. Pooling them quadruples n and shrinks the SE by ~1.5-1.9x
## here, which is pseudo-replication: it turned t = -3.2 into a real-looking
## result on V0.5_c1 when the genome-level value is -1.7.
d[, env := as.integer(sub(".*_env", "", cell))]
d[, gap := if (CS) PR_AUC_C_cs - PR_AUC_alpha_cs else PR_AUC_C - PR_AUC_alpha]
gen <- d[is.finite(gap), .(gap = mean(gap)), by = .(V, cc, env, tag)]
cat("\n=== C - alpha, SE over genomes (cell x env x tag) -- the independent unit ===\n")
print(gen[, .(n_genomes = .N, gap = round(mean(gap), 3), se = round(se(gap), 3),
              t = round(mean(gap) / se(gap), 1), C_wins = sum(gap > 0)),
          by = .(V, cc)][order(V, cc)])
cat("\n(row-level SEs, shown for comparison only -- these are pseudo-replicated)\n")
print(d[is.finite(gap), .(n_rows = .N, se_row = round(se(gap), 3)), by = .(V, cc)][order(V, cc)])

## ---- absolute vs relative advantage -----------------------------------
## Two questions that are easy to conflate, and that this design answers
## differently:
##   ABSOLUTE  how many PR-AUC points does C buy?
##   RELATIVE  how much of what is achievable does it buy?
## Baselines differ ~17-fold across cells (alpha AUC 0.03 to 0.54), so absolute
## differences are comparable WITHIN a cell and not across -- the same
## fold-vs-absolute problem as enrichment against a base rate that varies.
##
## Ratios are taken PER GENOME with engine x l_min averaged first, and only
## where the denominator is > 0: one genome has alpha AUC of exactly 0, which
## makes the ratio undefined and, worse, would enter a Spearman as +Inf and
## rank top. Medians are reported alongside means because the mean is inflated
## by the small-denominator tail (V1_c1.5: mean 64% vs median 40%).
## Runs on EVERY cell, including c2. c2 is dropped from the FIGURE because it
## varies two things at once, but excluding it from the test would flatter the
## absolute claim: c2 is the hardest cell and the one that breaks monotonicity
## at the low end, so dropping it moves absolute-gap-vs-power from rho -0.19
## (ns) to -0.38 (p = 0.003). The c2-excluded version is reported below as a
## sensitivity, not as the headline.
cs_ok <- all(c("PR_AUC_C_cs","PR_AUC_alpha_cs") %in% names(d_all)) && any(is.finite(d_all$PR_AUC_C_cs))
if (cs_ok) {
  gA <- d_all[is.finite(PR_AUC_C_cs) & is.finite(PR_AUC_alpha_cs),
              .(C = mean(PR_AUC_C_cs), A = mean(PR_AUC_alpha_cs)), by = .(V, cc, env, tag)]
} else {
  gA <- d_all[is.finite(PR_AUC_C) & is.finite(PR_AUC_alpha),
              .(C = mean(PR_AUC_C), A = mean(PR_AUC_alpha)), by = .(V, cc, env, tag)]
}
gA[, gap := C - A]
nz <- gA[A > 0][, rel := gap / A][]
cat(sprintf("\n=== absolute vs relative advantage (%d genomes; %d with alpha AUC == 0 excluded from ratios) ===\n",
            nrow(gA), sum(gA$A == 0)))
print(merge(
  gA[, .(power = round(mean(A), 3), abs_gap = round(mean(gap), 3), abs_se = round(se(gap), 3)),
     by = .(V, cc)],
  nz[, .(n_ratio = .N, rel_median = paste0(round(100 * median(rel), 1), "%"),
         rel_mean = paste0(round(100 * mean(rel), 1), "%")), by = .(V, cc)],
  by = c("V", "cc"))[order(-power)])

## Reported at TWO levels, deliberately, instead of one pooled p-value.
## Power varies mostly BETWEEN cells (~77% of its variance), so a pooled
## Spearman over ~79 genomes treats a four-point comparison as 79 independent
## observations and returns a p-value like 3e-06 that no reviewer should
## believe. That is the same pseudo-replication as the row-vs-genome problem
## above, one level up. The between-cell component has effective n = 4; the
## within-cell component is real but individually underpowered, and the ten
## environments inside a cell share burn-ins besides.
sp <- function(x, y) { ct <- suppressWarnings(stats::cor.test(x, y, method = "spearman"))
  c(rho = unname(ct$estimate), p = ct$p.value, n = length(x)) }
cat("\n=== does the advantage track power? reported at both levels ===\n")
vb <- gA[, .(m = mean(A)), by = .(V, cc)]
## proper sums-of-squares split; var() of the four cell means is NOT a
## decomposition (it ignores group sizes and can exceed the total)
gm  <- mean(gA$A)
ssb <- gA[, .N * (mean(A) - gm)^2, by = .(V, cc)][, sum(V1)]
sst <- sum((gA$A - gm)^2)
cat(sprintf("  power variance between cells: %.0f%%  within cells: %.0f%%  (SS split)\n",
            100 * ssb / sst, 100 * (1 - ssb / sst)))

cat("\n  BETWEEN cells (n = number of cells, the honest n for this component):\n")
cm <- nz[, .(power = median(A), rel = median(rel)), by = .(V, cc)][order(-power)]
print(cm[, .(V, cc, power = round(power, 3), rel_median = paste0(round(100 * rel, 1), "%"))])
if (nrow(cm) > 2) { r <- sp(cm$power, cm$rel)
  cat(sprintf("    Spearman on %d cell medians: rho %+.2f (p not meaningful at this n)\n", r["n"], r["rho"])) }

cat("\n  WITHIN cells (replication of the sign; each individually underpowered):\n")
wi <- nz[, as.list(sp(A, rel)), by = .(V, cc)]
print(wi[, .(V, cc, n, rho = round(rho, 3), p = round(p, 3))])
k <- sum(wi$rho < 0); n <- nrow(wi)
cat(sprintf("    %d of %d cells negative -- sign test p = %.3f\n", k, n,
            stats::binom.test(k, n, 0.5)$p.value))

## permutation keeping cell structure: shuffles rel WITHIN each cell, so the
## between-cell ordering cannot contribute. This is the defensible single
## number for the within-cell component.
set.seed(1)
obs <- suppressWarnings(stats::cor(nz$A, nz$rel, method = "spearman"))
perm <- replicate(2000, {
  z <- copy(nz)[, rel := rel[sample.int(.N)], by = .(V, cc)]
  suppressWarnings(stats::cor(z$A, z$rel, method = "spearman")) })
## The null preserves each cell's mean rel and mean A, so the BETWEEN-cell
## ordering survives shuffling and only the WITHIN-cell association is tested.
pp <- (1 + sum(perm <= obs)) / (1 + length(perm))
cat(sprintf("    within-cell permutation (2000x, shuffled inside cells): observed rho %+.3f, p %s\n",
            obs, if (pp <= 1/(1 + length(perm))) sprintf("< %.4f", 1/(1 + length(perm))) else sprintf("= %.3f", pp)))

cat("\n  sensitivity of the relative result to unstable denominators:\n")
for (thr in c(0.02, 0.05)) { h <- nz[A >= thr]
  if (nrow(h) > 5) { r <- sp(h$A, h$rel)
    cat(sprintf("    alpha AUC >= %.2f: rho %+.3f (n=%d, pooled -- see caveat above)\n",
                thr, r["rho"], r["n"])) } }
cat("  sensitivity to dropping c2 (the FIGURE's exclusion, not the test's):\n")
g2 <- gA[cc != 2]
if (nrow(g2) > 5) { r <- sp(g2$A, g2$gap)
  cat(sprintf("    absolute gap vs power, c2 excluded: rho %+.3f (n=%d); with c2: rho %+.3f (n=%d)\n",
              r["rho"], r["n"], sp(gA$A, gA$gap)["rho"], nrow(gA))) }
cat("\n  Read them as different questions: absolute is largest at intermediate\n")
cat("  power and vanishes where detection is impossible; relative rises as power\n")
cat("  falls, then plateaus. Neither is the other's correction.\n")
