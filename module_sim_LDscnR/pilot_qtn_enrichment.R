## =====================================================================
## module_sim_LDscnR / pilot_qtn_enrichment.R
##
## Fold enrichment for TAGGING a detectable QTN, as a function of how many
## markers you select. Three selectors on one axis, at equal budget:
##   ld_w      genotype-only; never sees the phenotype
##   C-score   ld_w integrated against the p-values
##   p-value   the association scan alone
##
## Equal budget is the point. Each rule has its own natural threshold and they
## select wildly different numbers of markers, so comparing them at their own
## operating points confounds the rule with the budget. Here k is held equal.
##
## THE C-SCORE LINE STOPS at k = #{C > 0}. C is exactly zero for >98% of
## markers, so "top k by C" is undefined beyond its support and extending the
## line would be inventing a ranking that does not exist. This is also why a
## genome-wide AUC of C is meaningless -- it is pinned near 0.5 by ties
## whatever the selection quality, which is a trap worth not falling into twice.
##
## Budgets start at 25. Below that a single marker moves enrichment by >10%
## and several datasets have top-25 sets containing no tagging marker at all,
## which is a real failure but rests on too few markers to plot honestly.
##
## Truth is flag_true_qtns(), the same detectable set the PR-AUC scoring uses;
## "tagging" means r2 >= score_thresholds(rho_r2 = 0.75)$r2min, per dataset.
##
## Read the ORDERING between selectors, which is consistent across the ten
## panels; do not read individual curve heights. n = 2 burn-in replicates,
## one environment, one cell.
##
## THE ENRICHMENT IS NOT GEOMETRY. Tested against a rotation null (ld_w shifted
## circularly along marker order, preserving its autocorrelation and the tagging
## pattern while breaking their alignment, 999 rotations): at the top 10%, 18 of
## 18 datasets significant at the rotation floor, null mean 1.00 against observed
## 1.40-6.23. Worth stating explicitly because a containment-style statistic over
## CLUSTERS does collapse against size-matched controls (99.8x -> 6.6x elsewhere
## in the project); a per-marker statistic cannot have that failure mode, and the
## spatial analogue was checked rather than assumed.
##
## THE SUB-1 CURVES IN pilot_26 / bgs ARE ONE DATASET, NOT A PATTERN. The
## commit message of b24e993 called them a reproducible blind spot on the
## strength of "both replicates, same direction". Hypergeometric test on the
## same comparison: rep2 is genuine depletion (0 observed against 7.3 expected,
## p = 5.3e-04, survives Bonferroni over all 16 datasets), rep1 is a coin flip
## (18 against 22.4, p = 0.18). The blind-spot INTERPRETATION -- a real cold spot
## need not contain a causal variant -- may well be right, but one significant
## dataset does not establish it. Note also that a fold of 0.00 from a 5.3% base
## rate expecting 7.3 markers is a degenerate estimator, the same class of error
## as scoring a mostly-zero C-score with a genome-wide AUC.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/pilot_qtn_enrichment.R
## Env: SIM_ROOT, REF_CELL (default V1_c1.5_env2), OUT
## =====================================================================
suppressMessages({ library(data.table); library(ggplot2); library(LDscnR) })
ROOT <- Sys.getenv("SIM_ROOT", "/Volumes/Nemo/Nemo_sim")
REF  <- Sys.getenv("REF_CELL", "V1_c1.5_env2")
OUT  <- Sys.getenv("OUT", "module_sim_LDscnR/figures")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
KS   <- unique(round(10^seq(log10(25), log10(5000), length.out = 34)))
QSTAR <- seq(0, 0.95, by = 0.05); ALPHA <- 0.05

SETS <- c(bgs5 = "regen_sim_data_bgs5", pilot_26 = "regen_sim_data_pilot_26",
          pilotC_26 = "regen_sim_data_pilotC_26", pilot_53 = "regen_sim_data_pilot_53",
          pilotC_53 = "regen_sim_data_pilotC_53")

curve_for <- function(f, nm) {
  x <- readRDS(f); m <- flag_true_qtns(as.data.table(x$map))
  if (!sum(m$true_pos_QTN %in% TRUE)) return(NULL)
  th  <- score_thresholds(as.data.table(x$LD_decay$decay_sum),
                          rho_r2 = 0.75, rho_d = 0.95, dmax_cap = 5e5)
  pos <- is.finite(m$max_LD_with_QTN) & m$max_LD_with_QTN >= th$r2min
  base <- mean(pos); if (!base) return(NULL)
  arm <- if (grepl("_nobgs_", basename(f))) "nobgs" else "bgs"
  out <- list()
  for (eng in c("emmax", "lfmm")) {
    p <- m[[if (eng == "emmax") "emx_p" else "lfmm_p"]]
    if (is.null(p) || all(is.na(p))) next
    C  <- ld_cscore(p, x$ld_ws, alpha = ALPHA, rho = colnames(x$ld_ws), qstar = QSTAR)
    Cv <- rep(0, nrow(m)); names(Cv) <- m$marker; Cv[names(C)] <- C
    nC <- sum(Cv > 0)
    ord <- list(`ld[w]` = order(-m$ld_w_095), `C-score` = order(-Cv),
                `p-value` = order(p, na.last = NA))
    for (sel in names(ord)) {
      ks <- KS[KS <= if (sel == "C-score") nC else length(ord[[sel]])]
      if (!length(ks)) next
      out[[paste(eng, sel)]] <- data.table(
        dataset = nm, arm = arm, engine = eng, selector = sel, k = ks,
        enrich = vapply(ks, function(k) mean(pos[ord[[sel]][seq_len(k)]]) / base, numeric(1)))
    }
  }
  rbindlist(out)
}

d <- rbindlist(lapply(names(SETS), function(nm) {
  ff <- list.files(file.path(ROOT, SETS[[nm]]), pattern = "[.]rds$", full.names = TRUE)
  if (nm == "bgs5") ff <- grep(sprintf("_%s[.]rds$", gsub("\\.", "[.]", REF)), ff, value = TRUE)
  rbindlist(lapply(ff, curve_for, nm = nm))
}), fill = TRUE)
if (!nrow(d)) stop("no curves built -- check SIM_ROOT and REF_CELL")

a <- d[, .(enrich = mean(enrich), n_rep = .N), by = .(dataset, arm, engine, selector, k)]
a[, dataset  := factor(dataset, levels = names(SETS))]
a[, selector := factor(selector, levels = c("p-value", "ld[w]", "C-score"))]

p <- ggplot(a, aes(k, enrich, colour = selector, linetype = arm)) +
  geom_hline(yintercept = 1, colour = "grey60", linewidth = .3) +
  geom_line(linewidth = .65) +
  facet_grid(engine ~ dataset) +
  scale_x_log10(breaks = c(25, 100, 1000), labels = c("25", "100", "1k")) +
  scale_y_log10() +
  scale_colour_manual(values = c(`p-value` = "#9E4630", `ld[w]` = "#6B7C84", `C-score` = "#1F6F8B"),
                      labels = c("p-value", "ld_w", "C-score"), name = NULL) +
  scale_linetype_manual(values = c(nobgs = "solid", bgs = "22"), name = NULL) +
  labs(x = "markers selected (log scale)", y = "fold enrichment for tagging a QTN",
       title = "How well does each rule pick markers that tag a causal variant?",
       subtitle = paste(sprintf("%s, two burn-in replicates averaged. Enrichment over each dataset's own base rate; 1 = chance.", REF),
                        "The C-score line ends where C stops being non-zero -- beyond that it has no ranking to offer.",
                        sep = "\n")) +
  theme_bw(base_size = 10) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "top", plot.subtitle = element_text(size = 8, colour = "grey35"))

f <- file.path(OUT, "pilot_qtn_enrichment.png")
ggsave(f, p, width = 13, height = 6, dpi = 150)
cat("  wrote", f, "\n\n")

## enrichment at the nearest budget to 100, averaged over arms and replicates
k100 <- a[, .SD[which.min(abs(k - 100))], by = .(dataset, arm, engine, selector)]
cat("=== fold enrichment at ~100 markers selected ===\n")
print(dcast(k100, dataset ~ selector, value.var = "enrich",
            fun.aggregate = function(z) round(mean(z), 1)))
