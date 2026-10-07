## =============================================================================
## final_analysis/R/25_tier_vs_region_size.R
##
## Independent audit (2026-10-07): does the floor-stability tier predict QTN
## overlap beyond region span and marker count? 3/3 regions are much larger
## (median span ~118 kb vs ~6 kb at 1/3), so the raw precision gradient could
## be a size/coverage effect. Logistic model TP ~ tier + log1p(span_kb) +
## log(n_markers) + method + cell + BGS tag, likelihood-ratio test for tier,
## rep-clustered bootstrap (1,000 draws over the 10 replicate maps) of the tier
## odds ratios, and precision by tier within span tertiles.
## Reads results/floor_stability_density_half_double_region_audit.tsv (R/23).
## =============================================================================
suppressMessages(library(data.table))
a <- fread(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis/results/floor_stability_density_half_double_region_audit.tsv"))
a[, `:=`(span_kb = (to - from) / 1e3, tier = factor(tier_effective), rep = as.integer(sub(".*_rep(\\d+)_.*", "\\1", combo_id)),
         cell = sub("^[a-z]+_(V[0-9.]+_c[0-9.]+)_.*", "\\1", combo_id), tag = sub("_.*", "", combo_id))]
size_by_tier <- a[, .(n = .N, median_span_kb = median(span_kb), median_markers = as.numeric(median(n_markers)), precision = mean(TP)), by = .(method, tier)][order(method, tier)]; print(size_by_tier)
m0 <- glm(TP ~ log1p(span_kb) + log(n_markers) + method + cell + tag, binomial, a)
m1 <- update(m0, . ~ . + tier)
print(summary(m1)$coefficients[c("tier2", "tier3", "log1p(span_kb)", "log(n_markers)"), ])
cat("LRT tier | span+markers: p =", anova(m0, m1, test = "LRT")$`Pr(>Chi)`[2], "\n")
## within span tertiles
a[, span_bin := cut(span_kb, quantile(span_kb, 0:3 / 3), include.lowest = TRUE, labels = c("short", "mid", "long"))]
print(dcast(a[, .(p = sprintf("%.2f (%d)", mean(TP), .N)), by = .(span_bin, tier)], span_bin ~ tier, value.var = "p"))
## rep-clustered bootstrap of tier ORs
set.seed(1); reps <- sort(unique(a$rep))
bt <- t(replicate(1000, { r <- sample(reps, replace = TRUE); x <- rbindlist(lapply(r, function(k) a[rep == k]))
  coef(glm(TP ~ log1p(span_kb) + log(n_markers) + method + cell + tag + tier, binomial, x))[c("tier2", "tier3")] }))
print(round(exp(rbind(est = coef(m1)[c("tier2", "tier3")], apply(bt, 2, quantile, c(.025, .975)))), 2))
or <- exp(rbind(est = coef(m1)[c("tier2", "tier3")], apply(bt, 2, quantile, c(.025, .975))))
R <- function(x) file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "results", x)
fwrite(size_by_tier, R("floor_stability_tier_size_by_tier.tsv"), sep = "\t")
fwrite(data.table(term = c("tier2", "tier3"), odds_ratio = or["est", ], ci_lo = or["2.5%", ], ci_hi = or["97.5%", ],
                  lrt_p_tier = anova(m0, m1, test = "LRT")$`Pr(>Chi)`[2]),
       R("floor_stability_tier_size_model.tsv"), sep = "\t")
fwrite(a[, .(n = .N, precision = mean(TP)), by = .(span_bin, tier)][order(span_bin, tier)],
       R("floor_stability_tier_by_span_tertile.tsv"), sep = "\t")
