## module_sim/08_sweep_aggregate.R
## Aggregate the legacy-vs-ours comparison across env replicates for one or more
## (V,c) conditions. Reads the consensus_V*_c*_env*.rds written by 06b (run 06a+
## 06b for every env first, e.g. via the loop in the header of this file), then
## reports mean +/- SE of Precision / Recall / PR / F1 per (condition, filter,
## method) and draws a precision-recall summary.
##
## Compute first (from LDscnR-paper/):
##   for cond in "2 1" "1 2"; do set -- $cond
##     for e in 1 2 3 4 5; do
##       [ -f module_sim/cache_V${1}_c${2}_env${e}.rds ] || Rscript module_sim/06a_run_caller.R $1 $2 $e
##       Rscript module_sim/06b_score.R $1 $2 $e 2
##     done; done
## Then:  Rscript module_sim/08_sweep_aggregate.R
## Output (git-ignored): module_sim/sweep_summary.rds, module_sim/sweep_PR.png

suppressMessages({ library(data.table); library(ggplot2) })
mod   <- "/Users/petrikem/gitlab/LDscnR-paper/module_sim"
dir_data <- file.path(mod, "data"); dir_fig <- file.path(mod, "figures")
for (d in c(dir_data, dir_fig)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
CONDS <- list(c("2", "1"), c("1", "2"))          # (V,c): moderate, hard-structure
ENV   <- 1:5

all <- rbindlist(lapply(CONDS, function(vc) {
  rbindlist(lapply(ENV, function(e) {
    f <- file.path(dir_data, sprintf("consensus_V%s_c%s_env%s.rds", vc[1], vc[2], e))
    if (!file.exists(f)) return(NULL)
    C <- readRDS(f)
    rbindlist(list(C$res, C$res_snpN, C$res_multi))[
      , `:=`(cond = sprintf("V%s_c%s", vc[1], vc[2]), env = e)]
  }), fill = TRUE)
}), fill = TRUE)
if (!nrow(all)) stop("no consensus_*.rds found -- run 06a+06b for the swept conditions first")

se <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) NA_real_ else sd(x) / sqrt(length(x)) }
agg <- all[, .(nrep = .N,
               Precision = mean(Precision, na.rm = TRUE), Precision_se = se(Precision),
               Recall    = mean(Recall,    na.rm = TRUE), Recall_se    = se(Recall),
               PR = mean(PR, na.rm = TRUE), PR_se = se(PR),
               F1 = mean(F1, na.rm = TRUE), F1_se = se(F1)),
           by = .(cond, filter, method)]
setorder(agg, cond, filter, -F1)
cat(sprintf("=== sweep: %d conditions x %d env (rows = mean over env) ===\n",
            length(CONDS), length(ENV)))
print(agg[, .(cond, filter, method, nrep,
              Precision = round(Precision, 3), Recall = round(Recall, 3),
              PR = round(PR, 3), F1 = round(F1, 3))])

## precision-recall summary: facet = condition x filter; point per method +/- SE
lvl <- c("all", ">=2 SNPs", ">=2 methods")
agg[, filter := factor(filter, levels = lvl)]
MCOL <- c("EMMAX single" = "#8C8C8C", "LFMM single" = "#3B9AB2",
          "EMMAX ld_w" = "#E1AF00", "LFMM ld_w" = "#F21A00")
p <- ggplot(agg, aes(Recall, Precision, color = method)) +
  geom_errorbar(aes(ymin = Precision - Precision_se, ymax = Precision + Precision_se), width = 0, alpha = 0.5) +
  geom_errorbarh(aes(xmin = Recall - Recall_se, xmax = Recall + Recall_se), height = 0, alpha = 0.5) +
  geom_point(size = 3) +
  facet_grid(cond ~ filter) +
  scale_color_manual(values = MCOL, name = NULL) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = sprintf("Legacy single-SNP vs ld_w+null (mean +/- SE over %d env)", length(ENV)),
       subtitle = "rows = condition; cols = region filter; up-and-right = better",
       x = "Recall", y = "Precision") +
  theme_bw(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(dir_fig, "sweep_PR.png"), p, width = 11, height = 7, dpi = 150)
saveRDS(list(agg = agg, all = all, conds = CONDS, env = ENV), file.path(dir_data, "sweep_summary.rds"))
cat("wrote sweep_PR.png + sweep_summary.rds\n")
