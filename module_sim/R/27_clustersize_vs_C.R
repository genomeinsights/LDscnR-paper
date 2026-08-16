## module_sim/R/27_clustersize_vs_C.R
## Cluster SIZE vs C-score, colored by TP/FP (sim truth) -- visualizes the l_min
## mechanism directly: do high-C SNPs form small (FP) or large (TP) clusters?
## Cluster the C>0 SNPs (tau-independent gate); per cluster record size, max C, mean C,
## and TP/FP (diagnose_OR_classification). Pool env1-5, both methods.
## Run: Rscript module_sim/R/27_clustersize_vs_C.R [V c]

source("module_sim/R/21_estimate.R")
suppressMessages(library(ggplot2))
a  <- commandArgs(trailingOnly = TRUE)
V  <- if (length(a) >= 1) a[1] else "2"
cc <- if (length(a) >= 2) a[2] else "1"; ENV <- as.character(1:5)
ALPHA_SIM <- c(0.001, 0.01, 0.05, 0.1)

CL <- rbindlist(lapply(ENV, function(e) { ca <- load_cache(V, cc, e)
  rbindlist(lapply(c("emx_p", "lfmm_p"), function(p) {
    Cv <- cscore_count(ca$map[[p]], ca$LDW, ca$grid$RHO, ca$grid$QSTAR, ALPHA_SIM)
    gate <- names(Cv)[Cv > 0]; reg <- cluster_from_cache(gate, ca$edges)
    if (!length(reg)) return(data.table())
    dg <- as.data.table(diagnose_OR_classification(reg, ca$map, ca$qtab, ca$th$r2min, ca$th$dmax))
    data.table(env = e, method = sub("_p$", "", p),
               size = lengths(reg),
               maxC = vapply(reg, function(m) max(Cv[m]), numeric(1)),
               meanC = vapply(reg, function(m) mean(Cv[m]), numeric(1)),
               is_TP = dg$is_TP) })) }))
CL[, class := ifelse(is_TP, "TP", "FP")]
cat(sprintf("clusters pooled over %d env: %d (EMMAX %d, LFMM %d)\n", length(ENV), nrow(CL),
            CL[method=="emx", .N], CL[method=="lfmm", .N]))
## how TP fraction rises with size
cat("\n=== TP fraction by size bin (pooled) ===\n")
CL[, szbin := cut(size, c(0,1,2,4,9,19,Inf), labels=c("1","2","3-4","5-9","10-19","20+"))]
print(dcast(CL[, .(TP_frac = round(mean(is_TP),3), n = .N), by=.(method, szbin)], szbin ~ method, value.var="TP_frac"))

p <- ggplot(CL, aes(maxC, size, color = class)) +
  geom_hline(yintercept = c(2, 10), linetype = 3, color = "grey50") +
  geom_jitter(height = 0.06, width = 0.01, size = 1.1, alpha = 0.6) +
  facet_wrap(~method) +
  scale_y_log10() + scale_color_manual(values = c(TP = "#1D3557", FP = "#E63946")) +
  annotate("text", x = 0.05, y = 2.3, label = "l_min=2", size = 3, color = "grey40", hjust = 0) +
  annotate("text", x = 0.05, y = 11, label = "l_min=10", size = 3, color = "grey40", hjust = 0) +
  labs(title = sprintf("V%s_c%s: cluster size vs C-score, colored by TP/FP (sim truth, 5 env)", V, cc),
       subtitle = "each point = one C>0 cluster; high-C small clusters are FP, large clusters are TP -> l_min separates them",
       x = "cluster max C-score", y = "cluster size (SNPs, log10)") +
  theme_bw(base_size = 11)
ggsave(file.path(dir_fig, sprintf("clustersize_vs_C_V%s_c%s.png", V, cc)), p, width = 10, height = 5, dpi = 150)
saveRDS(CL, file.path(dir_data, sprintf("clustersize_vs_C_V%s_c%s.rds", V, cc)))
cat("\nwrote clustersize_vs_C figure\n")
