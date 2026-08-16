## module_sticklebacks/13_clustersize_vs_C_3sp.R
## Cluster SIZE vs C-score for the empirical 3sp data (EMMAX, LFMM) -- the empirical
## analog of module_sim/27. No truth to colour by, so Eda (Chr4) is highlighted and the
## l_min=10 line shown. Expectation vs the sim: 3sp clusters reach LARGE sizes at high C
## (dense SNPs, real loci span many SNPs) -> l_min=10 is the appropriate cut here.
## Precomputed emx_p/lfmm_p (poster engine), C at alpha=0.05. Cluster C>0 SNPs with the
## 3sp convention rho_ld=0.95/rho_d=0.99<=500kb.
## Run: Rscript module_sticklebacks/13_clustersize_vs_C_3sp.R

suppressMessages({ library(data.table); library(ggplot2) })
source("module_sim/R/_config.R")                       # build_edge_cache, cluster_from_cache
mod <- "/Users/petrikem/gitlab/LDscnR-paper/module_sticklebacks"
R2LINK <- 0.1; DC <- 5e5                                # direct r^2 + 0.5 Mb split (keeps Chr1 inversion whole, resolves Chr4 loci)
o  <- readRDS(file.path(mod, "cscore_obs.rds")); sr <- as.data.table(o$sr)
decs <- as.data.table(readRDS(file.path(mod, "decay_sum_3sp.rds")))
e <- new.env(); load("/Users/petrikem/gitlab/LDscnR-paper/LFMM_3sp/data/3sp_data.RData", envir = e)
GTs <- e$GTs_3sp; colnames(GTs) <- e$map_3sp$marker; GTs <- GTs[, sr$marker]
uni <- sr[C_emx > 0 | C_lfmm > 0, marker]
edges <- build_edge_cache(uni, sr[, .(marker, Chr, Pos)], GTs, decs, r2_link = R2LINK, dcap = DC); rm(GTs); gc()

CE <- setNames(sr$C_emx, sr$marker); CL <- setNames(sr$C_lfmm, sr$marker); CH <- setNames(sr$Chr, sr$marker)
cl_tab <- function(m) { Cv <- if (m == "EMMAX") CE else CL
  reg <- cluster_from_cache(names(Cv)[Cv > 0], edges); if (!length(reg)) return(data.table())
  data.table(method = m, size = lengths(reg),
             maxC = vapply(reg, function(mk) max(Cv[mk]), numeric(1)),
             Chr  = vapply(reg, function(mk) names(sort(table(CH[mk]), decreasing = TRUE))[1], character(1))) }
CLU <- rbindlist(lapply(c("EMMAX", "LFMM"), cl_tab))
CLU[, grp := ifelse(Chr == "Chr1", "Chr1 (inversion)", ifelse(Chr == "Chr4", "Chr4 (Eda)", "other"))]
CLU[, grp := factor(grp, levels = c("other", "Chr1 (inversion)", "Chr4 (Eda)"))]
CLU[, method := factor(method, levels = c("EMMAX", "LFMM"))]

cat(sprintf("clustering: r2>=%.2f, %.1f Mb split | 3sp C>0 clusters: EMMAX %d, LFMM %d\n",
            R2LINK, DC/1e6, CLU[method=="EMMAX",.N], CLU[method=="LFMM",.N]))
for (m in c("EMMAX","LFMM")) for (ch in c("Chr1","Chr4")) {
  d <- CLU[method==m & Chr==ch][order(-size)]
  cat(sprintf("[%s] %s: %d clusters, %d with >10 SNPs (sizes: %s)\n", m, ch, nrow(d), sum(d$size>10),
              paste(head(d$size,6), collapse=","))) }

p <- ggplot(CLU, aes(maxC, size, color = grp, size = grp)) +
  geom_hline(yintercept = c(2, 10), linetype = 3, color = "grey50") +
  geom_jitter(data = CLU[grp == "other"], height = 0.05, width = 0.005, alpha = 0.35) +
  geom_jitter(data = CLU[grp != "other"], height = 0.05, width = 0.005, alpha = 0.9) +
  facet_wrap(~method) + scale_y_log10() +
  scale_color_manual(values = c("other" = "grey60", "Chr1 (inversion)" = "#1D3557", "Chr4 (Eda)" = "#E63946"), name = NULL) +
  scale_size_manual(values = c("other" = 1, "Chr1 (inversion)" = 2.2, "Chr4 (Eda)" = 2.2), guide = "none") +
  annotate("text", x = 0.02, y = 11, label = "l_min=10", size = 3, color = "grey40", hjust = 0) +
  labs(title = sprintf("3sp: cluster size vs C-score (r2>=%.1f, %.1f Mb split); Chr1 vs Chr4 highlighted", R2LINK, DC/1e6),
       subtitle = "each point = one C>0 cluster; Chr1 = one big inversion cluster (top); Chr4 = several loci",
       x = "cluster max C-score", y = "cluster size (SNPs, log10)") +
  theme_bw(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(mod, "fig_clustersize_vs_C_3sp.png"), p, width = 10, height = 5.5, dpi = 150)
cat("\nwrote fig_clustersize_vs_C_3sp.png\n")
