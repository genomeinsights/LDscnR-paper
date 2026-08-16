## module_sticklebacks/12_manhattan_lmin10.R
## C-score Manhattan at l_min=10 (poster setting): cluster the C>0 SNPs and keep only
## >=10-SNP regions, so only large robust clusters remain (vs the per-SNP view in 11,
## which shows every C>0 SNP and looks nothing like the poster). Precomputed emx_p /
## lfmm_p (the poster/caller engine). Distance-restricted clustering rho_ld=0.95 /
## rho_d=0.99 <=500kb (the 3sp convention). Eda = Chr4.
## Run from LDscnR-paper/:  Rscript module_sticklebacks/12_manhattan_lmin10.R [l_min]
## Output: module_sticklebacks/fig_cscore_manhattan_3sp_lmin{L}.png

suppressMessages({ library(data.table); library(ggplot2) })
source("module_sim/R/_config.R")                       # build_edge_cache, cluster_from_cache
mod <- "/Users/petrikem/gitlab/LDscnR-paper/module_sticklebacks"
a <- commandArgs(trailingOnly = TRUE); LMIN <- if (length(a) >= 1) as.integer(a[1]) else 10L
R2LINK <- 0.1; DC <- 5e5      # direct r^2 + 0.5 Mb split: Chr1 inversion stays one cluster, Chr4 loci resolve

o  <- readRDS(file.path(mod, "cscore_obs.rds")); sr <- as.data.table(o$sr)
decs <- as.data.table(readRDS(file.path(mod, "decay_sum_3sp.rds")))
e <- new.env(); load("/Users/petrikem/gitlab/LDscnR-paper/LFMM_3sp/data/3sp_data.RData", envir = e)
GTs <- e$GTs_3sp; colnames(GTs) <- e$map_3sp$marker; GTs <- GTs[, sr$marker]
universe <- sr[C_emx > 0 | C_lfmm > 0, marker]
edges <- build_edge_cache(universe, sr[, .(marker, Chr, Pos)], GTs, decs, r2_link = R2LINK, dcap = DC)
rm(GTs); gc()

## per method: cluster C>0 SNPs, keep >=LMIN-SNP regions, build a member table explicitly
CH <- setNames(sr$Chr, sr$marker); PS <- setNames(sr$Pos, sr$marker)
CE <- setNames(sr$C_emx, sr$marker); CL <- setNames(sr$C_lfmm, sr$marker)
build_mem <- function(m) {
  Cvec <- if (m == "EMMAX") CE else CL
  reg <- cluster_from_cache(names(Cvec)[Cvec > 0], edges)
  reg <- reg[lengths(reg) >= LMIN]
  if (!length(reg)) return(data.table(method = character()))
  rbindlist(lapply(seq_along(reg), function(i) {
    mk <- reg[[i]]
    data.table(method = m, region = i, marker = mk, Chr = CH[mk], Pos = PS[mk], C = Cvec[mk]) }))
}
mem <- rbindlist(list(build_mem("EMMAX"), build_mem("LFMM")), fill = TRUE)

## genome coordinates
chr_lev <- paste0("Chr", c(1:18, 20, 21))
clen <- sr[, .(mx = max(Pos)), by = Chr][, Chr := factor(Chr, levels = chr_lev)][order(Chr)]
clen[, off := cumsum(shift(mx, fill = 0))]; xoff <- setNames(clen$off, as.character(clen$Chr))
sr[, gx := Pos + xoff[Chr]]; mem[, gx := Pos + xoff[Chr]]
mem[, `:=`(Chr = factor(Chr, levels = chr_lev), method = factor(method, levels = c("EMMAX", "LFMM")))]
ax <- clen[, .(center = off + mx / 2, Chr)]

cat(sprintf("l_min=%d clusters: EMMAX %d regions, LFMM %d regions\n",
            LMIN, mem[method == "EMMAX", uniqueN(region)], mem[method == "LFMM", uniqueN(region)]))
for (m in c("EMMAX", "LFMM")) { cat(sprintf("\n[%s] >=%d-SNP C-clusters by chromosome:\n", m, LMIN))
  print(mem[method == m, .(regions = uniqueN(region), SNPs = .N, maxC = round(max(C), 3)), by = Chr][order(-regions)]) }

p <- ggplot() +
  geom_point(data = sr, aes(gx, pmax(C_emx, C_lfmm)), color = "grey85", size = 0.3) +
  geom_point(data = mem, aes(gx, C, color = factor(region %% 12)), size = 0.9) +
  facet_wrap(~method, ncol = 1) +
  scale_color_viridis_d(guide = "none") +
  scale_x_continuous(breaks = ax$center, labels = sub("Chr", "", ax$Chr), expand = c(0.01, 0)) +
  labs(title = sprintf("3sp C-score Manhattan at l_min=%d: only >=%d-SNP clusters coloured (r2>=%.1f, %.1f Mb split)", LMIN, LMIN, R2LINK, DC/1e6),
       subtitle = "grey = all SNPs (max C over methods); coloured = members of clusters passing l_min (one colour per cluster)",
       x = "Chromosome", y = "C-score") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
ggsave(file.path(mod, sprintf("fig_cscore_manhattan_3sp_lmin%d.png", LMIN)), p, width = 12, height = 6.5, dpi = 150)
cat(sprintf("\nwrote fig_cscore_manhattan_3sp_lmin%d.png\n", LMIN))
