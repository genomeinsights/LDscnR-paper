## module_sticklebacks/19_or_overlap.R
## Overlap of LFMM outlier regions with EMMAX outlier regions (full run 16). Recompute
## each method's surviving >=10-SNP regions (EMMAX @ null tau_C; LFMM @ GC-mapped tau_C),
## take each region's genomic interval [min,max Pos] per chromosome, and test overlap
## (same Chr, intervals within PAD). Question: does LFMM (higher power) CONTAIN EMMAX's
## high-confidence set? Output: overlap counts + a genome-wide two-track segment figure.
## Run: Rscript module_sticklebacks/19_or_overlap.R

suppressMessages({ library(data.table); library(ggplot2) })
source("module_sim/R/_config.R")
mod <- "/Users/petrikem/gitlab/LDscnR-paper/module_sticklebacks"
R2LINK <- 0.1; DC <- 5e5; LMIN <- 10L; PAD <- 1e5   # 100 kb slack for "same locus"
x <- readRDS(file.path(mod, "sim_machinery_3sp.rds"))
C_emx <- x$C_emx; C_lfmm <- x$C_lfmm; tc <- x$tau_c; tc_lfmm <- x$tau_lfmm
sr <- readRDS(file.path(mod, "snp_stats_aligned.rds")); setDT(sr)
decs <- as.data.table(readRDS(file.path(mod, "decay_sum_3sp.rds")))
e <- new.env(); load("/Users/petrikem/gitlab/LDscnR-paper/LFMM_3sp/data/3sp_data.RData", envir = e)
GTs <- e$GTs_3sp; colnames(GTs) <- e$map_3sp$marker; GTs <- GTs[, sr$marker]
POS <- setNames(sr$Pos, sr$marker); CH <- setNames(as.character(sr$Chr), sr$marker)
uni <- unique(c(names(C_emx)[C_emx > 0], names(C_lfmm)[C_lfmm > 0]))
edges <- build_edge_cache(uni, sr[, .(marker, Chr, Pos)], GTs, decs, r2_link = R2LINK, dcap = DC); rm(GTs); gc()

regions <- function(Cv, tau, method) { mk <- names(Cv)[Cv >= tau]; reg <- cluster_from_cache(mk, edges)
  reg <- reg[lengths(reg) >= LMIN]
  rbindlist(lapply(seq_along(reg), function(i) { m <- reg[[i]]; p <- POS[m]
    data.table(method = method, id = i, Chr = names(sort(table(CH[m]), decreasing = TRUE))[1],
               start = min(p), end = max(p), size = length(m), maxC = max(Cv[m])) })) }
Remx <- regions(C_emx, tc, "EMMAX"); Rlf <- regions(C_lfmm, tc_lfmm, "LFMM")

## overlap: same Chr and intervals within PAD
ov <- function(a, B) B[Chr == a$Chr & start - PAD <= a$end & end + PAD >= a$start]
Remx[, lfmm_hits := sapply(seq_len(.N), function(i) nrow(ov(Remx[i], Rlf)))]
Rlf[,  emx_hit  := sapply(seq_len(.N), function(i) nrow(ov(Rlf[i],  Remx)) > 0)]
cat(sprintf("EMMAX ORs: %d | LFMM ORs: %d\n", nrow(Remx), nrow(Rlf)))
cat(sprintf("EMMAX ORs overlapped by an LFMM OR: %d / %d\n", Remx[lfmm_hits > 0, .N], nrow(Remx)))
cat(sprintf("LFMM ORs overlapping an EMMAX OR:   %d / %d  (LFMM-only: %d)\n", Rlf[emx_hit == TRUE, .N], nrow(Rlf), Rlf[emx_hit == FALSE, .N]))
cat("\n=== EMMAX ORs and their LFMM overlap ===\n")
print(Remx[order(Chr), .(Chr, start, end, size, maxC = round(maxC,2), lfmm_hits)])

saveRDS(list(Remx = Remx, Rlf = Rlf), file.path(mod, "or_overlap_3sp.rds"))
## two-track segment figure (cumulative genome x)
chr_lev <- paste0("Chr", c(1:18,20,21))
clen <- sr[, .(mx = max(Pos)), by = Chr][, Chr := factor(Chr, levels = chr_lev)][order(Chr)]
clen[, off := cumsum(shift(mx, fill = 0))]; xoff <- setNames(clen$off, as.character(clen$Chr))
seg <- rbind(Remx[, .(method, Chr, start, end, overlap = lfmm_hits > 0)],
             Rlf[,  .(method, Chr, start, end, overlap = emx_hit)])
seg[, `:=`(x1 = start + xoff[Chr], x2 = end + xoff[Chr], y = ifelse(method == "EMMAX", 2, 1))]
ax <- clen[, .(center = off + mx/2, Chr)]
p <- ggplot(seg) +
  geom_segment(aes(x = x1, xend = x2, y = y, yend = y, color = overlap), linewidth = 4) +
  geom_point(aes(x = (x1+x2)/2, y = y, color = overlap), size = 1.5) +
  scale_y_continuous(breaks = c(1,2), labels = c("LFMM", "EMMAX"), limits = c(0.5, 2.5)) +
  scale_x_continuous(breaks = ax$center, labels = sub("Chr","",ax$Chr), expand = c(0.01,0)) +
  scale_color_manual(values = c(`TRUE` = "#2A9D8F", `FALSE` = "#E63946"),
                     labels = c(`TRUE` = "shared (EMMAX∩LFMM)", `FALSE` = "method-only"), name = NULL) +
  labs(title = sprintf("3sp outlier regions: LFMM (%d) vs EMMAX (%d); shared = teal, method-only = red",
                       nrow(Rlf), nrow(Remx)),
       subtitle = sprintf("EMMAX @ tau_C=%.2f, LFMM @ GC-mapped tau_C=%.2f, l_min=10, overlap within %.0f kb",
                          tc, tc_lfmm, PAD/1e3), x = "Chromosome", y = NULL) +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), legend.position = "top")
ggsave(file.path(mod, "fig_or_overlap_3sp.png"), p, width = 13, height = 3.6, dpi = 150)
cat("\nwrote fig_or_overlap_3sp.png\n")
