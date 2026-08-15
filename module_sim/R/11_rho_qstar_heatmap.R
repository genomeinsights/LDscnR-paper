## module_sim/R/11_rho_qstar_heatmap.R
## Reproduce the poster heatmap: number of outlier regions across the LD-filtering
## parameter grid (rho x q*).
##   rho = LD window size relative to LD decay -> which precomputed ld_w column
##         (ld_ws has rho_0.05 .. rho_0.99). Small rho = near-random filtering.
##   q*  = local-LD filtering threshold (ld_w quantile kept as candidates).
## For each (rho, q*): candidates = ld_w[rho] >= quantile(q*), BH-FDR within
## candidates < 0.05 -> significant SNPs, cluster them (r2>r2_link, distance-
## capped), count regions with >= L_MIN loci. LD-informed filtering (rho>0) should
## reveal more regions than near-random; a sweet spot in (rho, q*) maximises them.
## Run from LDscnR-paper/:  Rscript module_sim/R/11_rho_qstar_heatmap.R [V c env]
## Output (git-ignored): data/rhoqstar_V*.rds, figures/rhoqstar_V*.png

source("module_sim/R/_config.R")
suppressMessages(library(ggplot2))
a   <- commandArgs(trailingOnly = TRUE)
V   <- if (length(a) >= 1) a[1] else "2"
cc  <- if (length(a) >= 2) a[2] else "1"
env <- if (length(a) >= 3) a[3] else "1"
QSTAR  <- seq(0, 0.95, by = 0.05)
R2LINK <- 0.4                                   # poster connects r^2 > 0.4
FDR    <- 0.05

## pool GTs + map + the FULL ld_ws matrix (all rho columns)
files <- list.files(SIM_DATA, full.names = TRUE,
                    pattern = sprintf("^adapt_bgs_chr[0-9]+_V%s_c%s_env%s\\.rds$", V, cc, env))
reps <- as.integer(sub(".*_chr([0-9]+)_.*", "\\1", basename(files))); files <- files[order(reps)]
maps <- gts <- ldws <- decs <- vector("list", length(files))
for (i in seq_along(files)) {
  d <- readRDS(files[i]); m <- as.data.table(d$map)
  m[, `:=`(Chr = paste0("R", i, "_", Chr), marker = paste0("R", i, "_", marker))]
  G <- d$GTs; colnames(G) <- paste0("R", i, "_", colnames(G))
  lw <- d$ld_ws; rownames(lw) <- m$marker
  ds <- as.data.table(d$LD_decay$decay_sum); ds[, Chr := paste0("R", i, "_", Chr)]
  maps[[i]] <- m; gts[[i]] <- G; ldws[[i]] <- lw; decs[[i]] <- ds
}
map <- flag_true_positive_QTNs(rbindlist(maps, fill = TRUE))
GTs <- do.call(cbind, gts); LDW <- do.call(rbind, ldws)[map$marker, ]
decay_sum <- rbindlist(decs, fill = TRUE)
dcap <- d_from_rho(median(decay_sum$a_pred), 0.95)
RHO  <- colnames(LDW)
cat(sprintf("V%s_c%s_env%s: %d SNPs, %d chr | rho cols=%d, q* grid=%d\n",
    V, cc, env, nrow(map), uniqueN(map$Chr), length(RHO), length(QSTAR)))

heat <- function(pcol, meth) {
  pv <- map[[pcol]]
  rbindlist(lapply(RHO, function(rc) {
    lw <- LDW[, rc]
    rbindlist(lapply(QSTAR, function(q) {
      thr  <- stats::quantile(lw, q, na.rm = TRUE); cand <- which(lw >= thr)
      qv   <- stats::p.adjust(pv[cand], "BH"); sig <- map$marker[cand[qv < FDR]]
      sz   <- if (length(sig)) lengths(cluster_regions(sig, map, GTs, r2_link = R2LINK, dcap = dcap)) else integer(0)
      data.table(method = meth, rho = as.numeric(sub("rho_", "", rc)), qstar = q,
                 n_sig = length(sig), n_ge2 = sum(sz >= 2), n_ge5 = sum(sz >= 5), n_ge10 = sum(sz >= 10))
    }))
  }))
}
H <- rbind(heat("emx_p", "EMMAX"), heat("lfmm_p", "LFMM"))
saveRDS(H, file.path(dir_data, sprintf("rhoqstar_V%s_c%s_env%s.rds", V, cc, env)))

plot_one <- function(fill, lab) {
  ggplot(H, aes(factor(rho), factor(qstar), fill = .data[[fill]])) +
    geom_tile() + facet_wrap(~method) +
    scale_fill_viridis_c(option = "turbo", name = "outlier\nregions") +
    labs(x = expression(LD~window~size~relative~to~LD~decay~(rho)),
         y = "local LD filtering threshold (q*)",
         title = sprintf("V%s_c%s_env%s: outlier regions across rho x q*  (%s)", V, cc, env, lab)) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7), panel.grid = element_blank())
}
ggsave(file.path(dir_fig, sprintf("rhoqstar_V%s_c%s_env%s.png", V, cc, env)),
       plot_one("n_ge10", "regions with >=10 loci"), width = 12, height = 5, dpi = 150)
ggsave(file.path(dir_fig, sprintf("rhoqstar_ge2_V%s_c%s_env%s.png", V, cc, env)),
       plot_one("n_ge2", "regions with >=2 loci"), width = 12, height = 5, dpi = 150)
cat(sprintf("max regions >=10 loci: EMMAX=%d LFMM=%d | >=2 loci: EMMAX=%d LFMM=%d\n",
    H[method=="EMMAX", max(n_ge10)], H[method=="LFMM", max(n_ge10)],
    H[method=="EMMAX", max(n_ge2)], H[method=="LFMM", max(n_ge2)]))
cat("wrote rhoqstar figures\n")
