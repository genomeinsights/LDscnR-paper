## module_sim/R/09_rmsc.R
## RMSC (rejection-maximising selection curve) for every dataset and both
## methods: over a grid of ld_w quantile thresholds, keep SNPs above each, BH-FDR
## their p-values, count discoveries. A clear interior PEAK (q* > 0) means the
## ld_w filter raises power (the local-LD prioritisation applies); a monotone-
## decreasing / flat curve (no_elbow, q* = 0) means filtering only loses power ->
## nothing to gain from ld_w on that dataset/method. Needs only ld_w + p-values
## (no genotypes), so it is cheap.
## Run from LDscnR-paper/:  Rscript module_sim/R/09_rmsc.R
## Output (git-ignored): figures/rmsc_all.png

source("module_sim/R/_config.R")               # SIM_DATA, dirs, LDscnR::rmsc_threshold
suppressMessages(library(ggplot2))
CONDS <- list(c("2", "1"), c("1", "2"))
ENV   <- 1:5
GRID  <- seq(0, 0.98, by = 0.01)

pool_ldw_p <- function(V, cc, env) {              # light: ld_w + p only, no genotypes
  files <- list.files(SIM_DATA, full.names = TRUE,
                      pattern = sprintf("^adapt_bgs_chr[0-9]+_V%s_c%s_env%s\\.rds$", V, cc, env))
  reps <- as.integer(sub(".*_chr([0-9]+)_.*", "\\1", basename(files))); files <- files[order(reps)]
  rbindlist(lapply(files, function(f) {
    d <- readRDS(f); m <- as.data.table(d$map)
    lwc <- if ("rho_0.95" %in% colnames(d$ld_ws)) "rho_0.95" else "0.95"
    data.table(ld_w = d$ld_ws[, lwc], emx_p = m$emx_p, lfmm_p = m$lfmm_p)
  }))
}

res <- rbindlist(lapply(CONDS, function(vc) rbindlist(lapply(ENV, function(e) {
  D <- pool_ldw_p(vc[1], vc[2], e)
  one <- function(pcol, meth) { s <- rmsc_threshold(D[[pcol]], D$ld_w, grid = GRID)
    data.table(cond = sprintf("V%s_c%s", vc[1], vc[2]), env = e, method = meth,
               q = s$grid, rej = s$rejections, qstar = s$q_star, no_elbow = s$no_elbow) }
  rbind(one("emx_p", "EMMAX"), one("lfmm_p", "LFMM"))
}))))

stars <- unique(res[, .(cond, env, method, qstar, no_elbow,
                        peak = max(rej)), by = .(cond, env, method)][
                        , .(cond, env, method, qstar, no_elbow)])
stars <- res[, .SD[which.max(rej)], by = .(cond, env, method)]     # point at the peak
MCOL <- c(EMMAX = "#E1AF00", LFMM = "#3B9AB2")
p <- ggplot(res, aes(q, rej, color = method)) +
  geom_line(linewidth = 0.5) +
  geom_point(data = stars, aes(qstar, rej), size = 2) +
  geom_vline(data = stars, aes(xintercept = qstar, color = method), linetype = 3, linewidth = 0.3) +
  facet_grid(cond ~ env, scales = "free_y") +
  scale_color_manual(values = MCOL, name = NULL) +
  labs(title = "RMSC per dataset: discoveries vs ld_w quantile threshold (point = q*)",
       subtitle = "interior peak (q*>0) => ld_w filter adds power; flat/monotone (q*=0) => nothing to gain",
       x = "ld_w quantile threshold q", y = "discoveries (BH-FDR < 0.05)") +
  theme_bw(base_size = 10) + theme(legend.position = "top")
ggsave(file.path(dir_fig, "rmsc_all.png"), p, width = 13, height = 6, dpi = 150)
saveRDS(res, file.path(dir_data, "rmsc_all.rds"))

## summary: q* and whether an elbow exists, per dataset/method
cat("=== RMSC q* (peak threshold) and no_elbow flag ===\n")
print(dcast(unique(res[, .(cond, env, method, qstar, no_elbow)]),
            cond + env ~ method, value.var = c("qstar", "no_elbow")))
cat("wrote figures/rmsc_all.png\n")
