## module_sim/R/18e_null_fdr.R
## Extract the structured-null FDR(tau_C, l_min) from a null bundle (18d) -- a pure
## lookup, no genotypes, no EMMAX. Shows how the calibrated tau_C shifts with l_min
## (the region-size POST-filter): higher l_min removes scattered singletons -> lower
## null count -> lower tau_C for the same FDR. l_min counts REGIONS with >= l_min SNPs;
## the "per-SNP" curve (count gated SNPs) is shown as the l_min-free reference.
## Run: Rscript module_sim/R/18e_null_fdr.R [V c env]   (env may be "pooled": pool env1-5)

source("module_sim/R/_config.R")
suppressMessages(library(ggplot2))
a   <- commandArgs(trailingOnly = TRUE)
V   <- if (length(a) >= 1) a[1] else "2"
cc  <- if (length(a) >= 2) a[2] else "1"
env <- if (length(a) >= 3) a[3] else "1"
TAU <- seq(0.05, 1, by = 0.05); LMINS <- c(1L, 2L, 5L, 10L)

load_bundle <- function(e) readRDS(file.path(dir_data, sprintf("null_bundle_V%s_c%s_env%s.rds", V, cc, e)))
bundles <- if (env == "pooled") lapply(as.character(1:5), load_bundle) else list(load_bundle(env))

## count REGIONS with >= lmin SNPs among the C>=tau gate (per-SNP = count gated SNPs)
n_regions <- function(Csp, edges, tau, lmin) { mk <- names(Csp)[Csp >= tau]
  if (!length(mk)) return(0L); sum(lengths(cluster_from_cache(mk, edges)) >= lmin) }
n_snps <- function(Csp, tau) sum(Csp >= tau)

res <- rbindlist(lapply(c(as.list(LMINS), "perSNP"), function(lm) {
  rbindlist(lapply(TAU, function(t) {
    no <- nn <- 0
    for (bd in bundles) {
      if (identical(lm, "perSNP")) { no <- no + n_snps(bd$C_obs, t)
        nn <- nn + mean(vapply(bd$C_surr, n_snps, numeric(1), tau = t)) }
      else { no <- no + n_regions(bd$C_obs, bd$edges, t, lm)
        nn <- nn + mean(vapply(bd$C_surr, n_regions, numeric(1), edges = bd$edges, tau = t, lmin = lm)) } }
    data.table(l_min = if (identical(lm, "perSNP")) "perSNP" else as.character(lm),
               tau = t, n_obs = no, n_null = nn, FDR = pmin(1, nn / max(no, 1))) })) }))
res[, l_min := factor(l_min, levels = c("perSNP", as.character(LMINS)))]

tau_at <- function(d, q){ ok <- d[FDR <= q & n_obs > 0, tau]; if (length(ok)) min(ok) else NA_real_ }
cat(sprintf("=== structured-null tau_C by l_min (V%s_c%s_env%s) from ONE bundle ===\n", V, cc, env))
for (lv in levels(res$l_min)) { d <- res[l_min == lv]
  cat(sprintf("  l_min=%-6s  FDR<=0.10 -> tau_C=%.2f | FDR<=0.05 -> tau_C=%.2f\n", lv, tau_at(d,0.10), tau_at(d,0.05))) }

p <- ggplot(res, aes(tau, FDR, color = l_min)) + geom_line(linewidth = 0.8) +
  geom_hline(yintercept = c(0.05, 0.10), linetype = 2, color = "grey50") +
  scale_color_viridis_d(end = 0.9) + coord_cartesian(ylim = c(0, 1)) +
  labs(title = sprintf("V%s_c%s_env%s: structured-null FDR vs tau_C by l_min (one bundle)", V, cc, env),
       subtitle = "l_min applied only at extraction; higher l_min -> lower tau_C for same FDR", x = "tau_C", y = "empirical FDR") +
  theme_bw(base_size = 11)
ggsave(file.path(dir_fig, sprintf("null_fdr_by_lmin_V%s_c%s_env%s.png", V, cc, env)), p, width = 8, height = 5, dpi = 150)
saveRDS(res, file.path(dir_data, sprintf("null_fdr_by_lmin_V%s_c%s_env%s.rds", V, cc, env)))
cat("wrote null_fdr_by_lmin figure\n")
