## module_sim/R/30_rhold_reconcile.R
## Pin ONE canonical clustering rho_ld: sweep rho_ld, rebuild the r^2-edge cache at each
## (fixed 500 kb cap for distance, TP-match r^2 held fixed), and measure sim PR-AUC
## (EMMAX + LFMM, l_min=2, mean env1-5). 3sp is robust to rho_ld (Chr1=1 inversion,
## Chr4 EMMAX=2 >10 at all values) so the sim decides. Pick the rho_ld that holds PR-AUC.
## Run: Rscript module_sim/R/30_rhold_reconcile.R [V c]

source("module_sim/R/21_estimate.R")
suppressMessages(library(ggplot2))
a  <- commandArgs(trailingOnly = TRUE)
V  <- if (length(a) >= 1) a[1] else "2"
cc <- if (length(a) >= 2) a[2] else "1"; ENV <- as.character(1:5)
RHO_LD_SET <- c(0.75, 0.85, 0.90, 0.95); ALPHA_SIM <- c(0.001, 0.01, 0.05, 0.1)
QSTAR <- seq(0, 0.95, by = 0.05); TAUC <- seq(0.05, 1, by = 0.05); LMIN <- 2L; DCAP <- 1e5

prauc <- function(Cv, edges, cache) { pc <- rbindlist(lapply(TAUC, function(t) {
    mk <- names(Cv)[Cv >= t]; reg <- if (length(mk)) cluster_from_cache(mk, edges) else list()
    ev <- evaluate_ORs_qtn(reg[lengths(reg) >= LMIN], cache$map, cache$qtab, cache$th$r2min, cache$th$dmax)
    data.table(recall = ev$Recall, precision = if (is.na(ev$Precision)) NA_real_ else ev$Precision) }))
  pr_auc(pc$recall, pc$precision) }

R <- rbindlist(lapply(ENV, function(env) { ca <- load_cache(V, cc, env)
  ## re-pool GTs to rebuild edges at each rho_ld (only genotype-touching step)
  files <- list.files(SIM_DATA, full.names = TRUE, pattern = sprintf("^adapt_bgs_chr[0-9]+_V%s_c%s_env%s\\.rds$", V, cc, env))
  reps <- as.integer(sub(".*_chr([0-9]+)_.*", "\\1", basename(files))); files <- files[order(reps)]
  gts <- lapply(seq_along(files), function(i) { G <- readRDS(files[i])$GTs; colnames(G) <- paste0("R", i, "_", colnames(G)); G })
  GTs <- do.call(cbind, gts)[, ca$map$marker]
  C_emx  <- cscore_count(ca$map$emx_p,  ca$LDW, ca$grid$RHO, QSTAR, ALPHA_SIM)
  C_lfmm <- cscore_count(ca$map$lfmm_p, ca$LDW, ca$grid$RHO, QSTAR, ALPHA_SIM)
  rbindlist(lapply(RHO_LD_SET, function(rl) {
    edges <- build_edge_cache(ca$universe, ca$map, GTs, ca$decay_sum, rho_ld = rl, dcap = DCAP)
    data.table(env = env, rho_ld = rl,
               emx  = prauc(setNames(C_emx,  rownames(ca$LDW)), edges, ca),
               lfmm = prauc(setNames(C_lfmm, rownames(ca$LDW)), edges, ca)) })) }))
se <- function(x){ x <- x[is.finite(x)]; if (length(x) < 2) NA_real_ else sd(x)/sqrt(length(x)) }
agg <- R[, .(EMMAX = mean(emx), EMMAX_se = se(emx), LFMM = mean(lfmm), LFMM_se = se(lfmm)), by = rho_ld][order(rho_ld)]
cat("=== sim PR-AUC vs clustering rho_ld (mean 5 env, l_min=2, 500kb cap) ===\n")
print(agg[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)])
cat(sprintf("\nbest rho_ld: EMMAX -> %.2f | LFMM -> %.2f\n", agg$rho_ld[which.max(agg$EMMAX)], agg$rho_ld[which.max(agg$LFMM)]))

saveRDS(list(R = R, agg = agg), file.path(dir_data, sprintf("rhold_reconcile_V%s_c%s.rds", V, cc)))
m <- melt(agg, id.vars = "rho_ld", measure.vars = c("EMMAX","LFMM"), variable.name = "method", value.name = "PR_AUC")
p <- ggplot(m, aes(rho_ld, PR_AUC, color = method)) + geom_line() + geom_point(size = 2) +
  scale_color_manual(values = c(EMMAX = "#1D3557", LFMM = "#E63946")) +
  labs(title = sprintf("V%s_c%s: sim PR-AUC vs clustering rho_ld (mean 5 env, l_min=2, 500kb cap)", V, cc),
       subtitle = "flat => rho_ld choice is not critical for the sim; pick to match the empirical r^2 regime", x = "rho_ld (r2-link)", y = "PR-AUC") +
  theme_bw(base_size = 11)
ggsave(file.path(dir_fig, sprintf("rhold_reconcile_V%s_c%s.png", V, cc)), p, width = 7, height = 5, dpi = 150)
cat("wrote rhold_reconcile figure\n")
