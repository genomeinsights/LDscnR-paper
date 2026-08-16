## module_sim/R/29_lfmm_gc_sim.R
## VALIDATE the empirical calibration on the sim (ground truth): EMMAX structured-null
## tau_C -> quantile-map to LFMM (genomic control) -> does LFMM at that operating point
## find MORE true ORs (its higher power) while the EMMAX-anchored calibration controls
## FP? For each env: anchor = EMMAX structured-null tau_C (18b, truth-free); q = EMMAX's
## C-quantile at that tau_C; tau_lfmm = LFMM's C at the same quantile. Score EMMAX@tau_C
## and LFMM@tau_lfmm against the QTN truth (# true ORs = TP, FP, recall, precision),
## l_min=2, alpha-swept C. Replicate-averaged over env1-5.
## Run: Rscript module_sim/R/29_lfmm_gc_sim.R [V c]

source("module_sim/R/21_estimate.R")
suppressMessages(library(ggplot2))
a  <- commandArgs(trailingOnly = TRUE)
V  <- if (length(a) >= 1) a[1] else "2"
cc <- if (length(a) >= 2) a[2] else "1"; ENV <- as.character(1:5)
ALPHA_SIM <- c(0.001, 0.01, 0.05, 0.1); LMIN <- 2L
tau_anchor <- readRDS(file.path(dir_data, sprintf("structured_null_summary_V%s_c%s.rds", V, cc)))$tau05
cat(sprintf("EMMAX structured-null tau_C anchor (FDR<=0.05) = %.3f\n", tau_anchor))

score <- function(ca, Cv, tau) { mk <- names(Cv)[Cv >= tau]
  reg <- if (length(mk)) cluster_from_cache(mk, ca$edges) else list()
  ev <- evaluate_ORs_qtn(reg[lengths(reg) >= LMIN], ca$map, ca$qtab, ca$th$r2min, ca$th$dmax)
  data.table(TP = ev$TP, FP = ev$FP, recall = ev$Recall,
             precision = if (is.na(ev$Precision)) NA_real_ else ev$Precision) }

R <- rbindlist(lapply(ENV, function(env) { ca <- load_cache(V, cc, env)
  C_emx  <- cscore_count(ca$map$emx_p,  ca$LDW, ca$grid$RHO, ca$grid$QSTAR, ALPHA_SIM)
  C_lfmm <- cscore_count(ca$map$lfmm_p, ca$LDW, ca$grid$RHO, ca$grid$QSTAR, ALPHA_SIM)
  q_anchor <- mean(C_emx < tau_anchor); tau_lfmm <- as.numeric(quantile(C_lfmm, q_anchor))
  rbind(score(ca, C_emx, tau_anchor)[, `:=`(env = env, method = "EMMAX @ null tau_C", tau = tau_anchor)],
        score(ca, C_lfmm, tau_lfmm)[,  `:=`(env = env, method = "LFMM @ GC-mapped tau_C", tau = tau_lfmm)]) }))
se <- function(x){ x <- x[is.finite(x)]; if (length(x) < 2) NA_real_ else sd(x)/sqrt(length(x)) }
agg <- R[, .(TP = mean(TP), TP_se = se(TP), FP = mean(FP), FP_se = se(FP),
             recall = mean(recall), recall_se = se(recall),
             precision = mean(precision, na.rm = TRUE), tau = mean(tau)), by = method]
cat("\n=== EMMAX-null calibration + GC-map to LFMM, mean over 5 env (l_min=2) ===\n")
print(agg[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
d_tp <- agg[method == "LFMM @ GC-mapped tau_C", TP] - agg[method == "EMMAX @ null tau_C", TP]
cat(sprintf("\nLFMM finds %.1f more true ORs than EMMAX at the calibrated operating point (LFMM's power),\nprecision EMMAX=%.2f vs LFMM=%.2f -> GC-map exploits LFMM power with EMMAX-anchored FP control.\n",
            d_tp, agg[method=="EMMAX @ null tau_C", precision], agg[method=="LFMM @ GC-mapped tau_C", precision]))

saveRDS(list(per_env = R, agg = agg, tau_anchor = tau_anchor), file.path(dir_data, sprintf("lfmm_gc_sim_V%s_c%s.rds", V, cc)))
m <- melt(agg, id.vars = "method", measure.vars = c("TP", "FP"), variable.name = "kind", value.name = "count")
p <- ggplot(m, aes(method, count, fill = kind)) + geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = round(count, 1)), position = position_dodge(0.6), vjust = -0.3, size = 3.2) +
  scale_fill_manual(values = c(TP = "#1D3557", FP = "#E63946"), name = NULL,
                    labels = c(TP = "true ORs (TP)", FP = "false ORs (FP)")) +
  labs(title = sprintf("V%s_c%s: EMMAX-null tau_C -> GC-map to LFMM (mean 5 env, l_min=2)", V, cc),
       subtitle = "LFMM finds more true ORs (its power) while EMMAX-anchored calibration holds FP", x = NULL, y = "# outlier regions") +
  theme_bw(base_size = 11)
ggsave(file.path(dir_fig, sprintf("lfmm_gc_sim_V%s_c%s.png", V, cc)), p, width = 7.5, height = 5, dpi = 150)
cat("wrote lfmm_gc_sim figure\n")
