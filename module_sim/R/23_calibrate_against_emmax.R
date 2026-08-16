## module_sim/R/23_calibrate_against_emmax.R
## CALIBRATE OTHER METHODS AGAINST EMMAX -- genomic-control on the C-score.
## We can only run a cheap structured null for EMMAX (fast EMMAX; and EMMAX's F is
## well-behaved, median ~ null). So make EMMAX the ANCHOR: its structured null gives a
## tau_C at a common FDR target. Bring any other method onto EMMAX's C scale by
## quantile-matching its genome-wide C distribution to EMMAX's (a monotone QQ map --
## the C distribution is dominated by null/structure behaviour, few true QTN, so this
## matches the nulls, exactly like GC). Then the other method's threshold is the C at
## the SAME upper-tail quantile EMMAX's anchor sits at. VALIDATE on the sim (truth):
## does the EMMAX-referenced LFMM tau_C land at LFMM's own PR-optimum?
## Run from LDscnR-paper/:  Rscript module_sim/R/23_calibrate_against_emmax.R [V c]

source("module_sim/R/21_estimate.R")
a  <- commandArgs(trailingOnly = TRUE)
V  <- if (length(a) >= 1) a[1] else "2"
cc <- if (length(a) >= 2) a[2] else "1"; ENV <- as.character(1:5)
ALPHA_SIM <- c(0.001, 0.01, 0.05, 0.1); LMIN <- 2L; TAUC <- seq(0.05, 0.95, by = 0.025)

## EMMAX anchor tau_C from the structured null (FDR<=0.05, pooled over env, alpha-swept)
sn <- readRDS(file.path(dir_data, sprintf("structured_null_summary_V%s_c%s.rds", V, cc)))
tau_anchor <- sn$tau05
cat(sprintf("EMMAX anchor tau_C (structured-null FDR<=0.05, pooled) = %.3f\n", tau_anchor))

CA <- lapply(ENV, function(e) load_cache(V, cc, e))
Cs <- lapply(CA, function(ca) list(
  emx  = cscore_count(ca$map$emx_p,  ca$LDW, ca$grid$RHO, ca$grid$QSTAR, ALPHA_SIM),
  lfmm = cscore_count(ca$map$lfmm_p, ca$LDW, ca$grid$RHO, ca$grid$QSTAR, ALPHA_SIM)))

## pooled genome-wide C per method -> stable ECDF for the QQ transfer
allE <- unlist(lapply(Cs, `[[`, "emx")); allL <- unlist(lapply(Cs, `[[`, "lfmm"))
q_anchor <- mean(allE < tau_anchor)                       # upper-tail position of the anchor
tau_lfmm <- as.numeric(quantile(allL, q_anchor))          # LFMM C at the SAME quantile = GC transfer
cat(sprintf("anchor quantile position = %.5f  ->  EMMAX-referenced LFMM tau_C = %.3f\n", q_anchor, tau_lfmm))

## truth-based PR at a given per-method tau, one env
pr_at <- function(ca, pcol, Cv, tau) { mk <- names(Cv)[Cv >= tau]
  reg <- if (length(mk)) cluster_from_cache(mk, ca$edges) else list()
  ev <- evaluate_ORs_qtn(reg[lengths(reg) >= LMIN], ca$map, ca$qtab, ca$th$r2min, ca$th$dmax)
  data.table(PR = if (is.na(ev$PR)) 0 else ev$PR, recall = ev$Recall,
             precision = if (is.na(ev$Precision)) NA_real_ else ev$Precision) }

## per-env: EMMAX at anchor; LFMM at transferred tau vs LFMM's OWN PR-optimum
res <- rbindlist(lapply(seq_along(CA), function(i) { ca <- CA[[i]]
  emx  <- pr_at(ca, "emx_p",  Cs[[i]]$emx,  tau_anchor)[, `:=`(env = i, method = "EMMAX@anchor",       tau = tau_anchor)]
  lfmT <- pr_at(ca, "lfmm_p", Cs[[i]]$lfmm, tau_lfmm)[,   `:=`(env = i, method = "LFMM@EMMAX-ref",     tau = tau_lfmm)]
  own  <- TAUC[which.max(vapply(TAUC, function(t) pr_at(ca, "lfmm_p", Cs[[i]]$lfmm, t)$PR, numeric(1)))]
  lfmO <- pr_at(ca, "lfmm_p", Cs[[i]]$lfmm, own)[,        `:=`(env = i, method = "LFMM@own-optimum",   tau = own)]
  rbind(emx, lfmT, lfmO) }))
se <- function(x){ x <- x[is.finite(x)]; if (length(x) < 2) NA_real_ else sd(x)/sqrt(length(x)) }
agg <- res[, .(tau = mean(tau), tau_se = se(tau), PR = mean(PR), PR_se = se(PR),
               recall = mean(recall), precision = mean(precision, na.rm = TRUE)), by = method]
cat("\n=== EMMAX-referenced (genomic-control) calibration, mean over 5 env, l_min=2 ===\n")
print(agg[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)])
d <- agg[method == "LFMM@EMMAX-ref", PR] - agg[method == "LFMM@own-optimum", PR]
cat(sprintf("\nLFMM: EMMAX-referenced PR - own-optimum PR = %.4f (near 0 => calibrating LFMM against EMMAX WORKS)\n", d))
saveRDS(list(agg = agg, per_env = res, tau_anchor = tau_anchor, tau_lfmm = tau_lfmm, q_anchor = q_anchor),
        file.path(dir_data, sprintf("calibrate_against_emmax_V%s_c%s.rds", V, cc)))
cat("wrote calibrate_against_emmax result\n")
