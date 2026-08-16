## module_sticklebacks/09_cscore.R
## OBSERVED C-score (consistency across the rho x q* x alpha grid) on the empirical
## 3sp data, EMMAX and LFMM. This is the module_sim C-score machinery (13/14) applied
## genome-wide to real data: for each SNP, C = fraction of (rho, q*, alpha) cells where
## it is a candidate (ld_w >= quantile(q*)) AND BH-FDR<alpha among candidates. No truth
## here (empirical) -> the deliverable is the C landscape + which known loci (Eda/Chr4)
## light up, and n_obs vs tau_C to hand to the structured-null calibration (10).
## Run from LDscnR-paper/:  Rscript module_sticklebacks/09_cscore.R
## Output (git-ignored): module_sticklebacks/cscore_obs.rds

suppressMessages({ library(data.table) })
mod <- "/Users/petrikem/gitlab/LDscnR-paper/module_sticklebacks"
sr  <- readRDS(file.path(mod, "snp_stats_aligned.rds")); setDT(sr)
LDW <- readRDS("/Users/petrikem/gitlab/LDscnR-paper/3sp_data/ld_ws_3sp_MAF01.rds")[sr$marker, ]
RHO <- colnames(LDW); QSTAR <- seq(0, 0.95, by = 0.05); ALPHA_C <- 0.05   # fixed alpha (see _config cscore_count)
ncell <- length(RHO) * length(QSTAR) * length(ALPHA_C); TAU <- seq(0.02, 1, by = 0.02)

## per-marker consistency counter (same fast counter as module_sim/18)
Cscore <- function(pv) { cnt <- integer(nrow(sr))
  for (rc in RHO) { lw <- LDW[, rc]
    for (q in QSTAR) { thr <- stats::quantile(lw, q, na.rm = TRUE); cand <- which(lw >= thr)
      if (!length(cand)) next; qv <- stats::p.adjust(pv[cand], "BH")
      for (al in ALPHA_C) { h <- cand[qv < al]; if (length(h)) cnt[h] <- cnt[h] + 1L } } }
  cnt / ncell }

cat(sprintf("3sp: %d SNPs, %d chr; grid %d rho x %d q* x %d alpha = %d cells\n",
            nrow(sr), uniqueN(sr$Chr), length(RHO), length(QSTAR), length(ALPHA_C), ncell))
C_emx  <- Cscore(sr$emx_p)
C_lfmm <- Cscore(sr$lfmm_p)
sr[, `:=`(C_emx = C_emx, C_lfmm = C_lfmm)]

## n high-C SNPs vs tau_C, both methods
tab <- data.table(tau = TAU,
                  n_emx  = vapply(TAU, function(t) sum(C_emx  >= t), numeric(1)),
                  n_lfmm = vapply(TAU, function(t) sum(C_lfmm >= t), numeric(1)))
cat("\n=== n high-C SNPs vs tau_C ===\n"); print(tab[tau %in% seq(0.1, 1, 0.1)])

## which chromosomes carry the high-C signal (tau_C >= 0.5); Eda = Chr4
for (m in c("emx", "lfmm")) { Cc <- if (m=="emx") C_emx else C_lfmm
  hi <- sr[Cc >= 0.5, .N, by = Chr][order(-N)]
  cat(sprintf("\n[%s] chromosomes with C>=0.5 (n SNPs):\n", toupper(m))); print(head(hi, 8)) }

## Eda/Chr4 focus: max C and how many high-C SNPs there
cat("\n=== Chr4 (Eda region) ===\n")
cat(sprintf("EMMAX : max C = %.3f, n(C>=0.5) = %d\n", max(sr[Chr=="Chr4", C_emx]),  sr[Chr=="Chr4" & C_emx  >= 0.5, .N]))
cat(sprintf("LFMM  : max C = %.3f, n(C>=0.5) = %d\n", max(sr[Chr=="Chr4", C_lfmm]), sr[Chr=="Chr4" & C_lfmm >= 0.5, .N]))

saveRDS(list(sr = sr[, .(marker, Chr, Pos, maf, ld_w, emx_p, lfmm_p, C_emx, C_lfmm)],
             tab = tab, RHO = RHO, QSTAR = QSTAR, ALPHA_C = ALPHA_C, ncell = ncell),
        file.path(mod, "cscore_obs.rds"))
cat("\nwrote cscore_obs.rds\n")
