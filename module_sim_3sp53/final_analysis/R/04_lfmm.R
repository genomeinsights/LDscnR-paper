## =============================================================================
## final_analysis/R/04_lfmm.R
##
## Secondary association-engine check (CLAUDE_REANALYSIS_INSTRUCTIONS.md):
## "Retain LFMM only as a portability analysis... unrestricted marker-wise
## LFMM2 plus BH; the same marker p-values combined within phenotype-blind
## Stage-1 units by Simes plus BH; no LFMM consensus-dosage analysis."
##
## Split into its own stage, deliberately independent of R/03_emmax.R (reads
## the SAME R/02_build_ld_units.R bundle, but never emmax.rds) so the
## primary EMMAX analysis can finish and be audited without rerunning LFMM,
## per instructions.
##
## K=5 is now a JUSTIFIED choice, not a legacy default -- see 00_config.R's
## LFMM_K comment: the 80 sampled populations fall into 5 discrete spatial
## groups (4 grid corners + 1 centre) BY DESIGN, verifiable from population
## coordinates alone, independent of genotypes or phenotype (PK, 2026-09-09).
##
## Call pattern reused verbatim from module_sim_3sp53/R/03_scan.R's own
## .lfmm_scan() (LEA::write.lfmm/write.env -> lfmm2(K) -> lfmm2.test) --
## reusable implementation material per the instructions' source-of-truth
## list, not re-derived. pv$pvalues from lfmm2.test(genomic.control=TRUE) is
## already GC-corrected; pv$gif is the PRE-correction inflation factor, kept
## as a diagnostic, not re-applied to the p-values that are actually used.
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(LEA)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))

.lfmm_scan <- function(geno_mat, y, K) {
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  gf <- file.path(tmp, "geno.lfmm"); ef <- file.path(tmp, "grad.env")
  write.lfmm(geno_mat, gf); write.env(y, ef)
  proj <- lfmm2(gf, ef, K = K)
  pv <- suppressWarnings(lfmm2.test(proj, gf, ef, genomic.control = TRUE, full = TRUE))
  list(p = pv$pvalues, gif_precorrection = pv$gif)
}

run_lfmm <- function(tag, cell, rep, env, force = FALSE) {
  STAGE <- "04_lfmm"   ## LOCAL -- see 02_build_ld_units.R's STAGE-clobbering comment
  combo_id <- sprintf("%s_%s_rep%d_env%d", tag, cell, rep, env)
  ld_units_file <- file.path(stage_dir("02_build_ld_units", combo_id), "ld_units.rds")
  if (!file.exists(ld_units_file)) stop("R/02_build_ld_units.R has not produced: ", ld_units_file)

  INPUTS <- unname(ld_units_file)
  PARAMS <- list(alpha = ALPHA, size_floor = SIZE_FLOOR, lfmm_k = LFMM_K)
  if (!force && !stage_stale(STAGE, INPUTS, PARAMS, target = combo_id)) {
    return(readRDS(file.path(stage_dir(STAGE, combo_id), "lfmm.rds")))
  }

  say("=== %s: %s/%s/rep%d/env%d ===\n\n", STAGE, tag, cell, rep, env)
  b <- readRDS(ld_units_file)
  GTs <- b$GTs; map <- b$map; env_dt <- b$env; stage1 <- b$stage1
  y <- env_dt$env
  stopifnot("GTs/map correspondence" = identical(colnames(GTs), map$marker))

  say("[1] LFMM2 (K=%d, %s markers)\n", LFMM_K, format(ncol(GTs), big.mark = ","))
  lf <- .lfmm_scan(GTs, y, LFMM_K)
  say("    pre-correction gif=%.3f (pvalues already GC-corrected)\n", lf$gif_precorrection)

  q_snp <- stats::p.adjust(lf$p, method = "BH")
  sig_snp <- !is.na(q_snp) & q_snp <= ALPHA
  bh_crit_snp <- if (any(sig_snp)) max(lf$p[sig_snp]) else NA_real_
  say("    lfmm_snp: %d/%d significant (BH alpha=%.2f)\n", sum(sig_snp), length(lf$p), ALPHA)

  say("\n[2] lfmm_simes (Stage-1 units, Simes combination)\n")
  test_sim <- ld_outlier_test(stage1, map, lf$p, statistic = "simes", size_floor = SIZE_FLOOR,
                              alpha = ALPHA, assembly = "stage2_discovered", GTs = GTs,
                              LD_decay = b$LD_decay, score_threshold = REGION_ASSEMBLY$score_threshold,
                              distance_threshold = REGION_ASSEMBLY$distance_threshold)
  say("    %d/%d tested units significant\n", sum(test_sim$units$significant), nrow(test_sim$units))

  out <- list(
    marker = list(Chr = map$Chr, Pos = map$Pos, marker = map$marker,
                  p = lf$p, q = q_snp, significant_snp = sig_snp,
                  bh_crit_p_snp = bh_crit_snp, n_tested = length(lf$p), gif_precorrection = lf$gif_precorrection),
    lfmm_simes = test_sim$units,
    settings = list(tag = tag, cell = cell, rep = rep, env = env, params = PARAMS))

  OUT_DIR <- stage_dir(STAGE, combo_id)
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  OUT <- file.path(OUT_DIR, "lfmm.rds")
  saveRDS(out, OUT)
  write_receipt(STAGE, inputs = INPUTS, params = PARAMS, outputs = OUT, target = combo_id)
  say("\n[3] wrote %s\n", OUT)
  out
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 4) stop("Usage: Rscript R/04_lfmm.R <tag> <cell> <rep> <env>")
  invisible(run_lfmm(args[1], args[2], as.integer(args[3]), as.integer(args[4])))
}
