## =============================================================================
## final_analysis/R/03_emmax.R
##
## Primary EMMAX comparison (CLAUDE_REANALYSIS_INSTRUCTIONS.md): all three
## methods on the SAME MAF-filtered genotypes and the SAME Stage-1-derived
## GCTA relationship matrix from R/02_build_ld_units.R.
##
##   emmax_snp:       unrestricted marker-wise EMMAX, BH across ALL tested
##                     markers. Each significant marker is its own hypothesis.
##   emmax_simes:      marker p-values combined within phenotype-blind Stage-1
##                     units by Simes, BH across tested units.
##   emmax_consensus:  one consensus-dosage variable per eligible Stage-1
##                     unit, tested by EMMAX, BH across units.
##   emmax_snp_nonsingleton (diagnostic, NOT *_snp_clustered): the SAME
##     unrestricted marker p-values/BH decision as emmax_snp, but excludes
##     significant markers whose Stage-1 cluster size is 1 from the scored
##     set. Still scores every RETAINED marker as its own hypothesis (no
##     truth credit borrowed from neighbours) -- isolates "does dropping
##     singletons help" without emmax_snp_clustered's conflation of a marker
##     hit with discovery of its whole enclosing unit.
##
## Stage 2 (region assembly) is NOT run here -- instructions: "Stage 2 must
## not alter p-values, BH correction, or the rejection set... not needed for
## the primary simulation score." ld_outlier_test() is called with
## assembly="stage2_discovered" only because the package function requires
## an assembly argument; only $units (pre-assembly) is used downstream. A
## later, explicitly separate stage may add region-level reporting.
## =============================================================================
suppressMessages({library(data.table); library(LDscnR)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"), "R", "00_config.R"))

run_emmax <- function(tag, cell, rep, env, force = FALSE) {
  STAGE <- "03_emmax"   ## LOCAL -- see 02_build_ld_units.R's comment on this exact bug
  combo_id <- sprintf("%s_%s_rep%d_env%d", tag, cell, rep, env)
  ld_units_file <- file.path(stage_dir("02_build_ld_units", combo_id), "ld_units.rds")
  if (!file.exists(ld_units_file)) stop("R/02_build_ld_units.R has not produced: ", ld_units_file)

  INPUTS <- unname(ld_units_file)
  PARAMS <- list(alpha = ALPHA, size_floor = SIZE_FLOOR, unit_repr = UNIT_REPR)
  if (!force && !stage_stale(STAGE, INPUTS, PARAMS, target = combo_id)) {
    return(readRDS(file.path(stage_dir(STAGE, combo_id), "emmax.rds")))
  }

  say("=== %s: %s/%s/rep%d/env%d ===\n\n", STAGE, tag, cell, rep, env)
  b <- readRDS(ld_units_file)
  GTs <- b$GTs; map <- b$map; env_dt <- b$env; stage1 <- b$stage1; GRM <- b$GRM
  y <- env_dt$env
  stopifnot("GTs/map correspondence" = identical(colnames(GTs), map$marker),
            "GTs/env correspondence" = nrow(GTs) == nrow(env_dt))

  ## ---- per-marker EMMAX (basis for emmax_snp, emmax_snp_nonsingleton, emmax_simes) ----
  say("[1] per-marker EMMAX (%s markers)\n", format(ncol(GTs), big.mark = ","))
  Pm <- emmax_setup(GTs, GRM)
  pm_obs <- emmax_fast(Pm, y)
  stopifnot("per-marker p-vector must align to map" = length(pm_obs) == nrow(map))

  ## emmax_snp: unrestricted, BH over all tested markers
  q_snp <- stats::p.adjust(pm_obs, method = "BH")
  sig_snp <- !is.na(q_snp) & q_snp <= ALPHA
  bh_crit_snp <- if (any(sig_snp)) max(pm_obs[sig_snp]) else NA_real_
  say("    emmax_snp: %d/%d significant (BH alpha=%.2f)\n", sum(sig_snp), length(pm_obs), ALPHA)

  ## emmax_snp_nonsingleton: SAME p/q/BH decision, retained set excludes
  ## significant markers whose Stage-1 cluster size is 1.
  marker_cluster_size <- stats::setNames(stage1$map_snp$n_loci, stage1$map_snp$marker)[map$marker]
  sig_snp_nonsingleton <- sig_snp & !is.na(marker_cluster_size) & marker_cluster_size > 1
  say("    emmax_snp_nonsingleton: %d/%d significant (same BH decision, singletons excluded from the scored set)\n",
      sum(sig_snp_nonsingleton), length(pm_obs))

  ## ---- emmax_simes: marker p-values combined within Stage-1 units ------------
  say("\n[2] emmax_simes (Stage-1 units, Simes combination)\n")
  test_sim <- ld_outlier_test(stage1, map, pm_obs, statistic = "simes", size_floor = SIZE_FLOOR,
                              alpha = ALPHA, assembly = "stage2_discovered", GTs = GTs,
                              LD_decay = b$LD_decay, score_threshold = REGION_ASSEMBLY$score_threshold,
                              distance_threshold = REGION_ASSEMBLY$distance_threshold)
  say("    %d/%d tested units significant\n", sum(test_sim$units$significant), nrow(test_sim$units))

  ## ---- emmax_consensus: one EMMAX test per eligible Stage-1 unit -------------
  say("\n[3] emmax_consensus (Stage-1 units, %s)\n", UNIT_REPR)
  um <- ld_unit_matrix(GTs, stage1, map, size_floor = SIZE_FLOOR, repr = "consensus")
  Pu <- emmax_setup(um, GRM)
  pu_obs <- emmax_fast(Pu, y)
  test_con <- ld_outlier_test(stage1, map, pu_obs, statistic = "unit", size_floor = SIZE_FLOOR,
                              alpha = ALPHA, assembly = "stage2_discovered", GTs = GTs,
                              LD_decay = b$LD_decay, score_threshold = REGION_ASSEMBLY$score_threshold,
                              distance_threshold = REGION_ASSEMBLY$distance_threshold)
  say("    %d/%d tested units significant\n", sum(test_con$units$significant), nrow(test_con$units))

  out <- list(
    marker = list(Chr = map$Chr, Pos = map$Pos, marker = map$marker,
                  p = pm_obs, q = q_snp,
                  significant_snp = sig_snp, significant_snp_nonsingleton = sig_snp_nonsingleton,
                  bh_crit_p_snp = bh_crit_snp, n_tested = length(pm_obs)),
    emmax_simes = test_sim$units,
    emmax_consensus = test_con$units,
    settings = list(tag = tag, cell = cell, rep = rep, env = env, params = PARAMS))

  OUT_DIR <- stage_dir(STAGE, combo_id)
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  OUT <- file.path(OUT_DIR, "emmax.rds")
  saveRDS(out, OUT)
  write_receipt(STAGE, inputs = INPUTS, params = PARAMS, outputs = OUT, target = combo_id)
  say("\n[4] wrote %s\n", OUT)
  out
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 4) stop("Usage: Rscript R/03_emmax.R <tag> <cell> <rep> <env>")
  run_emmax(args[1], args[2], as.integer(args[3]), as.integer(args[4]))
}
