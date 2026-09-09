## module_sim_3sp53/R/14_random_removal_control.R
##
## PK, 2026-09-09, on the size-floor sweep (R/12_structured_null.R): "if we
## randomly remove the same number of outlier clusters (but independently
## of size), would we also see a decline in FPs, and is the decline in
## fig_structured_null_sizesweep relative to that?"
##
## A matched-cardinality random-subset control for the `group`-null,
## emmax_consensus OR emmax_simes arm (pass "simes" as the 3rd arg -- see
## the 2026-09-09 update below). At each SIZE_FLOOR_GRID floor f: instead of
## restricting the candidate pool to Stage-1 units with n_markers>=f (what
## R/12 does), draw N_RANDOM random subsets of THE SAME CARDINALITY from
## the full candidate pool, ignoring size, and compute realised_fdr the
## same way for each. If size-based restriction is doing something SIZE-
## SPECIFIC (concentrating true signal), it should beat this baseline. If
## it's merely benefiting from a smaller candidate pool the way ANY
## same-sized restriction would, the two should track each other.
##
## Requires full per-unit significance tracking across surrogate draws (a
## B x N_units logical matrix), not just sizes-of-significant-units the
## way R/12's .perm_sizes() keeps them -- units_base's unit_id is stable
## across every draw (built once, only p/q/significant change per draw),
## so re-slicing by an arbitrary random subset of unit_id afterwards is a
## legitimate, cheap re-tabulation (no refitting).
##
## [!] The per-combo random-arm ratio is computed PER DRAW
## (fdr_r = n_surr_r / max(n_obs_r, 1), kept only where n_obs_r > 0) and
## THEN averaged -- never mean(n_obs_r) first and max(.,1) after. The
## latter reintroduces the exact degenerate-denominator artefact fixed in
## R/12 (verified directly: 88%/100% of combos had a per-combo MEAN random
## n_obs below 1 at floor=20/50, which would silently floor the
## denominator almost every time and manufacture a spurious "random looks
## better" result on its own). Caught and fixed before any conclusion was
## drawn from it -- see module_sim_3sp53/README.md's 2026-09-09 write-up.
##
## SCOPE: full grid as of 2026-09-09 -- all 7 cells x both tags (bgs,
## nobgs), full 10x10 rep x env grid per cell (100 combos each, 1400
## total), group-null / emmax_consensus only. Result (see README.md for
## the full write-up): H0 ("size-based restriction beats a matched-
## cardinality random one") is REJECTED, identically in bgs and nobgs --
## which is itself informative: it rules out a recombination/BGS-specific
## structure confound as the driver and points instead at a consensus_
## dosage-specific noise/power effect (see ld_unit_matrix.R's own
## docstring: bigger units are a better-estimated, lower-noise summary
## variable -- more power for BOTH real signal and residual structure).
## [!] UPDATED 2026-09-09 -- PK: "run it on emmax_simes too." Added a 3rd
## arg (default "consensus"). For "simes": per-marker EMMAX on GTs
## directly (not the unit matrix), unit p-values via `.simes()`'s min(n*
## p_(i)/i) over each unit's member markers -- EXACTLY R/12's EMMAX-Simes
## section (Pm/pm_obs/statistic="simes"), a genuine multiple-comparisons
## combination across a unit's markers, unlike consensus_dosage's single
## averaged variable. N_PERM drops to 100 for simes, matching R/12's own
## N_PERM_SIMES (vs N_PERM_CONSENSUS=200) -- Simes is markedly slower
## per-draw (per-marker EMMAX + a .simes() pass per unit vs. one EMMAX
## call on the small unit matrix), so this keeps wall-time comparable.
## Output/pooling: consensus keeps the original rrc_<tag>_<cell>.rds name
## (no `arm` column in those older rows -- add one when pooling, it's
## implicitly "emmax_consensus"); simes writes rrc_<tag>_<cell>_simes.rds
## and DOES carry an explicit `arm` column.
##
## Usage: Rscript R/14_random_removal_control.R <tag> <cell> [consensus|simes]
## Writes out/14_random_removal_control/rrc_<tag>_<cell>[_simes].rds. Pool
## all per-cell files into results/random_removal_control_summary.rds with
## rbindlist(lapply(Sys.glob("out/14_random_removal_control/rrc_*.rds"),
## readRDS), fill=TRUE) + saveRDS -- fill=TRUE matters, an early batch of
## consensus files was written with 2 extra now-dropped intermediate
## columns, and simes files carry the extra `arm` column consensus ones
## don't.
suppressMessages({library(data.table); library(LDscnR)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53"), "R", "00_config.R"))

args <- commandArgs(trailingOnly = TRUE)
TAG <- args[1]; CELL <- args[2]
ARM <- if (length(args) >= 3 && !is.na(args[3])) args[3] else "consensus"
if (is.na(TAG) || is.na(CELL)) stop("Usage: Rscript R/14_random_removal_control.R <tag> <cell> [consensus|simes]")
if (!ARM %in% c("consensus", "simes")) stop("3rd arg must be \"consensus\" or \"simes\", got: ", ARM)
SIZE_FLOOR <- 2L; ALPHA <- 0.05
SIZE_FLOOR_GRID <- c(2, 3, 5, 10, 20, 50)
N_PERM <- if (ARM == "consensus") 200L else 100L  ## matches R/12's N_PERM_CONSENSUS / N_PERM_SIMES
N_RANDOM <- 100L    ## random subsets drawn per (combo, floor); vectorised rowSums makes this cheap
STATISTIC <- if (ARM == "consensus") "unit" else "simes"

say("=== R/14_random_removal_control: %s/%s, group null, emmax_%s ===\n\n", TAG, CELL, ARM)

all_rows <- list()
for (REP in 1:10) for (ENVN in 1:10) {
  bf <- file.path(PATHS$out, "02_bundle", sprintf("bundle_%s_rep%d_%s_env%d.rds", TAG, REP, CELL, ENVN))
  if (!file.exists(bf)) next
  bd <- readRDS(bf)
  GTs <- bd$GTs; map <- bd$map; stage1 <- bd$stage1; GRM <- bd$GRM; env <- bd$env
  y <- env$env; n_ind <- nrow(GTs)

  ## same group-null recipe as R/12_structured_null.R section 1
  pos <- unique(env[, .(pop, x, y)])
  set.seed(1L); km <- kmeans(pos[, .(x, y)], centers = 5, nstart = 10)
  pos[, pop_group := km$cluster]
  POPT <- merge(unique(env[, .(pop, env)]), pos[, .(pop, pop_group)], by = "pop")
  perm_group <- function() { pt <- copy(POPT)[, ep := sample(env), by = pop_group]; pt$ep[match(env$pop, pt$pop)] }

  units_base <- LDscnR:::.ld_outlier_units(stage1, map, SIZE_FLOOR)
  N_UNITS <- nrow(units_base)

  if (ARM == "consensus") {
    um <- ld_unit_matrix(GTs, stage1, map, size_floor = SIZE_FLOOR, repr = "consensus")
    Pu <- emmax_setup(um, GRM)
    p_obs <- emmax_fast(Pu, y)
  } else {
    Pu <- emmax_setup(GTs, GRM)
    p_obs <- emmax_fast(Pu, y)
  }
  obs_units <- LDscnR:::.ld_outlier_tested_units(stage1, map, p_obs, STATISTIC, SIZE_FLOOR, ALPHA, units = units_base)

  sig_matrix <- matrix(FALSE, nrow = N_PERM, ncol = N_UNITS)
  for (b in seq_len(N_PERM)) {
    set.seed(b)
    p_perm <- emmax_fast(Pu, perm_group())
    u <- LDscnR:::.ld_outlier_tested_units(stage1, map, p_perm, STATISTIC, SIZE_FLOOR, ALPHA, units = units_base)
    sig_matrix[b, ] <- u$significant
  }

  set.seed(42)
  for (f in SIZE_FLOOR_GRID) {
    size_ids <- which(units_base$n_markers >= f)
    m <- length(size_ids)
    n_obs_size <- sum(obs_units$significant[size_ids])
    n_surr_size <- mean(rowSums(sig_matrix[, size_ids, drop = FALSE]))

    rnd <- rbindlist(lapply(seq_len(N_RANDOM), function(r) {
      rid <- sample.int(N_UNITS, m)
      n_obs_r <- sum(obs_units$significant[rid])
      n_surr_r <- mean(rowSums(sig_matrix[, rid, drop = FALSE]))
      data.table(n_obs_r = n_obs_r, n_surr_r = n_surr_r, fdr_r = n_surr_r / max(n_obs_r, 1))
    }))

    all_rows[[length(all_rows) + 1]] <- data.table(
      tag = TAG, cell = CELL, arm = sprintf("emmax_%s", ARM), rep = REP, env = ENVN, size_floor = f, pool_m = m,
      n_obs_size = n_obs_size, n_surr_size = n_surr_size,
      n_random_draws_pos = sum(rnd$n_obs_r > 0),
      frac_random_draws_pos = mean(rnd$n_obs_r > 0),
      fdr_random_conditional = mean(rnd$fdr_r[rnd$n_obs_r > 0]))
  }
  say("  rep%d env%d done\n", REP, ENVN)
}

out <- rbindlist(all_rows)
out[, fdr_size := n_surr_size / max(n_obs_size, 1), by = seq_len(nrow(out))]

say("\n=== pooled (n_obs_size>0 only): size-based vs random-matched realised_fdr by floor ===\n")
print(out[n_obs_size > 0, .(n_combo = .N, mean_n_obs_size = mean(n_obs_size), mean_fdr_size = mean(fdr_size),
                             n_combo_any_random_pos = sum(n_random_draws_pos > 0),
                             mean_fdr_random = mean(fdr_random_conditional, na.rm = TRUE)),
          by = size_floor][order(size_floor)])

OUT_DIR <- file.path(PATHS$out, "14_random_removal_control")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(OUT_DIR, if (ARM == "consensus") sprintf("rrc_%s_%s.rds", TAG, CELL)
                          else sprintf("rrc_%s_%s_simes.rds", TAG, CELL))
saveRDS(out, OUT)
say("\nwrote %s (%d rows)\n", OUT, nrow(out))
