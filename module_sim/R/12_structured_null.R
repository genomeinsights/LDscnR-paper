## module_sim/R/12_structured_null.R
##
## Structured-null calibration check for the EMMAX consensus/Simes arms,
## matching module_3sp/module_9sp's manuscript methodology (materials_and_
## methods.tex, "Association models and empirical calibration" /
## "Nine-spined stickleback permutation stress test") and their actual code
## (module_3sp/R/03_EMMAX.R's perm_regional() + ld_outlier_perm()) as
## closely as the simulated design allows. PK: "run the structured nulls...
## pull the same numbers as we have done in module_3sp and module_9sp."
##
## THREE NULL SCHEMES (group, mvn, spatial -- see below), all built on the
## SAME machinery module_3sp/9sp's ld_outlier_perm() uses at level = "units"
## (the fast path -- no region assembly per surrogate): .ld_outlier_tested_
## units(), called directly here (via .run_null()/.perm_sizes(), not
## ld_outlier_perm() itself) so each surrogate's full table of significant
## units -- not just a count -- is available for the minimum-size sweep in
## section "minimum stage-1-unit-size sweep" below (PK: "can we say that
## larger stage-2 clusters are more likely to be true positives? We should
## sweep for minimum cluster size"). Free given what ld_outlier_perm()
## would compute anyway -- no extra emmax_fast() fits, just more of each
## call's own output kept.
##
##   group   Population-level ENV VALUE permutation, stratified within a
##           discrete grouping factor -- the direct analog of module_3sp's
##           regional-locality / module_9sp's lineage stratification.
##           Simulated design has no explicit locality/lineage field, but
##           PK: "there are 5 groups of populations (hence K=5 for LFMM)...
##           find it in the env file and population coordinates. The grid
##           is 48x48 but we have sampled in the corners and in the
##           middle." Confirmed directly (2026-09-08): 80 sampled
##           populations fall into five tight, well-separated 4x4 blocks
##           (four corners + centre) -- trivially recovered by
##           kmeans(k=5) on population (x,y). Permuting the population-
##           level environment value WITHIN each of these five blocks,
##           propagated to individuals, keeps kinship/relatedness and the
##           within-block spatial structure intact, exactly as
##           perm_regional() does for a discrete label.
##
##   mvn     s ~ MVN(0, GRM), orthogonalised against the observed
##           phenotype by Gram-Schmidt (resid(lm(s ~ y))) -- PK: "a regular
##           MVN shows that the method itself does not cause any false
##           positives (only structure)... would be a form of negative
##           control (if something doesn't work there, then something is
##           seriously wrong)." Reuses the exact draw recipe from LDscnR's
##           own structured_null() (eigen(K) -> Vk %*% sqrt(Lv) %*%
##           rnorm(n)) without its C-score machinery -- that function is
##           built around the superseded ld_cscore()/ld_regions() pipeline
##           (module_sim_LDscnR's domain), not the current
##           ld_outlier_test()/ld_outlier_perm() API this project uses
##           everywhere else. If `mvn`'s realised_fdr is also inflated,
##           that points at the method/code, not at how well the GRM
##           captures structure -- module_9sp's own "sensitivity analysis
##           of the null construction, not... a preferred replacement"
##           framing for this same comparison.
##
## SCOPE: one (tag,cell,rep,env) combo per call, same per-combo architecture
## as every other stage -- run across a REPRESENTATIVE 14-combo subset (1
## rep x env per cell x tag; see run_structured_null_sample.sh), not the
## full 1400-combo grid. module_3sp and module_9sp each ran this ONCE, on
## one dataset; matching that scale -- enough to see where calibration
## holds or breaks across the design space -- rather than an exhaustive
## resweep that would need refitting EMMAX hundreds of times per surrogate
## across 1400 combos.
##
## LFMM IS NOT CALIBRATED HERE, matching module_3sp/9sp's own choice:
## "LFMM phenotype permutations were not performed because each valid
## surrogate would require refitting the complete genome-wide latent-
## factor model."
##
## B = 200 (consensus, cheap) / 100 (Simes, rescans every marker per draw)
## -- module_3sp's own 1000/200 split, scaled down for a 14-combo sweep
## rather than the single-dataset case it was set for.
suppressMessages({library(data.table); library(LDscnR)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim"), "R", "00_config.R"))
STAGE <- "12_structured_null"
say("=== %s ===\n\n", STAGE)
invisible(check_ldscnr())

TARGET_TAG  <- Sys.getenv("SIM_TAG",  TAGS[1])
TARGET_CELL <- Sys.getenv("SIM_CELL", CELLS[1])
TARGET_ENV  <- as.integer(Sys.getenv("SIM_ENV",  ENVS[1]))
TARGET_REP  <- as.integer(Sys.getenv("SIM_REP",  REPS[1]))
combo_id <- sprintf("%s_%s_rep%d_env%d", TARGET_TAG, TARGET_CELL, TARGET_REP, TARGET_ENV)

N_PERM_CONSENSUS <- 200L
N_PERM_SIMES     <- 100L

BUNDLE_PATH <- file.path(PATHS$out, "02_bundle",
  sprintf("bundle_%s_rep%d_%s_env%d.rds", TARGET_TAG, TARGET_REP, TARGET_CELL, TARGET_ENV))
if (!file.exists(BUNDLE_PATH)) stop("R/02_bundle.R has not produced: ", basename(BUNDLE_PATH))
bd <- readRDS(BUNDLE_PATH)
GTs <- bd$GTs; map <- bd$map; stage1 <- bd$stage1; GRM <- bd$GRM; LD_decay <- bd$LD_decay; env <- bd$env
y <- env$env
say("[0] bundle: %d individuals x %s markers\n", nrow(GTs), format(ncol(GTs), big.mark = ","))

INPUTS <- BUNDLE_PATH
PARAMS <- list(size_floor = SIZE_FLOOR, alpha = ALPHA, unit_repr = UNIT_REPR,
               region_assembly = REGION_ASSEMBLY, n_perm_consensus = N_PERM_CONSENSUS,
               n_perm_simes = N_PERM_SIMES, n_groups = 5L,
               ## [!] bumped 2026-09-08: added the `spatial` null scheme --
               ## forces a rerun.
               spatial_null = TRUE,
               ## [!] bumped 2026-09-08: added the minimum-unit-size sweep
               ## -- forces a rerun. Literal, not a reference to
               ## SIZE_FLOOR_GRID (defined later in the file) -- keep the
               ## two in sync by hand if the grid ever changes.
               size_floor_grid = c(2, 3, 5, 10, 20, 50))
if (!stage_stale(STAGE, INPUTS, PARAMS, target = combo_id) && !nzchar(Sys.getenv("FORCE"))) {
  say("\nNothing to do. Set FORCE=1 to rerun anyway.\n"); quit(save = "no")
}
dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)

## ---- 1. the 5 spatial population groups ----------------------------------------
say("[1] population groups: kmeans(k=5) on sampled-population (x,y)\n")
pos <- unique(env[, .(pop, x, y)])
set.seed(1L)
km <- kmeans(pos[, .(x, y)], centers = 5, nstart = 10)
pos[, pop_group := km$cluster]
say("    %d populations in 5 groups, sizes: %s\n", nrow(pos), paste(table(pos$pop_group), collapse = ","))
POPT <- merge(unique(env[, .(pop, env)]), pos[, .(pop, pop_group)], by = "pop")

perm_group <- function() {
  pt <- copy(POPT)[, ep := sample(env), by = pop_group]
  pt$ep[match(env$pop, pt$pop)]
}

## ---- 2. MVN(0, GRM) negative control, orthogonalised against observed y --------
say("[2] MVN(0, GRM) negative control (structured_null()'s own draw recipe)\n")
eK <- eigen(GRM, symmetric = TRUE)
Lv <- pmax(eK$values, 0); Vk <- eK$vectors
n_ind <- nrow(GTs)
gen_mvn <- function() {
  s <- as.numeric(Vk %*% (sqrt(Lv) * stats::rnorm(n_ind)))
  as.numeric(stats::resid(stats::lm(s ~ y)))
}

## ---- 2b. spatial MVN, orthogonalised against observed y -------------------------
## ADDED 2026-09-08 (PK, on `group` being "very conservative": "the
## environmental gradient is very slow, hence the strong correlation
## between structure and env... is there a way to make the env orthogonal
## to the observed and test that as a third alternative?"). Confirmed
## directly (section 4 below / R/06_popgen_summary.R's new env_icc): ~88%
## of env variance sits BETWEEN the 5 population groups, so `group`'s
## within-group permutation barely perturbs the phenotype -- a weak null by
## construction, not a property of the population-structure hypothesis
## itself. `spatial` draws from genuine spatial autocorrelation (a Gaussian
## kernel over individual (x,y), independent of the population-group
## binning) and orthogonalises the SAME way as `mvn` -- structured_null()'s
## own basis = "spatial" recipe, reused directly rather than reinvented.
say("[2b] spatial MVN negative control (structured_null()'s basis = \"spatial\" recipe)\n")
coords <- as.matrix(env[, .(x, y)])
Dm <- as.matrix(stats::dist(coords)); l_bw <- stats::median(Dm[lower.tri(Dm)])
eK_sp <- eigen(exp(-0.5 * (Dm / l_bw)^2), symmetric = TRUE)
Lv_sp <- pmax(eK_sp$values, 0); Vk_sp <- eK_sp$vectors
gen_spatial <- function() {
  s <- as.numeric(Vk_sp %*% (sqrt(Lv_sp) * stats::rnorm(n_ind)))
  as.numeric(stats::resid(stats::lm(s ~ y)))
}

ROWS <- list()

## ---- minimum stage-1-unit-size sweep: does raising the floor improve calibration? --
## PK: "The permutation is an assertion independent of the truth, telling
## that probably a large portion of the outlier regions are false
## positives (but we know not all of them are). Can we say that larger
## stage-2 clusters are more likely to be true positives? We should sweep
## for minimum cluster size to see if the figure above improves?"
##
## Cluster size here is n_markers, the SAME size every FP-by-size figure
## this session has used (R/04_score.R's n_loci -- the tested stage-1 unit,
## not a further stage-2 grouping of units; this project's own cluster_
## detail/fp_by_size analyses have always scored at the stage-1-unit level,
## never re-aggregated to stage-2 regions, so this sweep matches what's
## already being reported everywhere else rather than introducing a new
## unit). FREE given what's already being computed above: no extra
## emmax_fast() calls -- .ld_outlier_tested_units() is cheap, and each
## surrogate's full table of significant units (with their sizes) is kept
## instead of collapsed to a single count.
SIZE_FLOOR_GRID <- c(2, 3, 5, 10, 20, 50)
## [!] LDscnR:::, not a bare call -- these are internal (dot-prefixed,
## unexported) functions; library(LDscnR) (unlike devtools::load_all(),
## which some other scripts in this project rely on and which attaches
## everything) does not put them on the search path.
units_base <- LDscnR:::.ld_outlier_units(stage1, map, SIZE_FLOOR)

.perm_sizes <- function(p_perm, B, statistic) {
  lapply(seq_len(B), function(b) {
    u <- LDscnR:::.ld_outlier_tested_units(stage1, map, p_perm(b), statistic, SIZE_FLOOR, ALPHA, units = units_base)
    u$n_markers[u$significant]
  })
}
## [!] CAVEAT (added 2026-09-09, ported from module_sim_3sp53 after PK
## caught it there: "are you accounting for the fact that the number of
## false positives also drop by randomly removing clusters from the
## outlier list?") -- max(n_obs, 1) avoids a divide-by-zero when the
## OBSERVED test has no significant unit at some floor, but the eligible
## Stage-1-unit pool shrinks fast with floor, so n_obs=0 becomes the
## MAJORITY case at large floors (confirmed in module_sim_3sp53's full
## grid: 36.5% of combos at floor=2, 95.2% at floor=50). For those rows
## realised_fdr is NOT an FDR estimate -- it's just mean_surrogate, the
## rate pure noise clears an emptying pool. This script is DORMANT
## (module_sim/ isn't the live pipeline -- see module_sim_3sp53/, whose
## R/12_structured_null.R and README.md carry the actual fix and the
## corrected numbers); n_obs is kept in the output so it CAN be
## conditioned on if this ever gets rerun for real, but don't trust the
## unconditional realised_fdr from this copy at face value.
.fdr_by_floor <- function(obs_units, surr_sizes, scheme, arm) {
  obs_sizes <- obs_units$n_markers[obs_units$significant]
  rbindlist(lapply(SIZE_FLOOR_GRID, function(f) {
    n_obs  <- sum(obs_sizes >= f)
    n_surr <- mean(vapply(surr_sizes, function(s) sum(s >= f), integer(1)))
    data.table(tag = TARGET_TAG, cell = TARGET_CELL, rep = TARGET_REP, env = TARGET_ENV,
              arm = arm, scheme = scheme, size_floor = f, n_obs = n_obs,
              mean_surrogate = n_surr, realised_fdr = n_surr / max(n_obs, 1))
  }))
}
SIZE_ROWS <- list()

## Runs the surrogate loop ONCE (via .perm_sizes()) and derives BOTH the
## headline realised_fdr (at SIZE_FLOOR, matching the original summary
## format) and the full floor-sweep table from the SAME B draws -- calling
## ld_outlier_perm() separately would redo every emmax_fast() fit a second
## time for no reason.
.run_null <- function(p_perm, B, statistic, obs_units, scheme, arm) {
  say("    [%s]: %d surrogates\n", scheme, B)
  surr_sizes <- .perm_sizes(p_perm, B, statistic)
  sweep <- .fdr_by_floor(obs_units, surr_sizes, scheme, arm)
  base <- sweep[size_floor == SIZE_FLOOR]
  say("    realised_fdr (floor=%d) = %.3f\n", SIZE_FLOOR, base$realised_fdr)
  list(summary = data.table(tag = TARGET_TAG, cell = TARGET_CELL, rep = TARGET_REP, env = TARGET_ENV,
                            arm = arm, scheme = scheme, B = B, observed = base$n_obs,
                            mean_surrogate = base$mean_surrogate, realised_fdr = base$realised_fdr),
       sweep = sweep)
}

## ---- 3. EMMAX consensus ---------------------------------------------------------
say("\n[3] EMMAX consensus\n")
um <- ld_unit_matrix(GTs, stage1, map, size_floor = SIZE_FLOOR, repr = UNIT_REPR)
Pu <- emmax_setup(um, GRM)
pu_obs <- emmax_fast(Pu, y)
test_con <- ld_outlier_test(stage1, map, pu_obs, statistic = "unit", size_floor = SIZE_FLOOR,
                            alpha = ALPHA, assembly = "stage2_discovered", GTs = GTs,
                            LD_decay = LD_decay, score_threshold = REGION_ASSEMBLY$score_threshold,
                            distance_threshold = REGION_ASSEMBLY$distance_threshold)
say("    observed: %d significant units\n", sum(test_con$units$significant))

p_perm_con_group <- function(bb) { set.seed(bb); emmax_fast(Pu, perm_group()) }
res_con_group <- .run_null(p_perm_con_group, N_PERM_CONSENSUS, "unit", test_con$units, "group", "emmax_consensus")
ROWS$con_group <- res_con_group$summary; SIZE_ROWS$con_group <- res_con_group$sweep

p_perm_con_mvn <- function(bb) { set.seed(bb + 1e6); emmax_fast(Pu, gen_mvn()) }
res_con_mvn <- .run_null(p_perm_con_mvn, N_PERM_CONSENSUS, "unit", test_con$units, "mvn", "emmax_consensus")
ROWS$con_mvn <- res_con_mvn$summary; SIZE_ROWS$con_mvn <- res_con_mvn$sweep

p_perm_con_spatial <- function(bb) { set.seed(bb + 2e6); emmax_fast(Pu, gen_spatial()) }
res_con_spatial <- .run_null(p_perm_con_spatial, N_PERM_CONSENSUS, "unit", test_con$units, "spatial", "emmax_consensus")
ROWS$con_spatial <- res_con_spatial$summary; SIZE_ROWS$con_spatial <- res_con_spatial$sweep

## ---- 4. EMMAX Simes ---------------------------------------------------------
say("\n[4] EMMAX Simes\n")
Pm <- emmax_setup(GTs, GRM)
pm_obs <- emmax_fast(Pm, y)
test_sim <- ld_outlier_test(stage1, map, pm_obs, statistic = "simes", size_floor = SIZE_FLOOR,
                            alpha = ALPHA, assembly = "stage2_discovered", GTs = GTs,
                            LD_decay = LD_decay, score_threshold = REGION_ASSEMBLY$score_threshold,
                            distance_threshold = REGION_ASSEMBLY$distance_threshold)
say("    observed: %d significant units\n", sum(test_sim$units$significant))

p_perm_sim_group <- function(bb) { set.seed(bb); emmax_fast(Pm, perm_group()) }
res_sim_group <- .run_null(p_perm_sim_group, N_PERM_SIMES, "simes", test_sim$units, "group", "emmax_simes")
ROWS$sim_group <- res_sim_group$summary; SIZE_ROWS$sim_group <- res_sim_group$sweep

p_perm_sim_mvn <- function(bb) { set.seed(bb + 1e6); emmax_fast(Pm, gen_mvn()) }
res_sim_mvn <- .run_null(p_perm_sim_mvn, N_PERM_SIMES, "simes", test_sim$units, "mvn", "emmax_simes")
ROWS$sim_mvn <- res_sim_mvn$summary; SIZE_ROWS$sim_mvn <- res_sim_mvn$sweep

p_perm_sim_spatial <- function(bb) { set.seed(bb + 2e6); emmax_fast(Pm, gen_spatial()) }
res_sim_spatial <- .run_null(p_perm_sim_spatial, N_PERM_SIMES, "simes", test_sim$units, "spatial", "emmax_simes")
ROWS$sim_spatial <- res_sim_spatial$summary; SIZE_ROWS$sim_spatial <- res_sim_spatial$sweep

## ---- 5. save ------------------------------------------------------------------
summary_row <- rbindlist(ROWS)
size_sweep <- rbindlist(SIZE_ROWS)
print(summary_row)
print(size_sweep)
OUT <- file.path(stage_dir(STAGE), sprintf("structnull_%s_rep%d_%s_env%d.rds",
                                           TARGET_TAG, TARGET_REP, TARGET_CELL, TARGET_ENV))
dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(summary = summary_row, size_sweep = size_sweep, pop_groups = pos,
            test_con = test_con, test_sim = test_sim), OUT)
write_receipt(STAGE, inputs = INPUTS, params = PARAMS, outputs = OUT, target = combo_id)
say("\n[5] wrote %s\n    receipt: %s\n", OUT, receipt_path(STAGE, combo_id))
