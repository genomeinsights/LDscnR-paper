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
## TWO NULL SCHEMES, BOTH FED THROUGH THE SAME ld_outlier_perm() MACHINERY
## module_3sp/9sp use (level = "units", the fast path -- no region assembly
## per surrogate):
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
               n_perm_simes = N_PERM_SIMES, n_groups = 5L)
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

.summarise_perm <- function(null, scheme, arm) {
  data.table(tag = TARGET_TAG, cell = TARGET_CELL, rep = TARGET_REP, env = TARGET_ENV,
            arm = arm, scheme = scheme, B = null$params$B, observed = null$observed,
            mean_surrogate = mean(null$surrogates), p = null$p, realised_fdr = null$realised_fdr)
}
ROWS <- list()

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

say("    [group] null: %d population-group surrogates\n", N_PERM_CONSENSUS)
p_perm_con_group <- function(bb) { set.seed(bb); emmax_fast(Pu, perm_group()) }
null_con_group <- ld_outlier_perm(test_con, stage1, map, p_perm_con_group, GTs = GTs, LD_decay = LD_decay,
                                  B = N_PERM_CONSENSUS, level = "units", verbose = TRUE)
say("    realised_fdr = %.3f (p = %.4f)\n", null_con_group$realised_fdr, null_con_group$p)
ROWS$con_group <- .summarise_perm(null_con_group, "group", "emmax_consensus")

say("    [mvn] negative control: %d surrogates\n", N_PERM_CONSENSUS)
p_perm_con_mvn <- function(bb) { set.seed(bb + 1e6); emmax_fast(Pu, gen_mvn()) }
null_con_mvn <- ld_outlier_perm(test_con, stage1, map, p_perm_con_mvn, GTs = GTs, LD_decay = LD_decay,
                                B = N_PERM_CONSENSUS, level = "units", verbose = TRUE)
say("    realised_fdr = %.3f (p = %.4f)\n", null_con_mvn$realised_fdr, null_con_mvn$p)
ROWS$con_mvn <- .summarise_perm(null_con_mvn, "mvn", "emmax_consensus")

## ---- 4. EMMAX Simes ---------------------------------------------------------
say("\n[4] EMMAX Simes\n")
Pm <- emmax_setup(GTs, GRM)
pm_obs <- emmax_fast(Pm, y)
test_sim <- ld_outlier_test(stage1, map, pm_obs, statistic = "simes", size_floor = SIZE_FLOOR,
                            alpha = ALPHA, assembly = "stage2_discovered", GTs = GTs,
                            LD_decay = LD_decay, score_threshold = REGION_ASSEMBLY$score_threshold,
                            distance_threshold = REGION_ASSEMBLY$distance_threshold)
say("    observed: %d significant units\n", sum(test_sim$units$significant))

say("    [group] null: %d population-group surrogates\n", N_PERM_SIMES)
p_perm_sim_group <- function(bb) { set.seed(bb); emmax_fast(Pm, perm_group()) }
null_sim_group <- ld_outlier_perm(test_sim, stage1, map, p_perm_sim_group, GTs = GTs, LD_decay = LD_decay,
                                  B = N_PERM_SIMES, level = "units", verbose = TRUE)
say("    realised_fdr = %.3f (p = %.4f)\n", null_sim_group$realised_fdr, null_sim_group$p)
ROWS$sim_group <- .summarise_perm(null_sim_group, "group", "emmax_simes")

say("    [mvn] negative control: %d surrogates\n", N_PERM_SIMES)
p_perm_sim_mvn <- function(bb) { set.seed(bb + 1e6); emmax_fast(Pm, gen_mvn()) }
null_sim_mvn <- ld_outlier_perm(test_sim, stage1, map, p_perm_sim_mvn, GTs = GTs, LD_decay = LD_decay,
                                B = N_PERM_SIMES, level = "units", verbose = TRUE)
say("    realised_fdr = %.3f (p = %.4f)\n", null_sim_mvn$realised_fdr, null_sim_mvn$p)
ROWS$sim_mvn <- .summarise_perm(null_sim_mvn, "mvn", "emmax_simes")

## ---- 5. save ------------------------------------------------------------------
summary_row <- rbindlist(ROWS)
print(summary_row)
OUT <- file.path(stage_dir(STAGE), sprintf("structnull_%s_rep%d_%s_env%d.rds",
                                           TARGET_TAG, TARGET_REP, TARGET_CELL, TARGET_ENV))
dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(summary = summary_row, pop_groups = pos,
            test_con = test_con, test_sim = test_sim,
            null_con_group = null_con_group, null_con_mvn = null_con_mvn,
            null_sim_group = null_sim_group, null_sim_mvn = null_sim_mvn), OUT)
write_receipt(STAGE, inputs = INPUTS, params = PARAMS, outputs = OUT, target = combo_id)
say("\n[5] wrote %s\n    receipt: %s\n", OUT, receipt_path(STAGE, combo_id))
