## =============================================================================
## Does background selection remove some isolated marker-wise noise?
##
## This downstream sensitivity analysis uses the already-scored primary c=1
## simulations. BGS and no-BGS runs sharing a map/burn-in identifier are kept
## paired. All three selection strengths and all ten environmental
## continuations remain together when map/burn-in identifiers are resampled.
## =============================================================================
suppressMessages(library(data.table))
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis"),
                 "R", "00_config.R"))
say("=== 18_bgs_marker_noise ===\n\n")

C1_CELLS <- c("V0.5_c1", "V1_c1", "V2_c1")
N_BOOT <- 2000L
BOOT_SEED <- 180921L

red <- fread(file.path(MODULE_ROOT, "results", "simulation_marker_qtn_redundancy_raw.tsv"))[
  cell %chin% C1_CELLS & method == "emmax_snp"]
sing <- fread(file.path(MODULE_ROOT, "results", "simulation_marker_qtn_singleton_raw.tsv"))
perf <- fread(file.path(MODULE_ROOT, "results", "simulation_performance.tsv"))[
  cell %chin% C1_CELLS & method == "emmax_snp_region"]
score_raw <- readRDS(file.path(MODULE_ROOT, "results", "simulation_summary_full.rds"))$raw[
  cell %chin% C1_CELLS &
    method %chin% c("emmax_snp_region", "emmax_simes_region", "emmax_consensus_region")]

expected <- CJ(tag = TAGS_ALL, cell = C1_CELLS, rep = REPS_ALL, env = ENVS_ALL)
stopifnot(nrow(red) == nrow(expected), nrow(sing) == nrow(expected),
          uniqueN(red[, .(tag, cell, rep, env)]) == nrow(expected),
          uniqueN(sing[, .(tag, cell, rep, env)]) == nrow(expected),
          nrow(perf) == length(TAGS_ALL) * length(C1_CELLS),
          nrow(score_raw) == nrow(expected) * 3L)

## Pool within map/burn-in identity before resampling. This keeps the three
## selection strengths, ten environmental continuations and paired BGS
## treatments together in every bootstrap draw.
map_counts <- merge(
  red[, .(marker_calls = sum(n_significant),
          qtn_linked_calls = sum(n_qtn_linked_naive),
          unmatched_calls = sum(n_false_positive)), by = .(rep, tag)],
  sing[, .(singleton_calls = sum(n_singleton),
           singleton_unmatched = sum(unmatched_singleton)), by = .(rep, tag)],
  by = c("rep", "tag"))

wide <- dcast(map_counts, rep ~ tag,
              value.var = c("marker_calls", "qtn_linked_calls", "unmatched_calls",
                            "singleton_calls", "singleton_unmatched"))
stopifnot(nrow(wide) == length(REPS_ALL), !anyNA(wide))

summarise_draw <- function(z) {
  ratio_reduction <- function(stem)
    100 * (1 - sum(z[[paste0(stem, "_bgs")]]) / sum(z[[paste0(stem, "_nobgs")]]))
  frac <- function(stem, tag)
    sum(z[[paste0(stem, "_", tag)]]) / sum(z[[paste0("singleton_calls_", tag)]])
  f_no <- frac("singleton_unmatched", "nobgs")
  f_bg <- frac("singleton_unmatched", "bgs")
  c(marker_call_reduction_pct = ratio_reduction("marker_calls"),
    qtn_linked_call_reduction_pct = ratio_reduction("qtn_linked_calls"),
    unmatched_call_reduction_pct = ratio_reduction("unmatched_calls"),
    singleton_call_reduction_pct = ratio_reduction("singleton_calls"),
    singleton_unmatched_call_reduction_pct = ratio_reduction("singleton_unmatched"),
    singleton_unmatched_nobgs_pct = 100 * f_no,
    singleton_unmatched_bgs_pct = 100 * f_bg,
    singleton_unmatched_decrease_pp = 100 * (f_no - f_bg))
}

point <- summarise_draw(wide)
set.seed(BOOT_SEED)
draws <- replicate(N_BOOT, summarise_draw(wide[sample(.N, .N, replace = TRUE)]))
out <- data.table(
  metric = names(point),
  estimate = as.numeric(point),
  lo = apply(draws, 1L, quantile, probs = 0.025, na.rm = TRUE),
  hi = apply(draws, 1L, quantile, probs = 0.975, na.rm = TRUE),
  n_boot = N_BOOT)

## The number of common markers tested is an exact analysis-output count, not
## a bootstrap estimate. It checks whether fewer significant calls under BGS
## merely reflect fewer SNPs entering the association scan.
tests <- perf[, .(marker_tests = sum(n_tested)), by = tag]
test_change <- 100 * (tests[tag == "bgs", marker_tests] /
                        tests[tag == "nobgs", marker_tests] - 1)
out <- rbind(out, data.table(metric = "post_maf_marker_test_increase_pct",
                             estimate = test_change, lo = NA_real_, hi = NA_real_, n_boot = 0L))

counts <- merge(
  red[, .(marker_calls = sum(n_significant),
          qtn_linked_calls = sum(n_qtn_linked_naive),
          unmatched_calls = sum(n_false_positive)), by = tag],
  sing[, .(singleton_calls = sum(n_singleton),
           singleton_unmatched = sum(unmatched_singleton)), by = tag],
  by = "tag")
counts <- merge(counts, tests, by = "tag")
counts[, singleton_unmatched_fraction := singleton_unmatched / singleton_calls]

## Compare the aggregation-minus-marker precision and recall contrasts between
## matched BGS and no-BGS simulations within each bootstrap draw.
method_order <- c("emmax_snp_region", "emmax_simes_region", "emmax_consensus_region")
score_map <- score_raw[, .(TP = sum(TP), FP = sum(FP),
                           recovered = sum(n_recovered), opportunities = sum(n_detectable_qtn)),
                       by = .(rep, tag, method)]
score_map[, key := factor(paste(tag, method),
                          levels = as.vector(outer(TAGS_ALL, method_order, paste)))]
setorder(score_map, rep, key)
score_mats <- lapply(c("TP", "FP", "recovered", "opportunities"), function(v)
  xtabs(as.formula(paste(v, "~ rep + key")), score_map))
names(score_mats) <- c("TP", "FP", "recovered", "opportunities")
keys <- colnames(score_mats$TP)

aggregation_contrasts <- function(weights) {
  z <- lapply(score_mats, function(m) colSums(m * weights))
  precision <- z$TP / (z$TP + z$FP)
  recall <- z$recovered / z$opportunities
  names(precision) <- names(recall) <- keys
  val <- function(tag, method, x) unname(x[paste(tag, method)])
  contrast <- function(tag, method, x)
    val(tag, method, x) - val(tag, "emmax_snp_region", x)
  c(
    simes_precision_gain_nobgs = contrast("nobgs", "emmax_simes_region", precision),
    simes_precision_gain_bgs = contrast("bgs", "emmax_simes_region", precision),
    simes_precision_bgs_interaction =
      contrast("bgs", "emmax_simes_region", precision) -
      contrast("nobgs", "emmax_simes_region", precision),
    consensus_precision_gain_nobgs = contrast("nobgs", "emmax_consensus_region", precision),
    consensus_precision_gain_bgs = contrast("bgs", "emmax_consensus_region", precision),
    consensus_precision_bgs_interaction =
      contrast("bgs", "emmax_consensus_region", precision) -
      contrast("nobgs", "emmax_consensus_region", precision),
    simes_recall_change_nobgs = contrast("nobgs", "emmax_simes_region", recall),
    simes_recall_change_bgs = contrast("bgs", "emmax_simes_region", recall),
    simes_recall_bgs_interaction =
      contrast("bgs", "emmax_simes_region", recall) -
      contrast("nobgs", "emmax_simes_region", recall),
    consensus_recall_change_nobgs = contrast("nobgs", "emmax_consensus_region", recall),
    consensus_recall_change_bgs = contrast("bgs", "emmax_consensus_region", recall),
    consensus_recall_bgs_interaction =
      contrast("bgs", "emmax_consensus_region", recall) -
      contrast("nobgs", "emmax_consensus_region", recall))
}

interaction_point <- aggregation_contrasts(rep(1, length(REPS_ALL)))
set.seed(BOOT_SEED + 1L)
interaction_draws <- replicate(N_BOOT, {
  weights <- tabulate(sample(seq_along(REPS_ALL), length(REPS_ALL), replace = TRUE),
                      nbins = length(REPS_ALL))
  aggregation_contrasts(weights)
})
interaction_out <- data.table(
  metric = names(interaction_point), estimate = as.numeric(interaction_point),
  lo = apply(interaction_draws, 1L, quantile, probs = 0.025, na.rm = TRUE),
  hi = apply(interaction_draws, 1L, quantile, probs = 0.975, na.rm = TRUE),
  n_boot = N_BOOT)

fwrite(out, file.path(MODULE_ROOT, "results", "simulation_bgs_marker_noise.tsv"), sep = "\t")
fwrite(counts, file.path(MODULE_ROOT, "results", "simulation_bgs_marker_noise_counts.tsv"), sep = "\t")
fwrite(interaction_out,
       file.path(MODULE_ROOT, "results", "simulation_bgs_aggregation_interaction.tsv"), sep = "\t")
print(out)
print(counts)
print(interaction_out)

write_receipt("18_bgs_marker_noise",
              inputs = c("results/simulation_marker_qtn_redundancy_raw.tsv",
                         "results/simulation_marker_qtn_singleton_raw.tsv",
                         "results/simulation_performance.tsv",
                         "results/simulation_summary_full.rds"),
              params = list(c1_cells = C1_CELLS, n_bootstrap = N_BOOT,
                            bootstrap_seed = BOOT_SEED,
                            bootstrap_unit = "paired map/burn-in identifier"),
              outputs = c("results/simulation_bgs_marker_noise.tsv",
                          "results/simulation_bgs_marker_noise_counts.tsv",
                          "results/simulation_bgs_aggregation_interaction.tsv"))
cat("BGS_MARKER_NOISE_DONE\n")
