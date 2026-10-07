## =============================================================================
## final_analysis/R/23_floor_stability_tiers_and_summary.R
##
## EXPLORATORY, additive analysis (PK, 2026-09-27) -- part 3 of 3 (see R/21's
## header for the full question/scope). Does NOT touch R/12_floor_
## decomposition.R, results/simulation_stage2_region_details.rds, or any
## manuscript value -- reads the latter only as a SEPARATE reference point
## (existing floor=2 canonical analysis), never blended with the stability
## filter's own numbers, per PK's explicit instruction that changing the
## canonical floor and imposing a stability filter are two different
## interventions.
##
## Canonical (99.7%-floor) regions are the base set. A canonical region is
## marked "recovered" at a side floor (99.5% or 99.9%) using
## match_region_sets() (R/helpers_floor_stability.R): >=1 shared core_snp AND
## compatible chromosome/physical overlap with some region at that side
## floor. `tier_naive` credits every recovered side-check; `tier_effective`
## does NOT credit a side whose floor is IDENTICAL to the canonical floor for
## that combo (a collapsed floor -- recovering against yourself is not an
## independent confirmation). Headline reporting uses `tier_effective`;
## `tier_naive` is kept in the per-region audit table for transparency.
##
## Pooling: counts summed across the 600 c=1 combos before computing
## precision/recall (never means of per-combo ratios), per-combo unique-QTN
## counts summed when pooling recall's numerator (QTN identifiers are local
## to each simulated dataset, never shared across combos -- same convention
## R/12 uses). Uncertainty: map/burn-in cluster bootstrap, clustered by `rep`
## alone (the shared-map unit established this session for R/13/R/14/R/15/
## R/20 -- the c=1 cells' full 10 reps are all present here, unlike the
## "balanced" 5-per-cell subsets used elsewhere), keeping envs/methods/tags
## paired within a resampled rep.
## =============================================================================
suppressMessages({library(data.table); library(ggplot2)})
MODULE_ROOT <- path.expand("~/gitlab/LDscnR-paper/module_sim_3sp53/final_analysis")
source(file.path(MODULE_ROOT, "R", "00_config.R"))
source(file.path(MODULE_ROOT, "R", "helpers_stage2_truth.R"))   ## bootstrap_rep_matrix(), ci_quantile()
source(file.path(MODULE_ROOT, "R", "helpers_floor_stability.R"))
say("=== 23_floor_stability_tiers_and_summary ===\n\n")

C1_CELLS <- c("V0.5_c1", "V1_c1", "V2_c1")
PRIMARY_METHODS <- c("emmax_consensus", "emmax_simes")

## ---- 1. load ------------------------------------------------------------------
sel_wide <- fread(file.path(MODULE_ROOT, "results", fs_name("floor_selection_wide.tsv")))
region_dt <- readRDS(file.path(MODULE_ROOT, "results", fs_name("region_details.rds")))
combo_meta <- fread(file.path(MODULE_ROOT, "results", fs_name("combo_meta.tsv")))
say("[1] loaded %d combos' floor selection, %s region rows, %d combo-meta rows\n",
    nrow(sel_wide), format(nrow(region_dt), big.mark = ","), nrow(combo_meta))
methods_present <- intersect(c(PRIMARY_METHODS, "lfmm_simes"), unique(region_dt$method))

## ---- 2. per-region cross-floor matching + tier assignment ---------------------
## One row per canonical (floor_997) region, for one (combo, method).
.tier_one_combo_method <- function(cid, m, s_row) {
  canon <- region_dt[combo_id == cid & method == m & floor == s_row$floor_997]
  if (!nrow(canon)) return(NULL)
  side995 <- region_dt[combo_id == cid & method == m & floor == s_row$floor_995]
  side999 <- region_dt[combo_id == cid & method == m & floor == s_row$floor_999]

  m995 <- match_region_sets(canon, side995, use_core_snp = TRUE)
  m999 <- match_region_sets(canon, side999, use_core_snp = TRUE)
  m995_ov <- match_region_sets(canon, side995, use_core_snp = FALSE)
  m999_ov <- match_region_sets(canon, side999, use_core_snp = FALSE)

  canon[, `:=`(
    combo_id = cid, method = m,
    recovered_995 = m995$recovered, n_matches_995 = m995$n_matches,
    recovered_999 = m999$recovered, n_matches_999 = m999$n_matches,
    recovered_995_overlap_only = m995_ov$recovered, recovered_999_overlap_only = m999_ov$recovered,
    collapsed_995_997 = s_row$collapsed_995_997, collapsed_997_999 = s_row$collapsed_997_999)]
  canon[, tier_naive := 1L + as.integer(recovered_995) + as.integer(recovered_999)]
  canon[, tier_effective := 1L + as.integer(recovered_995 & !collapsed_995_997) +
                             as.integer(recovered_999 & !collapsed_997_999)]
  canon[, tier_overlap_only := 1L + as.integer(recovered_995_overlap_only & !collapsed_995_997) +
                                as.integer(recovered_999_overlap_only & !collapsed_997_999)]
  canon
}

say("[2] cross-floor matching for %d combos x %d methods\n", nrow(sel_wide), length(methods_present))
tier_rows <- list()
for (i in seq_len(nrow(sel_wide))) {
  s <- sel_wide[i]
  for (m in methods_present) {
    r <- .tier_one_combo_method(s$combo_id, m, s)
    if (!is.null(r)) tier_rows[[length(tier_rows) + 1]] <- r
  }
}
tiered <- rbindlist(tier_rows, fill = TRUE)
say("    %s canonical regions tiered (%d had zero canonical regions and contribute none)\n",
    format(nrow(tiered), big.mark = ","),
    nrow(sel_wide) * length(methods_present) - uniqueN(tiered[, .(combo_id, method)]))

## ---- 3. ambiguous one-to-many match audit --------------------------------------
ambiguous <- tiered[n_matches_995 > 1 | n_matches_999 > 1]
say("[3] %d/%d canonical regions have an ambiguous (>1) side-floor match (audited, not silently resolved)\n",
    nrow(ambiguous), nrow(tiered))

## ---- 4. per-region audit table -------------------------------------------------
audit_cols <- tiered[, .(combo_id, method, Chr, from, to, n_markers, n_units,
                        member_markers = vapply(member_markers, function(x) paste(x, collapse = ";"), character(1)),
                        member_units = vapply(member_units, function(x) paste(x, collapse = ";"), character(1)),
                        linked_qtn = vapply(linked_qtn, function(x) paste(x, collapse = ";"), character(1)),
                        TP, recovered_995, recovered_999, n_matches_995, n_matches_999,
                        collapsed_995_997, collapsed_997_999, tier_naive, tier_effective, tier_overlap_only)]
audit_cols <- merge(audit_cols, sel_wide[, .(combo_id, floor_995, floor_997, floor_999)], by = "combo_id")
fwrite(audit_cols, file.path(MODULE_ROOT, "results", fs_name("region_audit.tsv")), sep = "\t")
say("[4] wrote results/floor_stability_region_audit.tsv (%d rows)\n", nrow(audit_cols))

## ---- 5. zero-call frequency (never silently dropped) ---------------------------
intended <- CJ(combo_id = combo_meta$combo_id, method = methods_present)
intended <- intended[!(method == "lfmm_simes" & !combo_id %chin% combo_meta[have_lfmm == TRUE, combo_id])]
called <- unique(tiered[, .(combo_id, method, has_call = TRUE)])
zero_call <- merge(intended, called, by = c("combo_id", "method"), all.x = TRUE)
zero_call[is.na(has_call), has_call := FALSE]
zero_call <- merge(zero_call, combo_meta[, .(combo_id, tag, cell, rep)], by = "combo_id")
zero_call[, V := sub("_c.*", "", cell)]
zero_summary <- zero_call[, .(n_combos = .N, n_zero_call = sum(!has_call), frac_zero_call = mean(!has_call)),
                          by = .(method, V, tag)]
fwrite(zero_summary, file.path(MODULE_ROOT, "results", fs_name("zero_call_summary.tsv")), sep = "\t")
say("[5] wrote results/floor_stability_zero_call_summary.tsv -- overall zero-call rate by method:\n")
print(zero_call[, .(n = .N, n_zero = sum(!has_call), frac_zero = mean(!has_call)), by = method])

## ---- 6. pooled precision/recall/FP-reduction, by call set ----------------------
## Call sets: ALL (tier_effective>=1, i.e. every canonical region), GE2
## (>=2/3), EQ3 (=3/3); incremental EQ1/EQ2/EQ3 partition ALL. Recall's
## denominator is summed detectable QTN across ALL 600 combos (combo_meta),
## fixed for every call set/method -- per-combo unique-recovered-QTN counts
## (never a re-uniqued global set -- QTN ids are local to each combo) summed
## for the numerator, restricted to whichever regions the call set retains.
n_dq_total <- sum(combo_meta$n_detectable_qtn)

.per_combo_stats <- function(d) {
  ## d: a subset of `tiered` (one call set x method). One row per (combo,rep,
  ## cell,tag) with TP/FP/n_recovered summed/deduped WITHIN that combo.
  if (!nrow(d)) return(data.table(combo_id = character(), TP = integer(), FP = integer(), n_recovered = integer()))
  d[, .(TP = sum(TP), FP = sum(!TP),
       n_recovered = length(unique(unlist(linked_qtn[TP])))), by = combo_id]
}

CALL_SETS <- list(
  all_ge1     = function(d) d,
  ge2         = function(d) d[tier_effective >= 2],
  eq3         = function(d) d[tier_effective == 3],
  incr_eq1    = function(d) d[tier_effective == 1],
  incr_eq2    = function(d) d[tier_effective == 2],
  incr_eq3    = function(d) d[tier_effective == 3]
)

## `n_dq_denom`: detectable-QTN denominator, summed ONLY over the combos in
## `rep_lookup`'s scope (the full 600 for the headline; a stratum's own
## combos for the by-V/by-tag breakdown) -- NEVER the global 600-combo total
## when reporting a stratum, since that would silently answer a different
## question ("this stratum's share of ALL combos' QTNs", not "recall within
## this stratum"). Passed explicitly, never re-derived from a hardcoded
## global inside this function, so a caller cannot forget to restrict it.
.pool_and_ci <- function(cs_dt, method_label, call_set_label, rep_lookup, n_dq_denom) {
  pc <- .per_combo_stats(cs_dt)
  pc <- merge(rep_lookup, pc, by = "combo_id", all.x = TRUE)
  pc[is.na(TP), `:=`(TP = 0L, FP = 0L, n_recovered = 0L)]     ## zero-call combos are REAL zeros, kept
  pooled_TP <- sum(pc$TP); pooled_FP <- sum(pc$FP); pooled_rec <- sum(pc$n_recovered)
  precision <- if (pooled_TP + pooled_FP > 0) pooled_TP / (pooled_TP + pooled_FP) else NA_real_
  recall <- pooled_rec / n_dq_denom

  rep_agg <- pc[, .(TP = sum(TP), FP = sum(FP), n_recovered = sum(n_recovered)), by = rep][match(REPS_ALL, rep)]
  rep_agg[is.na(TP), `:=`(TP = 0L, FP = 0L, n_recovered = 0L)]
  bs <- bootstrap_rep_matrix(as.matrix(rep_agg[, .(TP, FP, n_recovered)]), N_BOOTSTRAP, SEEDS[["bootstrap"]])
  prec_b <- bs[, "TP"] / pmax(bs[, "TP"] + bs[, "FP"], 1)
  rec_b <- bs[, "n_recovered"] / n_dq_denom
  data.table(method = method_label, call_set = call_set_label,
            n_regions = nrow(cs_dt), TP = pooled_TP, FP = pooled_FP,
            precision = precision, precision_ci_lo = ci_quantile(prec_b)[1], precision_ci_hi = ci_quantile(prec_b)[2],
            recall = recall, recall_ci_lo = ci_quantile(rec_b)[1], recall_ci_hi = ci_quantile(rec_b)[2],
            .bs_TP = list(bs[, "TP"]), .bs_FP = list(bs[, "FP"]), .bs_rec = list(bs[, "n_recovered"]))
}

say("\n[6] pooled precision/recall by call set (map-cluster bootstrap, clustered by rep)\n")
rep_lookup_full <- combo_meta[, .(combo_id, rep)]
summary_rows <- list()
boot_store <- list()   ## keep bootstrap vectors per (method,call_set) for paired contrasts below
for (m in methods_present) {
  d_m <- tiered[method == m]
  for (cs_name in names(CALL_SETS)) {
    row <- .pool_and_ci(CALL_SETS[[cs_name]](d_m), m, cs_name, rep_lookup_full, n_dq_total)
    boot_store[[paste(m, cs_name)]] <- list(TP = row$.bs_TP[[1]], FP = row$.bs_FP[[1]], rec = row$.bs_rec[[1]])
    row[, c(".bs_TP", ".bs_FP", ".bs_rec") := NULL]
    summary_rows[[length(summary_rows) + 1]] <- row
  }
}
summary_dt <- rbindlist(summary_rows)
## FP reduction relative to all_ge1, same method -- computed from the pooled
## point estimates already in summary_dt (a simple relative comparison, not
## itself bootstrapped -- see section 7 for the paired-uncertainty version).
summary_dt <- merge(summary_dt, summary_dt[call_set == "all_ge1", .(method, FP_all = FP)], by = "method")
summary_dt[, fp_reduction_vs_all := ifelse(FP_all > 0, 1 - FP / FP_all, NA_real_)]
summary_dt[, FP_all := NULL]
fwrite(summary_dt, file.path(MODULE_ROOT, "results", fs_name("summary.tsv")), sep = "\t")
say("[7] wrote results/floor_stability_summary.tsv\n")
print(summary_dt[call_set %in% c("all_ge1", "ge2", "eq3"),
                 .(method, call_set, n_regions, TP, FP, precision, recall, fp_reduction_vs_all)])
say("\n    incremental groups (do precision/recall move monotonically with tier?):\n")
print(summary_dt[call_set %in% c("incr_eq1", "incr_eq2", "incr_eq3"),
                 .(method, call_set, n_regions, TP, FP, precision, recall)])

## ---- 7. paired contrasts: does the stability filter's precision gain survive? --
say("\n[8] paired bootstrap contrasts: ge2 - all_ge1, eq3 - all_ge1 (precision, recall)\n")
contrast_rows <- list()
for (m in methods_present) for (cs in c("ge2", "eq3")) {
  base <- boot_store[[paste(m, "all_ge1")]]; comp <- boot_store[[paste(m, cs)]]
  prec_base <- base$TP / pmax(base$TP + base$FP, 1); prec_comp <- comp$TP / pmax(comp$TP + comp$FP, 1)
  rec_base <- base$rec / n_dq_total; rec_comp <- comp$rec / n_dq_total
  d_prec <- prec_comp - prec_base; d_rec <- rec_comp - rec_base
  s_base <- summary_dt[method == m & call_set == "all_ge1"]; s_comp <- summary_dt[method == m & call_set == cs]
  contrast_rows[[length(contrast_rows) + 1]] <- data.table(
    method = m, comparison = paste0(cs, "_vs_all_ge1"),
    diff_precision = s_comp$precision - s_base$precision,
    diff_precision_ci_lo = ci_quantile(d_prec)[1], diff_precision_ci_hi = ci_quantile(d_prec)[2],
    diff_recall = s_comp$recall - s_base$recall,
    diff_recall_ci_lo = ci_quantile(d_rec)[1], diff_recall_ci_hi = ci_quantile(d_rec)[2])
}
contrast_dt <- rbindlist(contrast_rows)
fwrite(contrast_dt, file.path(MODULE_ROOT, "results", fs_name("contrasts.tsv")), sep = "\t")
say("[9] wrote results/floor_stability_contrasts.tsv\n")
print(contrast_dt)

## ---- 8. by-V and by-tag breakdowns ---------------------------------------------
say("\n[10] by-V and by-BGS(tag) breakdowns, all_ge1/ge2/eq3\n")
strat_rows <- list()
combo_strat <- combo_meta[, .(combo_id, rep, V = sub("_c.*", "", cell), tag)]
for (m in methods_present) {
  d_m <- tiered[method == m]
  for (strat_var in c("V", "tag")) {
    for (lev in unique(combo_strat[[strat_var]])) {
      combos_lev <- combo_strat[get(strat_var) == lev, combo_id]
      rep_lookup_lev <- combo_strat[get(strat_var) == lev, .(combo_id, rep)]
      d_lev <- d_m[combo_id %chin% combos_lev]
      n_dq_lev <- sum(combo_meta[combo_id %chin% combos_lev, n_detectable_qtn])
      for (cs_name in c("all_ge1", "ge2", "eq3")) {
        row <- .pool_and_ci(CALL_SETS[[cs_name]](d_lev), m, cs_name, rep_lookup_lev, n_dq_lev)
        row[, c(".bs_TP", ".bs_FP", ".bs_rec") := NULL]
        row[, `:=`(stratum_var = strat_var, stratum_level = lev)]
        strat_rows[[length(strat_rows) + 1]] <- row
      }
    }
  }
}
strat_dt <- rbindlist(strat_rows)
fwrite(strat_dt, file.path(MODULE_ROOT, "results", fs_name("summary_by_strata.tsv")), sep = "\t")
say("[11] wrote results/floor_stability_summary_by_strata.tsv\n")

## ---- 9. reference comparison: existing floor=2 canonical (SEPARATE, never blended) --
say("\n[12] reference only -- existing floor=2 canonical analysis, same c=1 combos (NOT the stability filter)\n")
ref_file <- file.path(MODULE_ROOT, "results", "simulation_stage2_region_details.rds")
if (file.exists(ref_file)) {
  ref <- readRDS(ref_file)
  ## [!] EXTENDED (PK follow-up, 2026-09-27): originally EMMAX-only, but
  ## `tiered` tracks lfmm_simes throughout this script too -- omitting it
  ## here left the reference comparison, and the new cross-check below,
  ## silently incomplete for one of the three primary methods.
  ref_method_map <- c(emmax_consensus_region = "emmax_consensus", emmax_simes_region = "emmax_simes",
                      lfmm_simes_region = "lfmm_simes")
  ref_c1 <- ref[cell %chin% C1_CELLS & method %chin% names(ref_method_map)]
  ref_c1[, method := ref_method_map[method]]
  ref_c1[, combo_id := sprintf("%s_%s_rep%d_env%d", tag, cell, rep, env)]
  ## recall's numerator: per-combo unique recovered QTN, summed (same
  ## convention as .per_combo_stats() above) -- the denominator (n_dq_total)
  ## is the SAME 600-combo detectable-QTN total already computed from
  ## combo_meta, since detectable status is a property of the combo's
  ## genome, not of which floor/stage reported the region.
  ref_recall <- ref_c1[, .(n_recovered = length(unique(unlist(linked_qtn[TP])))), by = .(method, combo_id)][
    , .(n_recovered_total = sum(n_recovered)), by = method]
  ref_summary <- ref_c1[, .(n_regions = .N, TP = sum(TP), FP = sum(!TP)), by = method]
  ref_summary[, precision := TP / pmax(TP + FP, 1)]
  ref_summary <- merge(ref_summary, ref_recall, by = "method")
  ref_summary[, recall := n_recovered_total / n_dq_total]
  fwrite(ref_summary, file.path(MODULE_ROOT, "results", fs_name("reference_floor2.tsv")), sep = "\t")
  say("    wrote results/floor_stability_reference_floor2.tsv (floor=2 canonical, c=1 cells only, EMMAX methods)\n")
  print(ref_summary)
} else {
  say("    %s not found on this host -- reference comparison skipped, not fabricated\n", ref_file)
}

## ---- 10. cross-check vs. floor=2 canonical: are higher-floor regions genuinely NEW? --
## PK's follow-up question (2026-09-27), after seeing the tier results: are
## there regions found at the higher floors that floor=2 never reports? A
## region can appear only at a higher floor for two very different reasons --
## it is a genuinely novel candidate the looser floor=2 test missed, or it
## only clears BH there because a much stricter floor has far FEWER
## competing tests (a smaller candidate set, an easier per-test threshold) --
## and these two explanations make very different predictions about that
## region's precision. Same match_region_sets() criterion as the cross-floor
## tiering above (core_snp + overlap), just matched against `ref_c1`
## (section 9) instead of another floor_stability floor.
if (exists("ref_c1")) {
  say("\n[14] cross-check vs. floor=2 canonical: are higher-floor regions genuinely new?\n")
  .found_at_floor2 <- function(region_source) {
    rows <- list()
    for (m in unique(region_source$method)) {
      side <- ref_c1[method == m]
      combos_m <- unique(region_source[method == m, combo_id])
      for (cid in combos_m) {
        canon <- copy(region_source[method == m & combo_id == cid])
        s <- side[combo_id == cid]
        mm <- match_region_sets(canon, s, use_core_snp = TRUE)
        canon[, found_at_floor2 := mm$recovered]
        rows[[length(rows) + 1]] <- canon
      }
    }
    rbindlist(rows, fill = TRUE)
  }

  ## (a) across ALL three higher floors pooled (995/997/999) -- "are there
  ## regions at ANY higher floor floor=2 never reports, and how good are they?"
  all_higher <- .found_at_floor2(region_dt)
  new_vs_known <- all_higher[, .(n = .N, n_new = sum(!found_at_floor2), frac_new = mean(!found_at_floor2),
                                 precision_new = mean(TP[!found_at_floor2]),
                                 precision_known = mean(TP[found_at_floor2])), by = method]
  say("    across all three higher floors pooled (995/997/999):\n")
  print(new_vs_known)
  fwrite(all_higher[, .(combo_id, method, floor, Chr, from, to, n_markers, TP, found_at_floor2)],
        file.path(MODULE_ROOT, "results", fs_name("vs_floor2_all_floors.tsv")), sep = "\t")

  ## (b) restricted to the CANONICAL (99.7%-floor) regions, cross-tabulated
  ## against stability tier -- "does the tiering precision gain just track
  ## re-discovering floor=2's own hits, or is it doing something more?"
  ## [!] member_units MUST be kept -- match_region_sets() needs it for the
  ## core_snp criterion. Dropping it (an earlier version of this line did)
  ## does not error: cr$member_units resolves to NULL, `mu %chin% NULL` is
  ## FALSE for every side region, so every region silently comes back
  ## "not recovered" instead of failing loudly. Caught by checking the
  ## actual output (frac_found_at_floor2 == 0 for every tier, implausible on
  ## its face and contradicted by the ad hoc version of this same check run
  ## just before writing this section), not assumed correct from a clean exit code.
  canon_only <- .found_at_floor2(tiered[, .(combo_id, method, Chr, from, to, n_markers, member_units, TP, tier_effective)])
  by_tier <- canon_only[, .(n = .N, frac_found_at_floor2 = mean(found_at_floor2), precision = mean(TP)),
                        by = tier_effective][order(tier_effective)]
  by_tier_found <- canon_only[, .(n = .N, precision = mean(TP)), by = .(found_at_floor2, tier_effective)][
    order(found_at_floor2, tier_effective)]
  say("\n    canonical (99.7%%-floor) regions: found-at-floor=2 rate and precision BY TIER:\n")
  print(by_tier)
  say("\n    ... split further by found_at_floor2 x tier (does stability filtering just recover floor=2's hits?):\n")
  print(by_tier_found)
  fwrite(canon_only[, .(combo_id, method, Chr, from, to, n_markers, TP, tier_effective, found_at_floor2)],
        file.path(MODULE_ROOT, "results", fs_name("vs_floor2_canonical_by_tier.tsv")), sep = "\t")
  fwrite(by_tier_found, file.path(MODULE_ROOT, "results", fs_name("vs_floor2_by_tier_summary.tsv")), sep = "\t")
  say("[15] wrote results/floor_stability_vs_floor2_all_floors.tsv, _canonical_by_tier.tsv, _by_tier_summary.tsv\n")
} else {
  say("\n[14] ref_c1 not available (floor=2 reference file missing) -- cross-check vs floor=2 skipped\n")
}

cat("TIERS_AND_SUMMARY_DONE\n")
