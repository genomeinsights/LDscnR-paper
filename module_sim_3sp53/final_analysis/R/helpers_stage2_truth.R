## =============================================================================
## final_analysis/R/helpers_stage2_truth.R
##
## Shared Stage-2 assembly + truth-scoring + bootstrap primitives, refactored
## out of 05_score_truth.R's private closures so the canonical pipeline, the
## floor-decomposition sweep (R/12) and the null-truth-calibration analysis
## (R/13) provably run the SAME code (ADDITIONAL_ANALYSES_AUDIT.md's brief:
## "The floor sweep, null analysis and existing canonical analysis must use
## the same implementation").
##
## assemble_stage2()/score_stage2_regions() are a strict behavioural superset
## of 05_score_truth.R's former .run_stage2()/.assemble_regions_from_units()/
## .assemble_regions_from_markers(): same ld_prune_and_eMLG() call, same
## seeding logic, same $groups$members provenance (never a physical-bounds
## sweep) -- now also returning per-region structural detail (Chr/from/to/
## n_markers/n_units/member ids), not just a compact members-list.
##
## [!] STABLE UNIT IDENTITY, found while implementing the floor sweep (R/12):
## LDscnR:::.ld_outlier_units()/ld_outlier_test()$units assign `unit_id` as
## seq_len() over whichever Stage-1 clusters clear THAT CALL's size_floor --
## it is NOT a stable per-cluster identity across floors. A cluster that is
## unit_id=7 at floor=2 is a different integer at floor=5 once smaller
## clusters drop out of the numbering. "Compute raw Stage-1 p-values once,
## subset by unit_id at each floor" (the natural reading of the spec's own
## optimisation) would therefore silently join the wrong cluster's p-value
## to the wrong unit at every floor but the one the p-values were computed
## at. Fixed here by keying every Stage-1-unit-level object in this file on
## `core_snp` (the cluster's representative marker -- one per cluster, by
## construction, and completely floor-independent) rather than unit_id.
## stage1$clusters itself (unfiltered by any floor) is the one stable table
## everything below is ultimately keyed against.
## =============================================================================
suppressMessages({library(data.table)})

## ---- which Stage-1 clusters get fed into Stage-2, two seeding routes --------
## Both take `stage1` (unfiltered -- stage1$clusters covers every cluster
## regardless of floor) and return a `cl_sub` data.table (subset of
## stage1$clusters) ready for assemble_stage2().

## [!] ORDER MATTERS AND IS NOW UNIFIED (PK, 2026-09-18 review). ld_prune_
## and_eMLG()'s distance-restricted dynamic cut merges ADJACENT input
## clusters, so its output partition depends on the ORDER clusters are
## handed to it, not just the SET. An earlier version of this file gave the
## two seeding routes two DIFFERENT orderings specifically to reproduce two
## different pre-refactor historical behaviours bit-for-bit (unit-seeded:
## genomic-span sorted; marker-seeded: stage1$clusters' own native,
## essentially arbitrary discovery-order row order, confirmed empirically to
## differ from genomic order at 30-99% of cluster positions across a sample
## of combos -- i.e. not "already genomic in practice"). PK's review pointed
## out that "reproducing historical output" is the wrong target once the
## SAME reporting rule is claimed to apply to every method: if arms are fed
## to Stage 2 in different orders, some of their FP/region-count
## differences are an artefact of that inconsistency, not of the methods
## themselves. Both routes are now genomic-span sorted. This DOES change the
## canonical marker-wise numbers slightly (quantified on the full 1,400-
## combo grid before making this change: emmax_snp_region 6,658->6,623
## regions, precision 0.1236->0.1228; lfmm_snp_region 11,569->11,581,
## precision 0.1104->0.1086 -- small in aggregate, though 8.5%/16.3% of
## individual EMMAX/LFMM combos get a different region count under the two
## orderings). 05_score_truth.R was regenerated for the full grid after this
## change (see ADDITIONAL_ANALYSES_AUDIT.md), and values_simulation.tex's
## affected macros were updated to match. Neither ordering is "more
## correct" in any principled sense -- ld_prune_and_eMLG()'s order-
## sensitivity is a pre-existing property of the package function, not
## introduced or fixed here -- but a single, consistent choice is required
## for "the same Stage-2 rule" to mean what it says.
## [!] ONE match() call total, not one per cluster (PK, flagged 2026-09-17):
## the original `vapply(cl_sub$members, function(mm) ...match(mm, mp$marker))`
## called match() once per cluster in cl_sub, each call re-hashing the FULL
## marker universe (mp$marker, ~15,000 markers per combo) from scratch --
## exactly the anti-pattern LDscnR:::.marker_positions()'s own header comment
## warns about (a named-vector-style lookup repeated per group costs whole
## seconds where one hash-once match() costs milliseconds). This function is
## called on every Stage-2 assembly, including every null replicate in
## Analysis 1 -- likely the dominant cost behind that run's slower-than-
## smoke-test-extrapolated wall time. Fixed by hashing mp$marker once (via
## the package's own .marker_positions()) and aggregating per-cluster minima
## with a single grouped data.table op, mirroring .ld_outlier_units()'s own
## span-construction idiom rather than inventing a new one.
.genomic_sort_clusters <- function(cl_sub, map) {
  if (!nrow(cl_sub)) return(cl_sub)
  mp <- data.table::as.data.table(map)
  flat_marker <- unlist(cl_sub$members, use.names = FALSE)
  flat_cluster <- rep.int(seq_len(nrow(cl_sub)), lengths(cl_sub$members))
  flat_idx <- LDscnR:::.marker_positions(flat_marker, mp$marker)
  span <- data.table::data.table(cl = flat_cluster, Pos = mp$Pos[flat_idx])[, .(from = min(Pos)), by = cl]
  data.table::setkey(span, cl)
  from_pos <- span[.(seq_len(nrow(cl_sub)))]$from
  cl_sub[order(cl_sub$Chr, from_pos)]
}

## Unit-seeded (Simes/consensus): `sig_core_snps` = core_snp of every
## Stage-1 unit that cleared BH at whichever floor the caller tested.
## Genomic-span sorted -- see header comment above.
stage2_seed_from_units <- function(stage1, map, sig_core_snps) {
  cl <- data.table::as.data.table(stage1$clusters)
  if (!length(sig_core_snps)) return(cl[0L])
  .genomic_sort_clusters(cl[cl$core_snp %chin% sig_core_snps], map)
}

## Marker-seeded (unrestricted comparator): every phenotype-blind Stage-1
## cluster -- INCLUDING SINGLETONS, unfiltered by any floor -- containing
## >=1 significant marker is "discovered". Genomic-span sorted, same as the
## unit-seeded route -- see header comment above.
stage2_seed_from_markers <- function(stage1, map, sig_markers) {
  cl <- data.table::as.data.table(stage1$clusters)
  if (!length(sig_markers)) return(cl[0L])
  disc <- vapply(cl$members, function(mm) any(mm %chin% sig_markers), logical(1))
  .genomic_sort_clusters(cl[disc], map)
}

## ---- Stage-2 assembly: cl_sub -> detailed reported-region table -------------
## Identical mechanics to the pre-refactor .run_stage2(): re-examine the
## already-significant clusters directly from genotypes (no window
## restriction), consolidate via ld_prune_and_eMLG()'s distance-restricted,
## quality-gated dynamic cut. Returns constituent MARKERS from
## pr$groups$members verbatim -- never a [Chr,from,to] bounds sweep.
##
## Beyond the old compact members-list, also returns: n_units (how many of
## the INPUT Stage-1 clusters in cl_sub ended up inside this region -- looked
## up via each input cluster's first member marker, since ld_prune_and_eMLG()
## only merges/splits whole input clusters, never scrambles individual
## markers across an input cluster's boundary) and member_units (their
## core_snp ids, i.e. the stable identity described above).
##
## @return data.table: region_id, Chr, from, to, n_markers, n_units,
##   member_markers (list-col), member_units (list-col, core_snp values).
##   Zero rows (with correct columns) when cl_sub is empty.
assemble_stage2 <- function(stage1, map, GTs, LD_decay, cl_sub,
                            score_threshold, distance_threshold) {
  empty <- data.table::data.table(
    region_id = integer(), Chr = character(), from = numeric(), to = numeric(),
    n_markers = integer(), n_units = integer(),
    member_markers = list(), member_units = list())
  if (!nrow(cl_sub)) return(empty)

  mk_sub <- unlist(cl_sub$members, use.names = FALSE)
  ms_sub <- data.table::as.data.table(stage1$map_snp)[marker %chin% mk_sub]
  sub <- structure(list(map_snp = ms_sub, clusters = cl_sub, pruned = cl_sub$core_snp),
                   class = "ld_complexity_reduction")
  pr <- LDscnR::ld_prune_and_eMLG(
    GTs = GTs[, mk_sub, drop = FALSE], stage1 = sub,
    ld_w_col = "ld_w_095", ld_w_threshold = 0,
    LD_decay = LD_decay, min_r2_rho = stage1$params$rho,
    score_threshold = score_threshold, distance_threshold = distance_threshold,
    compute_unflagged_eMLG = FALSE, min_n_loci_eMLG = 1, min_n_loci_flag = 1, cores = 1)
  g <- data.table::as.data.table(pr$groups)
  if (!nrow(g)) return(empty)

  g_marker <- unlist(g$members, use.names = FALSE)
  g_region <- rep.int(seq_len(nrow(g)), lengths(g$members))
  mp <- data.table::as.data.table(map)
  g_idx <- LDscnR:::.marker_positions(g_marker, mp$marker)

  spans <- data.table::data.table(region_id = g_region, Chr = as.character(mp$Chr[g_idx]),
                                  Pos = mp$Pos[g_idx])[
    , .(Chr = Chr[1], from = min(Pos), to = max(Pos), n_markers = .N), by = region_id]

  ## Each INPUT cluster (row of cl_sub) belongs to exactly one output region --
  ## look it up via any one of its own member markers (the first).
  first_marker <- vapply(cl_sub$members, `[`, character(1), 1L)
  marker_to_region <- stats::setNames(g_region, g_marker)
  cl_region <- unname(marker_to_region[first_marker])
  n_units <- data.table::data.table(region_id = cl_region)[!is.na(region_id), .(n_units = .N), by = region_id]

  out <- merge(spans, n_units, by = "region_id", all.x = TRUE)
  out[is.na(n_units), n_units := 0L]   # defensive; should not occur if cl_sub fed pr's own input
  out[, member_markers := split(g_marker, g_region)[as.character(region_id)]]
  out[, member_units := split(cl_sub$core_snp, cl_region)[as.character(region_id)]]
  data.table::setorder(out, Chr, from)
  out[]
}

## ---- attach truth status to an assembled-region table -----------------------
## `qtn_lut_match` = the combo's already-computed marker-QTN match table
## (LD/distance-pass, restricted to detectable QTN) -- built once per combo
## exactly as 05_score_truth.R already does; never recomputed here. A region
## is TP iff >=1 of its ACTUAL constituent markers (member_markers, never
## markers merely inside [from,to]) matches a detectable QTN.
##
## @return `regions` with two new columns: TP (logical), linked_qtn (list-col
##   of detectable-QTN marker ids linked to that region, possibly empty).
score_stage2_regions <- function(regions, qtn_lut_match) {
  if (!nrow(regions)) {
    regions[, TP := logical(0)]
    regions[, linked_qtn := list()]
    return(regions[])
  }
  regions[, linked_qtn := lapply(member_markers, function(mm)
    unique(qtn_lut_match[marker %in% mm, qtn_marker]))]
  regions[, TP := lengths(linked_qtn) > 0]
  regions[]
}

## ---- map/burn-in cluster bootstrap, generalised -------------------------------
## The n_rep x B multiplicity-matrix / crossprod() idiom 08_bgs_validation.R
## and 09_truth_sensitivity.R each hand-rolled inline, generalised to an
## arbitrary n_rep x k matrix of summable per-rep statistics (TP, FP,
## n_detectable_qtn, n_recovered, null-region counts, ...) instead of being
## specific to precision/recall. Resamples the `n_rep` map/burn-in rows WITH
## REPLACEMENT, B times; every column of `rep_stat_matrix` is resampled by
## the SAME draw per replicate b, so contrasts computed from the same call's
## output are paired (methods, BGS treatments, whatever the caller put in
## separate columns).
##
## @param rep_stat_matrix Numeric matrix, one row per map/burn-in rep, one
##   column per statistic to bootstrap (already summed across that rep's
##   environmental continuations, and across whatever else is being pooled).
## @param B,seed As elsewhere in this pipeline (N_BOOTSTRAP, SEEDS["bootstrap"]).
## @return B x ncol(rep_stat_matrix) matrix of bootstrap sums (same column
##   names as the input). Divide/ratio these yourself for the statistic you
##   actually want (precision, a null/obs ratio, ...) -- this function only
##   does the resampling-and-summing.
bootstrap_rep_matrix <- function(rep_stat_matrix, B, seed) {
  rep_stat_matrix <- as.matrix(rep_stat_matrix)
  n_rep <- nrow(rep_stat_matrix)
  set.seed(seed)
  draws <- matrix(sample.int(n_rep, size = n_rep * B, replace = TRUE), nrow = n_rep, ncol = B)
  mult <- apply(draws, 2, tabulate, nbins = n_rep)
  out <- crossprod(mult, rep_stat_matrix)
  colnames(out) <- colnames(rep_stat_matrix)
  out
}

ci_quantile <- function(x, probs = c(0.025, 0.975)) stats::quantile(x, probs, na.rm = TRUE, names = FALSE)
