## =============================================================================
## final_analysis/R/helpers_floor_stability.R
##
## Shared helpers for the floor-stability analysis (R/21-R/23). New and
## additive -- nothing here touches R/12_floor_decomposition.R or
## R/helpers_stage2_truth.R. Two pieces of genuinely new logic this analysis
## needed that no existing script provided:
##   1. select_stability_floors(): phenotype-blind selection of the integer
##      size floor closest to a target % reduction in test count, per combo.
##   2. match_region_sets(): cross-floor region matching (core_snp-sharing +
##      physical overlap, or physical-overlap-only as a sensitivity variant),
##      used to decide whether a canonical (99.7%-floor) region is "recovered"
##      at a side floor (99.5% or 99.9%).
## Everything else this analysis needs (assemble_stage2(), score_stage2_
## regions(), stage2_seed_from_units()/_markers(), bootstrap_rep_matrix(),
## ci_quantile()) is reused UNCHANGED from R/helpers_stage2_truth.R.
## =============================================================================
suppressMessages({library(data.table)})

## ---- 1. phenotype-blind floor selection --------------------------------------
## `units1`: LDscnR:::.ld_outlier_units(stage1, map, 1L)'s output for ONE combo
## (every Stage-1 cluster, singletons included -- the same stable floor=1 base
## R/12 uses). `n_total_markers`: nrow(map) for that combo (the "all assayed
## markers" baseline the % reduction is relative to). `targets`: reduction
## fractions, e.g. c(0.995, 0.997, 0.999).
##
## "Test count" at floor f = number of ELIGIBLE STAGE-1 UNITS with
## n_markers >= f (the unit becomes one test once floor-filtered testing is
## used) -- reduction(f) = 1 - n_eligible(f) / n_total_markers. n_eligible(f)
## is a non-increasing step function of f; candidate floors are every integer
## from 1 to max(units1$n_markers) (past the largest cluster's size,
## n_eligible is 0 for every larger floor too, never closer to any target
## below 100%, so the range need not extend further).
##
## Ties (two floors equally close to a target) resolve to the SMALLER floor --
## an arbitrary but fixed, cluster-size-only rule, decided before any
## association result is examined, consistent with the phenotype-blind
## requirement (`which.min()` returns the first, i.e. smallest-floor, index
## on a tie).
##
## @return data.table: target, floor, achieved_reduction, n_eligible_units.
select_stability_floors <- function(units1, n_total_markers, targets = c(0.995, 0.997, 0.999)) {
  sizes <- units1$n_markers
  max_size <- if (length(sizes)) max(sizes) else 1L
  floors <- seq_len(max_size)
  ## n_eligible(f) for every candidate floor in one pass: a reverse
  ## cumulative count of cluster sizes, not a per-floor rescan of `sizes`.
  size_tab <- tabulate(sizes, nbins = max_size)               # count at EXACTLY size s
  n_eligible <- rev(cumsum(rev(size_tab)))                    # count >= f, for f = 1..max_size
  reduction <- 1 - n_eligible / n_total_markers

  rows <- lapply(targets, function(tgt) {
    best <- which.min(abs(reduction - tgt))                   # ties -> first (smallest) floor
    data.table(target = tgt, floor = floors[best],
              achieved_reduction = reduction[best], n_eligible_units = n_eligible[best])
  })
  rbindlist(rows)
}

## ---- 2. cross-floor region matching -------------------------------------------
## `canonical`/`side`: assemble_stage2()+score_stage2_regions() output
## (region_id, Chr, from, to, member_markers, member_units [core_snp list],
## TP, linked_qtn, ...) for the SAME combo+method, at the canonical (99.7%)
## floor and one side floor (99.5% or 99.9%) respectively.
##
## Primary criterion (use_core_snp = TRUE): >=1 shared core_snp (member_units)
## AND overlapping [Chr,from,to]. Sensitivity criterion (use_core_snp = FALSE):
## physical overlap alone, ignoring shared-unit membership entirely -- run
## both and compare, never silently pick one (PK's explicit request).
##
## @return one row per canonical region_id: recovered (logical), n_matches
##   (count of side regions satisfying the criterion -- >1 flags an ambiguous
##   one-to-many match), matched_side_region_ids (list-col, possibly empty).
match_region_sets <- function(canonical, side, use_core_snp = TRUE) {
  if (!nrow(canonical))
    return(data.table(region_id = integer(), recovered = logical(), n_matches = integer(),
                      matched_side_region_ids = list()))
  if (!nrow(side))
    return(data.table(region_id = canonical$region_id, recovered = FALSE, n_matches = 0L,
                      matched_side_region_ids = replicate(nrow(canonical), integer(0), simplify = FALSE)))

  rows <- vector("list", nrow(canonical))
  for (i in seq_len(nrow(canonical))) {
    cr <- canonical[i]
    overlap <- side$Chr == cr$Chr & side$from <= cr$to & cr$from <= side$to
    if (use_core_snp) {
      cr_units <- cr$member_units[[1]]
      shares_unit <- vapply(side$member_units, function(mu) any(mu %chin% cr_units), logical(1))
      hit <- overlap & shares_unit
    } else {
      hit <- overlap
    }
    matched_ids <- side$region_id[hit]
    rows[[i]] <- data.table(region_id = cr$region_id, recovered = length(matched_ids) > 0,
                            n_matches = length(matched_ids),
                            matched_side_region_ids = list(matched_ids))
  }
  rbindlist(rows)
}

## ---- 3. floor scheme (2026-10-07) -------------------------------------------
## FLOOR_SCHEME selects how the three floors are chosen. "pct" (default) is
## the original 99.5/99.7/99.9% test-reduction profile and keeps every output
## name unchanged. The density schemes make the CANONICAL floor the
## simulation analysis's fixed floor of 2 (PK 2026-10-07: the marker-density
## rule, max(2, round(markers per Mb / 250)), gives 2 in 567/600 c=1 runs and
## 3 in the 33 on the shortest map, so 2 is used throughout; the density value
## is still recorded per combo in `markers_per_mb`), and add two
## STRICTER side floors (a canonical floor of 2 has no admissible looser
## side): "density_x2" uses 2f and 4f, "density_x1.5" uses round(1.5f) and 2f.
## "density_half_double" (PK's primary design) is symmetric instead: f/2 and
## 2f, i.e. floors 1, 2 and 4 here. Floor 1 (singletons included) is only a
## comparison floor, never the reported one.
## The canonical floor sits in the floor_997 column and the two sides in
## floor_995 (less strict side) and floor_999 (stricter side), so R/22-R/23's
## tier logic is reused unchanged. Density-scheme outputs are written under
## floor_stability_<scheme>_* names and never overwrite the pct results.
FLOOR_SCHEME <- Sys.getenv("FLOOR_SCHEME", "pct")
stopifnot(FLOOR_SCHEME %in% c("pct", "density_half_double", "density_x2", "density_x1.5"))
DENSITY_MARKERS_PER_MB <- 250
fs_name <- function(x) {
  if (FLOOR_SCHEME == "pct") paste0("floor_stability_", x) else paste0("floor_stability_", FLOOR_SCHEME, "_", x)
}
fs_stage <- function(x) if (FLOOR_SCHEME == "pct") x else paste0(x, "_", FLOOR_SCHEME)

## Density-rule canonical floor plus two stricter sides, in the same long
## shape select_stability_floors() returns (`target` 0.995/0.997/0.999 is
## used here only as a role code: less-strict side / canonical / stricter side).
select_density_floors <- function(units1, map, scheme = FLOOR_SCHEME, canonical = 2L) {
  m <- data.table::as.data.table(map)
  span_mb <- m[, .(mb = (max(Pos) - min(Pos)) / 1e6), by = Chr][, sum(mb)]
  dens <- nrow(m) / span_mb
  f <- canonical
  sides <- switch(scheme, density_half_double = c(max(1L, as.integer(round(f / 2))), 2L * f),
                          density_x2 = c(2L * f, 4L * f),
                          `density_x1.5` = c(as.integer(round(1.5 * f)), 2L * f))
  fl <- c(sides[1], f, sides[2])
  sizes <- units1$n_markers
  data.table(target = c(0.995, 0.997, 0.999), floor = fl,
             achieved_reduction = vapply(fl, function(x) 1 - sum(sizes >= x) / nrow(m), numeric(1)),
             n_eligible_units = vapply(fl, function(x) sum(sizes >= x), integer(1)),
             markers_per_mb = dens)
}
