## =============================================================================
## module_3sp/R/helpers_floor_tiers.R
##
## Floor-stability tiers for an empirical panel under the marker-density floor
## rule (PK, 2026-10-07): the reported floor f = round(markers per Mb / 250) is
## compared with 0.5f and 2f, OBSERVED ONLY (the permutation null stays at f,
## in 03_EMMAX.R). Sourced by module_3sp/R/14_floor_tiers.R and
## module_9sp/R/04_floor_tiers.R.
##
## Cost: each statistic is fitted ONCE. Consensus p-values come from one
## emmax_fast() over the unit matrix at the LOWEST floor; a higher floor's
## units are a subset of those columns and each unit's EMMAX p-value does not
## depend on which other units are tested, so only BH and Stage-2 assembly
## are redone per floor. Simes uses one marker-level scan for every floor.
##
## Tier rule (mirrors the simulation analysis, final_analysis/R/23): a
## reported region is recovered at a side floor when at least one of its
## constituent significant Stage-1 units (matched by core_snp, which is
## floor-independent; unit_id is not) is also significant there AND the
## region physically overlaps a side-floor region. Tier = 1 + number of side
## floors recovering it (1/3, 2/3, 3/3).
## =============================================================================

## floor = round(markers per Mb / 250), minimum 2; density over the summed
## per-chromosome span (max - min Pos).
density_floor <- function(map, markers_per_mb = 250) {
  m <- data.table::as.data.table(map)
  span_mb <- m[, .(mb = (max(Pos) - min(Pos)) / 1e6), by = Chr][, sum(mb)]
  dens <- nrow(m) / span_mb
  list(markers_per_mb = dens, floor = max(2L, as.integer(round(dens / markers_per_mb))))
}

floor_tiers <- function(GTs, map, stage1, GRM, LD_decay, y, f, unit_repr, region_assembly, alpha) {
  floors <- c(lo = max(1L, as.integer(round(f / 2))), canon = as.integer(f), hi = 2L * as.integer(f))
  units_lo <- LDscnR:::.ld_outlier_units(stage1, map, floors[["lo"]])
  um <- ld_unit_matrix(GTs, stage1, map, size_floor = floors[["lo"]], repr = unit_repr)
  stopifnot("ld_unit_matrix columns must follow .ld_outlier_units() order" =
              identical(colnames(um), as.character(units_lo$unit_id)))
  con_p <- stats::setNames(emmax_fast(emmax_setup(um, GRM), y), units_lo$core_snp)
  pm <- emmax_fast(emmax_setup(GTs, GRM), y)

  run_one <- function(stat, fl) {
    uf <- LDscnR:::.ld_outlier_units(stage1, map, fl)
    p <- if (stat == "unit") unname(con_p[uf$core_snp]) else pm
    tt <- ld_outlier_test(stage1, map, p, statistic = stat, size_floor = fl, alpha = alpha,
                          assembly = "stage2_discovered", GTs = GTs, LD_decay = LD_decay,
                          score_threshold = region_assembly$score_threshold,
                          distance_threshold = region_assembly$distance_threshold)
    u <- data.table::copy(tt$units)
    stopifnot(identical(u$unit_id, uf$unit_id))
    u[, core_snp := uf$core_snp]
    list(test = tt, sig = u[significant == TRUE], regions = data.table::as.data.table(tt$regions))
  }

  out <- list()
  for (stat in c("unit", "simes")) {
    res <- lapply(floors, function(fl) run_one(stat, fl)); names(res) <- names(floors)
    reg <- data.table::copy(res$canon$regions)
    if (!nrow(reg)) next
    reg[, region := sprintf("%s:%d-%d", Chr, as.integer(from), as.integer(to))]
    recovered <- function(side) vapply(seq_len(nrow(reg)), function(i) {
      r <- reg[i]
      cores <- res$canon$sig[Chr == r$Chr & to >= r$from & from <= r$to, core_snp]
      sr <- res[[side]]$regions
      any(cores %in% res[[side]]$sig$core_snp) &&
        nrow(sr[Chr == r$Chr & to >= r$from & from <= r$to]) > 0
    }, logical(1))
    reg[, `:=`(statistic = stat, floor_lo = floors[["lo"]], floor = floors[["canon"]], floor_hi = floors[["hi"]],
               found_lo = recovered("lo"), found_hi = recovered("hi"))]
    reg[, tier := 1L + found_lo + found_hi]
    ## every floor's regions tiered against the OTHER two floors (same rule):
    ## which floor do single-floor survivors come from?
    all_fl <- rbindlist(lapply(names(floors), function(base) {
      rb <- data.table::copy(res[[base]]$regions); if (!nrow(rb)) return(NULL)
      n_other <- vapply(seq_len(nrow(rb)), function(i) {
        r <- rb[i]; cores <- res[[base]]$sig[Chr == r$Chr & to >= r$from & from <= r$to, core_snp]
        sum(vapply(setdiff(names(floors), base), function(side) {
          sr <- res[[side]]$regions
          any(cores %in% res[[side]]$sig$core_snp) && nrow(sr[Chr == r$Chr & to >= r$from & from <= r$to]) > 0
        }, logical(1)))
      }, numeric(1))
      rb[, `:=`(statistic = stat, base = base, floor = floors[[base]], n_floors = 1L + as.integer(n_other))]
    }), fill = TRUE)
    out[[stat]] <- list(regions = reg, all_floors = all_fl, tests = lapply(res, `[[`, "test"))
  }
  out$floors <- floors
  out
}
