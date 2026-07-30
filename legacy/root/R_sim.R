######################################################
## OR-level TP/FP/FN classification for the pooled
## simulated genome, following the focal-QTN assignment
## protocol (Va-based true positives, LD/distance-gated
## focal QTN assignment, duplicate resolution, satellite-
## cluster / over-merging behaviour handled implicitly).
##
## Requires: outlier_regions_sim.R (precompute_LD_edges,
## LD_igraph_components, get_potential_outliers_sim,
## run_one_grid_sim) to already be sourced/loaded.
######################################################

library(data.table)

#----------------------------------------------------------
# 1) True-positive QTN definition
#    p_Va,l = beta_l^2 * p_l*(1-p_l) / sum(...), among
#    QTNs with MAF > 0.1 on the SAME QTN chromosome.
#    TP QTN if p_Va >= 0.05.
#----------------------------------------------------------
flag_true_positive_QTNs <- function(map, va_col = "Va", maf_col = "MAF",
                                    maf_min = 0.1, p_va_min = 0.05) {
  map <- copy(map)
  map[, true_pos_QTN := FALSE]

  qtn_ok <- map$type == "QTN" & map[[maf_col]] > maf_min

  map[qtn_ok, sum_Va_ok := sum(get(va_col)), by = Chr_9sp]
  map[qtn_ok, p_Va_ok := get(va_col) / sum_Va_ok]
  map[qtn_ok & p_Va_ok >= p_va_min, true_pos_QTN := TRUE]

  map
}

#----------------------------------------------------------
# 2) Precompute LD (r^2) and distance between every
#    candidate outlier marker and every QTN on the same
#    chromosome (restricted to markers that ever appear as
#    outliers, for efficiency; max_bp bounds the search
#    window generously above d_{rho=0.95}).
#----------------------------------------------------------
precompute_QTN_LD <- function(GTs, map, candidate_markers, max_bp = 2e6, cores = 1) {

  map_sub <- map[marker %in% candidate_markers]
  chr_levels <- unique(map_sub$Chr_9sp)

  out <- parallel::mclapply(chr_levels, function(ch) {

    chr_map      <- map[Chr_9sp == ch]
    qtn_markers  <- chr_map[type == "QTN", marker]
    cand_markers <- chr_map[marker %in% candidate_markers, marker]

    if (length(qtn_markers) == 0 || length(cand_markers) == 0) return(NULL)

    gts_qtn <- as.matrix(GTs[, qtn_markers, drop = FALSE])
    gts_cnd <- as.matrix(GTs[, cand_markers, drop = FALSE])
    storage.mode(gts_qtn) <- "double"
    storage.mode(gts_cnd) <- "double"

    R2 <- cor(gts_cnd, gts_qtn, use = "pairwise.complete.obs")^2  # cand x qtn

    pos_cnd <- setNames(chr_map$Pos[match(cand_markers, chr_map$marker)], cand_markers)
    pos_qtn <- setNames(chr_map$Pos[match(qtn_markers, chr_map$marker)], qtn_markers)

    dt <- data.table(
      marker     = rep(cand_markers, times = length(qtn_markers)),
      qtn_marker = rep(qtn_markers, each = length(cand_markers)),
      r2         = as.vector(R2),
      dist_bp    = abs(rep(pos_cnd, times = length(qtn_markers)) -
                         rep(pos_qtn, each = length(cand_markers)))
    )
    dt[is.finite(max_bp)][dist_bp <= max_bp]

  }, mc.cores = cores)

  rbindlist(out[!vapply(out, is.null, logical(1))], use.names = TRUE)
}

#----------------------------------------------------------
# 3) Assign a focal QTN to one OR (cluster of marker names),
#    resolving multi-candidate ties via LD with the OR's
#    non-QTN ("neutral" outlier) members.
#    Returns the assigned qtn_marker (or NA) and the r2
#    evidence used for the assignment (for duplicate
#    resolution across ORs afterwards).
#----------------------------------------------------------
assign_focal_qtn_one_OR <- function(cluster_markers, qtn_ld_table,
                                    qtn_marker_set, r2_min_focal, d_max_focal) {

  rows <- qtn_ld_table[marker %in% cluster_markers &
                         r2 > r2_min_focal & dist_bp < d_max_focal]

  if (nrow(rows) == 0) return(list(qtn = NA_character_, evidence = NA_real_))

  ## Stage 1: candidate QTNs = any QTN meeting r2/distance criteria
  ## with ANY SNP in the OR
  best_per_qtn <- rows[, .(max_r2 = max(r2)), by = qtn_marker]

  if (nrow(best_per_qtn) == 1) {
    return(list(qtn = best_per_qtn$qtn_marker[1], evidence = best_per_qtn$max_r2[1]))
  }

  ## Stage 2: tie-break using LD with non-QTN members of the OR only
  neutral_markers <- setdiff(cluster_markers, qtn_marker_set)

  if (length(neutral_markers) > 0) {
    rows_neutral <- rows[marker %in% neutral_markers]
    if (nrow(rows_neutral) > 0) {
      best_neutral <- rows_neutral[, .(max_r2 = max(r2)), by = qtn_marker]
      setorder(best_neutral, -max_r2)
      winner <- best_neutral$qtn_marker[1]
      return(list(qtn = winner, evidence = best_per_qtn[qtn_marker == winner, max_r2]))
    }
  }

  ## Fall back: highest overall LD among the tied candidates
  setorder(best_per_qtn, -max_r2)
  list(qtn = best_per_qtn$qtn_marker[1], evidence = best_per_qtn$max_r2[1])
}

#----------------------------------------------------------
# 4) Evaluate one set of ORs (one grid point x one p_name):
#    assign focal QTNs, resolve duplicates (best-evidence OR
#    keeps the QTN, others become unassigned -> FP), then
#    compute TP / FP / FN / Precision / Recall / PR.
#----------------------------------------------------------
evaluate_ORs <- function(clusters, map, qtn_ld_table,
                         r2_min_focal, d_max_focal) {

  qtn_marker_set   <- map[type == "QTN", marker]
  true_pos_markers <- map[true_pos_QTN == TRUE, marker]

  if (length(clusters) == 0) {
    return(list(
      TP = 0L, FP = 0L, FN = length(true_pos_markers),
      Precision = NA_real_, Recall = 0, PR = 0,
      assignments = data.table(CL_id = integer(0), qtn = character(0), evidence = numeric(0))
    ))
  }

  assign_list <- lapply(clusters, assign_focal_qtn_one_OR,
                        qtn_ld_table = qtn_ld_table,
                        qtn_marker_set = qtn_marker_set,
                        r2_min_focal = r2_min_focal, d_max_focal = d_max_focal)

  assignments <- data.table(
    CL_id    = seq_along(clusters),
    qtn      = vapply(assign_list, function(x) x$qtn, character(1)),
    evidence = vapply(assign_list, function(x) x$evidence, numeric(1))
  )

  ## resolve duplicates: same QTN assigned to >1 OR -> keep highest
  ## evidence OR, unassign (NA) the rest (they become FP as satellites)
  assigned <- assignments[!is.na(qtn)]
  if (nrow(assigned) > 0) {
    setorder(assigned, qtn, -evidence)
    keep_id <- assigned[, .SD[1], by = qtn]$CL_id
    drop_id <- setdiff(assigned$CL_id, keep_id)
    assignments[CL_id %in% drop_id, qtn := NA_character_]
  }

  TP_markers <- assignments[!is.na(qtn) & qtn %in% true_pos_markers, qtn]
  TP <- length(TP_markers)
  FP <- sum(is.na(assignments$qtn)) + sum(!is.na(assignments$qtn) & !(assignments$qtn %in% true_pos_markers))
  FN <- length(setdiff(true_pos_markers, TP_markers))

  Precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  Recall    <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  PR        <- if (!is.na(Precision) && !is.na(Recall)) Precision * Recall else NA_real_

  list(TP = TP, FP = FP, FN = FN, Precision = Precision, Recall = Recall, PR = PR,
       assignments = assignments)
}

#----------------------------------------------------------
# 5) Apply evaluate_ORs() across the full outliers_dt grid
#    (output of run_one_grid_sim, rbind'ed over the param grid)
#----------------------------------------------------------
score_outlier_grid <- function(outliers_dt, map, qtn_ld_table, p_names,
                               r2_min_focal, d_max_focal) {

  dt <- copy(outliers_dt)

  for (nm in p_names) {
    res <- lapply(dt[[nm]], evaluate_ORs,
                  map = map, qtn_ld_table = qtn_ld_table,
                  r2_min_focal = r2_min_focal, d_max_focal = d_max_focal)

    dt[, paste0("TP_", nm)        := vapply(res, `[[`, numeric(1), "TP")]
    dt[, paste0("FP_", nm)        := vapply(res, `[[`, numeric(1), "FP")]
    dt[, paste0("FN_", nm)        := vapply(res, `[[`, numeric(1), "FN")]
    dt[, paste0("Precision_", nm) := vapply(res, `[[`, numeric(1), "Precision")]
    dt[, paste0("Recall_", nm)    := vapply(res, `[[`, numeric(1), "Recall")]
    dt[, paste0("PR_", nm)        := vapply(res, `[[`, numeric(1), "PR")]
  }

  dt[]
}

#----------------------------------------------------------
# Example usage
#----------------------------------------------------------
# source("outlier_regions_sim.R")
#
# conc_data <- readRDS("./parsed_sim_data_genomes/adapt_bgs_V2_c1_env1_genome.rds")
# map   <- conc_data$map
# GTs   <- conc_data$GTs
# ld_ws <- conc_data$ld_ws
#
# map[, sim_id := ceiling(Chr_9sp / 2)]
# map <- flag_true_positive_QTNs(map)          # true_pos_QTN column added
#
# p_cols <- c(EMX = "emx_p")
# th_ldw_grid <- c(0, 0.5, 0.75, 0.9, 0.95)
# r2_grid     <- seq(0.6, 0.9, by = 0.1)
# lmin_grid   <- c(5, 10, 20)
# bp_th_cluster <- 5e5                          # fixed 500kb clustering distance (tau_d)
#
# potential_outliers <- get_potential_outliers_sim(
#   map, ld_ws, th_ldw_grid, p_cols, group_col = "sim_id", alpha = 0.05
# )
#
# el_potential <- precompute_LD_edges(
#   GTs = GTs[, potential_outliers, drop = FALSE],
#   map = map[marker %in% potential_outliers],
#   r2_min = 0.1, max_bp = bp_th_cluster, cores = 1
# )
#
# qtn_ld_table <- precompute_QTN_LD(
#   GTs = GTs, map = map, candidate_markers = potential_outliers,
#   max_bp = 2e6, cores = 1
# )
#
# ## focal-QTN assignment thresholds: r2 at rho=0.75, distance at rho=0.95
# ## -- EITHER from decay params (if you load LD_decay$decay_sum):
# # r2_min_focal <- median(ld_from_rho(b = decay_sum$b, c = decay_sum$c, rho = 0.75))
# # d_max_focal  <- median(d_from_rho(a = decay_sum$a, rho = 0.95))
# ## -- OR fixed placeholders (adjust as appropriate):
# r2_min_focal <- 0.2
# d_max_focal  <- 5e5
#
# param_grid <- CJ(rho = colnames(ld_ws), th_ldw = th_ldw_grid)
#
# outliers_persim <- rbindlist(lapply(seq_len(nrow(param_grid)), function(i) {
#   pars <- param_grid[i]
#   run_one_grid_sim(
#     map = map, el = el_potential, ld_ws = ld_ws,
#     rho = pars$rho, th_ldw = pars$th_ldw,
#     p_cols = p_cols, group_col = "sim_id", alpha = 0.05,
#     r2_grid = r2_grid, lmin_grid = lmin_grid,
#     bp_th = bp_th_cluster, cores = 1
#   )
# }), fill = TRUE)
#
# scored <- score_outlier_grid(
#   outliers_persim, map, qtn_ld_table, names(p_cols),
#   r2_min_focal = r2_min_focal, d_max_focal = d_max_focal
# )
#
# scored[r2_th == 0.7 & l_min == 10,
#        .(rho, th_ldw, TP_EMX, FP_EMX, FN_EMX, Precision_EMX, Recall_EMX, PR_EMX)]
