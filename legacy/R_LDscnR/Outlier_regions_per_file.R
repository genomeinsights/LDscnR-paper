######################################################
## Per-simulation-file outlier-region (OR) pipeline
##
## Works directly on the individual parsed_sim_data/*.rds
## files (GTs, map, env, LD_decay, ld_ws) -- no concatenation
## needed. Each file already represents one joint analysis
## (emx_p / lfmm_p computed once across both its chromosomes),
## so filtering + FDR correction is a single global step, and
## LD_decay is used directly to get real decay-based focal-QTN
## thresholds instead of placeholders.
######################################################

library(data.table)
library(igraph)
library(parallel)
library(ggplot2)
library(wesanderson)

#----------------------------------------------------------
# LD edge precomputation (within-chromosome pairs only)
#----------------------------------------------------------
precompute_LD_edges <- function(GTs, map, r2_min = 0.1, max_bp = Inf, cores = 1) {
  markers <- intersect(colnames(GTs), map$marker)
  if (length(markers) == 0) stop("No overlapping markers between GTs and map.")

  map_sub <- copy(map[marker %in% markers])
  setkey(map_sub, Chr, Pos)
  chr_levels <- unique(map_sub$Chr)

  out <- mclapply(chr_levels, function(ch) {
    chr_map <- map_sub[Chr == ch]
    chr_markers <- chr_map$marker

    if (length(chr_markers) < 2) {
      return(data.table(Chr = ch, marker1 = chr_markers, marker2 = chr_markers, r2 = 1, dist_bp = 0))
    }

    gts <- as.matrix(GTs[, chr_markers, drop = FALSE])
    storage.mode(gts) <- "double"

    R2 <- cor(gts, use = "pairwise.complete.obs")^2
    R2[is.na(R2)] <- 0
    diag(R2) <- 0

    idx <- which(R2 >= r2_min, arr.ind = TRUE)
    idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]

    if (nrow(idx) == 0) {
      return(data.table(Chr = ch, marker1 = chr_markers, marker2 = chr_markers, r2 = 1, dist_bp = 0))
    }

    pos <- chr_map$Pos
    dt <- data.table(
      Chr = ch,
      marker1 = chr_markers[idx[, 1]],
      marker2 = chr_markers[idx[, 2]],
      r2 = R2[idx],
      dist_bp = abs(pos[idx[, 1]] - pos[idx[, 2]])
    )
    if (is.finite(max_bp)) dt <- dt[dist_bp <= max_bp]
    dt
  }, mc.cores = cores)

  out <- rbindlist(out, use.names = TRUE, fill = TRUE)
  setkey(out, Chr, marker1, marker2)
  out
}

#----------------------------------------------------------
# LD clusters from edge list
#----------------------------------------------------------
LD_igraph_components <- function(el, markers, r2_th = 0.8, bp_th = Inf) {
  markers <- unique(markers)
  if (length(markers) == 0) return(data.table(marker = character(), CL_id = integer(), n_loci = integer()))
  if (length(markers) == 1) return(data.table(marker = markers, CL_id = 1L, n_loci = 1L))

  edges <- el[marker1 %in% markers & marker2 %in% markers & r2 >= r2_th]
  if (is.finite(bp_th)) edges <- edges[dist_bp <= bp_th]

  if (nrow(edges) == 0) return(data.table(marker = markers, CL_id = seq_along(markers), n_loci = 1L))

  g <- graph_from_data_frame(edges[, .(from = marker1, to = marker2)], directed = FALSE,
                             vertices = data.table(name = markers))
  comp <- components(g)
  clusters <- data.table(marker = names(comp$membership), CL_id = as.integer(comp$membership))
  clusters[, n_loci := comp$csize[CL_id]]
  clusters
}

#----------------------------------------------------------
# Potential outliers across the whole (rho x th_ldw) grid
# -- single global FDR correction (no replicate grouping needed)
#----------------------------------------------------------
get_potential_outliers <- function(map, ld_ws, th_ldw_grid, p_cols, alpha_grid = 0.05) {
  ## q < alpha is monotonic in alpha, so the set of candidates at the
  ## LOOSEST (largest) alpha in the grid is a strict superset of the
  ## candidates at every smaller alpha -- no need to loop alpha here,
  ## just use the max once. Individual alpha values are applied properly
  ## later, in run_one_grid()'s own outlier-calling step.
  alpha_max <- max(alpha_grid)

  common_markers <- intersect(map$marker, rownames(ld_ws))
  map_sub <- map[marker %in% common_markers]
  potential <- character()

  for (rho in colnames(ld_ws)) {
    message("processing rho = ", rho)
    ld_vec <- ld_ws[map_sub$marker, rho]

    for (th_ldw in th_ldw_grid) {
      if (th_ldw > 0) {
        keep <- ld_vec > quantile(ld_vec, th_ldw, na.rm = TRUE)
        keep[is.na(keep)] <- FALSE
      } else {
        keep <- rep(TRUE, length(ld_vec))
      }
      if (!any(keep)) next

      for (p_col in p_cols) {
        q <- p.adjust(map_sub[[p_col]][keep], method = "fdr")
        potential <- c(potential, map_sub$marker[keep][q < alpha_max])
      }
    }
  }
  unique(potential)
}

#----------------------------------------------------------
# Empty result helper
#----------------------------------------------------------
empty_result <- function(rho, th_ldw, alpha, p_names, r2_grid, lmin_grid) {
  out <- CJ(r2_th = r2_grid, l_min = lmin_grid)
  out[, `:=`(th_ldw = th_ldw, rho = rho, alpha = alpha)]
  for (nm in p_names) out[, (nm) := list(list(character()))]
  out[]
}

#----------------------------------------------------------
# One grid point (rho, th_ldw) -- single global FDR correction
#----------------------------------------------------------
run_one_grid <- function(map, el = NULL, ld_ws, rho, th_ldw, p_cols,
                         p_names = names(p_cols), alpha = 0.05,
                         r2_grid, lmin_grid, bp_th = Inf, cores = 1) {
  stopifnot(length(p_cols) == length(p_names))

  common_markers <- intersect(map$marker, rownames(ld_ws))
  map_sub <- copy(map[marker %in% common_markers])
  ld_sub  <- ld_ws[map_sub$marker, , drop = FALSE]

  if (th_ldw > 0) {
    keep <- ld_sub[, rho] > quantile(ld_sub[, rho], th_ldw, na.rm = TRUE)
    keep[is.na(keep)] <- FALSE
  } else {
    keep <- rep(TRUE, nrow(ld_sub))
  }

  if (!any(keep)) return(cbind(empty_result(rho, th_ldw, alpha, p_names, r2_grid, lmin_grid), n_loci = 0))

  markers_keep <- map_sub[keep, marker]
  outliers <- setNames(vector("list", length(p_cols)), p_names)

  for (i in seq_along(p_cols)) {
    p_col <- p_cols[i]
    nm <- p_names[i]
    q <- p.adjust(unlist(map_sub[keep, ..p_col]), method = "fdr")
    outliers[[nm]] <- markers_keep[q < alpha]
  }

  if (length(unique(unlist(outliers))) == 0) {
    return(cbind(empty_result(rho, th_ldw, alpha, p_names, r2_grid, lmin_grid), n_loci = length(which(keep))))
  }

  if (is.null(el)) {
    all_outliers <- unique(unlist(outliers))
    el <- precompute_LD_edges(GTs = GTs[, all_outliers, drop = FALSE],
                              map = map_sub[marker %in% all_outliers],
                              r2_min = 0.1, max_bp = 1e6, cores = 1)
  }

  out <- rbindlist(mclapply(r2_grid, function(r2_th) {
    clusters <- lapply(outliers, function(markers) {
      LD_igraph_components(el = el, markers = markers, r2_th = r2_th, bp_th = bp_th)
    })
    rbindlist(lapply(lmin_grid, function(l_min) {
      row <- data.table(r2_th = r2_th, l_min = l_min, th_ldw = th_ldw, rho = rho, alpha = alpha)
      for (nm in p_names) {
        tmp <- clusters[[nm]][n_loci >= l_min, ]
        cls <- split(tmp$marker, tmp$CL_id)
        if (length(cls) == 0) {
          row[, (nm) := list(list(character(0)))]
        } else if (length(cls) == 1) {
          ## data.table's `:=` unwraps a length-1 list when nrow==1==length(cls),
          ## storing the raw marker vector instead of a length-1 list containing
          ## it -- double-wrap to force it to stay nested as ONE cluster.
          row[, (nm) := list(list(cls))]
        } else {
          row[, (nm) := list(cls)]
        }
      }
      row
    }), fill = TRUE)
  }, mc.cores = cores), fill = TRUE)

  out[, n_loci := length(which(keep))]
}

#----------------------------------------------------------
# True-positive QTN definition: p_Va >= 0.05 among QTNs
# with MAF > 0.1 on the SAME chromosome
#----------------------------------------------------------
flag_true_positive_QTNs <- function(map, va_col = "Va", maf_col = "MAF",
                                    maf_min = 0.1, p_va_min = 0.05) {
  map <- copy(map)
  map[, true_pos_QTN := FALSE]
  qtn_ok <- map$type == "QTN" & map[[maf_col]] > maf_min

  map[qtn_ok, sum_Va_ok := sum(get(va_col)), by = Chr]
  map[qtn_ok, p_Va_ok := get(va_col) / sum_Va_ok]
  map[qtn_ok & p_Va_ok >= p_va_min, true_pos_QTN := TRUE]
  map
}

#----------------------------------------------------------
# Focal-QTN thresholds from a file's own LD_decay object:
# r2 at rho=0.75, distance at rho=0.95, for the QTN chromosome
#----------------------------------------------------------
get_focal_qtn_thresholds <- function(LD_decay, chr_qtn = "Chr1", rho_r2 = 0.75, rho_d = 0.95) {
  ds <- LD_decay$decay_sum[Chr == chr_qtn]
  if (nrow(ds) != 1) {
    stop("Expected exactly one decay_sum row for chr '", chr_qtn,
         "', found ", nrow(ds), ". Check LD_decay$decay_sum$Chr labels.")
  }
  list(
    r2_min_focal = ld_from_rho(b = ds$b, c = ds$c, rho = rho_r2),
    d_max_focal  = d_from_rho(a = ds$a, rho = rho_d)
  )
}

## these two come from compute_ld_structure.R -- re-declared here for
## standalone use if that script isn't sourced
if (!exists("ld_from_rho")) {
  ld_from_rho <- function(b, c = 1, rho) b + (c - b) * (1 - rho)
}
if (!exists("d_from_rho")) {
  d_from_rho <- function(a, rho) rho / (a * (1 - rho))
}

#----------------------------------------------------------
# Precompute LD (r^2) + distance between candidate outlier
# markers and every QTN on the same chromosome
#----------------------------------------------------------
precompute_QTN_LD <- function(GTs, map, candidate_markers, max_bp = 2e6, cores = 1) {
  chr_levels <- unique(map[marker %in% candidate_markers, Chr])

  out <- mclapply(chr_levels, function(ch) {
    chr_map      <- map[Chr == ch]
    qtn_markers  <- chr_map[type == "QTN", marker]
    cand_markers <- chr_map[marker %in% candidate_markers, marker]
    if (length(qtn_markers) == 0 || length(cand_markers) == 0) return(NULL)

    gts_qtn <- as.matrix(GTs[, qtn_markers, drop = FALSE])
    gts_cnd <- as.matrix(GTs[, cand_markers, drop = FALSE])
    storage.mode(gts_qtn) <- "double"
    storage.mode(gts_cnd) <- "double"

    R2 <- cor(gts_cnd, gts_qtn, use = "pairwise.complete.obs")^2

    pos_cnd <- setNames(chr_map$Pos[match(cand_markers, chr_map$marker)], cand_markers)
    pos_qtn <- setNames(chr_map$Pos[match(qtn_markers, chr_map$marker)], qtn_markers)

    dt <- data.table(
      marker     = rep(cand_markers, times = length(qtn_markers)),
      qtn_marker = rep(qtn_markers, each = length(cand_markers)),
      r2         = as.vector(R2),
      dist_bp    = abs(rep(pos_cnd, times = length(qtn_markers)) -
                         rep(pos_qtn, each = length(cand_markers)))
    )
    dt[dist_bp <= max_bp]
  }, mc.cores = cores)

  rbindlist(out[!vapply(out, is.null, logical(1))], use.names = TRUE)
}

#----------------------------------------------------------
# Assign a focal QTN to one OR, with the two-stage tie-break
#----------------------------------------------------------
assign_focal_qtn_one_OR <- function(cluster_markers, qtn_ld_table, qtn_marker_set,
                                    r2_min_focal, d_max_focal) {
  rows <- qtn_ld_table[marker %in% cluster_markers & r2 > r2_min_focal & dist_bp < d_max_focal]
  if (nrow(rows) == 0) return(list(qtn = NA_character_, evidence = NA_real_))

  best_per_qtn <- rows[, .(max_r2 = max(r2)), by = qtn_marker]
  if (nrow(best_per_qtn) == 1) {
    return(list(qtn = best_per_qtn$qtn_marker[1], evidence = best_per_qtn$max_r2[1]))
  }

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

  setorder(best_per_qtn, -max_r2)
  list(qtn = best_per_qtn$qtn_marker[1], evidence = best_per_qtn$max_r2[1])
}

#----------------------------------------------------------
# Classify each OR (cluster) individually: assign focal QTN,
# resolve duplicates, flag TP/FP -- one row per cluster.
# evaluate_ORs() below just summarizes this to counts; the
# per-cluster table itself is what the Manhattan plot needs.
#----------------------------------------------------------
classify_ORs <- function(clusters, map, qtn_ld_table, r2_min_focal, d_max_focal) {

  qtn_marker_set   <- map[type == "QTN", marker]
  true_pos_markers <- map[true_pos_QTN == TRUE, marker]

  if (length(clusters) == 0) {
    return(data.table(CL_id = integer(0), n_loci = integer(0),
                      qtn = character(0), evidence = numeric(0), is_TP = logical(0)))
  }

  assign_list <- lapply(clusters, assign_focal_qtn_one_OR,
                        qtn_ld_table = qtn_ld_table, qtn_marker_set = qtn_marker_set,
                        r2_min_focal = r2_min_focal, d_max_focal = d_max_focal)

  assignments <- data.table(
    CL_id    = seq_along(clusters),
    n_loci   = vapply(clusters, length, integer(1)),
    qtn      = vapply(assign_list, function(x) x$qtn, character(1)),
    evidence = vapply(assign_list, function(x) x$evidence, numeric(1))
  )

  assigned <- assignments[!is.na(qtn)]
  if (nrow(assigned) > 0) {
    setorder(assigned, qtn, -evidence)
    keep_id <- assigned[, .SD[1], by = qtn]$CL_id
    drop_id <- setdiff(assigned$CL_id, keep_id)
    assignments[CL_id %in% drop_id, qtn := NA_character_]
  }

  assignments[, is_TP := !is.na(qtn) & qtn %in% true_pos_markers]
  assignments[]
}

#----------------------------------------------------------
# DIAGNOSTIC: like classify_ORs(), but keeps the pre-dedup
# candidate QTN for every cluster (even ones that lost
# duplicate resolution), whether that candidate itself clears
# the Va-based true-positive threshold, and which cluster beat
# it and by how much. Use this to tell apart:
#   (a) "correct" satellite-cluster suppression -- candidate_qtn
#       is a real true-positive QTN, but this cluster lost to a
#       stronger nearby cluster pointing at the SAME QTN and
#       physically close by (over-clustering artifact)
#   (b) candidate_qtn is NOT a true-positive QTN at all (Va < 5%
#       threshold) -- correctly FP regardless of clustering
#   (c) two DIFFERENT, physically separated peaks both matched
#       to the same QTN and competing -- check dist_bp for the
#       winning vs losing evidence to see if d_max_focal is too
#       permissive (over-merging distinct signals as duplicates)
#----------------------------------------------------------
diagnose_OR_classification <- function(clusters, map, qtn_ld_table, r2_min_focal, d_max_focal) {

  qtn_marker_set   <- map[type == "QTN", marker]
  true_pos_lookup  <- setNames(map$true_pos_QTN, map$marker)
  pos_lookup       <- setNames(map$Pos, map$marker)

  if (length(clusters) == 0) return(data.table())

  assign_list <- lapply(clusters, assign_focal_qtn_one_OR,
                        qtn_ld_table = qtn_ld_table, qtn_marker_set = qtn_marker_set,
                        r2_min_focal = r2_min_focal, d_max_focal = d_max_focal)

  dt <- data.table(
    CL_id         = seq_along(clusters),
    n_loci        = vapply(clusters, length, integer(1)),
    cluster_pos   = vapply(clusters, function(cl) median(pos_lookup[cl], na.rm = TRUE), numeric(1)),
    candidate_qtn = vapply(assign_list, function(x) x$qtn, character(1)),
    evidence      = vapply(assign_list, function(x) x$evidence, numeric(1))
  )
  dt[, candidate_qtn_pos := pos_lookup[candidate_qtn]]
  dt[, candidate_qtn_is_true_positive := true_pos_lookup[candidate_qtn]]
  dt[, dist_to_candidate_qtn := abs(cluster_pos - candidate_qtn_pos)]

  dt[, final_qtn := candidate_qtn]
  dt[, dropped_by_dedup := FALSE]
  dt[, beaten_by_CL_id := NA_integer_]
  dt[, beaten_by_evidence := NA_real_]

  assigned <- dt[!is.na(candidate_qtn)]
  if (nrow(assigned) > 0) {
    setorder(assigned, candidate_qtn, -evidence)
    winners <- assigned[, .SD[1], by = candidate_qtn]
    keep_id <- winners$CL_id
    drop_id <- setdiff(assigned$CL_id, keep_id)

    dt[CL_id %in% drop_id, `:=`(final_qtn = NA_character_, dropped_by_dedup = TRUE)]

    winner_lookup <- winners[, .(candidate_qtn, beaten_by_CL_id = CL_id, beaten_by_evidence = evidence)]
    dt <- merge(dt, winner_lookup, by = "candidate_qtn", all.x = TRUE, suffixes = c("", ".win"))
    dt[dropped_by_dedup == TRUE, `:=`(beaten_by_CL_id = beaten_by_CL_id.win,
                                      beaten_by_evidence = beaten_by_evidence.win)]
    dt[, c("beaten_by_CL_id.win", "beaten_by_evidence.win") := NULL]
  }

  dt[, is_TP := !is.na(final_qtn) & true_pos_lookup[final_qtn]]
  setorder(dt, CL_id)
  dt[]
}

#----------------------------------------------------------
# Evaluate one set of ORs: compute TP / FP / FN / Precision /
# Recall / PR from classify_ORs()
#----------------------------------------------------------
evaluate_ORs <- function(clusters, map, qtn_ld_table, r2_min_focal, d_max_focal) {
  true_pos_markers <- map[true_pos_QTN == TRUE, marker]

  assignments <- classify_ORs(clusters, map, qtn_ld_table, r2_min_focal, d_max_focal)

  if (nrow(assignments) == 0) {
    return(list(TP = 0L, FP = 0L, FN = length(true_pos_markers),
                Precision = NA_real_, Recall = 0, PR = 0, assignments = assignments))
  }

  TP_markers <- assignments[is_TP == TRUE, qtn]
  TP <- length(TP_markers)
  FP <- sum(!assignments$is_TP)
  FN <- length(setdiff(true_pos_markers, TP_markers))

  Precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  Recall    <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  PR        <- if (!is.na(Precision) && !is.na(Recall)) Precision * Recall else NA_real_

  list(TP = TP, FP = FP, FN = FN, Precision = Precision, Recall = Recall, PR = PR,
       assignments = assignments)
}

#----------------------------------------------------------
# Score a full outliers_dt grid (rbind'ed over param_grid).
# Retains a per-grid-point "assignments_<nm>" list-column
# (cluster-level TP/FP labels) so plot_OR_manhattan() can
# pull directly from the scored table without recomputing.
#----------------------------------------------------------
score_outlier_grid <- function(outliers_dt, map, qtn_ld_table, p_names, r2_min_focal, d_max_focal) {
  dt <- copy(outliers_dt)
  for (nm in p_names) {
    res <- lapply(dt[[nm]], evaluate_ORs, map = map, qtn_ld_table = qtn_ld_table,
                  r2_min_focal = r2_min_focal, d_max_focal = d_max_focal)
    dt[, paste0("TP_", nm)        := vapply(res, `[[`, numeric(1), "TP")]
    dt[, paste0("FP_", nm)        := vapply(res, `[[`, numeric(1), "FP")]
    dt[, paste0("FN_", nm)        := vapply(res, `[[`, numeric(1), "FN")]
    dt[, paste0("Precision_", nm) := vapply(res, `[[`, numeric(1), "Precision")]
    dt[, paste0("Recall_", nm)    := vapply(res, `[[`, numeric(1), "Recall")]
    dt[, paste0("PR_", nm)        := vapply(res, `[[`, numeric(1), "PR")]
    dt[, paste0("assignments_", nm) := lapply(res, `[[`, "assignments")]
  }
  dt[]
}

#----------------------------------------------------------
# Annotate a map with call status (not called / FP OR / TP OR)
# for one grid point's clusters+assignments -- no plotting.
# Factored out so plot_OR_manhattan() and the multi-chromosome
# genome-wide version can share the same annotation logic.
#----------------------------------------------------------
annotate_OR_calls <- function(map, value_col, clusters, assignments, ld_ws = NULL, rho = NULL,
                              th_ldw = NULL, p_col = NULL) {

  map_sub <- copy(map)

  if (identical(value_col, "ld_w")) {
    if (is.null(ld_ws) || is.null(rho)) {
      stop("value_col = 'ld_w' requires both `ld_ws` and `rho` to be supplied.")
    }
    rho <- as.character(rho)
    if (!rho %in% colnames(ld_ws)) {
      stop("rho = '", rho, "' not found in colnames(ld_ws): ",
           paste(colnames(ld_ws), collapse = ", "))
    }
    ## match by marker name via match(), not direct rowname subsetting --
    ## robust to markers that aren't present in ld_ws (e.g. filtered out
    ## upstream, or after chromosome relabeling), which would otherwise
    ## throw a subscript-out-of-bounds error rather than giving NA
    row_idx <- match(map_sub$marker, rownames(ld_ws))
    y_val <- rep(NA_real_, length(row_idx))
    ok <- !is.na(row_idx)
    y_val[ok] <- ld_ws[row_idx[ok], rho]
    map_sub[, y_val := y_val]

  } else if (identical(value_col, "fdr_q")) {
    if (is.null(ld_ws) || is.null(rho) || is.null(th_ldw) || is.null(p_col)) {
      stop("value_col = 'fdr_q' requires `ld_ws`, `rho`, `th_ldw`, and `p_col` to be supplied.")
    }
    rho <- as.character(rho)
    if (!rho %in% colnames(ld_ws)) {
      stop("rho = '", rho, "' not found in colnames(ld_ws): ",
           paste(colnames(ld_ws), collapse = ", "))
    }

    ## replicate run_one_grid()'s filtering EXACTLY, so the q-values shown
    ## are the ones that actually produced the clusters being plotted:
    ## restrict to markers common to map & ld_ws, then threshold on ld_w
    ## at th_ldw (or keep everyone common if th_ldw <= 0), THEN FDR-correct
    ## only the kept subset -- everything else stays NA (not tested).
    row_idx   <- match(map_sub$marker, rownames(ld_ws))
    is_common <- !is.na(row_idx)

    ld_vec <- rep(NA_real_, nrow(map_sub))
    ld_vec[is_common] <- ld_ws[row_idx[is_common], rho]

    if (th_ldw > 0) {
      keep <- is_common & !is.na(ld_vec) & (ld_vec > quantile(ld_vec[is_common], th_ldw, na.rm = TRUE))
    } else {
      keep <- is_common
    }

    y_val <- rep(NA_real_, nrow(map_sub))
    if (any(keep)) {
      y_val[keep] <- -log10(p.adjust(map_sub[[p_col]][keep], method = "fdr"))
    }
    map_sub[, y_val := y_val]

  } else {
    map_sub[, y_val := get(value_col)]
  }

  cl_dt <- if (length(clusters) > 0) {
    rbindlist(lapply(seq_along(clusters), function(i)
      data.table(marker = clusters[[i]], CL_id = i)))
  } else {
    data.table(marker = character(0), CL_id = integer(0))
  }

  if (nrow(cl_dt) > 0 && nrow(assignments) > 0) {
    cl_dt <- merge(cl_dt, assignments[, .(CL_id, is_TP)], by = "CL_id", all.x = TRUE)
  } else {
    cl_dt[, is_TP := logical(0)]
  }

  map_sub <- merge(map_sub, cl_dt, by = "marker", all.x = TRUE)

  map_sub[, status := fifelse(is.na(CL_id), "not called",
                              fifelse(is_TP, "true positive OR", "false positive OR"))]
  map_sub[, status := factor(status, levels = c("not called", "false positive OR", "true positive OR"))]
  map_sub[]
}

#----------------------------------------------------------
# Shared rendering step, used by both the single-file and
# multi-chromosome ("genome") Manhattan plots
#----------------------------------------------------------
.render_OR_manhattan <- function(map_sub, title = NULL, chr_levels = NULL, nrow_facets = 1,
                                 value_label = "value", color_by = "status",
                                 wes_palette_name = "Zissou1") {

  if (!is.null(chr_levels)) map_sub[, Chr := factor(Chr, levels = chr_levels)]

  p <- ggplot(map_sub, aes(x = indx, y = y_val)) +
    geom_point(aes(shape = type == "QTN", size = type == "QTN"), alpha = 0.85) +
    scale_shape_manual(values = c(`TRUE` = 3, `FALSE` = 16), guide = "none") +
    scale_size_manual(values = c(`TRUE` = 3, `FALSE` = 1), guide = "none") +
    facet_wrap(~ Chr, scales = "free_x", nrow = nrow_facets) +
    labs(x = "Genomic index", y = value_label, color = NULL, title = title) +
    theme_minimal(base_size = 13) +
    theme(strip.text = element_text(face = "bold"), legend.position = "bottom",
          axis.text.x = element_blank(), axis.ticks.x = element_blank())

  if (identical(color_by, "status")) {
    ## discrete: not called / false positive OR / true positive OR
    cols <- c("not called" = "grey70", "false positive OR" = "steelblue3", "true positive OR" = "firebrick3")
    p <- p + aes(color = status) + scale_color_manual(values = cols)
  } else {
    ## continuous column (e.g. max_LD_with_QTN) via a wesanderson continuous palette
    if (!color_by %in% names(map_sub)) {
      stop("color_by = '", color_by, "' not found in the data. ",
           "Use 'status' for the TP/FP/not-called coloring, or a numeric column name.")
    }
    p <- p + aes(color = .data[[color_by]]) +
      scale_color_gradientn(colors = wes_palette(wes_palette_name, 100, type = "continuous"),
                            name = color_by, na.value = "grey85")
  }

  p
}

#----------------------------------------------------------
# Manhattan-style plot of one grid point's outlier calls:
# true-positive ORs, false-positive ORs, and un-called markers.
# Uses raw (unadjusted) -log10(p) as a stable visual background
# so the same plot is comparable across different filtering
# settings -- only the coloring of "called" points changes,
# not the underlying point cloud.
#----------------------------------------------------------
plot_OR_manhattan <- function(map, value_col, clusters, assignments, title = NULL,
                              value_label = value_col, ld_ws = NULL, rho = NULL,
                              th_ldw = NULL, p_col = NULL,
                              color_by = "status", wes_palette_name = "Zissou1") {
  map_sub <- annotate_OR_calls(map, value_col, clusters, assignments, ld_ws = ld_ws, rho = rho,
                               th_ldw = th_ldw, p_col = p_col)
  if (!"indx" %in% names(map_sub)) map_sub[, indx := .I]
  .render_OR_manhattan(map_sub, title = title, nrow_facets = 1, value_label = value_label,
                       color_by = color_by, wes_palette_name = wes_palette_name)
}

#----------------------------------------------------------
# Genome-wide view: run the SAME grid point independently on
# each of several replicate files (e.g. the 10 chr1..chr10
# simulations for one V/env), then stack the results side by
# side for one combined figure. Each file is still filtered/
# FDR-corrected/clustered entirely on its own -- chromosomes
# are only relabeled here, for display, exactly as in the
# (now-abandoned) concatenation scheme (Chr1 -> Chr(2i-1) QTN,
# Chr2 -> Chr(2i) neutral) so the panels read left-to-right as
# a 20-chromosome genome. This does NOT re-pool FDR/GRM/LD
# across files -- purely a plotting convenience.
#----------------------------------------------------------
#----------------------------------------------------------
# Expensive step: run the per-file OR pipeline ONCE for a
# given grid point, across all replicate files, and cache
# everything needed to plot it repeatedly (relabeled map,
# clusters, assignments, ld_ws) -- WITHOUT rendering. Re-run
# only when the grid point (rho/th_ldw/r2_th/l_min) or the
# file set changes; changing value_col/color_by afterwards is
# then a cheap plotting-only step via
# plot_OR_manhattan_genome_cached().
#----------------------------------------------------------
score_OR_genome <- function(file_paths, p_cols, p_name,
                            sel_rho, sel_th_ldw, sel_r2_th, sel_l_min, sel_alpha = 0.05,
                            bp_th_cluster = 5e5, rho_r2_focal = 0.25,
                            rho_d_focal = 0.99, cores = 1,
                            r2_grid_scale = c("r2", "rho")) {

  r2_grid_scale <- match.arg(r2_grid_scale)

  if (!p_name %in% names(p_cols)) {
    stop("p_name = '", p_name, "' not found in names(p_cols) = ",
         paste(names(p_cols), collapse = ", "),
         ". p_name must exactly match one of the names given to p_cols ",
         "(e.g. p_cols <- c(EMX = \"emx_p\") means p_name must be \"EMX\", not \"EMMAX\").")
  }

  ## order files by the chrN embedded in the filename (chr1..chr10, ...)
  chr_num <- as.integer(gsub("[^0-9]", "",
                             regmatches(basename(file_paths), regexpr("chr[0-9]+", basename(file_paths)))))
  if (any(is.na(chr_num))) stop("Could not parse chromosome number from one or more filenames.")
  file_paths <- file_paths[order(chr_num)]

  maps        <- vector("list", length(file_paths))
  clusters_l  <- vector("list", length(file_paths))
  assign_l    <- vector("list", length(file_paths))
  ld_ws_l     <- vector("list", length(file_paths))

  for (i in seq_along(file_paths)) {
    fp <- file_paths[i]
    message("Processing ", basename(fp))

    scored <- tryCatch(
      run_and_score_one_sim_file(
        fp, p_cols = p_cols, th_ldw_grid = sel_th_ldw,
        r2_grid = sel_r2_th, lmin_grid = sel_l_min,
        bp_th_cluster = bp_th_cluster, rho_r2_focal = rho_r2_focal,
        rho_d_focal = rho_d_focal, alpha_grid = sel_alpha, cores = cores,
        rho_grid = sel_rho, r2_grid_scale = r2_grid_scale
      ),
      error = function(e) { message("  FAILED: ", conditionMessage(e)); NULL }
    )

    d <- readRDS(fp)
    map_i   <- flag_true_positive_QTNs(d$map)
    ld_ws_i <- d$ld_ws   ## always cache ld_ws so value_col="ld_w" works later without a re-run

    if (is.null(scored) || nrow(scored) == 0) {
      ## no potential outliers anywhere in this file at this setting --
      ## still cache the chromosome, entirely un-called, rather than dropping it
      clusters    <- list()
      assignments <- data.table(CL_id = integer(0), n_loci = integer(0),
                                qtn = character(0), evidence = numeric(0), is_TP = logical(0))
    } else {
      match_col <- if (identical(r2_grid_scale, "rho")) "r2_grid_rho" else "r2_th"
      row <- scored[rho == sel_rho & th_ldw == sel_th_ldw &
                      get(match_col) == sel_r2_th & l_min == sel_l_min]
      clusters    <- row[[p_name]][[1]]
      assignments <- row[[paste0("assignments_", p_name)]][[1]]
    }

    ## relabel chromosomes/markers for DISPLAY only (Chr1->Chr(2i-1) QTN, Chr2->Chr(2i) neutral)
    old_marker <- map_i$marker
    new_chr    <- ifelse(map_i$Chr == "Chr1", paste0("Chr", 2 * i - 1), paste0("Chr", 2 * i))
    old_prefix <- paste0(map_i$Chr, ":")
    new_prefix <- paste0(new_chr, ":")
    new_marker <- mapply(function(m, o, n) sub(o, n, m, fixed = TRUE),
                         old_marker, old_prefix, new_prefix, USE.NAMES = FALSE)
    marker_map <- setNames(new_marker, old_marker)

    map_i[, Chr := new_chr]
    map_i[, marker := new_marker]

    matched <- match(rownames(ld_ws_i), names(marker_map))
    rownames(ld_ws_i) <- ifelse(is.na(matched), rownames(ld_ws_i), marker_map[matched])

    clusters_relabeled <- lapply(clusters, function(cl) unname(marker_map[cl]))

    maps[[i]]       <- map_i
    clusters_l[[i]] <- clusters_relabeled
    assign_l[[i]]   <- assignments
    ld_ws_l[[i]]    <- ld_ws_i
  }

  chr_order <- unlist(lapply(maps, function(m) unique(m$Chr)))
  chr_order <- chr_order[order(as.integer(gsub("[^0-9]", "", chr_order)))]

  or_per_chr <- vapply(seq_along(maps), function(i) {
    length(clusters_l[[i]])
  }, integer(1))
  message(sum(or_per_chr > 0), " of ", length(maps),
          " chromosomes have >=1 called OR at this setting")

  list(
    files = file_paths, maps = maps, clusters = clusters_l,
    assignments = assign_l, ld_ws = ld_ws_l, chr_order = chr_order,
    p_name = p_name, sel_rho = sel_rho, sel_th_ldw = sel_th_ldw,
    sel_r2_th = sel_r2_th, sel_l_min = sel_l_min, r2_grid_scale = r2_grid_scale
  )
}

#----------------------------------------------------------
# Cheap step: render the genome-wide plot from an object
# returned by score_OR_genome(), for any value_col/color_by
# combination -- no re-scoring, just annotate_OR_calls() +
# rbind + render per call.
#----------------------------------------------------------
plot_OR_manhattan_genome_cached <- function(scored_genome, value_col, value_label = value_col,
                                            title = NULL, nrow_facets = 1, p_col = NULL,
                                            color_by = "status", wes_palette_name = "Zissou1") {

  all_annot <- lapply(seq_along(scored_genome$maps), function(i) {
    annot_i <- annotate_OR_calls(scored_genome$maps[[i]], value_col,
                                 scored_genome$clusters[[i]], scored_genome$assignments[[i]],
                                 ld_ws = scored_genome$ld_ws[[i]], rho = scored_genome$sel_rho,
                                 th_ldw = scored_genome$sel_th_ldw, p_col = p_col)
    annot_i[, file := basename(scored_genome$files[i])]
    annot_i
  })

  combined <- rbindlist(all_annot, fill = TRUE)
  combined[, Chr := factor(Chr, levels = scored_genome$chr_order)]
  setorder(combined, Chr, Pos)
  combined[, indx := .I]

  if (is.null(title)) {
    r2_label <- if (identical(scored_genome$r2_grid_scale, "rho")) "r2_rho" else "r2_th"
    title <- sprintf("%s | rho=%s, th_ldw=%.2f, %s=%.2f, l_min=%d",
                     scored_genome$p_name, scored_genome$sel_rho, scored_genome$sel_th_ldw,
                     r2_label, scored_genome$sel_r2_th, scored_genome$sel_l_min)
  }

  .render_OR_manhattan(combined, title = title, chr_levels = scored_genome$chr_order,
                       nrow_facets = nrow_facets, value_label = value_label,
                       color_by = color_by, wes_palette_name = wes_palette_name)
}

#----------------------------------------------------------
# Convenience wrapper: score + plot in one call (same
# interface as before). For repeated plots of the same grid
# point with different value_col/color_by, call
# score_OR_genome() once and reuse plot_OR_manhattan_genome_cached()
# instead of calling this repeatedly.
#----------------------------------------------------------
plot_OR_manhattan_genome <- function(file_paths, p_cols, p_name, value_col,
                                     sel_rho, sel_th_ldw, sel_r2_th, sel_l_min, sel_alpha = 0.05,
                                     bp_th_cluster = 5e5, rho_r2_focal = 0.25,
                                     rho_d_focal = 0.99, cores = 1,
                                     title = NULL, value_label = value_col, nrow_facets = 1,
                                     p_col = NULL, color_by = "status", wes_palette_name = "Zissou1",
                                     r2_grid_scale = c("r2", "rho")) {

  scored_genome <- score_OR_genome(
    file_paths = file_paths, p_cols = p_cols, p_name = p_name,
    sel_rho = sel_rho, sel_th_ldw = sel_th_ldw, sel_r2_th = sel_r2_th, sel_l_min = sel_l_min,
    sel_alpha = sel_alpha, bp_th_cluster = bp_th_cluster, rho_r2_focal = rho_r2_focal,
    rho_d_focal = rho_d_focal, cores = cores, r2_grid_scale = r2_grid_scale
  )

  plot_OR_manhattan_genome_cached(
    scored_genome, value_col = value_col, value_label = value_label,
    title = title, nrow_facets = nrow_facets, p_col = p_col, color_by = color_by,
    wes_palette_name = wes_palette_name
  )
}


plot_OR_manhattan_from_results <- function(map, all_results, value_col, p_name,
                                           sel_rho, sel_th_ldw, sel_r2_th, sel_l_min,
                                           value_label = value_col, ld_ws = NULL, p_col = NULL,
                                           color_by = "status", wes_palette_name = "Zissou1",
                                           r2_grid_scale = c("r2", "rho")) {

  r2_grid_scale <- match.arg(r2_grid_scale)

  if (!p_name %in% names(all_results)) {
    known_cols <- setdiff(names(all_results), c("r2_th", "l_min", "th_ldw", "rho", "n_loci", "file"))
    known_cols <- known_cols[!grepl("^(TP_|FP_|FN_|Precision_|Recall_|PR_|assignments_)", known_cols)]
    stop("p_name = '", p_name, "' is not a column in all_results (its clusters would silently be ",
         "NULL -> nothing plotted as called). Available cluster columns: ",
         paste(known_cols, collapse = ", "))
  }

  match_col <- if (identical(r2_grid_scale, "rho")) "r2_grid_rho" else "r2_th"
  if (!match_col %in% names(all_results)) {
    stop("r2_grid_scale = 'rho' requires an `r2_grid_rho` column in all_results ",
         "(produced when run_and_score_all()/run_and_score_one_sim_file() was called ",
         "with r2_grid_scale = 'rho'). This table doesn't have one.")
  }

  row <- all_results[rho == sel_rho & th_ldw == sel_th_ldw &
                       get(match_col) == sel_r2_th & l_min == sel_l_min]
  if (nrow(row) != 1) {
    stop("Expected exactly 1 matching row, found ", nrow(row),
         ". Check sel_rho/sel_th_ldw/sel_r2_th/sel_l_min match values present in all_results.")
  }

  clusters    <- row[[p_name]][[1]]
  assignments <- row[[paste0("assignments_", p_name)]][[1]]

  r2_label <- if (identical(r2_grid_scale, "rho")) "r2_rho" else "r2_th"
  title <- sprintf("%s | rho=%s, th_ldw=%.2f, %s=%.2f, l_min=%d",
                   p_name, sel_rho, sel_th_ldw, r2_label, sel_r2_th, sel_l_min)

  plot_OR_manhattan(map, value_col, clusters, assignments, title = title,
                    value_label = value_label, ld_ws = ld_ws, rho = sel_rho,
                    th_ldw = sel_th_ldw, p_col = p_col,
                    color_by = color_by, wes_palette_name = wes_palette_name)
}



#----------------------------------------------------------
# Full pipeline for ONE simulation file
#----------------------------------------------------------
run_and_score_one_sim_file <- function(file_path, p_cols, th_ldw_grid, r2_grid, lmin_grid,
                                       bp_th_cluster = 5e5, rho_r2_focal = 0.75,
                                       rho_d_focal = 0.95, alpha_grid = 0.05, cores = 1,
                                       rho_grid = NULL, r2_grid_scale = c("r2", "rho"),
                                       return_intermediates = FALSE) {

  r2_grid_scale <- match.arg(r2_grid_scale)

  d <- readRDS(file_path)
  map <- flag_true_positive_QTNs(d$map)
  GTs <- d$GTs
  ld_ws <- d$ld_ws
  LD_decay <- d$LD_decay

  ## candidate pool built at the loosest (max) alpha -- see get_potential_outliers()
  potential_outliers <- get_potential_outliers(map, ld_ws, th_ldw_grid, p_cols, alpha_grid = alpha_grid)
  if (length(potential_outliers) == 0) {
    message("No potential outliers for ", basename(file_path))
    return(NULL)
  }

  el_potential <- precompute_LD_edges(
    GTs = GTs[, potential_outliers, drop = FALSE],
    map = map[marker %in% potential_outliers],
    r2_min = 0.1, max_bp = bp_th_cluster, cores = cores
  )

  qtn_ld_table <- precompute_QTN_LD(GTs = GTs, map = map, candidate_markers = potential_outliers,
                                    max_bp = 2e6, cores = cores)

  thr <- get_focal_qtn_thresholds(LD_decay, chr_qtn = "Chr1", rho_r2 = rho_r2_focal, rho_d = rho_d_focal)

  ## fetch the QTN chromosome's decay row once -- used both for the optional
  ## r2_grid rho conversion below, and always to report which rho the fixed
  ## bp_th_cluster distance corresponds to on this file's own decay curve
  ds <- LD_decay$decay_sum[Chr == "Chr1"]
  if (nrow(ds) != 1) {
    stop("Expected exactly one decay_sum row for chr 'Chr1', found ", nrow(ds),
         ". Check LD_decay$decay_sum$Chr labels.")
  }
  ## rho equivalent to the fixed clustering distance, i.e. the inverse of
  ## d_from_rho(): rho = a*d / (1 + a*d) -- same formula compute_LD_decay()
  ## already uses internally for rho_slide_raw
  bp_th_cluster_rho <- (ds$a * bp_th_cluster) / (1 + ds$a * bp_th_cluster)

  ## optionally treat r2_grid as rho values and convert to actual r^2
  ## clustering thresholds using the QTN chromosome's own fitted LD-decay
  ## curve -- same source (LD_decay$decay_sum) as the focal-QTN thresholds
  ## above, so "rho" means the same thing in both places for this file.
  r2_grid_rho <- NULL
  if (identical(r2_grid_scale, "rho")) {
    r2_grid_rho <- r2_grid                       ## keep the original rho values for traceability
    r2_grid <- ld_from_rho(b = ds$b, c = ds$c, rho = r2_grid)
    message("r2_grid (rho -> r2) for ", basename(file_path), ": ",
            paste(sprintf("%.3f->%.3f", r2_grid_rho, r2_grid), collapse = ", "))
  }

  rho_vals <- if (is.null(rho_grid)) colnames(ld_ws) else as.character(rho_grid)
  param_grid <- CJ(rho = rho_vals, th_ldw = th_ldw_grid, alpha = alpha_grid)

  outliers <- rbindlist(lapply(seq_len(nrow(param_grid)), function(i) {
    pars <- param_grid[i]
    run_one_grid(map = map, el = el_potential, ld_ws = ld_ws,
                 rho = pars$rho, th_ldw = pars$th_ldw,
                 p_cols = p_cols, alpha = pars$alpha,
                 r2_grid = r2_grid, lmin_grid = lmin_grid,
                 bp_th = bp_th_cluster, cores = cores)
  }), fill = TRUE)

  scored <- score_outlier_grid(outliers, map, qtn_ld_table, names(p_cols),
                               r2_min_focal = thr$r2_min_focal, d_max_focal = thr$d_max_focal)

  if (!is.null(r2_grid_rho)) {
    ## attach the rho that produced each row's r2_th, so callers can select
    ## grid points by rho without needing to know the converted r2 value
    rho_lookup <- setNames(r2_grid_rho, as.character(r2_grid))
    scored[, r2_grid_rho := rho_lookup[as.character(r2_th)]]
  }

  ## report the clustering distance threshold on both scales -- bp_th_cluster
  ## is fixed across files, but the rho it corresponds to differs by file
  ## since each has its own decay rate `a`
  scored[, bp_th_cluster     := bp_th_cluster]
  scored[, bp_th_cluster_rho := bp_th_cluster_rho]

  scored[, file := basename(file_path)]

  if (return_intermediates) {
    return(list(scored = scored[], map = map, el = el_potential, qtn_ld_table = qtn_ld_table,
                thr = thr, bp_th_cluster = bp_th_cluster, bp_th_cluster_rho = bp_th_cluster_rho))
  }
  scored[]
}

#----------------------------------------------------------
# Loop over all simulation files and aggregate
#----------------------------------------------------------
run_and_score_all <- function(parsed_folder, p_cols, th_ldw_grid, r2_grid, lmin_grid,
                              bp_th_cluster = 5e5, rho_r2_focal = 0.75, rho_d_focal = 0.95,
                              alpha_grid = 0.05, cores = 1, r2_grid_scale = c("r2", "rho")) {

  r2_grid_scale <- match.arg(r2_grid_scale)
  files <- list.files(parsed_folder, pattern = "\\.rds$", full.names = TRUE)

  all_scored <- rbindlist(lapply(files, function(f) {
    message("Processing ", basename(f))
    tryCatch(
      run_and_score_one_sim_file(f, p_cols, th_ldw_grid, r2_grid, lmin_grid,
                                 bp_th_cluster, rho_r2_focal, rho_d_focal, alpha_grid, cores,
                                 r2_grid_scale = r2_grid_scale),
      error = function(e) {
        message("  FAILED (", basename(f), "): ", conditionMessage(e))
        NULL
      }
    )
  }), fill = TRUE)

  ## attach filename-derived metadata (sim, bgs, Chr, V, c, env) for aggregation
  base  <- sub("\\.rds$", "", all_scored$file)
  parts <- tstrsplit(base, "_", fixed = TRUE)
  all_scored[, `:=`(sim = parts[[1]], bgs = parts[[2]], Chr = parts[[3]],
                    V = parts[[4]], c = parts[[5]], env = parts[[6]])]
  all_scored[]
}

#----------------------------------------------------------
# Pooled (micro-averaged) Precision/Recall/PR across replicates
#
# Sums TP/FP/FN across replicate files WITHIN each grid point
# (rho, th_ldw, r2_th, l_min[, group_vars]) before computing
# ratios -- this is what gives you higher-resolution PR estimates
# without needing to concatenate the underlying genotype data or
# share FDR corrections across replicates. Naively averaging each
# replicate's own Precision/Recall instead ("macro-averaging")
# breaks down badly when many replicates have TP=FP=0 (undefined
# ratio) or tiny counts (noisy 0/1 ratios dominating the mean).
#----------------------------------------------------------
aggregate_PR <- function(all_scored, p_names, group_vars = character(0)) {

  by_cols <- c(intersect(c("rho", "th_ldw", "r2_th", "l_min", "alpha"), names(all_scored)), group_vars)

  out <- list()
  for (nm in p_names) {
    tp_col <- paste0("TP_", nm); fp_col <- paste0("FP_", nm); fn_col <- paste0("FN_", nm)

    agg <- all_scored[, .(
      TP = sum(get(tp_col), na.rm = TRUE),
      FP = sum(get(fp_col), na.rm = TRUE),
      FN = sum(get(fn_col), na.rm = TRUE),
      n_files = .N
    ), by = by_cols]

    agg[, Precision := ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_)]
    agg[, Recall    := ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)]
    agg[, PR        := ifelse(!is.na(Precision) & !is.na(Recall), Precision * Recall, NA_real_)]
    agg[, method := nm]

    out[[nm]] <- agg
  }

  rbindlist(out, use.names = TRUE)
}



######################################################
## Consistency (C-score) and AUC-PR framework
##
## All four functions below pool over (rho x th_ldw x l_min)
## for a FIXED r2_th, i.e. r2_th is deliberately NOT integrated
## out -- it changes which markers get grouped into which OR in
## the first place, so it isn't just another significance-style
## threshold the way rho/th_ldw/l_min are.
##
## Requires `scored` = output of run_and_score_one_sim_file(...,
## return_intermediates = TRUE)$scored (or a row-subset/rbind of
## several such tables for the same file).
######################################################

#----------------------------------------------------------
# 1) Per-MARKER consistency score (works for empirical or
#    simulated data -- matches summarise_stability() in
#    ld_w_filtering_3sp.R, generalized to also integrate out l_min).
#    fixed_r2_th = NULL (default) pools r2_th in too, for a single
#    overall "how often was this marker ever called, across
#    everything" score -- e.g. for attaching to map for plotting.
#    Pass a specific value to restrict to one r2_th, matching the
#    formal C-score-sweep usage where r2_th is kept separate.
#----------------------------------------------------------
compute_marker_C_score <- function(scored, p_name, fixed_r2_th = NULL, filter_fun = NULL) {

  sub <- if (is.null(filter_fun)) scored else filter_fun(scored)
  sub <- if (is.null(fixed_r2_th)) sub else sub[r2_th == fixed_r2_th]
  n_grid <- nrow(sub)
  if (n_grid == 0) {
    msg <- if (is.null(fixed_r2_th)) "No rows in `scored` (after filter_fun, if supplied)." else paste0("No rows in `scored` at r2_th = ", fixed_r2_th, " (after filter_fun, if supplied).")
    stop(msg)
  }

  tab <- table(unlist(sub[[p_name]]))

  ## `p_name` may have ZERO calls across the ENTIRE grid (e.g. res was not
  ## NULL because SOME p_col had potential outliers somewhere, but this
  ## particular p_name never did on this file/chromosome). table(NULL)
  ## then has length 0, and names(tab) is NULL rather than character(0) --
  ## passing marker = NULL into data.table() silently DROPS the marker
  ## column instead of creating an empty one, which breaks any downstream
  ## join on "marker" (e.g. attach_C_score_to_map()). Guard explicitly so
  ## the returned table always has a well-typed, present marker column.
  if (length(tab) == 0) {
    return(data.table(marker = character(0), n_grid = n_grid,
                      n_called = integer(0), C_score = numeric(0)))
  }

  data.table(marker = names(tab), n_grid = n_grid,
             n_called = as.integer(tab), C_score = as.numeric(tab) / n_grid)
}

#----------------------------------------------------------
# 2) Per-QTN "recall stability" (SIMULATION ONLY -- needs
#    true_pos_QTN and the assignments_<p_name> list-column).
#    For every known true-positive QTN, including ones NEVER
#    recovered (C_score = 0), the fraction of grid points where
#    it was correctly assigned as SOME OR's focal QTN.
#    fixed_r2_th = NULL (default) pools r2_th in too, same as above.
#    filter_fun, if supplied, is applied to `scored` FIRST (e.g.
#    function(x) x[th_ldw == 0] to isolate the no-LD-filtering
#    C-score, or function(x) x[alpha == 0.01] as a difficulty proxy).
#----------------------------------------------------------
compute_QTN_C_score <- function(scored, p_name, fixed_r2_th = NULL, map, filter_fun = NULL) {

  sub <- if (is.null(filter_fun)) scored else filter_fun(scored)
  sub <- if (is.null(fixed_r2_th)) sub else sub[r2_th == fixed_r2_th]
  n_grid <- nrow(sub)
  if (n_grid == 0) {
    msg <- if (is.null(fixed_r2_th)) "No rows in `scored` (after filter_fun, if supplied)." else paste0("No rows in `scored` at r2_th = ", fixed_r2_th, " (after filter_fun, if supplied).")
    stop(msg)
  }

  true_pos_markers <- map[true_pos_QTN == TRUE, marker]

  recovered_list <- lapply(sub[[paste0("assignments_", p_name)]], function(a) {
    if (is.null(a) || nrow(a) == 0) return(character(0))
    unique(a[is_TP == TRUE, qtn])
  })

  tab <- table(unlist(recovered_list))

  out <- data.table(qtn = true_pos_markers, n_grid = n_grid)
  out[, n_recovered := as.integer(tab[qtn])]
  out[is.na(n_recovered), n_recovered := 0L]
  out[, C_score := n_recovered / n_grid]
  out[]
}

#----------------------------------------------------------
# Attach a compute_marker_C_score() result onto `map` as a new
# column, so it can be used directly as value_col/color_by in
# the Manhattan-plotting functions. Markers absent from
# marker_C_score (never called anywhere in the grid) get
# fill_value (0 by default) rather than NA, since "never called"
# IS a C_score of exactly 0, not a missing/unknown value.
#----------------------------------------------------------
attach_C_score_to_map <- function(map, marker_C_score, col_name = "C_score", fill_value = 0) {
  map <- copy(map)
  map[, (col_name) := fill_value]

  if (!"marker" %in% names(marker_C_score)) {
    stop("`marker_C_score` has no `marker` column -- likely produced by an older/unguarded ",
         "compute_marker_C_score() call when a p_name had zero calls across the whole grid. ",
         "Every marker in `map` will keep fill_value = ", fill_value, ".")
  }
  if (nrow(marker_C_score) == 0) return(map[])   ## nothing to join; every marker stays at fill_value

  map[marker_C_score, on = "marker", (col_name) := i.C_score]
  map[]
}

#----------------------------------------------------------
# 3) AUC of cummax(metric) under random search order -- generic
#    core, reused for both the raw grid (below) and the
#    C-score-threshold sweep (further down). Repeatedly shuffle
#    `values`, track the running best value as if trying them
#    one at a time in that random order, average across shuffles,
#    and take the area under the resulting curve (x = fraction
#    tried, y = expected best-so-far). Starts at (0,0) so AUC is
#    comparable across different-sized value sets.
#----------------------------------------------------------
auc_cummax <- function(values, n_perm = 200, seed = NULL) {

  values[is.na(values)] <- 0
  n <- length(values)
  if (n == 0) return(list(trajectory = data.table(), AUC = NA_real_))

  if (!is.null(seed)) set.seed(seed)

  cummax_mat <- replicate(n_perm, cummax(sample(values)))
  if (is.null(dim(cummax_mat))) cummax_mat <- matrix(cummax_mat, ncol = n_perm)
  mean_traj <- rowMeans(cummax_mat)

  x <- c(0, seq_len(n) / n)
  y <- c(0, mean_traj)
  AUC <- sum((y[-1] + y[-length(y)]) / 2 * diff(x))

  list(
    trajectory = data.table(n_tried = seq_len(n), frac_tried = seq_len(n) / n, mean_cummax = mean_traj),
    AUC = AUC
  )
}

#----------------------------------------------------------
# AUC of cummax(PR) over the raw (rho x th_ldw x l_min x r2_th x
# alpha) grid. POOLS r2_th by default (fixed_r2_th = NULL) -- see
# above. ALSO pools TP/FP/FN across whatever rows share the same
# grid point (rho/th_ldw/r2_th/l_min/alpha) BEFORE computing PR,
# same rationale as aggregate_PR(): a single replicate's grid
# point very often has TP=0 (few/no true positives findable on
# that one chromosome at that setting), which mechanically forces
# PR=0. Feeding those raw per-file zeros into the random-search
# AUC understates real performance -- pooling counts first, then
# computing PR, then AUC, fixes this. For a single-file `scored`
# table (one row per grid point already) this pooling is a
# mathematical no-op, so behaviour there is unchanged; it only
# actually matters when `scored` spans multiple replicate files.
#----------------------------------------------------------
auc_cummax_PR <- function(scored, p_name, fixed_r2_th = NULL, n_perm = 200, seed = NULL) {

  sub <- if (is.null(fixed_r2_th)) scored else scored[r2_th == fixed_r2_th]

  by_cols <- intersect(c("rho", "th_ldw", "r2_th", "l_min", "alpha"), names(sub))
  tp_col <- paste0("TP_", p_name); fp_col <- paste0("FP_", p_name); fn_col <- paste0("FN_", p_name)

  pooled <- sub[, .(
    TP = sum(get(tp_col), na.rm = TRUE),
    FP = sum(get(fp_col), na.rm = TRUE),
    FN = sum(get(fn_col), na.rm = TRUE)
  ), by = by_cols]

  pooled[, Precision := ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_)]
  pooled[, Recall    := ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)]
  pooled[, PR        := ifelse(!is.na(Precision) & !is.na(Recall), Precision * Recall, NA_real_)]

  res <- auc_cummax(pooled$PR, n_perm = n_perm, seed = seed)
  res$pooled_r2_th <- is.null(fixed_r2_th)
  res$n_grid_points <- nrow(pooled)
  res
}

#----------------------------------------------------------
# 4) C-score-threshold sweep: re-cluster and re-classify the
#    set of markers passing each C-score threshold (instead of
#    "significant at one specific grid point"), for one r2_th.
#    `el` and `qtn_ld_table` should come from the SAME file's
#    run_and_score_one_sim_file(..., return_intermediates = TRUE)
#    call that produced `scored` / the marker_C_score table --
#    reused as-is, not recomputed.
#----------------------------------------------------------
run_C_score_threshold_grid <- function(marker_C_score, map, el, qtn_ld_table,
                                       r2_min_focal, d_max_focal,
                                       C_score_grid, r2_th, l_min = 1, bp_th = Inf) {

  rbindlist(lapply(C_score_grid, function(C_th) {

    called_markers <- marker_C_score[C_score >= C_th, marker]

    if (length(called_markers) == 0) {
      return(data.table(C_score_threshold = C_th, r2_th = r2_th, l_min = l_min,
                        TP = 0L, FP = 0L, FN = NA_integer_,
                        Precision = NA_real_, Recall = 0, PR = 0))
    }

    clusters_dt <- LD_igraph_components(el = el, markers = called_markers, r2_th = r2_th, bp_th = bp_th)
    clusters_dt <- clusters_dt[n_loci >= l_min]
    cls <- split(clusters_dt$marker, clusters_dt$CL_id)

    res <- evaluate_ORs(cls, map, qtn_ld_table, r2_min_focal, d_max_focal)

    data.table(C_score_threshold = C_th, r2_th = r2_th, l_min = l_min,
               TP = res$TP, FP = res$FP, FN = res$FN,
               Precision = res$Precision, Recall = res$Recall, PR = res$PR)
  }))
}

#----------------------------------------------------------
# Full sweep: r2_th x C_score_threshold. Computes marker C-scores
# fresh at EACH r2_th (since OR membership -- and therefore
# C-score -- depends on r2_th), then re-clusters/classifies at
# each C-score threshold within that r2_th.
#----------------------------------------------------------
run_C_score_sweep <- function(scored, map, el, qtn_ld_table, p_name,
                              r2_min_focal, d_max_focal,
                              r2_th_grid, C_score_grid, l_min = 1, bp_th = Inf,
                              filter_fun = NULL) {

  rbindlist(lapply(r2_th_grid, function(r2_th) {
    mc <- compute_marker_C_score(scored, p_name, fixed_r2_th = r2_th, filter_fun = filter_fun)
    run_C_score_threshold_grid(mc, map, el, qtn_ld_table, r2_min_focal, d_max_focal,
                               C_score_grid = C_score_grid, r2_th = r2_th,
                               l_min = l_min, bp_th = bp_th)
  }))
}

#----------------------------------------------------------
# Multi-method convenience wrapper: run_C_score_sweep() once per
# p_name, tagging each with a `method` column, for direct
# side-by-side comparison (e.g. EMX vs LFMM).
#----------------------------------------------------------
run_C_score_sweep_multi <- function(scored, map, el, qtn_ld_table, p_names,
                                    r2_min_focal, d_max_focal,
                                    r2_th_grid, C_score_grid, l_min = 1, bp_th = Inf,
                                    filter_fun = NULL) {

  rbindlist(lapply(p_names, function(nm) {
    out <- run_C_score_sweep(scored, map, el, qtn_ld_table, p_name = nm,
                             r2_min_focal = r2_min_focal, d_max_focal = d_max_focal,
                             r2_th_grid = r2_th_grid, C_score_grid = C_score_grid,
                             l_min = l_min, bp_th = bp_th, filter_fun = filter_fun)
    out[, method := nm]
    out[]
  }))
}

#----------------------------------------------------------
# Run the C-score sweep across ALL replicate files (e.g. the 10
# chr1..chr10 simulations for one V/env), tagging each row with
# `file`. Each file needs its OWN el/qtn_ld_table/map (LD
# structure is chromosome-specific, unlike the raw PR values
# auc_cummax_PR() pools directly), so this reruns
# run_and_score_one_sim_file(..., return_intermediates = TRUE)
# once per file. The payoff: auc_cummax_PR_Cscore_sweep() then
# pools n_files x length(r2_th_grid) points per C_score_threshold
# instead of just length(r2_th_grid) -- e.g. 10x more data points
# from 10 replicate files, same logic as pooling the raw grid.
#----------------------------------------------------------
run_C_score_sweep_all_files <- function(file_paths, p_cols, p_names, th_ldw_grid, r2_grid, lmin_grid,
                                        r2_th_grid, C_score_grid, l_min_cluster = 1,
                                        bp_th_cluster = 5e5, rho_r2_focal = 0.75, rho_d_focal = 0.95,
                                        alpha_grid = 0.05, cores = 1, r2_grid_scale = c("r2", "rho"),
                                        bp_th = Inf, filter_fun = NULL) {

  r2_grid_scale <- match.arg(r2_grid_scale)

  out <- rbindlist(lapply(file_paths, function(fp) {
    message("Processing ", basename(fp))

    res <- tryCatch(
      run_and_score_one_sim_file(
        fp, p_cols = p_cols, th_ldw_grid = th_ldw_grid, r2_grid = r2_grid, lmin_grid = lmin_grid,
        bp_th_cluster = bp_th_cluster, rho_r2_focal = rho_r2_focal, rho_d_focal = rho_d_focal,
        alpha_grid = alpha_grid, cores = cores, r2_grid_scale = r2_grid_scale, return_intermediates = TRUE
      ),
      error = function(e) { message("  FAILED (", basename(fp), "): ", conditionMessage(e)); NULL }
    )
    if (is.null(res)) return(NULL)

    sweep <- run_C_score_sweep_multi(
      scored = res$scored, map = res$map, el = res$el, qtn_ld_table = res$qtn_ld_table,
      p_names = p_names, r2_min_focal = res$thr$r2_min_focal, d_max_focal = res$thr$d_max_focal,
      r2_th_grid = r2_th_grid, C_score_grid = C_score_grid, l_min = l_min_cluster, bp_th = bp_th,
      filter_fun = filter_fun
    )
    sweep[, file := basename(fp)]
    sweep[]
  }), fill = TRUE)

  out[]
}

#----------------------------------------------------------
# Collapse the C-score sweep's r2_th axis via the SAME random-
# search-AUC logic as auc_cummax_PR() -- for each C-score
# threshold (and `method` if present), pool PR across r2_th_grid
# and report one AUC. This is the final, doubly-robust summary:
# rho/th_ldw/l_min integrated out via the C-score itself, r2_th
# integrated out via this AUC step -- leaving C_score_threshold
# as the only axis you sweep and report, the same way you'd
# sweep alpha when comparing two testing methods.
#----------------------------------------------------------
auc_cummax_PR_Cscore_sweep <- function(C_sweep, n_perm = 200, seed = NULL) {

  group_cols <- intersect(c("C_score_threshold", "method"), names(C_sweep))
  pool_cols  <- intersect(c("r2_th", "l_min"), names(C_sweep))

  ## Pool TP/FP/FN across replicate files (a run_C_score_sweep_all_files()
  ## table has one row per file per r2_th; a single-file table already has
  ## just one row per r2_th, so this is a no-op there) WITHIN each r2_th,
  ## before computing PR -- same rationale as auc_cummax_PR(): per-file PR
  ## at one r2_th is very often TP=0, and feeding those raw zeros into the
  ## AUC understates real performance.
  pooled <- C_sweep[, .(
    TP = sum(TP, na.rm = TRUE), FP = sum(FP, na.rm = TRUE), FN = sum(FN, na.rm = TRUE)
  ), by = c(group_cols, pool_cols)]

  pooled[, Precision := ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_)]
  pooled[, Recall    := ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)]
  pooled[, PR        := ifelse(!is.na(Precision) & !is.na(Recall), Precision * Recall, NA_real_)]

  pooled[, {
    res <- auc_cummax(PR, n_perm = n_perm, seed = seed)
    list(AUC_PR = res$AUC, n_r2_th = .N)
  }, by = group_cols]
}

#----------------------------------------------------------
# ONE overall AUC per method, pooling BOTH C_score_threshold AND
# r2_th into a single random-search space -- the direct analogue
# of auc_cummax_PR() applied to the C-score framework instead of
# the raw (rho x th_ldw x l_min x r2_th) grid. Where
# auc_cummax_PR_Cscore_sweep() keeps C_score_threshold as a
# reportable axis (like sweeping alpha) and only integrates out
# r2_th, THIS function integrates out BOTH, giving one number per
# method that answers "how easily do you find a good (C, r2_th)
# combination by random search" -- directly comparable to
# auc_cummax_PR()'s "how easily do you find a good
# (rho, th_ldw, l_min, r2_th) combination", so you can compare
# e.g. a C-score built from th_ldw=0 only (via filter_fun when
# building C_sweep) against one built from the full grid, or
# compare across different alpha_grid values, using the same
# single number both times.
#----------------------------------------------------------
auc_cummax_PR_Cscore_overall <- function(C_sweep, n_perm = 200, seed = NULL) {

  group_cols <- intersect("method", names(C_sweep))
  pool_cols  <- intersect(c("C_score_threshold", "r2_th", "l_min"), names(C_sweep))

  ## Pool TP/FP/FN across replicate files (if present) WITHIN each
  ## (C_score_threshold, r2_th) combination, before computing PR --
  ## same rationale as auc_cummax_PR()/auc_cummax_PR_Cscore_sweep().
  pooled <- C_sweep[, .(
    TP = sum(TP, na.rm = TRUE), FP = sum(FP, na.rm = TRUE), FN = sum(FN, na.rm = TRUE)
  ), by = c(group_cols, pool_cols)]

  pooled[, Precision := ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_)]
  pooled[, Recall    := ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)]
  pooled[, PR        := ifelse(!is.na(Precision) & !is.na(Recall), Precision * Recall, NA_real_)]

  if (length(group_cols) == 0) {
    res <- auc_cummax(pooled$PR, n_perm = n_perm, seed = seed)
    return(data.table(AUC_PR = res$AUC, n_grid_points = nrow(pooled)))
  }

  pooled[, {
    res <- auc_cummax(PR, n_perm = n_perm, seed = seed)
    list(AUC_PR = res$AUC, n_grid_points = .N)
  }, by = group_cols]
}

#----------------------------------------------------------
# Comparison plot: AUC(PR, integrating over r2_th) vs C-score
# threshold, one line per method if a `method` column is present.
#----------------------------------------------------------
plot_Cscore_AUC_curve <- function(auc_by_Cscore, title = NULL) {

  has_method <- "method" %in% names(auc_by_Cscore)

  p <- ggplot(auc_by_Cscore, aes(x = C_score_threshold, y = AUC_PR))

  if (has_method) {
    p <- p + geom_line(aes(color = method), linewidth = 1) +
      geom_point(aes(color = method), size = 2)
  } else {
    p <- p + geom_line(linewidth = 1, color = "steelblue4") +
      geom_point(size = 2, color = "steelblue4")
  }

  p +
    labs(x = "C-score threshold", y = "AUC of cummax(PR), r2_th integrated out",
         color = NULL, title = title) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom")
}

#----------------------------------------------------------
# Comparison plot for the "no LD-filtering" vs "full LD-filtering"
# question across a range of alpha values (alpha as a cheap,
# ad-hoc proxy for detection difficulty, similar in effect to V
# but requiring no new simulations). Expects a data.table with
# columns alpha, AUC_PR, filter (e.g. "no LD-filtering" /
# "full LD-filtering"), and optionally method, built by looping
# auc_cummax_PR_Cscore_overall() over alpha values and the two
# filter_fun variants -- see example usage below.
#----------------------------------------------------------
plot_AUC_vs_alpha <- function(auc_vs_alpha, title = NULL) {

  has_method <- "method" %in% names(auc_vs_alpha)

  p <- ggplot(auc_vs_alpha, aes(x = alpha, y = AUC_PR, color = filter, linetype = filter)) +
    geom_line(linewidth = 1) + geom_point(size = 2) +
    scale_x_log10()

  if (has_method) p <- p + facet_wrap(~ method)

  p +
    labs(x = "alpha (log scale) -- stricter (harder) to the left",
         y = "AUC of cummax(PR), C_score_threshold & r2_th integrated out",
         color = NULL, linetype = NULL, title = title) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom")
}

#----------------------------------------------------------
# Example usage

#----------------------------------------------------------
# res <- run_and_score_one_sim_file(
#   "./parsed_sim_data/adapt_bgs_chr1_V2_c1_env1.rds",
#   p_cols = p_cols, th_ldw_grid = th_ldw_grid, r2_grid = r2_grid, lmin_grid = lmin_grid,
#   return_intermediates = TRUE
# )
#
# ## (1) per-marker C-score at one r2_th
# marker_C <- compute_marker_C_score(res$scored, p_name = "LFMM", fixed_r2_th = 0.6)
#
# ## (2) per-QTN recall stability (simulation only)
# qtn_C <- compute_QTN_C_score(res$scored, p_name = "LFMM", fixed_r2_th = 0.6, map = res$map)
#
# ## (3) AUC of cummax(PR), POOLING r2_th along with rho/th_ldw/l_min (new default)
# auc_res <- auc_cummax_PR(res$scored, p_name = "LFMM", n_perm = 200, seed = 1)
# auc_res$AUC
# ## old behaviour (one fixed r2_th) still available if wanted:
# auc_res_fixed <- auc_cummax_PR(res$scored, p_name = "LFMM", fixed_r2_th = 0.6, n_perm = 200, seed = 1)
#
# ## (4) C-score-threshold sweep across r2_th, for BOTH methods at once
# C_sweep <- run_C_score_sweep_multi(
#   scored = res$scored, map = res$map, el = res$el, qtn_ld_table = res$qtn_ld_table,
#   p_names = c("EMX", "LFMM"), r2_min_focal = res$thr$r2_min_focal, d_max_focal = res$thr$d_max_focal,
#   r2_th_grid = r2_grid, C_score_grid = seq(0, 1, by = 0.1), l_min = 1
# )
#
# ## (5) integrate r2_th out of the C-score sweep via the same AUC logic,
# ## leaving ONE curve per method over C_score_threshold -- the actual
# ## comparison plot, doubly robust to both the raw grid AND r2_th
# auc_by_Cscore <- auc_cummax_PR_Cscore_sweep(C_sweep, n_perm = 200, seed = 1)
# plot_Cscore_AUC_curve(auc_by_Cscore)
#
# ## (6) ONE overall AUC pooling BOTH C_score_threshold and r2_th --
# ## this is the actual headline number to compare "no LD-filtering"
# ## vs "full LD-filtering" C-scores, and across alpha as a difficulty
# ## proxy. filter_fun is applied to each file's scored table before
# ## the C-score/re-clustering step, so th_ldw == 0 isolates the
# ## no-LD-filtering variant. run_C_score_sweep_all_files() reruns
# ## the per-file pipeline once per (alpha, filter) combination below --
# ## expensive, but each combination needs its own el/qtn_ld_table
# ## since the C-score-threshold step re-clusters markers.
# alpha_vals <- c(0.001, 0.01, 0.05, 0.1, 0.2)
# p_names    <- names(p_cols)   # e.g. c("EMX", "LFMM")
#
# auc_vs_alpha <- rbindlist(lapply(alpha_vals, function(a) {
#   rbindlist(lapply(c("no LD-filtering", "full LD-filtering"), function(filt) {
#
#     ff <- if (filt == "no LD-filtering") function(x) x[th_ldw == 0] else NULL
#
#     sweep <- run_C_score_sweep_all_files(
#       file_paths = files_v2_env1, p_cols = p_cols, p_names = p_names,
#       th_ldw_grid = th_ldw_grid, r2_grid = r2_grid, lmin_grid = lmin_grid,
#       r2_th_grid = r2_grid, C_score_grid = seq(0, 1, by = 0.1),
#       alpha_grid = a, filter_fun = ff
#     )
#
#     auc <- auc_cummax_PR_Cscore_overall(sweep, n_perm = 200, seed = 1)
#     auc[, `:=`(alpha = a, filter = filt)]
#     auc[]
#   }))
# }))
#
# plot_AUC_vs_alpha(auc_vs_alpha)

######################################################
## Genome-wide C-score / AUC across all replicate files
##
## Split into three parts, mirroring score_OR_genome() /
## plot_OR_manhattan_genome_cached():
##   1) fetch_C_score_data_all_files()      -- expensive, run once
##   2) compute_C_score_genome_from_data()  -- cheap, rerun freely
##      with different fixed_r2_th / filter_fun / p_name
##   3) plot_C_score_genome()               -- cheap, replot freely
## compute_C_score_genome() is kept as a thin wrapper of (1)+(2)
## for one-shot convenience when you don't need to iterate.
######################################################

#----------------------------------------------------------
# 1) EXPENSIVE: run the per-file (rho x th_ldw x alpha x r2_th x
#    l_min) grid for ALL replicate files ONCE, and relabel
#    chromosomes/markers for display immediately (Chr1->Chr(2i-1)
#    QTN, Chr2->Chr(2i) neutral) -- both the map AND the cluster/
#    assignment marker names inside `scored`, so the returned
#    per-file `scored` tables are already genome-labeled and
#    directly combinable. No C-score/AUC computed yet, and no
#    filtering applied yet -- both happen cheaply afterwards.
#----------------------------------------------------------
fetch_C_score_data_all_files <- function(file_paths, p_cols, th_ldw_grid, r2_grid, lmin_grid,
                                         alpha_grid = 0.05, bp_th_cluster = 5e5,
                                         rho_r2_focal = 0.75, rho_d_focal = 0.95,
                                         cores = 1, r2_grid_scale = c("r2", "rho")) {

  r2_grid_scale <- match.arg(r2_grid_scale)

  chr_num <- as.integer(gsub("[^0-9]", "",
                             regmatches(basename(file_paths), regexpr("chr[0-9]+", basename(file_paths)))))
  if (any(is.na(chr_num))) stop("Could not parse chromosome number from one or more filenames.")
  file_paths <- file_paths[order(chr_num)]

  maps           <- vector("list", length(file_paths))
  scored_by_file <- vector("list", length(file_paths))

  for (i in seq_along(file_paths)) {
    fp <- file_paths[i]
    message("Processing ", basename(fp))

    res <- tryCatch(
      run_and_score_one_sim_file(
        fp, p_cols = p_cols, th_ldw_grid = th_ldw_grid, r2_grid = r2_grid, lmin_grid = lmin_grid,
        bp_th_cluster = bp_th_cluster, rho_r2_focal = rho_r2_focal, rho_d_focal = rho_d_focal,
        alpha_grid = alpha_grid, cores = cores, r2_grid_scale = r2_grid_scale, return_intermediates = TRUE
      ),
      error = function(e) { message("  FAILED: ", conditionMessage(e)); NULL }
    )

    if (is.null(res)) {
      ## no potential outliers anywhere in this file at ANY setting --
      ## keep the map (so the chromosome still shows up), scored = NULL
      d     <- readRDS(fp)
      map_i <- flag_true_positive_QTNs(d$map)
      scored_i <- NULL
    } else {
      map_i    <- res$map
      scored_i <- res$scored
    }

    ## relabel chromosomes/markers for DISPLAY only
    old_marker <- map_i$marker
    new_chr    <- ifelse(map_i$Chr == "Chr1", paste0("Chr", 2 * i - 1), paste0("Chr", 2 * i))
    old_prefix <- paste0(map_i$Chr, ":")
    new_prefix <- paste0(new_chr, ":")
    new_marker <- mapply(function(m, o, n) sub(o, n, m, fixed = TRUE),
                         old_marker, old_prefix, new_prefix, USE.NAMES = FALSE)
    marker_map <- setNames(new_marker, old_marker)

    map_i[, Chr := new_chr]
    map_i[, marker := new_marker]

    if (!is.null(scored_i)) {
      ## relabel marker names inside the cluster list-columns AND the
      ## assignments_<p_name> qtn field, so scored_i is self-consistent
      ## with the relabeled map without needing to redo this later
      p_names_here <- intersect(names(p_cols), names(scored_i))
      for (nm in p_names_here) {
        scored_i[[nm]] <- lapply(scored_i[[nm]], function(clusters_at_row) {
          lapply(clusters_at_row, function(cl) unname(marker_map[cl]))
        })

        assign_col <- paste0("assignments_", nm)
        if (assign_col %in% names(scored_i)) {
          scored_i[[assign_col]] <- lapply(scored_i[[assign_col]], function(a) {
            if (is.null(a) || nrow(a) == 0) return(a)
            a <- copy(a)
            a[!is.na(qtn), qtn := unname(marker_map[qtn])]
            a[]
          })
        }
      }
    }

    maps[[i]] <- map_i
    scored_by_file[[i]] <- scored_i
  }

  names(scored_by_file) <- basename(file_paths)
  list(maps = maps, scored_by_file = scored_by_file, files = file_paths)
}

#----------------------------------------------------------
# 2) CHEAP: compute C-score + AUC and combine into one genome-
#    wide table, from data already fetched above. `filter_fun`,
#    if supplied, is applied to EACH file's scored table BEFORE
#    computing C-score/AUC -- e.g.
#      filter_fun = function(x) x[alpha <= 0.1]
#    to restrict the alpha (or th_ldw/rho/r2_th/l_min) range
#    post-hoc, without re-running the expensive fetch step.
#----------------------------------------------------------
compute_C_score_genome_from_data <- function(fetched, p_name, fixed_r2_th = NULL,
                                             filter_fun = NULL,
                                             auc_n_perm = 200, auc_seed = NULL) {

  maps     <- fetched$maps
  auc_list <- vector("list", length(maps))

  for (i in seq_along(maps)) {
    map_i    <- copy(maps[[i]])
    scored_i <- fetched$scored_by_file[[i]]

    if (!is.null(filter_fun) && !is.null(scored_i)) {
      scored_i <- filter_fun(scored_i)
    }

    if (is.null(scored_i) || nrow(scored_i) == 0 || !p_name %in% names(scored_i)) {
      ## nothing left for this p_name (either never had outliers, or
      ## filter_fun() filtered every row out) -- C_score = 0 everywhere,
      ## AUC undefined, rather than erroring
      map_i[, C_score := 0]
      auc_val <- NA_real_
    } else {
      marker_C <- compute_marker_C_score(scored_i, p_name = p_name, fixed_r2_th = fixed_r2_th)
      map_i    <- attach_C_score_to_map(map_i, marker_C)
      auc_res  <- auc_cummax_PR(scored_i, p_name = p_name, fixed_r2_th = fixed_r2_th,
                                n_perm = auc_n_perm, seed = auc_seed)
      auc_val  <- auc_res$AUC
    }

    maps[[i]] <- map_i
    ## AUC is a per-FILE quantity -- tag it to BOTH of this file's
    ## display chromosomes so it can be shown in either facet's strip label
    auc_list[[i]] <- data.table(Chr = unique(map_i$Chr), file = basename(fetched$files[i]), AUC = auc_val)
  }

  combined_map <- rbindlist(maps, fill = TRUE)
  chr_order <- unique(combined_map$Chr)
  chr_order <- chr_order[order(as.integer(gsub("[^0-9]", "", chr_order)))]
  combined_map[, Chr := factor(Chr, levels = chr_order)]
  setorder(combined_map, Chr, Pos)
  combined_map[, indx := .I]

  list(map = combined_map, chr_order = chr_order, auc_by_file = rbindlist(auc_list, fill = TRUE))
}

#----------------------------------------------------------
# 3) Thin wrapper: fetch + compute in one call, for one-shot use
#    when you don't need to try several fixed_r2_th/filter_fun
#    combinations. For iterating, call fetch_C_score_data_all_files()
#    once and compute_C_score_genome_from_data() repeatedly instead.
#----------------------------------------------------------
compute_C_score_genome <- function(file_paths, p_cols, p_name, th_ldw_grid, r2_grid, lmin_grid,
                                   fixed_r2_th = NULL, filter_fun = NULL, bp_th_cluster = 5e5,
                                   rho_r2_focal = 0.75, rho_d_focal = 0.95,
                                   alpha_grid = 0.05, cores = 1, r2_grid_scale = c("r2", "rho"),
                                   auc_n_perm = 200, auc_seed = NULL) {

  fetched <- fetch_C_score_data_all_files(
    file_paths = file_paths, p_cols = p_cols, th_ldw_grid = th_ldw_grid,
    r2_grid = r2_grid, lmin_grid = lmin_grid, alpha_grid = alpha_grid,
    bp_th_cluster = bp_th_cluster, rho_r2_focal = rho_r2_focal, rho_d_focal = rho_d_focal,
    cores = cores, r2_grid_scale = r2_grid_scale
  )

  compute_C_score_genome_from_data(
    fetched, p_name = p_name, fixed_r2_th = fixed_r2_th, filter_fun = filter_fun,
    auc_n_perm = auc_n_perm, auc_seed = auc_seed
  )
}

#----------------------------------------------------------
# Cheap replot step: any value_col/color_by combination, no
# recomputation. Defaults to plotting AND coloring by C_score.
# AUC is shown in each facet's strip label by default (same
# value repeated for a file's paired QTN/neutral chromosomes).
#----------------------------------------------------------
plot_C_score_genome <- function(C_score_genome, value_col = "C_score", value_label = "C-score",
                                color_by = "C_score", wes_palette_name = "Zissou1",
                                show_auc_in_strip = TRUE, title = NULL, nrow_facets = 2) {

  map_plot <- copy(C_score_genome$map)
  map_plot[, y_val := get(value_col)]
  map_plot[, status := "not called"]   ## unused unless color_by == "status"
  if (!"indx" %in% names(map_plot)) map_plot[, indx := .I]

  chr_levels <- as.character(C_score_genome$chr_order)

  if (show_auc_in_strip) {
    auc_lookup <- setNames(round(C_score_genome$auc_by_file$AUC, 3), C_score_genome$auc_by_file$Chr)
    strip_labels <- sprintf("%s (AUC=%s)", chr_levels,
                            ifelse(is.na(auc_lookup[chr_levels]), "NA", auc_lookup[chr_levels]))
    map_plot[, Chr := factor(Chr, levels = chr_levels, labels = strip_labels)]
    chr_levels <- strip_labels
  }

  .render_OR_manhattan(map_plot, title = title, chr_levels = chr_levels,
                       nrow_facets = nrow_facets, value_label = value_label,
                       color_by = color_by, wes_palette_name = wes_palette_name)
}

#----------------------------------------------------------
# Example usage
#----------------------------------------------------------
# ## one-shot convenience (fetch + compute in one call):
# C_score_genome <- compute_C_score_genome(
#   file_paths = files_v2_env1, p_cols = p_cols, p_name = "LFMM",
#   th_ldw_grid = th_ldw_grid, r2_grid = r2_grid, lmin_grid = lmin_grid
#   ## fixed_r2_th = NULL by default -- pools r2_th into the C-score/AUC too
# )
#
# ## OR, to try several fixed_r2_th / filters without re-running the
# ## expensive per-file pipeline each time:
# fetched <- fetch_C_score_data_all_files(
#   file_paths = files_v2_env1, p_cols = p_cols,
#   th_ldw_grid = th_ldw_grid, r2_grid = r2_grid, lmin_grid = lmin_grid,
#   alpha_grid = c(0.01, 0.05, 0.1, 0.2)
# )
#
# C_score_genome_all      <- compute_C_score_genome_from_data(fetched, p_name = "LFMM")
# C_score_genome_strict_a <- compute_C_score_genome_from_data(
#   fetched, p_name = "LFMM", filter_fun = function(x) x[alpha <= 0.05]
# )
# C_score_genome_r06      <- compute_C_score_genome_from_data(
#   fetched, p_name = "LFMM", fixed_r2_th = 0.6
# )
#
# plot_C_score_genome(C_score_genome_all)                       # C_score as both y and color
# plot_C_score_genome(C_score_genome_all, value_col = "lfmm_F", value_label = "LFMM F",
#                     color_by = "C_score")                     # F-stat y-axis, C_score coloring


