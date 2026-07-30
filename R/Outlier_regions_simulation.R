######################################################
## SIMULATION-SPECIFIC outlier-region (OR) pipeline
##
## Requires outlier_regions_generic.R to be sourced/loaded first --
## everything here builds on its clustering/C-score engine, adding
## the layer that only simulated data can provide: known true QTNs,
## and therefore TP/FP/FN, Precision/Recall/PR, and AUC.
##
## source("outlier_regions_generic.R")   # if not already loaded
######################################################

## these two come from compute_ld_structure.R -- re-declared here for
## standalone use if that script isn't sourced
if (!exists("ld_from_rho")) {
  ld_from_rho <- function(b, c = 1, rho) b + (c - b) * (1 - rho)
}
if (!exists("d_from_rho")) {
  d_from_rho <- function(a, rho) rho / (a * (1 - rho))
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
# MARKER-LEVEL (no outlier-region clustering) true-positive flag.
# A marker is a true positive if it's within the SAME r2_min_focal/
# d_max_focal thresholds (rho=0.75/0.95, via qtn_ld_table) of ANY
# true-positive QTN (true_pos_QTN == TRUE) -- i.e. the identical
# "close enough to count" rule already used for OR-level focal-QTN
# assignment (assign_focal_qtn_one_OR()), just applied per-marker
# instead of per-cluster. This keeps the marker-level comparison
# (C-score vs raw alpha, no clustering) and the OR-level comparison
# (LD-filtering vs l_min-only, WITH clustering) using one consistent
# ground-truth definition, differing only in whether clustering/
# focal-QTN-duplicate-resolution sits between "called" and "true".
#----------------------------------------------------------
flag_marker_true_positive <- function(map, qtn_ld_table, r2_min_focal, d_max_focal) {
  true_pos_qtns <- map[true_pos_QTN == TRUE, marker]

  close_to_true_qtn <- qtn_ld_table[
    qtn_marker %in% true_pos_qtns & r2 > r2_min_focal & dist_bp < d_max_focal,
    unique(marker)
  ]

  map <- copy(map)
  map[, true_pos_marker := marker %in% close_to_true_qtn]
  map[]
}

#----------------------------------------------------------
# MARKER-LEVEL (no outlier-region clustering, no focal-QTN
# duplicate resolution) TP/FP/FN/Precision/Recall/PR. Each called
# marker is independently TP or FP based on flag_marker_true_positive()'s
# proximity flag -- use this (with pr_curve_auc()) for the
# C-score-vs-raw-alpha comparison. Keep evaluate_ORs() (OR-based)
# for the LD-filtering-vs-l_min-only comparison, which DOES still
# consider outlier regions.
#----------------------------------------------------------
evaluate_markers <- function(called_markers, map) {
  true_pos_markers <- map[true_pos_marker == TRUE, marker]

  TP <- length(intersect(called_markers, true_pos_markers))
  FP <- length(setdiff(called_markers, true_pos_markers))
  FN <- length(setdiff(true_pos_markers, called_markers))

  Precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  Recall    <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  PR        <- if (!is.na(Precision) && !is.na(Recall)) Precision * Recall else NA_real_

  list(TP = TP, FP = FP, FN = FN, Precision = Precision, Recall = Recall, PR = PR)
}

#----------------------------------------------------------
# Marker-level C-score-threshold curve for ONE file: NO clustering,
# NO focal-QTN duplicate resolution at THIS step -- each marker with
# C_score >= C_th is independently evaluated via evaluate_markers().
# This function itself takes a precomputed marker_C_score table and
# just thresholds it, so it has no dependence on r2_th/l_min here --
# but how marker_C_score was BUILT matters a great deal: r2_th and
# l_min still shape the underlying C-score whenever l_min > 1 is part
# of the pooled grid (a marker only survives an l_min > 1 filter if
# it clustered with enough neighbors at that r2_th), which is exactly
# what gives C-score its expected advantage over raw alpha -- isolated
# significant SNPs get systematically lower C-scores. See
# run_marker_C_score_curve_from_results() for how the grid is pooled.
#----------------------------------------------------------
run_marker_C_score_curve <- function(marker_C_score, map, C_score_grid) {
  rbindlist(lapply(C_score_grid, function(C_th) {
    called <- marker_C_score[C_score >= C_th, marker]
    res <- evaluate_markers(called, map)
    data.table(C_score_threshold = C_th, TP = res$TP, FP = res$FP, FN = res$FN,
               Precision = res$Precision, Recall = res$Recall, PR = res$PR)
  }))
}

#----------------------------------------------------------
# Marker-level raw-alpha baseline curve for ONE file: extracts, for
# each alpha value already swept in `scored` (at th_ldw = 0, l_min = 1
# -- so cluster membership can't affect which markers are "called"),
# the plain FDR-significant marker set, evaluated marker-by-marker via
# evaluate_markers() -- no clustering involved, matching
# run_marker_C_score_curve() above for a fair, same-footing comparison.
# l_min = 1 here is DELIBERATE, not incidental: this is meant to be
# the true floor -- raw significance with NEITHER LD-filtering NOR
# any l_min-based clustering cleanup -- so it's the correct baseline
# against C-score, which (via run_marker_C_score_curve_from_results(),
# pooling the full l_min_grid) DOES benefit from l_min. r2_th is fixed
# to a single arbitrary value present in `scored`, since (at l_min = 1
# specifically) every r2_th gives an identical marker set -- this
# does NOT generalize to l_min > 1, where r2_th matters again.
#----------------------------------------------------------
extract_marker_alpha_baseline <- function(scored, map, p_name, sel_th_ldw = 0) {
  sub <- scored[th_ldw == sel_th_ldw & l_min == 1]
  if (nrow(sub) == 0) stop("No rows in `scored` at th_ldw = ", sel_th_ldw, " & l_min = 1.")
  sub <- sub[r2_th == sub$r2_th[1]]   ## arbitrary but fixed -- inert at l_min = 1 only

  rbindlist(lapply(seq_len(nrow(sub)), function(i) {
    called <- unlist(sub[[p_name]][[i]])
    res <- evaluate_markers(called, map)
    data.table(alpha = sub$alpha[i], TP = res$TP, FP = res$FP, FN = res$FN,
               Precision = res$Precision, Recall = res$Recall, PR = res$PR)
  }))
}

#----------------------------------------------------------
# Multi-file convenience wrappers: build the marker-level C-score
# curve / alpha baseline for every file in a results list (e.g.
# from run_and_score_all(..., return_intermediates = TRUE)), tagged
# with `file`. flag_marker_true_positive() is cheap (a lookup on
# already-cached qtn_ld_table/thr, no clustering/re-computation), so
# it's computed fresh per file here rather than needing to be cached
# separately. Pool across files afterward with a simple
#   curves[, .(TP = sum(TP), FP = sum(FP), FN = sum(FN)), by = <threshold col>]
# before pr_curve_auc()/plot_pr_curve() -- same pooling principle as
# aggregate_PR(), just at the marker level instead of the OR level.
#
# run_marker_C_score_curve_from_results() pools the FULL grid by
# default (fixed_r2_th = NULL, filter_fun = NULL) -- deliberately
# including the full l_min_grid, not just l_min = 1, since that's
# exactly the mechanism giving C-score its expected edge over raw
# alpha: an isolated significant SNP only survives the l_min = 1
# grid points, never l_min > 1 (it has no neighbors to cluster
# with), so it accumulates a low C-score, while a genuinely
# clustered true-signal SNP survives l_min = 1, 5, 10, 20 alike and
# gets a high one. Use filter_fun to isolate sub-questions, e.g.
# function(x) x[th_ldw == 0] to ask "does l_min/clustering alone,
# with NO LD-filtering, already beat raw alpha" -- the natural
# 3-way ladder is: raw alpha (l_min = 1 only, via
# extract_marker_alpha_baseline(), the true floor) < C-score with
# th_ldw = 0 only < C-score with the full th_ldw range.
#----------------------------------------------------------
run_marker_C_score_curve_from_results <- function(results_list, p_name, C_score_grid,
                                                  fixed_r2_th = NULL, filter_fun = NULL) {

  rbindlist(lapply(names(results_list), function(nm) {
    res <- results_list[[nm]]
    if (is.null(res)) return(NULL)

    map_i <- flag_marker_true_positive(res$map, res$qtn_ld_table, res$thr$r2_min_focal, res$thr$d_max_focal)
    marker_C <- compute_marker_C_score(res$scored, p_name = p_name, fixed_r2_th = fixed_r2_th,
                                       filter_fun = filter_fun)

    curve <- run_marker_C_score_curve(marker_C, map_i, C_score_grid)
    curve[, file := nm]
    curve[]
  }), fill = TRUE)
}

extract_marker_alpha_baseline_from_results <- function(results_list, p_name, sel_th_ldw = 0) {
  rbindlist(lapply(names(results_list), function(nm) {
    res <- results_list[[nm]]
    if (is.null(res)) return(NULL)

    map_i <- flag_marker_true_positive(res$map, res$qtn_ld_table, res$thr$r2_min_focal, res$thr$d_max_focal)
    curve <- extract_marker_alpha_baseline(res$scored, map_i, p_name = p_name, sel_th_ldw = sel_th_ldw)
    curve[, file := nm]
    curve[]
  }), fill = TRUE)
}

#----------------------------------------------------------
# MARKER-LEVEL score+label extraction, for use with a STANDARD
# PR-AUC package (PRROC::pr.curve()) instead of the custom
# pr_curve_auc()/threshold-table approach above. Since the
# marker-level comparison already reduces to "one continuous score
# + one binary true/false label per marker" (via
# flag_marker_true_positive()), this is exactly the classic
# per-instance classification setup PRROC expects -- no threshold
# table needed, no manual trapezoidal integration.
#
# Pooling across files here means concatenating (score, label)
# PAIRS across all files, NOT matching markers by name -- marker
# name collisions across files (the "Chr1:12345" problem discussed
# earlier) are a non-issue here, since each row is one independent
# instance regardless of which file it came from; we're never
# asking "what is THIS marker's score", only "does a high score
# correlate with true-positive status" in aggregate.
#----------------------------------------------------------

#' C-score as the per-marker score. Includes ALL markers (not just
#' ones ever called -- C_score = 0 for never-called markers, via
#' attach_C_score_to_map()), since PRROC needs every instance
#' represented, not just positive calls.
get_marker_Cscore_labels <- function(res, p_name, fixed_r2_th = NULL, filter_fun = NULL) {
  map_i <- flag_marker_true_positive(res$map, res$qtn_ld_table, res$thr$r2_min_focal, res$thr$d_max_focal)
  marker_C <- compute_marker_C_score(res$scored, p_name = p_name, fixed_r2_th = fixed_r2_th, filter_fun = filter_fun)
  map_i <- attach_C_score_to_map(map_i, marker_C)
  map_i[, .(marker, score = C_score, label = true_pos_marker)]
}

#' Raw significance (-log10(p)) as the per-marker score, straight
#' from the ORIGINAL p-value column on map -- the "raw alpha"
#' comparator, at the finest possible resolution (every observed
#' p-value is a candidate threshold, not just a discrete alpha_grid).
get_marker_alpha_labels <- function(res, p_col) {
  map_i <- flag_marker_true_positive(res$map, res$qtn_ld_table, res$thr$r2_min_focal, res$thr$d_max_focal)
  map_i[, .(marker, score = -log10(get(p_col)), label = true_pos_marker)]
}

get_marker_Cscore_labels_from_results <- function(results_list, p_name, fixed_r2_th = NULL, filter_fun = NULL) {
  rbindlist(lapply(names(results_list), function(nm) {
    res <- results_list[[nm]]
    if (is.null(res)) return(NULL)
    tab <- get_marker_Cscore_labels(res, p_name, fixed_r2_th = fixed_r2_th, filter_fun = filter_fun)
    tab[, file := nm]
    tab[]
  }), fill = TRUE)
}

get_marker_alpha_labels_from_results <- function(results_list, p_col) {
  rbindlist(lapply(names(results_list), function(nm) {
    res <- results_list[[nm]]
    if (is.null(res)) return(NULL)
    tab <- get_marker_alpha_labels(res, p_col)
    tab[, file := nm]
    tab[]
  }), fill = TRUE)
}

#----------------------------------------------------------
# Standard AUC-PR via PRROC::pr.curve(), given a (score, label)
# table from any of the functions above. PRROC's convention:
# scores.class0 = scores for the POSITIVE class, scores.class1 =
# scores for the NEGATIVE class.
#
# NOTE: PRROC returns BOTH $auc.integral (linear interpolation
# between observed points) and $auc.davis.goadrich (the
# interpolation-corrected version from Davis & Goadrich 2006).
# Linear interpolation is known to be misleading specifically for
# PR curves (unlike ROC), and the difference is largest exactly
# under severe class imbalance -- true QTNs vs. background markers
# here being a clear case of that. Prefer $auc.davis.goadrich
# unless you have a specific reason to use $auc.integral.
#----------------------------------------------------------
auc_pr_PRROC <- function(score_label_table, curve = FALSE) {
  if (!requireNamespace("PRROC", quietly = TRUE)) {
    stop("This function needs the PRROC package: install.packages('PRROC').")
  }
  pos <- score_label_table[label == TRUE, score]
  neg <- score_label_table[label == FALSE, score]
  if (length(pos) == 0 || length(neg) == 0) {
    stop("Need at least one TRUE- and one FALSE-labeled marker to compute a PR curve.")
  }
  PRROC::pr.curve(scores.class0 = pos, scores.class1 = neg, curve = curve)
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
# Score a full outliers_dt grid (rbind'ed over param_grid, from
# run_and_cluster()'s generic $outliers). Retains a per-grid-point
# "assignments_<nm>" list-column (cluster-level TP/FP labels) so
# plot_OR_manhattan() can pull directly from the scored table
# without recomputing.
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
# Per-QTN "recall stability" (SIMULATION ONLY -- needs
# true_pos_QTN and the assignments_<p_name> list-column).
# For every known true-positive QTN, including ones NEVER
# recovered (C_score = 0), the fraction of grid points where
# it was correctly assigned as SOME OR's focal QTN.
# fixed_r2_th = NULL (default) pools r2_th in too, matching
# compute_marker_C_score()'s default. filter_fun, if supplied, is
# applied to `scored` FIRST (e.g. function(x) x[th_ldw == 0] to
# isolate the no-LD-filtering C-score, or function(x) x[alpha ==
# 0.01] as a difficulty proxy).
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
# AUC of cummax(metric) under random search order -- generic MATH
# (doesn't need truth), but kept here since its only real use in
# this pipeline is summarizing PR (a simulation-only metric).
# Repeatedly shuffle `values`, track the running best value as if
# trying them one at a time in that random order, average across
# shuffles, and take the area under the resulting curve (x =
# fraction tried, y = expected best-so-far). Starts at (0,0) so
# AUC is comparable across different-sized value sets.
#
# Use this for UNORDERED, multi-dimensional nuisance parameter
# combinations (rho x th_ldw x l_min x alpha x r2_th) where there
# is no single "more to less stringent" ordering to sweep instead.
# For a single, naturally-ordered threshold (e.g. C-score alone,
# or alpha alone, at fixed r2_th), use a standard monotonic
# PR-curve AUC instead (order by threshold, trapezoidal rule) --
# NOT YET IMPLEMENTED here; needed for the C-score-vs-raw-alpha
# baseline comparison agreed on separately.
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
# alpha) grid -- the "modified/random-search AUC-PR" for the raw
# LD-filtering method, integrating across ALL of its own nuisance
# parameters at once (they have no natural ordering as a group).
# POOLS r2_th by default (fixed_r2_th = NULL); ALSO pools TP/FP/FN
# across whatever rows share the same grid point BEFORE computing
# PR, same rationale as aggregate_PR(): a single replicate's grid
# point very often has TP=0, which mechanically forces PR=0.
# Feeding those raw per-file zeros into the random-search AUC
# understates real performance -- pooling counts first, then
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
# Pooled (micro-averaged) Precision/Recall/PR across replicates
#
# Sums TP/FP/FN across replicate files WITHIN each grid point
# (rho, th_ldw, r2_th, l_min, alpha[, group_vars]) before computing
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

#----------------------------------------------------------
# C-score-threshold sweep: builds on the GENERIC clustering engine
# (run_ORs_C_score_sweep() / get_ORs_at_C_score_threshold() in
# outlier_regions_generic.R) and layers truth-scoring (evaluate_ORs)
# on top, at ONE r2_th. `el` and `qtn_ld_table` should come from the
# SAME file's run_and_score_one_sim_file(..., return_intermediates =
# TRUE) call that produced `scored` / the marker_C_score table --
# reused as-is, not recomputed.
#
# r2_th_linear: the actual raw r^2 used for clustering, if r2_th is
# on the rho scale -- defaults to r2_th itself (correct only if
# r2_th is ALREADY raw r^2). When calling this directly (not via
# run_C_score_sweep(), which looks this up automatically from
# `scored`), pass r2_th_linear explicitly if you're working on the
# rho scale, e.g. r2_th_linear = unique(scored[r2_th == sel_r2_th, r2_th_linear]).
#----------------------------------------------------------
run_C_score_threshold_grid <- function(marker_C_score, map, el, qtn_ld_table,
                                       r2_min_focal, d_max_focal,
                                       C_score_grid, r2_th, r2_th_linear = r2_th,
                                       l_min = 1, bp_th = Inf) {

  or_table <- get_ORs_at_C_score_threshold(marker_C_score, el, C_score_grid, r2_th = r2_th,
                                           r2_th_linear = r2_th_linear, l_min = l_min, bp_th = bp_th)

  rbindlist(lapply(seq_len(nrow(or_table)), function(i) {
    cls <- or_table$clusters[[i]]
    res <- evaluate_ORs(cls, map, qtn_ld_table, r2_min_focal, d_max_focal)
    data.table(C_score_threshold = or_table$C_score_threshold[i], r2_th = or_table$r2_th[i],
               l_min = or_table$l_min[i],
               TP = res$TP, FP = res$FP, FN = res$FN,
               Precision = res$Precision, Recall = res$Recall, PR = res$PR)
  }))
}

#----------------------------------------------------------
# Full sweep: r2_th x C_score_threshold, with truth-scoring, built
# on the generic run_ORs_C_score_sweep() engine.
#----------------------------------------------------------
run_C_score_sweep <- function(scored, map, el, qtn_ld_table, p_name,
                              r2_min_focal, d_max_focal,
                              r2_th_grid, C_score_grid, l_min = 1, bp_th = Inf,
                              filter_fun = NULL) {

  or_table <- run_ORs_C_score_sweep(scored, el, p_name, r2_th_grid, C_score_grid,
                                    l_min = l_min, bp_th = bp_th, filter_fun = filter_fun)

  rbindlist(lapply(seq_len(nrow(or_table)), function(i) {
    cls <- or_table$clusters[[i]]
    res <- evaluate_ORs(cls, map, qtn_ld_table, r2_min_focal, d_max_focal)
    data.table(C_score_threshold = or_table$C_score_threshold[i], r2_th = or_table$r2_th[i],
               l_min = or_table$l_min[i],
               TP = res$TP, FP = res$FP, FN = res$FN,
               Precision = res$Precision, Recall = res$Recall, PR = res$PR)
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
# Full pipeline for ONE simulation file: loads the file, calls the
# GENERIC run_and_cluster() for FDR+clustering, then layers
# truth-scoring (QTN-LD, focal-QTN thresholds, score_outlier_grid)
# on top. This is the simulation-specific "file orchestrator" --
# for empirical data, load map/GTs/ld_ws your own way and call
# run_and_cluster() directly (no truth layer needed/available).
#----------------------------------------------------------
run_and_score_one_sim_file <- function(file_path, p_cols, th_ldw_grid, r2_grid, lmin_grid,
                                       bp_th_cluster = 5e5, rho_r2_focal = 0.75,
                                       rho_d_focal = 0.95, alpha_grid = 0.05, cores = 1,
                                       rho_grid = NULL, r2_grid_scale = c("rho", "r2"),
                                       return_intermediates = FALSE) {

  r2_grid_scale <- match.arg(r2_grid_scale)

  d <- readRDS(file_path)
  map <- flag_true_positive_QTNs(d$map)
  GTs <- d$GTs
  ld_ws <- d$ld_ws
  ## LDscnR::compute_ld_w() returns ld_ws with dimnames = list(NULL, "rho_*"),
  ## i.e. no rownames, but the OR engine indexes ld_ws by marker name. Rows are
  ## row-aligned with map (same order, same nrow -- verified: map$ld_w_095 ==
  ## ld_ws[,"rho_0.95"] exactly), so recover the marker rownames from map.
  if (is.null(rownames(ld_ws))) rownames(ld_ws) <- d$map$marker
  LD_decay <- d$LD_decay

  ds <- LD_decay$decay_sum[Chr == "Chr1"]
  if (nrow(ds) != 1) {
    stop("Expected exactly one decay_sum row for chr 'Chr1', found ", nrow(ds),
         ". Check LD_decay$decay_sum$Chr labels.")
  }

  generic_res <- run_and_cluster(
    map = map, GTs = GTs, ld_ws = ld_ws, p_cols = p_cols,
    th_ldw_grid = th_ldw_grid, r2_grid = r2_grid, lmin_grid = lmin_grid,
    bp_th_cluster = bp_th_cluster, alpha_grid = alpha_grid, cores = cores,
    rho_grid = rho_grid, r2_grid_scale = r2_grid_scale,
    decay_abc = list(a = ds$a, b = ds$b, c = ds$c)
  )

  if (is.null(generic_res)) {
    message("No potential outliers for ", basename(file_path))
    return(NULL)
  }

  qtn_ld_table <- precompute_QTN_LD(GTs = GTs, map = map, candidate_markers = generic_res$potential_outliers,
                                    max_bp = 2e6, cores = cores)

  thr <- get_focal_qtn_thresholds(LD_decay, chr_qtn = "Chr1", rho_r2 = rho_r2_focal, rho_d = rho_d_focal)

  scored <- score_outlier_grid(generic_res$outliers, map, qtn_ld_table, names(p_cols),
                               r2_min_focal = thr$r2_min_focal, d_max_focal = thr$d_max_focal)
  scored[, file := basename(file_path)]

  if (return_intermediates) {
    return(list(scored = scored[], map = map, el = generic_res$el, qtn_ld_table = qtn_ld_table,
                thr = thr, bp_th_cluster = generic_res$bp_th_cluster,
                bp_th_cluster_rho = generic_res$bp_th_cluster_rho))
  }
  scored[]
}

#----------------------------------------------------------
# Attach filename-derived metadata (sim, bgs, Chr, V, c, env) to
# any table with a `file` column, matching extraxt_params()'s
# convention (sim_bgs_Chr_V_c_env.rds). Shared helper -- used by
# run_and_score_all() below, and directly usable on
# run_C_score_sweep_from_results()'s output too.
#----------------------------------------------------------
attach_sim_metadata <- function(dt, file_col = "file") {
  base  <- sub("\\.rds$", "", dt[[file_col]])
  parts <- tstrsplit(base, "_", fixed = TRUE)
  dt[, `:=`(sim = parts[[1]], bgs = parts[[2]], Chr = parts[[3]],
            V = parts[[4]], c = parts[[5]], env = parts[[6]])]
  dt[]
}

#----------------------------------------------------------
# Loop over all simulation files. With return_intermediates = TRUE
# (the point of this option), returns a NAMED LIST (one element per
# file) of the FULL per-file results (scored/map/el/qtn_ld_table/
# thr/...) instead of a single stacked `scored` table -- so you can
# later run compute_marker_C_score()/run_C_score_sweep_multi() (via
# run_C_score_sweep_from_results(), below) with ANY r2_th_grid/
# C_score_grid you like, as many times as you like, without
# re-running the expensive per-file pipeline at all. Each file's
# $scored table gets sim/bgs/Chr/V/c/env metadata attached either way.
#----------------------------------------------------------
run_and_score_all <- function(parsed_folder, p_cols, th_ldw_grid, r2_grid, lmin_grid,
                              bp_th_cluster = 5e5, rho_r2_focal = 0.75, rho_d_focal = 0.95,
                              alpha_grid = 0.05, cores = 1, r2_grid_scale = c("rho", "r2"),
                              return_intermediates = FALSE,pattern="\\.rds$") {

  r2_grid_scale <- match.arg(r2_grid_scale)
  files <- list.files(parsed_folder, pattern = pattern, full.names = TRUE)

  results <- lapply(files, function(f) {
    message("Processing ", basename(f))
    tryCatch(
      run_and_score_one_sim_file(f, p_cols, th_ldw_grid, r2_grid, lmin_grid,
                                 bp_th_cluster, rho_r2_focal, rho_d_focal, alpha_grid, cores,
                                 r2_grid_scale = r2_grid_scale,
                                 return_intermediates = return_intermediates),
      error = function(e) {
        message("  FAILED (", basename(f), "): ", conditionMessage(e))
        NULL
      }
    )
  })
  names(results) <- basename(files)

  if (return_intermediates) {
    ## attach metadata onto each file's $scored table individually --
    ## nothing else gets touched/stacked, so you pay nothing for files
    ## or analyses you don't end up using
    for (nm in names(results)) {
      if (!is.null(results[[nm]])) results[[nm]]$scored <- attach_sim_metadata(results[[nm]]$scored)
    }
    return(results)
  }

  all_scored <- rbindlist(results, fill = TRUE)
  attach_sim_metadata(all_scored)
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
                                        alpha_grid = 0.05, cores = 1, r2_grid_scale = c("rho", "r2"),
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
# Extract and stack just the $scored tables from a results list
# (e.g. from run_and_score_all(..., return_intermediates = TRUE)),
# for raw-grid analyses (aggregate_PR(), auc_cummax_PR()) that
# don't need el/qtn_ld_table/map -- cheap, no recomputation.
#----------------------------------------------------------
combine_scored_list <- function(results_list) {
  scored_list <- lapply(results_list, function(r) if (is.null(r)) NULL else r$scored)
  rbindlist(scored_list, fill = TRUE)
}

#----------------------------------------------------------
# Run the C-score-threshold sweep across all files in a results
# list (e.g. from run_and_score_all(..., return_intermediates =
# TRUE)), WITHOUT rerunning the expensive per-file pipeline --
# uses each file's already-computed map/el/qtn_ld_table/thr
# directly. Call this as many times as you like with DIFFERENT
# r2_th_grid/C_score_grid values to explore ad hoc, all cheaply,
# from the SAME one-time expensive fetch -- this is the actual
# "run the expensive analyses only once" pattern: commit to
# th_ldw_grid/r2_grid/lmin_grid/alpha_grid up front (needed for
# the per-file FDR/clustering step itself), but decide the
# C-score-specific r2_th_grid/C_score_grid later, as many times
# as needed, for free.
#----------------------------------------------------------
run_C_score_sweep_from_results <- function(results_list, p_names, r2_th_grid, C_score_grid,
                                           l_min = 1, bp_th = Inf, filter_fun = NULL) {

  rbindlist(lapply(names(results_list), function(nm) {
    res <- results_list[[nm]]
    if (is.null(res)) return(NULL)

    sweep <- run_C_score_sweep_multi(
      scored = res$scored, map = res$map, el = res$el, qtn_ld_table = res$qtn_ld_table,
      p_names = p_names, r2_min_focal = res$thr$r2_min_focal, d_max_focal = res$thr$d_max_focal,
      r2_th_grid = r2_th_grid, C_score_grid = C_score_grid, l_min = l_min, bp_th = bp_th,
      filter_fun = filter_fun
    )
    sweep[, file := nm]
    sweep[]
  }), fill = TRUE)
}

#----------------------------------------------------------
# Collapse the C-score sweep's r2_th axis via the SAME random-
# search-AUC logic as auc_cummax_PR() -- for each C-score
# threshold (and `method` if present), pool PR across r2_th_grid
# and report one AUC. Keeps C_score_threshold as a reportable
# axis (like sweeping alpha) and only integrates out r2_th.
#----------------------------------------------------------
auc_cummax_PR_Cscore_sweep <- function(C_sweep, n_perm = 200, seed = NULL) {

  group_cols <- intersect(c("C_score_threshold", "method"), names(C_sweep))
  pool_cols  <- intersect(c("r2_th", "l_min"), names(C_sweep))

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
# r2_th into a single random-search space. NOTE: per the most
# recent design discussion, this is likely NOT the right tool for
# the final C-score-vs-raw-alpha comparison (C-score, given fixed
# r2_th/clustering thresholds, is a single ordered threshold and
# should get a standard monotonic PR-curve AUC instead) -- kept
# here for the "explore th_ldw=0 vs full range" exploratory use
# discussed, not necessarily for the manuscript's headline figure.
#----------------------------------------------------------
auc_cummax_PR_Cscore_overall <- function(C_sweep, n_perm = 200, seed = NULL) {

  group_cols <- intersect("method", names(C_sweep))
  pool_cols  <- intersect(c("C_score_threshold", "r2_th", "l_min"), names(C_sweep))

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
# Standard monotonic precision-recall-curve AUC. Unlike
# auc_cummax()/auc_cummax_PR() (which handle an UNORDERED,
# multi-dimensional set of nuisance parameters via random-search
# shuffling), this is for a SINGLE, naturally ordered threshold --
# e.g. C-score at a fixed r2_th/l_min, or raw alpha at that same
# fixed r2_th/l_min -- where "more/less stringent" has a real
# meaning and a classical PR curve applies directly: order by
# Recall, trapezoidal rule, no shuffling needed. This is the tool
# for the C-score-vs-raw-alpha baseline comparison, computed
# identically for both so they're on equal footing.
#
# `dt` needs TP/FP/FN columns (default names TP/FP/FN; override
# via TP_col/FP_col/FN_col if needed), one row per threshold value
# already evaluated (e.g. one row per C_score_threshold from
# run_C_score_threshold_grid() at one r2_th, or one row per alpha
# from a matching alpha-only sweep). Ties in Recall are collapsed
# to the single highest Precision seen at that Recall (standard PR-
# curve convention -- a curve should never show precision
# decreasing at constant or increasing recall purely because of
# threshold-ordering noise). No extrapolation beyond the observed
# range -- the AUC reflects only the span of thresholds actually
# evaluated, not an assumed point at Recall=0 or Recall=1; if you
# want the two curves (C-score vs alpha) comparable in absolute
# terms, make sure both sweeps span a similar Recall range.
#----------------------------------------------------------
pr_curve_auc <- function(dt, TP_col = "TP", FP_col = "FP", FN_col = "FN") {

  TP <- dt[[TP_col]]; FP <- dt[[FP_col]]; FN <- dt[[FN_col]]

  Precision <- ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_)
  Recall    <- ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)

  curve <- data.table(Recall = Recall, Precision = Precision)
  curve <- curve[!is.na(Recall) & !is.na(Precision)]

  if (nrow(curve) < 2) {
    return(list(curve = curve[], AUC = NA_real_, n_thresholds = nrow(curve)))
  }

  ## collapse ties in Recall to the highest Precision observed at that Recall
  curve <- curve[, .(Precision = max(Precision)), by = Recall]
  setorder(curve, Recall)

  AUC <- sum(diff(curve$Recall) * (curve$Precision[-1] + curve$Precision[-nrow(curve)]) / 2)

  list(curve = curve[], AUC = AUC, n_thresholds = nrow(curve))
}

#----------------------------------------------------------
# Convenience plot: the classical PR curve itself (Recall on x,
# Precision on y), one line per group if `group_col` is supplied
# (e.g. "method", or a column distinguishing "C-score" vs "raw
# alpha" baseline curves). Useful alongside pr_curve_auc() to see
# *where* along the curve two methods differ, not just the AUC.
#----------------------------------------------------------
plot_pr_curve <- function(dt, group_col = NULL, TP_col = "TP", FP_col = "FP", FN_col = "FN", title = NULL) {

  if (is.null(group_col)) {
    res <- pr_curve_auc(dt, TP_col = TP_col, FP_col = FP_col, FN_col = FN_col)
    curve <- res$curve
    p <- ggplot(curve, aes(x = Recall, y = Precision)) +
      geom_line(linewidth = 1, color = "steelblue4") + geom_point(size = 2, color = "steelblue4")
  } else {
    curves <- rbindlist(lapply(unique(dt[[group_col]]), function(g) {
      res <- pr_curve_auc(dt[get(group_col) == g], TP_col = TP_col, FP_col = FP_col, FN_col = FN_col)
      res$curve[, (group_col) := g]
      res$curve[]
    }))
    p <- ggplot(curves, aes(x = Recall, y = Precision, color = .data[[group_col]])) +
      geom_line(linewidth = 1) + geom_point(size = 2)
  }

  p + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(color = NULL, title = title) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom")
}

#----------------------------------------------------------
# Plot the random-search AUC trajectories (mean best-PR-found-so-
# far vs. fraction of the grid tried) from one or more
# auc_cummax_PR() results, with each method's AUC value embedded
# directly in its legend label (e.g. "LFMM (AUC=0.231)").
#
# Takes a NAMED LIST of already-computed auc_cummax_PR() results,
# so nothing gets recomputed -- e.g.:
#   plot_auc_trajectories(list(
#     EMMAX = auc_cummax_PR(all_scored, p_name = "EMX",  n_perm = 1000, seed = 1),
#     LFMM  = auc_cummax_PR(all_scored, p_name = "LFMM", n_perm = 1000, seed = 1)
#   ))
#----------------------------------------------------------
plot_auc_trajectories <- function(auc_results, title = NULL) {

  if (is.null(names(auc_results)) || any(names(auc_results) == "")) {
    stop("`auc_results` must be a NAMED list, e.g. list(EMMAX = auc_cummax_PR(...), LFMM = auc_cummax_PR(...)).")
  }

  traj <- rbindlist(lapply(names(auc_results), function(nm) {
    tr <- copy(auc_results[[nm]]$trajectory)
    if (nrow(tr) == 0) return(NULL)
    tr[, method := nm]
    tr[]
  }), fill = TRUE)

  auc_lookup <- vapply(auc_results, function(x) x$AUC, numeric(1))
  method_labels <- sprintf("%s (AUC=%.3f)", names(auc_lookup), auc_lookup)
  label_lookup <- setNames(method_labels, names(auc_lookup))

  traj[, method_label := factor(label_lookup[method], levels = method_labels)]

  ggplot(traj, aes(x = frac_tried, y = mean_cummax, color = method_label)) +
    geom_line(linewidth = 1) +
    labs(x = "Fraction of grid tried (random search)", y = "Mean best PR found so far",
         color = NULL, title = title) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom")
}

#----------------------------------------------------------
# Extract the raw-alpha baseline curve from an already-computed
# grid `scored` table (from run_and_score_one_sim_file()/
# run_and_score_all()) -- no new pipeline needed, since alpha was
# already swept there. Filters to th_ldw = 0 (no LD-filtering, no
# C-score) and ONE fixed r2_th/l_min, matching whatever you used
# for the C-score comparison, so pr_curve_auc() on this and on
# run_C_score_threshold_grid()'s output are directly comparable --
# same clustering rule, only the calling criterion differs.
#----------------------------------------------------------
extract_alpha_baseline <- function(scored, p_name, sel_r2_th, sel_l_min, sel_th_ldw = 0) {
  sub <- scored[th_ldw == sel_th_ldw & r2_th == sel_r2_th & l_min == sel_l_min]
  data.table(
    alpha = sub$alpha,
    TP = sub[[paste0("TP_", p_name)]],
    FP = sub[[paste0("FP_", p_name)]],
    FN = sub[[paste0("FN_", p_name)]]
  )
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
# "full LD-filtering"), and optionally method.
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
# GENOME-WIDE display across all replicate files (e.g. the 10
# chr1..chr10 simulations for one V/env). This relabeling scheme
# (Chr1->Chr(2i-1) QTN, Chr2->Chr(2i) neutral) is SPECIFIC to the
# simulation design (10 replicate chromosome-pairs stacked into a
# display genome) and has no empirical-data equivalent -- real
# data already has its own real chromosomes, no relabeling needed.
#----------------------------------------------------------

#----------------------------------------------------------
# Expensive step: run the per-file OR pipeline ONCE for a given
# grid point, across all replicate files, and cache everything
# needed to plot it repeatedly -- WITHOUT rendering.
#----------------------------------------------------------
score_OR_genome <- function(file_paths, p_cols, p_name,
                            sel_rho, sel_th_ldw, sel_r2_th, sel_l_min, sel_alpha = 0.05,
                            bp_th_cluster = 5e5, rho_r2_focal = 0.25,
                            rho_d_focal = 0.99, cores = 1,
                            r2_grid_scale = c("rho", "r2")) {

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
      row <- scored[rho == sel_rho & th_ldw == sel_th_ldw &
                      r2_th == sel_r2_th & l_min == sel_l_min]
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
# Convenience wrapper: score + plot in one call.
#----------------------------------------------------------
plot_OR_manhattan_genome <- function(file_paths, p_cols, p_name, value_col,
                                     sel_rho, sel_th_ldw, sel_r2_th, sel_l_min, sel_alpha = 0.05,
                                     bp_th_cluster = 5e5, rho_r2_focal = 0.25,
                                     rho_d_focal = 0.99, cores = 1,
                                     title = NULL, value_label = value_col, nrow_facets = 1,
                                     p_col = NULL, color_by = "status", wes_palette_name = "Zissou1",
                                     r2_grid_scale = c("rho", "r2")) {

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

#----------------------------------------------------------
# Plot from an already-computed all_results/scored table (e.g.
# from run_and_score_one_sim_file()), for one specific grid point.
#----------------------------------------------------------
plot_OR_manhattan_from_results <- function(map, all_results, value_col, p_name,
                                           sel_rho, sel_th_ldw, sel_r2_th, sel_l_min,
                                           value_label = value_col, ld_ws = NULL, p_col = NULL,
                                           color_by = "status", wes_palette_name = "Zissou1",
                                           r2_grid_scale = c("rho", "r2")) {

  r2_grid_scale <- match.arg(r2_grid_scale)

  if (!p_name %in% names(all_results)) {
    known_cols <- setdiff(names(all_results), c("r2_th", "l_min", "th_ldw", "rho", "n_loci", "file"))
    known_cols <- known_cols[!grepl("^(TP_|FP_|FN_|Precision_|Recall_|PR_|assignments_)", known_cols)]
    stop("p_name = '", p_name, "' is not a column in all_results (its clusters would silently be ",
         "NULL -> nothing plotted as called). Available cluster columns: ",
         paste(known_cols, collapse = ", "))
  }

  row <- all_results[rho == sel_rho & th_ldw == sel_th_ldw &
                       r2_th == sel_r2_th & l_min == sel_l_min]
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
# Genome-wide C-score / AUC across all replicate files. Same
# "expensive fetch once, cheap recompute" split as score_OR_genome().
#----------------------------------------------------------

fetch_C_score_data_all_files <- function(file_paths, p_cols, th_ldw_grid, r2_grid, lmin_grid,
                                         alpha_grid = 0.05, bp_th_cluster = 5e5,
                                         rho_r2_focal = 0.75, rho_d_focal = 0.95,
                                         cores = 1, r2_grid_scale = c("rho", "r2")) {

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
# Cheap: build the genome-wide C-score display object directly
# from an already-cached results list (run_and_score_all(...,
# return_intermediates = TRUE)) -- no recomputation, no re-running
# the per-file pipeline.
#
# `results_list` should be ONE coherent group of files that share
# a genome-wide display scheme (e.g. the 10 chr1..chr10 replicates
# for one V/env) -- subset all_results to that group before
# calling this, e.g.:
#   sel <- names(all_results)[grepl("_V2_c1_env1\\.rds$", names(all_results))]
#   compute_C_score_genome_from_results(all_results[sel], ...)
#
# IMPORTANT: C-score is computed SEPARATELY per file (each file's
# own scored table, own marker set) -- markers/chromosomes are
# relabeled for DISPLAY only AFTER that, exactly as in
# fetch_C_score_data_all_files(). Do NOT compute C-score on a
# pre-pooled multi-file table instead (e.g. via
# compute_marker_C_score(combine_scored_list(all_results), ...)) --
# every file internally reuses the SAME marker labels (e.g.
# "Chr1:12345"), so pooling the raw grid across files BEFORE
# computing C-score silently conflates unrelated markers from
# different chromosomes/files into one (inflated, meaningless)
# C-score. This function avoids that by construction.
#----------------------------------------------------------
compute_C_score_genome_from_results <- function(results_list, p_name, fixed_r2_th = NULL,
                                                filter_fun = NULL, auc_n_perm = 200, auc_seed = NULL) {

  file_names <- names(results_list)
  if (is.null(file_names) || any(file_names == "")) {
    stop("`results_list` must be a NAMED list keyed by filename (as returned by ",
         "run_and_score_all(..., return_intermediates = TRUE)).")
  }
  if (any(vapply(results_list, is.null, logical(1)))) {
    bad <- file_names[vapply(results_list, is.null, logical(1))]
    stop("No cached result for: ", paste(bad, collapse = ", "), ". These files had no potential ",
         "outliers anywhere, so run_and_score_all() stored NULL for them -- there's no map to plot. ",
         "Drop them from results_list before calling this, or re-load their map directly via readRDS() ",
         "if you want them shown as an all-C_score=0 chromosome.")
  }

  chr_num <- as.integer(gsub("[^0-9]", "", regmatches(file_names, regexpr("chr[0-9]+", file_names))))
  if (any(is.na(chr_num))) stop("Could not parse chromosome number from one or more result names.")
  file_names <- file_names[order(chr_num)]

  maps     <- vector("list", length(file_names))
  auc_list <- vector("list", length(file_names))

  for (i in seq_along(file_names)) {
    nm  <- file_names[i]
    res <- results_list[[nm]]

    map_i    <- copy(res$map)
    scored_i <- if (is.null(filter_fun)) res$scored else filter_fun(res$scored)

    if (nrow(scored_i) == 0 || !p_name %in% names(scored_i)) {
      map_i[, C_score := 0]
      auc_val <- NA_real_
    } else {
      marker_C <- compute_marker_C_score(scored_i, p_name = p_name, fixed_r2_th = fixed_r2_th)
      map_i    <- attach_C_score_to_map(map_i, marker_C)
      auc_res  <- auc_cummax_PR(scored_i, p_name = p_name, fixed_r2_th = fixed_r2_th,
                                n_perm = auc_n_perm, seed = auc_seed)
      auc_val  <- auc_res$AUC
    }

    ## relabel chromosomes/markers for DISPLAY only, AFTER C-score is
    ## already computed on this file's own (unambiguous) marker set
    old_marker <- map_i$marker
    new_chr    <- ifelse(map_i$Chr == "Chr1", paste0("Chr", 2 * i - 1), paste0("Chr", 2 * i))
    old_prefix <- paste0(map_i$Chr, ":")
    new_prefix <- paste0(new_chr, ":")
    new_marker <- mapply(function(m, o, n) sub(o, n, m, fixed = TRUE),
                         old_marker, old_prefix, new_prefix, USE.NAMES = FALSE)
    map_i[, Chr := new_chr]
    map_i[, marker := new_marker]

    maps[[i]] <- map_i
    auc_list[[i]] <- data.table(Chr = unique(new_chr), file = nm, AUC = auc_val)
  }

  combined_map <- rbindlist(maps, fill = TRUE)
  chr_order <- unique(combined_map$Chr)
  chr_order <- chr_order[order(as.integer(gsub("[^0-9]", "", chr_order)))]
  combined_map[, Chr := factor(Chr, levels = chr_order)]
  setorder(combined_map, Chr, Pos)
  combined_map[, indx := .I]

  list(map = combined_map, chr_order = chr_order, auc_by_file = rbindlist(auc_list, fill = TRUE))
}

compute_C_score_genome <- function(file_paths, p_cols, p_name, th_ldw_grid, r2_grid, lmin_grid,
                                   fixed_r2_th = NULL, filter_fun = NULL, bp_th_cluster = 5e5,
                                   rho_r2_focal = 0.75, rho_d_focal = 0.95,
                                   alpha_grid = 0.05, cores = 1, r2_grid_scale = c("rho", "r2"),
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
# res <- run_and_score_one_sim_file(
#   "./parsed_sim_data/adapt_bgs_chr1_V2_c1_env1.rds",
#   p_cols = p_cols, th_ldw_grid = th_ldw_grid, r2_grid = r2_grid, lmin_grid = lmin_grid,
#   return_intermediates = TRUE
# )
#
# marker_C <- compute_marker_C_score(res$scored, p_name = "LFMM", fixed_r2_th = 0.6)
# qtn_C    <- compute_QTN_C_score(res$scored, p_name = "LFMM", fixed_r2_th = 0.6, map = res$map)
# auc_res  <- auc_cummax_PR(res$scored, p_name = "LFMM", n_perm = 200, seed = 1)
#
# C_sweep <- run_C_score_sweep_multi(
#   scored = res$scored, map = res$map, el = res$el, qtn_ld_table = res$qtn_ld_table,
#   p_names = c("EMX", "LFMM"), r2_min_focal = res$thr$r2_min_focal, d_max_focal = res$thr$d_max_focal,
#   r2_th_grid = r2_grid, C_score_grid = seq(0, 1, by = 0.1), l_min = 1
# )
# auc_by_Cscore <- auc_cummax_PR_Cscore_sweep(C_sweep, n_perm = 200, seed = 1)
# plot_Cscore_AUC_curve(auc_by_Cscore)
#
# ## Standard PR-curve comparison: C-score vs raw alpha, at the SAME fixed
# ## r2_th/l_min -- this isolates "does C-score beat raw significance" from
# ## any clustering-rule difference.
# sel_r2_th <- 0.7; sel_l_min <- 1
#
# C_curve <- run_C_score_threshold_grid(
#   marker_C_score = compute_marker_C_score(res$scored, p_name = "LFMM", fixed_r2_th = sel_r2_th),
#   map = res$map, el = res$el, qtn_ld_table = res$qtn_ld_table,
#   r2_min_focal = res$thr$r2_min_focal, d_max_focal = res$thr$d_max_focal,
#   C_score_grid = seq(0, 1, by = 0.05), r2_th = sel_r2_th, l_min = sel_l_min
# )
# alpha_curve <- extract_alpha_baseline(res$scored, p_name = "LFMM",
#                                       sel_r2_th = sel_r2_th, sel_l_min = sel_l_min)
#
# pr_curve_auc(C_curve)$AUC
# pr_curve_auc(alpha_curve)$AUC
#
# C_curve[, source := "C-score"]; alpha_curve[, source := "raw alpha"]
# plot_pr_curve(rbind(C_curve, alpha_curve, fill = TRUE), group_col = "source")

