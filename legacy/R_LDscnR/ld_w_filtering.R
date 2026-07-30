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
get_potential_outliers_sim <- function(map, ld_ws, th_ldw_grid, p_cols, alpha = 0.05) {
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
        potential <- c(potential, map_sub$marker[keep][q < alpha])
      }
    }
  }
  unique(potential)
}

#----------------------------------------------------------
# Empty result helper
#----------------------------------------------------------
empty_result <- function(rho, th_ldw, p_names, r2_grid, lmin_grid) {
  out <- CJ(r2_th = r2_grid, l_min = lmin_grid)
  out[, `:=`(th_ldw = th_ldw, rho = rho)]
  for (nm in p_names) out[, (nm) := list(list(character()))]
  out[]
}

#----------------------------------------------------------
# One grid point (rho, th_ldw) -- single global FDR correction
#----------------------------------------------------------
run_one_grid_sim <- function(map, el = NULL, ld_ws, rho, th_ldw, p_cols,
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

  if (!any(keep)) return(cbind(empty_result(rho, th_ldw, p_names, r2_grid, lmin_grid), n_loci = 0))

  markers_keep <- map_sub[keep, marker]
  outliers <- setNames(vector("list", length(p_cols)), p_names)

  for (i in seq_along(p_cols)) {
    p_col <- p_cols[i]
    nm <- p_names[i]
    q <- p.adjust(unlist(map_sub[keep, ..p_col]), method = "fdr")
    outliers[[nm]] <- markers_keep[q < alpha]
  }

  if (length(unique(unlist(outliers))) == 0) {
    return(cbind(empty_result(rho, th_ldw, p_names, r2_grid, lmin_grid), n_loci = length(which(keep))))
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
      row <- data.table(r2_th = r2_th, l_min = l_min, th_ldw = th_ldw, rho = rho)
      for (nm in p_names) {
        tmp <- clusters[[nm]][n_loci >= l_min, ]
        cls <- split(tmp$marker, tmp$CL_id)
        if (length(cls) >= 1) row[, (nm) := list(cls)] else row[, (nm) := list(list(character(0)))]
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
# Evaluate one set of ORs: assign, resolve duplicates,
# compute TP / FP / FN / Precision / Recall / PR
#----------------------------------------------------------
evaluate_ORs <- function(clusters, map, qtn_ld_table, r2_min_focal, d_max_focal) {
  qtn_marker_set   <- map[type == "QTN", marker]
  true_pos_markers <- map[true_pos_QTN == TRUE, marker]

  if (length(clusters) == 0) {
    return(list(TP = 0L, FP = 0L, FN = length(true_pos_markers),
                Precision = NA_real_, Recall = 0, PR = 0))
  }

  assign_list <- lapply(clusters, assign_focal_qtn_one_OR,
                        qtn_ld_table = qtn_ld_table, qtn_marker_set = qtn_marker_set,
                        r2_min_focal = r2_min_focal, d_max_focal = d_max_focal)

  assignments <- data.table(
    CL_id    = seq_along(clusters),
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

  TP_markers <- assignments[!is.na(qtn) & qtn %in% true_pos_markers, qtn]
  TP <- length(TP_markers)
  FP <- sum(is.na(assignments$qtn)) + sum(!is.na(assignments$qtn) & !(assignments$qtn %in% true_pos_markers))
  FN <- length(setdiff(true_pos_markers, TP_markers))

  Precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  Recall    <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  PR        <- if (!is.na(Precision) && !is.na(Recall)) Precision * Recall else NA_real_

  list(TP = TP, FP = FP, FN = FN, Precision = Precision, Recall = Recall, PR = PR)
}

#----------------------------------------------------------
# Score a full outliers_dt grid (rbind'ed over param_grid)
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
  }
  dt[]
}

#----------------------------------------------------------
# Full pipeline for ONE simulation file
#----------------------------------------------------------
run_and_score_one_sim_file <- function(file_path, p_cols, th_ldw_grid, r2_grid, lmin_grid,
                                       bp_th_cluster = 5e5, rho_r2_focal = 0.75,
                                       rho_d_focal = 0.95, alpha = 0.05, cores = 1) {

  d <- readRDS(file_path)
  map <- flag_true_positive_QTNs(d$map)
  GTs <- d$GTs
  ld_ws <- d$ld_ws
  LD_decay <- d$LD_decay

  potential_outliers <- get_potential_outliers_sim(map, ld_ws, th_ldw_grid, p_cols, alpha = alpha)
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

  param_grid <- CJ(rho = colnames(ld_ws), th_ldw = th_ldw_grid)

  outliers <- rbindlist(lapply(seq_len(nrow(param_grid)), function(i) {
    pars <- param_grid[i]
    run_one_grid_sim(map = map, el = el_potential, ld_ws = ld_ws,
                     rho = pars$rho, th_ldw = pars$th_ldw,
                     p_cols = p_cols, alpha = alpha,
                     r2_grid = r2_grid, lmin_grid = lmin_grid,
                     bp_th = bp_th_cluster, cores = cores)
  }), fill = TRUE)

  scored <- score_outlier_grid(outliers, map, qtn_ld_table, names(p_cols),
                               r2_min_focal = thr$r2_min_focal, d_max_focal = thr$d_max_focal)
  scored[, file := basename(file_path)]
  scored[]
}

#----------------------------------------------------------
# Loop over all simulation files and aggregate
#----------------------------------------------------------
run_and_score_all <- function(parsed_folder, p_cols, th_ldw_grid, r2_grid, lmin_grid,
                              bp_th_cluster = 5e5, rho_r2_focal = 0.75, rho_d_focal = 0.95,
                              alpha = 0.05, cores = 1) {

  files <- list.files(parsed_folder, pattern = "\\.rds$", full.names = TRUE)

  all_scored <- rbindlist(lapply(files, function(f) {
    message("Processing ", basename(f))
    tryCatch(
      run_and_score_one_sim_file(f, p_cols, th_ldw_grid, r2_grid, lmin_grid,
                                 bp_th_cluster, rho_r2_focal, rho_d_focal, alpha, cores),
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
# Example usage
#----------------------------------------------------------
p_cols <- c(EMX = "emx_p")   # add LFMM = "lfmm_p" for both
th_ldw_grid <- c(0, 0.5, 0.75, 0.9, 0.95)
r2_grid     <- seq(0.6, 0.9, by = 0.1)
lmin_grid   <- c(5, 10, 20)

all_results <- run_and_score_all(
  parsed_folder = "./parsed_sim_data",
  p_cols = p_cols, th_ldw_grid = th_ldw_grid,
  r2_grid = r2_grid, lmin_grid = lmin_grid,
  bp_th_cluster = 5e5, rho_r2_focal = 0.75, rho_d_focal = 0.95
)

## average Precision/Recall/PR across the 10 chr replicates, per (V, rho, th_ldw, r2_th, l_min)
summary_by_V <- all_results[r2_th == 0.7 & l_min == 10,
  .(mean_Precision = mean(Precision_EMX, na.rm = TRUE),
    mean_Recall    = mean(Recall_EMX, na.rm = TRUE),
    mean_PR        = mean(PR_EMX, na.rm = TRUE)),
  by = .(V, rho = factor(rho, levels = unique(rho)), th_ldw)]
