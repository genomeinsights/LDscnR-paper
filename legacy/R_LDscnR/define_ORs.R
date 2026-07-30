######################################################
## GENERIC outlier-region (OR) pipeline
##
## Everything in this file works on empirical OR simulated data --
## nothing here requires knowing true causal loci. This is the
## actual empirical-data pipeline: FDR + th_ldw filtering ->
## LD-clustering into ORs -> consistency (C-score) across a grid
## of (rho, th_ldw, l_min, alpha) -> final OR calls from
## thresholding the C-score, at a fixed (r2_th, distance) pair
## chosen by convenience/robustness rather than by optimizing
## against ground truth (which empirical data doesn't have).
##
## Simulation-only code (anything needing known true QTNs --
## truth-scoring, TP/FP/FN, Precision/Recall/PR, AUC) lives in
## outlier_regions_simulation.R, which sources/depends on this file.
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
# One grid point (rho, th_ldw, alpha) -- single global FDR
# correction, clustered across the full (r2_th x l_min) grid.
# Returns raw OR clusters ONLY -- no ground truth involved.
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
# Generic per-dataset orchestrator: takes already-loaded map/GTs/
# ld_ws (NOT a file path -- empirical data will load these its own
# way, e.g. via GDS, and simulated data via readRDS; both funnel
# into this same function). Loops the (rho x th_ldw x alpha) grid
# via run_one_grid(), returning raw OR clusters -- no ground truth.
#
# decay_abc, if supplied, is a list(a=, b=, c=) from a fitted
# LD-decay curve, used ONLY to (i) convert r2_grid from the rho
# scale to actual r^2 thresholds if r2_grid_scale = "rho", and
# (ii) report which rho the fixed bp_th_cluster distance
# corresponds to (bp_th_cluster_rho), for reference. Leave NULL
# if you don't have/want a decay-curve-relative scale.
#----------------------------------------------------------
run_and_cluster <- function(map, GTs, ld_ws, p_cols, th_ldw_grid, r2_grid, lmin_grid,
                            bp_th_cluster = 5e5, alpha_grid = 0.05, cores = 1,
                            rho_grid = NULL, r2_grid_scale = c("r2", "rho"),
                            decay_abc = NULL) {

  r2_grid_scale <- match.arg(r2_grid_scale)

  potential_outliers <- get_potential_outliers(map, ld_ws, th_ldw_grid, p_cols, alpha_grid = alpha_grid)
  if (length(potential_outliers) == 0) {
    message("No potential outliers found.")
    return(NULL)
  }

  el <- precompute_LD_edges(
    GTs = GTs[, potential_outliers, drop = FALSE],
    map = map[marker %in% potential_outliers],
    r2_min = 0.1, max_bp = bp_th_cluster, cores = cores
  )

  bp_th_cluster_rho <- NA_real_
  r2_grid_rho <- NULL

  if (!is.null(decay_abc)) {
    ## rho equivalent to the fixed clustering distance, i.e. the inverse of
    ## d_from_rho(): rho = a*d / (1 + a*d)
    bp_th_cluster_rho <- (decay_abc$a * bp_th_cluster) / (1 + decay_abc$a * bp_th_cluster)

    if (identical(r2_grid_scale, "rho")) {
      r2_grid_rho <- r2_grid                       ## keep original rho values for traceability
      r2_grid <- ld_from_rho(b = decay_abc$b, c = decay_abc$c, rho = r2_grid)
      message("r2_grid (rho -> r2): ", paste(sprintf("%.3f->%.3f", r2_grid_rho, r2_grid), collapse = ", "))
    }
  } else if (identical(r2_grid_scale, "rho")) {
    stop("r2_grid_scale = 'rho' requires decay_abc = list(a=, b=, c=) to convert rho -> r2.")
  }

  rho_vals <- if (is.null(rho_grid)) colnames(ld_ws) else as.character(rho_grid)
  param_grid <- CJ(rho = rho_vals, th_ldw = th_ldw_grid, alpha = alpha_grid)

  outliers <- rbindlist(lapply(seq_len(nrow(param_grid)), function(i) {
    pars <- param_grid[i]
    run_one_grid(map = map, el = el, ld_ws = ld_ws,
                 rho = pars$rho, th_ldw = pars$th_ldw,
                 p_cols = p_cols, alpha = pars$alpha,
                 r2_grid = r2_grid, lmin_grid = lmin_grid,
                 bp_th = bp_th_cluster, cores = cores)
  }), fill = TRUE)

  if (!is.null(r2_grid_rho)) {
    rho_lookup <- setNames(r2_grid_rho, as.character(r2_grid))
    outliers[, r2_grid_rho := rho_lookup[as.character(r2_th)]]
  }
  outliers[, bp_th_cluster := bp_th_cluster]
  outliers[, bp_th_cluster_rho := bp_th_cluster_rho]

  list(outliers = outliers, el = el, potential_outliers = potential_outliers,
       bp_th_cluster = bp_th_cluster, bp_th_cluster_rho = bp_th_cluster_rho)
}

#----------------------------------------------------------
# Small shared helper: auto-detect whether to filter a fixed r2
# threshold on the linear scale (r2_th) or the rho scale
# (r2_grid_rho, present when run_and_cluster()/
# run_and_score_one_sim_file() was called with r2_grid_scale =
# "rho") -- prefers r2_grid_rho whenever it's present, since a
# fixed_r2_th value supplied by the caller almost always means
# "the rho I actually chose", not "the raw r^2 it happened to map
# to on this particular file's decay curve".
#----------------------------------------------------------
.r2_col <- function(dt) if ("r2_grid_rho" %in% names(dt)) "r2_grid_rho" else "r2_th"

#----------------------------------------------------------
# Per-MARKER consistency score (C-score). GENERIC -- pools rows
# of `outliers`/`scored` sharing the same r2_th (or all of them,
# if fixed_r2_th = NULL) into "fraction of grid points this
# marker was called in". This is what final empirical outlier
# calls are based on: a marker (or OR built from called markers)
# passes if C_score >= some chosen threshold tau_C.
#----------------------------------------------------------
compute_marker_C_score <- function(scored, p_name, fixed_r2_th = NULL, filter_fun = NULL) {

  sub <- if (is.null(filter_fun)) scored else filter_fun(scored)
  sub <- if (is.null(fixed_r2_th)) sub else sub[get(.r2_col(sub)) == fixed_r2_th]
  n_grid <- nrow(sub)
  if (n_grid == 0) {
    msg <- if (is.null(fixed_r2_th)) "No rows in `scored` (after filter_fun, if supplied)." else paste0("No rows in `scored` at r2_th = ", fixed_r2_th, " (after filter_fun, if supplied).")
    stop(msg)
  }

  tab <- table(unlist(sub[[p_name]]))

  ## `p_name` may have ZERO calls across the ENTIRE grid. table(NULL) then
  ## has length 0, and names(tab) is NULL rather than character(0) --
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
# Attach a compute_marker_C_score() result onto `map` as a new
# column, so it can be used directly as value_col/color_by in
# the Manhattan-plotting functions below. Markers absent from
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
# GENERIC engine: given a marker_C_score table (from
# compute_marker_C_score(), at one r2_th), re-cluster the set of
# markers passing each C-score threshold in C_score_grid, at that
# r2_th. Returns raw OR clusters ONLY -- no ground truth. This
# (and run_ORs_C_score_sweep() below) is the actual final-outlier-
# calling engine for empirical data: build C-score once, then call
# ORs by thresholding it at a single, convenience-chosen tau_C and
# r2_th, with no truth-based optimization needed or possible.
#----------------------------------------------------------
get_ORs_at_C_score_threshold <- function(marker_C_score, el, C_score_grid, r2_th, l_min = 1, bp_th = Inf) {

  rbindlist(lapply(C_score_grid, function(C_th) {

    called_markers <- marker_C_score[C_score >= C_th, marker]
    row <- data.table(C_score_threshold = C_th, r2_th = r2_th, l_min = l_min)

    if (length(called_markers) == 0) {
      row[, clusters := list(list(character(0)))]
      return(row[])
    }

    clusters_dt <- LD_igraph_components(el = el, markers = called_markers, r2_th = r2_th, bp_th = bp_th)
    clusters_dt <- clusters_dt[n_loci >= l_min]
    cls <- split(clusters_dt$marker, clusters_dt$CL_id)

    ## same length-1 list-column subtlety as run_one_grid() -- see comment there
    if (length(cls) == 0) {
      row[, clusters := list(list(character(0)))]
    } else if (length(cls) == 1) {
      row[, clusters := list(list(cls))]
    } else {
      row[, clusters := list(cls)]
    }
    row[]
  }))
}

#----------------------------------------------------------
# GENERIC: sweep BOTH r2_th_grid and C_score_grid, computing
# marker C-scores fresh at each r2_th (since OR membership --
# and therefore C-score -- depends on r2_th), then re-clustering
# at each C-score threshold within that r2_th. Returns raw OR
# clusters per (r2_th, C_score_threshold) combination -- no
# ground truth. This is the direct empirical-data entry point:
# the final call is "take get_ORs_at_C_score_threshold()'s
# clusters at your chosen (r2_th, tau_C)".
#----------------------------------------------------------
run_ORs_C_score_sweep <- function(scored, el, p_name, r2_th_grid, C_score_grid,
                                  l_min = 1, bp_th = Inf, filter_fun = NULL) {

  rbindlist(lapply(r2_th_grid, function(r2_th) {
    mc <- compute_marker_C_score(scored, p_name, fixed_r2_th = r2_th, filter_fun = filter_fun)
    get_ORs_at_C_score_threshold(mc, el, C_score_grid, r2_th, l_min = l_min, bp_th = bp_th)
  }))
}

#----------------------------------------------------------
# Annotate a map with call status for one grid point's clusters
# -- no plotting. GENERIC: if `assignments` is NULL (the empirical/
# no-ground-truth case), status collapses to 2 levels ("called (OR)"
# / "not called"). If `assignments` is supplied (simulation, from
# classify_ORs() in outlier_regions_simulation.R), status keeps the
# full 3-level TP/FP/not-called distinction.
#----------------------------------------------------------
annotate_OR_calls <- function(map, value_col, clusters, assignments = NULL, ld_ws = NULL, rho = NULL,
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

  if (!is.null(assignments)) {
    ## SIMULATION-STYLE coloring (assignments come from classify_ORs() in
    ## outlier_regions_simulation.R) -- full TP/FP/not-called distinction
    if (nrow(cl_dt) > 0 && nrow(assignments) > 0) {
      cl_dt <- merge(cl_dt, assignments[, .(CL_id, is_TP)], by = "CL_id", all.x = TRUE)
    } else {
      cl_dt[, is_TP := logical(0)]
    }
    map_sub <- merge(map_sub, cl_dt, by = "marker", all.x = TRUE)
    map_sub[, status := fifelse(is.na(CL_id), "not called",
                                fifelse(is_TP, "true positive OR", "false positive OR"))]
    map_sub[, status := factor(status, levels = c("not called", "false positive OR", "true positive OR"))]
  } else {
    ## GENERIC coloring (no ground truth) -- called vs not called only
    map_sub <- merge(map_sub, cl_dt, by = "marker", all.x = TRUE)
    map_sub[, status := fifelse(is.na(CL_id), "not called", "called (OR)")]
    map_sub[, status := factor(status, levels = c("not called", "called (OR)"))]
  }
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
    status_levels <- levels(map_sub$status)
    cols <- if (identical(status_levels, c("not called", "called (OR)"))) {
      c("not called" = "grey70", "called (OR)" = "firebrick3")
    } else {
      c("not called" = "grey70", "false positive OR" = "steelblue3", "true positive OR" = "firebrick3")
    }
    p <- p + aes(color = status) + scale_color_manual(values = cols)
  } else {
    ## continuous column (e.g. C_score, max_LD_with_QTN) via a wesanderson continuous palette
    if (!color_by %in% names(map_sub)) {
      stop("color_by = '", color_by, "' not found in the data. ",
           "Use 'status' for the called/not-called (or TP/FP) coloring, or a numeric column name.")
    }
    p <- p + aes(color = .data[[color_by]]) +
      scale_color_gradientn(colors = wes_palette(wes_palette_name, 100, type = "continuous"),
                            name = color_by, na.value = "grey85")
  }

  p
}

#----------------------------------------------------------
# Manhattan-style plot of one grid point's outlier calls.
# `assignments = NULL` (default) gives the GENERIC 2-level
# called/not-called view -- pass assignments from classify_ORs()
# (simulation only) for the 3-level TP/FP/not-called view.
# Uses raw (unadjusted) -log10(p) as a stable visual background
# so the same plot is comparable across different filtering
# settings -- only the coloring of "called" points changes,
# not the underlying point cloud.
#----------------------------------------------------------
plot_OR_manhattan <- function(map, value_col, clusters, assignments = NULL, title = NULL,
                              value_label = value_col, ld_ws = NULL, rho = NULL,
                              th_ldw = NULL, p_col = NULL,
                              color_by = "status", wes_palette_name = "Zissou1") {
  map_sub <- annotate_OR_calls(map, value_col, clusters, assignments = assignments, ld_ws = ld_ws, rho = rho,
                               th_ldw = th_ldw, p_col = p_col)
  if (!"indx" %in% names(map_sub)) map_sub[, indx := .I]
  .render_OR_manhattan(map_sub, title = title, nrow_facets = 1, value_label = value_label,
                       color_by = color_by, wes_palette_name = wes_palette_name)
}

#----------------------------------------------------------
# Example usage (empirical or simulated -- no ground truth needed)
#----------------------------------------------------------
# ## p_cols: named vector mapping method label -> p-value column name on `map`
# p_cols <- c(EMX = "emx_p", LFMM = "lfmm_p")
#
# ## `map`, `GTs`, `ld_ws` loaded however is appropriate for your data
# ## (readRDS for simulated files; your own GDS-based loader for empirical)
# res <- run_and_cluster(
#   map = map, GTs = GTs, ld_ws = ld_ws, p_cols = p_cols,
#   th_ldw_grid = c(0, 0.5, 0.75, 0.9, 0.95),
#   r2_grid = seq(0.6, 0.9, by = 0.1), lmin_grid = c(1, 5, 10, 20),
#   alpha_grid = 0.05
# )
#
# ## per-marker C-score, pooled across the whole (rho, th_ldw, l_min, alpha)
# ## grid at one r2_th -- this is the actual final calling criterion
# marker_C <- compute_marker_C_score(res$outliers, p_name = "LFMM", fixed_r2_th = 0.7)
# map_C <- attach_C_score_to_map(map, marker_C)
#
# ## final OR calls: threshold the C-score at your chosen tau_C, at your
# ## chosen (fixed) r2_th -- the actual empirical-data outlier-calling step
# ORs <- get_ORs_at_C_score_threshold(
#   marker_C_score = marker_C, el = res$el, C_score_grid = 0.2, r2_th = 0.7, l_min = 1
# )
# final_clusters <- ORs$clusters[[1]]   # list of character vectors, one per OR
#
# ## Manhattan plot of the C-score itself, generic (no truth) coloring:
# plot_OR_manhattan(map_C, value_col = "C_score", clusters = final_clusters,
#                   value_label = "C-score", color_by = "C_score")
