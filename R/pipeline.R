library(data.table)
library(igraph)
library(parallel)

get_candidate_markers <- function(map, ld_ws, p_col, qt_grid, rho_grid, alpha_grid) {
  alpha_min <- min(alpha_grid)
  #rho = "0.5"
  cand <- rbindlist(lapply(rho_grid, function(rho) {

    ld_vec <- ld_ws[,rho]

    rbindlist(lapply(qt_grid, function(qt) {

      keep <- ld_vec > quantile(ld_vec, qt, na.rm = TRUE)

      if (!any(keep, na.rm = TRUE)) {
        return(data.table(marker = character()))
      }

      tmp <- copy(map[keep])
      tmp[, p_tmp := get(p_col)]
      tmp[, q_tmp := p.adjust(p_tmp, method = "fdr")]
      tmp[, logq_tmp := -log10(q_tmp)]

      tmp[logq_tmp > alpha_min, .(marker)]

    }))

  }))

  unique(cand$marker)
}

precompute_LD_edges <- function(
    GTs,
    map,
    r2_min = 0.1,
    max_bp = Inf,
    cores = 1
) {

  library(data.table)

  markers <- colnames(GTs)

  if (length(markers) == 0) {
    stop("No markers found in GTs.")
  }

  map_sub <- copy(map[marker %in% markers])

  setkey(map_sub, Chr, Pos)

  chr_levels <- unique(map_sub$Chr)

  out <- parallel::mclapply(chr_levels, function(ch) {

    chr_map <- map_sub[Chr == ch]

    chr_markers <- chr_map$marker

    if (length(chr_markers) < 2) {
      return(NULL)
    }

    gts <- as.matrix(GTs[, chr_markers, drop = FALSE])

    storage.mode(gts) <- "double"

    # pairwise LD
    R2 <- cor(gts, use = "pairwise.complete.obs")^2

    R2[is.na(R2)] <- 0
    diag(R2) <- 0

    idx <- which(R2 >= r2_min, arr.ind = TRUE)

    idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]

    if (nrow(idx) == 0) {
      return(NULL)
    }

    pos <- chr_map$Pos

    dt <- data.table(
      marker1 = chr_markers[idx[, 1]],
      marker2 = chr_markers[idx[, 2]],
      r2 = R2[idx],
      pos1 = pos[idx[, 1]],
      pos2 = pos[idx[, 2]]
    )

    dt[, dist_bp := abs(pos1 - pos2)]

    if (is.finite(max_bp)) {
      dt <- dt[dist_bp <= max_bp]
    }

    dt[, Chr := ch]

    dt[, c("pos1", "pos2") := NULL]

    setcolorder(dt, c(
      "Chr",
      "marker1",
      "marker2",
      "r2",
      "dist_bp"
    ))

    dt

  }, mc.cores = cores)

  out <- rbindlist(out, use.names = TRUE, fill = TRUE)

  setkey(out, Chr, marker1, marker2)

  out
}

LD_igraph_components <- function(
    GTs,
    map,
    el,
    markers,
    r2_th = 0.8,
    bp_th = Inf
) {
  markers <- intersect(markers, colnames(GTs))

  if (length(markers) == 0) {
    return(data.table(marker = character(),
                      CL_id = integer(),
                      n_loci = integer()))
  }

  if (length(markers) == 1) {
    return(data.table(marker = markers,
                      CL_id = 1L,
                      n_loci = 1L))
  }

  # subset precomputed edge list
  edges <- el[
    marker1 %in% markers &
      marker2 %in% markers &
      r2 >= r2_th
  ]

  if (is.finite(bp_th)) {
    edges <- edges[dist_bp <= bp_th]
  }

  # if no edges, every SNP is singleton
  if (nrow(edges) == 0) {
    return(data.table(
      marker = markers,
      CL_id = seq_along(markers),
      n_loci = 1L
    ))
  }

  g <- graph_from_data_frame(
    edges[, .(from = marker1, to = marker2)],
    directed = FALSE,
    vertices = data.table(name = markers)
  )

  comp <- components(g)

  clusters <- data.table(
    marker = names(comp$membership),
    CL_id = as.integer(comp$membership)
  )

  clusters[, n_loci := comp$csize[CL_id]]

  clusters
}

make_param_grid <- function(
    qt_grid,
    alpha_grid,
    rho_grid,
    rho_ld_grid,
    l_min_grid
) {
  CJ(
    qt = qt_grid,
    alpha = alpha_grid,
    rho = rho_grid,
    rho_ld = rho_ld_grid,
    l_min = l_min_grid
  )
}

run_OR_scan_one <- function(
    map,
    GTs,
    el,
    ld_ws,
    p_col,
    n_markers,
    qt,
    alpha,
    rho,
    rho_ld,
    l_min_grid,
    max_gap,
    decay_sum = NULL,
    ld_from_rho = NULL,
    default_bp_th = Inf,
    use_decay_threshold = FALSE
) {

  empty_out <- function() {
    data.table(
      CL_uid = character(),
      OR_size = integer(),
      phy_uid = character(),
      qt = numeric(),
      alpha = numeric(),
      rho = character(),
      rho_ld = numeric(),
      l_min = numeric()
    )
  }

  # p_col <- "lfmm_p"
  if (is.na(qt) || rho == "none") {
    dt <- copy(map)
    dt[, p_method := get(p_col)]
    dt[, q_method := p.adjust(p_method, method = "fdr")]
    dt[, logq_method := -log10(q_method)]


  } else {
    dt <- map[ld_ws[, rho]>quantile(ld_ws[, rho], qt, na.rm = TRUE)]
    if (nrow(dt) == 0) return(empty_out())

    dt[, p_method := get(p_col)]
    dt[, q_method := p.adjust(p_method, method = "fdr")]
    dt[, logq_method := -log10(q_method)]

  }

  peaks <- dt[logq_method > alpha & marker %in% colnames(GTs), .(Chr, marker, Pos)]

  #peaks <- dt[marker %in% colnames(GTs)]

  if (nrow(peaks) == 0) return(empty_out())

  setorder(peaks, Chr, Pos)

  peaks[, phy_cluster :=
          cumsum(is.na(Pos - shift(Pos)) |
                   (Pos - shift(Pos)) > max_gap),
        by = Chr]

  peaks[, phy_uid := paste0(Chr, "_phy", phy_cluster)]

  phy_cls <- split(peaks$marker, peaks$phy_uid)
  # i <- 1
  sl_cls <- rbindlist(lapply(seq_along(phy_cls), function(i) {
    markers_i <- phy_cls[[i]]
    phy_uid_i <- names(phy_cls)[i]

    chr_i <- map[marker %in% markers_i, unique(Chr)][1]

    bp_th_i <- default_bp_th

    if (use_decay_threshold) {
      if (is.null(decay_sum) || is.null(ld_from_rho)) {
        stop("If use_decay_threshold=TRUE, provide ld_decay and ld_from_rho.")
      }

      c_chr <- decay_sum[Chr == chr_i, c_pred][1]
      b_chr <- decay_sum[Chr == chr_i, b][1]
      ld_th_i <- ld_from_rho(b_chr,c_chr, rho = rho_ld)
    }

    out <- LD_igraph_components(
      GTs = GTs,
      map = map,
      el = el,
      markers = markers_i,
      r2_th = ld_th_i,
      bp_th = bp_th_i
    )

    out[, phy_uid := phy_uid_i]
    out[, CL_uid := paste0(phy_uid, "_ld", CL_id)]

    out
  }), fill = TRUE)

  sl_cls <- map[,.(marker,Chr,Pos)][sl_cls,on = "marker"]

  if (nrow(sl_cls) == 0) return(empty_out())

  out <- rbindlist(lapply(l_min_grid, function(l_min) {

    tmp <- sl_cls[n_loci > l_min]

    if (nrow(tmp) == 0) return(NULL)

    ors <- tmp[, {
      if (.N == 0 || all(is.na(Pos))) {
        NULL
      } else {
        .(
          Chr = Chr[1],
          start = min(Pos, na.rm = TRUE),
          end = max(Pos, na.rm = TRUE),
          width = max(Pos, na.rm = TRUE) -
            min(Pos, na.rm = TRUE),
          OR_size = .N,
          OR = list(marker)
        )
      }
    }, by = .(phy_uid, CL_uid)]

    if (nrow(ors) == 0) return(NULL)

    ors[, `:=`(
      qt = qt,
      alpha = alpha,
      rho = rho,
      rho_ld = rho_ld,
      l_min = l_min
    )]

    ors

  }), fill = TRUE)


}

compute_C_score <- function(
    map,
    GTs,
    ld_ws,
    p_col,
    param_grid,
    l_min_grid,
    max_gap = 1e6,
    decay_sum = NULL,
    ld_from_rho = NULL,
    default_bp_th = Inf,
    use_decay_threshold = FALSE,
    cores = 1
) {


  cand <- get_candidate_markers(
    map = map,
    ld_ws = ld_ws,
    p_col = p_col,
    qt_grid = qt_grid,
    rho_grid = rho_grid,
    alpha_grid = alpha_grid
  )

  cand_unfiltered <- {
    tmp <- copy(map)
    tmp[, p_tmp := get(p_col)]
    tmp[, q_tmp := p.adjust(p_tmp, method = "fdr")]
    tmp[, logq_tmp := -log10(q_tmp)]
    tmp[logq_tmp > min(alpha_grid), marker]
  }

  cand <- unique(c(cand, cand_unfiltered))

  el <- precompute_LD_edges(
    GTs = GTs[, cand, drop = FALSE],
    map = map[marker %in% cand],
    r2_min = min(rho_ld_grid),
    max_bp = max_gap
  )

  #i <- 1
  #pars$qt = NA
  #pars$alpha = 1.31
  #pars$rho_ld = 0.9
  #pars$rho = "none"
  # tmp <- out[l_min==0]
  # p_col = "emx_p"
  # tmp[, is_true_positive := sapply(OR, function(mrk) {
  #   any(mrk %in% true_positives)
  # })]
  # tmp[,table(is_true_positive)]
  #map[emx_q<0.05]
  #use_decay_threshold = TRUE
  # i <- 200
  res <- rbindlist(
    mclapply(seq_len(nrow(param_grid)), function(i) {
      pars <- param_grid[i]

      out <- run_OR_scan_one(
        map = map,
        GTs = GTs_sub,
        ld_ws = ld_ws,
        el = el,
        n_markers = n_markers,
        p_col = p_col,
        qt = pars$qt,
        alpha = pars$alpha,
        rho = pars$rho,
        rho_ld = pars$rho_ld,
        l_min_grid = l_min_grid,
        max_gap = max_gap,
        decay_sum = decay_sum,
        ld_from_rho = ld_from_rho,
        default_bp_th = default_bp_th,
        use_decay_threshold = TRUE
      )
    }, mc.cores = cores),
    fill = TRUE
  )

}

ld_ws <- precalculate_ld_w(c(seq(0.1,0.95,by=0.05),0.99),ld_decay)

rho_ld_grid <- c(seq(0.1, 0.95, by = 0.1), 0.99)

qt_grid <- c(seq(0.1, 0.95, by = 0.05), 0.99)

alpha_grid <- -log10(c(
  0.05,
  0.01,
  0.005,
  0.001,
  0.0005,
  0.0001
))

l_min_grid <- seq(0, 10, by = 2)

rho_grid <- colnames(ld_ws)

param_grid <- rbind(
  ## filtered
  CJ(
    qt = qt_grid,
    alpha = alpha_grid,
    rho = rho_grid,
    rho_ld = rho_ld_grid
  ),
  ## non-filtered
  CJ(
    qt = NA,
    alpha = alpha_grid,
    rho = "none",
    rho_ld = rho_ld_grid
  ))

res <- compute_C_score( map,
                        GTs,
                        ld_ws,
                        p_col = "lfmm_p",
                        param_grid,
                        l_min_grid,
                        max_gap = 1e6,
                        decay_sum = ld_decay$decay_sum,
                        ld_from_rho = ld_from_rho,
                        default_bp_th = Inf,
                        use_decay_threshold = TRUE,
                        cores = 8)



