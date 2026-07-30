## ======================================================================
## Reusable functions from ld_w_filtering_3sp.R
## (LD-edge precomputation, LD-clustering into outlier regions,
##  sequential ld_w-filtering + FDR grid search)
## Logic unchanged from the original 3sp analysis script - only pulled
## out from the exploratory/plotting code so it can be sourced cleanly.
## ======================================================================

library(data.table)
library(igraph)
library(parallel)

# ----------------------------
# LD edge precomputation
# ----------------------------
precompute_LD_edges <- function(
    GTs,
    map,
    r2_min = 0.1,
    max_bp = Inf,
    cores = 1
) {
  markers <- intersect(colnames(GTs), map$marker)
  if (length(markers) == 0) stop("No overlapping markers between GTs and map.")

  map_sub <- copy(map[marker %in% markers])
  setkey(map_sub, Chr, Pos)

  chr_levels <- unique(map_sub$Chr)

  out <- mclapply(chr_levels, function(ch) {
    chr_map <- map_sub[Chr == ch]
    chr_markers <- chr_map$marker

    if (length(chr_markers) < 2) {
      return(data.table(
        Chr = ch,
        marker1 = chr_markers,
        marker2 = chr_markers,
        r2 = 1,
        dist_bp = 0
      ))
    }

    gts <- as.matrix(GTs[, chr_markers, drop = FALSE])
    storage.mode(gts) <- "double"

    R2 <- cor(gts, use = "pairwise.complete.obs")^2
    R2[is.na(R2)] <- 0
    diag(R2) <- 0

    idx <- which(R2 >= r2_min, arr.ind = TRUE)
    idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]

    if (nrow(idx) == 0) {
      return(data.table(
        Chr = ch,
        marker1 = chr_markers,
        marker2 = chr_markers,
        r2 = 1,
        dist_bp = 0
      ))
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


# ----------------------------
# LD clusters from edge list
# ----------------------------
LD_igraph_components <- function(
    el,
    markers,
    r2_th = 0.8,
    bp_th = Inf
) {
  markers <- unique(markers)

  if (length(markers) == 0) {
    return(data.table(marker = character(), CL_id = integer(), n_loci = integer()))
  }

  if (length(markers) == 1) {
    return(data.table(marker = markers, CL_id = 1L, n_loci = 1L))
  }

  edges <- el[
    marker1 %in% markers &
      marker2 %in% markers &
      r2 >= r2_th
  ]

  if (is.finite(bp_th)) edges <- edges[dist_bp <= bp_th]

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


# ----------------------------
# Empty result, generalized
# ----------------------------
empty_result <- function(
    rho,
    th_ldw,
    p_names,
    r2_grid,
    lmin_grid
) {
  out <- CJ(r2_th = r2_grid, l_min = lmin_grid)

  out[, `:=`(
    th_ldw = th_ldw,
    rho = rho
  )]

  for (nm in p_names) {
    out[, (nm) := list(list(character()))]
  }

  out[]
}


# ----------------------------
# One grid point, generalized
# ----------------------------
run_one_grid <- function(
    map,
    el=NULL,
    ld_ws,
    rho,
    th_ldw,
    p_cols,
    p_names = names(p_cols),
    alpha = 0.05,
    r2_grid,
    lmin_grid,
    bp_th = Inf,
    cores
) {
  stopifnot(length(p_cols) == length(p_names))

  # Important: align ld_ws and map by marker before filtering
  common_markers <- intersect(map$marker, rownames(ld_ws))


  map_sub <- copy(map[marker %in% common_markers])
  ld_sub <- ld_ws[map_sub$marker, , drop = FALSE]

  if(th_ldw>0){
    keep <- ld_sub[, rho] > quantile(ld_sub[, rho], th_ldw, na.rm = TRUE)
    keep[is.na(keep)] <- FALSE
  }else{
    keep <- rep(TRUE,nrow(ld_sub))
  }
  #table(keep)

  if (!any(keep)) {
    return(cbind(empty_result(rho, th_ldw, p_names, r2_grid, lmin_grid),n_loci=0))
  }

  markers_keep <- map_sub[keep, marker]

  outliers <- setNames(vector("list", length(p_cols)), p_names)

  #i <- 1
  for (i in seq_along(p_cols)) {
    p_col <- p_cols[i]
    nm <- p_names[i]

    q <- p.adjust(unlist(map_sub[keep, ..p_col]), method = "fdr")
    outliers[[nm]] <- markers_keep[q < alpha]
  }

  if (length(unique(unlist(outliers))) == 0) {
    return(cbind(empty_result(rho, th_ldw, p_names, r2_grid, lmin_grid),n_loci=length(which(keep))))
  }


  if(is.null(el)){
    all_outliers <- unique(unlist(outliers))
    el <- precompute_LD_edges(
      GTs = GTs[, all_outliers, drop = FALSE],
      map = map_sub[marker %in% all_outliers],
      r2_min = 0.1,
      max_bp = 1e6,
      cores = 1
    )
  }
  #r2_th = 0.8
  out <- rbindlist(mclapply(r2_grid, function(r2_th) {
    clusters <- lapply(outliers, function(markers) {
      LD_igraph_components(
        el = el,
        markers = markers,
        r2_th = r2_th,
        bp_th = bp_th
      )
    })

    out <- rbindlist(lapply(lmin_grid, function(l_min) {
      row <- data.table(
        r2_th = r2_th,
        l_min = l_min,
        th_ldw = th_ldw,
        rho = rho
      )
      #l_min = 1
      #nm = p_names[1]
      for (nm in p_names) {
        tmp <- clusters[[nm]][n_loci >= l_min, ]
        cls <- split(tmp$marker,tmp$CL_id)

        if(length(cls)>1){
          row[, (nm) := list(cls)]
        }else{
          row[, (nm) := list(list(cls))]
        }
      }
      row

    }), fill = TRUE)
  },mc.cores=cores), fill = TRUE)

  out[,n_loci:=length(which(keep))]

}

# ----------------------------
# Potential outliers, generalized
# ----------------------------
get_potential_outliers <- function(
    map,
    ld_ws,
    th_ldw_grid,
    p_cols,
    alpha = 0.05
) {
  common_markers <- intersect(map$marker, rownames(ld_ws))

  potential <- character()

  for (rho in colnames(ld_ws)) {
    message("processing rho = ",rho)
    ld_vec <- ld_ws[, rho]

    potential <- unique(potential)
    for (th_ldw in th_ldw_grid) {
      keep <- ld_vec > quantile(ld_vec, th_ldw, na.rm = TRUE)

      keep[is.na(keep)] <- FALSE

      if (!any(keep)) next
      ##p_col <- p_cols[1]
      for (p_col in p_cols) {
        q <- p.adjust(unlist(map[keep, ..p_col]), method = "fdr")
        potential <- c(potential, map[keep, marker][q < alpha])
      }
    }
  }

  return(unique(potential))

}

# ----------------------------
# Summarize stability (proportion of param combos a marker is called in)
# ----------------------------
summarise_stability <- function(outliers, map, p_names) {
  map_C <- copy(map)
  for (nm in p_names) {
    C <- outliers[, table(unlist(get(nm))) / .N]
    if(length(C)>0){
      C <- data.table(C)
      setnames(C, c("V1", "N"), c("marker", paste0("C_", nm)))
      map_C <- C[map_C, on = "marker"]
      map_C[is.na(get(paste0("C_", nm))), (paste0("C_", nm)) := 0]
    }
  }

  map_C[]
}
