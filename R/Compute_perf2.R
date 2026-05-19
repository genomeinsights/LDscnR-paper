library(data.table)
library(igraph)
library(parallel)
ld_from_rho <- function(b, c = 1, rho){
  b + (c - b) * (1 - rho)
}
get_candidate_markers <- function(map, ld_ws, p_col, qt_grid, rho_grid, alpha_grid) {
  alpha_max <- max(alpha_grid)
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

      tmp[q_tmp < alpha_max, .(marker)]

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
      rho_ld = numeric()

    )
  }

  # p_col <- "lfmm_p"
  if (is.na(qt) || rho == "none") {
    dt <- copy(map)
    dt[, p_method := get(p_col)]
    dt[, q_method := p.adjust(p_method, method = "fdr")]

    #rho <- "0.5"
    #qt <- 0.8
  } else {
    dt <- map[ld_ws[, rho]>quantile(ld_ws[, rho], qt, na.rm = TRUE)]
    if (nrow(dt) == 0) return(empty_out())

    dt[, p_method := get(p_col)]
    dt[, q_method := p.adjust(p_method, method = "fdr")]

  }

  peaks <- dt[q_method < alpha & marker %in% colnames(GTs), .(Chr, marker, Pos)]


  if (nrow(peaks) == 0) return(empty_out())

  setorder(peaks, Chr, Pos)

  peaks[, phy_cluster :=
          cumsum(is.na(Pos - shift(Pos)) |
                   (Pos - shift(Pos)) > max_gap),
        by = Chr]

  peaks[, phy_uid := paste0(Chr, "_phy", phy_cluster)]

  phy_cls <- split(peaks$marker, peaks$phy_uid)
  # i <- 1
  # rho_ld = 0.5
  sl_cls <- rbindlist(lapply(seq_along(phy_cls), function(i) {
    markers_i <- phy_cls[[i]]
    phy_uid_i <- names(phy_cls)[i]

    chr_i <- map[marker %in% markers_i, unique(Chr)][1]

    bp_th_i <- default_bp_th

    if (use_decay_threshold) {
      if (is.null(decay_sum)) {
        stop("If use_decay_threshold=TRUE, provide ld_decay ")
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

  ors <- sl_cls[, {
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
    rho_ld = rho_ld
  )]

  ors


}

OR_scan <- function(
    map,
    GTs,
    ld_ws,
    p_col,
    param_grid,
    alpha_grid,
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
    tmp[q_tmp < max(alpha_grid), marker]
  }

  cand <- unique(c(cand, cand_unfiltered))

  if(length(cand)>=5){

      el <- precompute_LD_edges(
        GTs = GTs[, cand, drop = FALSE],
        map = map[marker %in% cand],
        r2_min = min(rho_ld_grid),
        max_bp = max_gap)


    # i <- 100
    res <- rbindlist(
      mclapply(seq_len(nrow(param_grid)), function(i) {
        pars <- param_grid[i]

        out <- run_OR_scan_one(
          map = map,
          GTs = GTs[,cand],
          ld_ws = ld_ws,
          el = el,
          p_col = p_col,
          qt = pars$qt,
          alpha = pars$alpha,
          rho = pars$rho,
          rho_ld = pars$rho_ld,
          max_gap = max_gap,
          decay_sum = decay_sum,
          default_bp_th = default_bp_th,
          use_decay_threshold = TRUE
        )
      }, mc.cores = cores),
      fill = TRUE
    )
  }

}

compute_performance <- function(res,map,GTs){
  res[, param_id := paste(alpha, qt, rho, rho_ld, sep = "_")]
  res_split <- split(res, res$param_id)


  true_positives <- map[true_pos == TRUE, marker]
  n_true <- uniqueN(map[true_QTN == TRUE, marker])

  perf <- rbindlist(mclapply(res_split, function(par_comb) {

    OR <- par_comb$OR

    true_pos_OR <- sapply(OR, function(mrk) {
      any(mrk %in% true_positives)
    })

    QTN_f <- sapply(OR, function(mrk) {
      tmp <- map[marker %in% mrk & !is.na(focal_QTN)]
      if (nrow(tmp) == 0) return(NA_character_)
      tmp[which.max(max_LD_with_QTN), focal_QTN]
    })

    max_ld <- sapply(OR, function(mrk) {
      tmp <- map[marker %in% mrk]
      max(tmp$max_LD_with_QTN, na.rm = TRUE)
    })

    eval_dt <- data.table(
      QTN_f = QTN_f,
      true_pos_OR = true_pos_OR,
      max_ld = max_ld
    )

    # keep only the best OR per focal QTN
    setorder(eval_dt, QTN_f, -max_ld)
    eval_dt[, keep_TP := !duplicated(QTN_f) & !is.na(QTN_f)]

    TP <- eval_dt[, sum(true_pos_OR & keep_TP, na.rm = TRUE)]
    FP <- eval_dt[, sum(!true_pos_OR | !keep_TP, na.rm = TRUE)]

    captured_QTNs <- unique(eval_dt[true_pos_OR & keep_TP, QTN_f])
    FN <- n_true - length(captured_QTNs)

    out <- data.table(
      par_comb[1, .(qt, alpha, rho, rho_ld)],
      TP = TP,
      FP = FP,
      FN = FN,
      precision = TP / (TP + FP),
      recall = TP / (TP + FN)
    )

    out[, PR := precision * recall]
    out

  }, mc.cores = cores))
}


rho_ld_grid <- c(seq(0.1, 0.95, by = 0.1), 0.99)

qt_grid <- c(seq(0.1, 0.95, by = 0.05), 0.99)

alpha_grid <- c(
  0.05,
  0.01,
  0.005,
  0.001,
  0.0005,
  0.0001
)

cores = 8
folder <- "./single_SNP_results/"
out_folder <- "./OR_performance/"
files <- list.files(folder,full.names = TRUE)
done <- list.files(out_folder,full.names = TRUE)

files <- files[!basename(files) %in% basename(done)]
#file <- files[1]


for(file in files){
  message("Processing ",basename(file))

  data <- readRDS(file)

  rho_grid <- colnames(data$ld_ws)

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

  param_grid_R2_ld_w <- CJ(
    qt = qt_grid,
    rho = rho_grid
  )


  #i <- 1
  r2_F_ld <- rbindlist(mclapply(seq_len(nrow(param_grid_R2_ld_w)), function(i) {

    pars <- param_grid_R2_ld_w[i]

    data$map[,ld_w := data$ld_ws[, param_grid_R2_ld_w[i]$rho]]

    dt <- data$map[ld_w>quantile(ld_w, param_grid_R2_ld_w[i]$qt, na.rm = TRUE)]

    data.table(R2_F_ld=dt[,cor(lfmm_F,ld_w)^2],n_loci=nrow(dt))

  },mc.cores=cores))

  r2_F_ld[,E:=R2_F_ld*(n_loci/nrow(data$map))^2]
  r2_F_ld <- cbind(param_grid_R2_ld_w,r2_F_ld)


  decay_sum <- rbindlist(lapply(data$ld_decay_chr_9sp,function(x) x$decay_sum))

  message("Scanning ORs ")
  res_lfmm <- OR_scan(map=data$map,
                 GTs=data$GTs,
                 ld_ws=data$ld_ws,
                 p_col = "lfmm_p",
                 param_grid=param_grid,
                 alpha_grid = alpha_grid,
                 max_gap = 1e6,
                 decay_sum = decay_sum,
                 default_bp_th = Inf,
                 use_decay_threshold = TRUE,
                 cores = 8)

  res_emx <- OR_scan(map=data$map,
                        GTs=data$GTs,
                        ld_ws=data$ld_ws,
                        p_col = "emx_p",
                        param_grid,
                        max_gap = 1e6,
                        alpha_grid = alpha_grid,
                        decay_sum = decay_sum,
                        default_bp_th = Inf,
                        use_decay_threshold = TRUE,
                        cores = 8)

  message("Analysing performanc ")
  if(!is.null(res_lfmm)){
    perf_lfmm <- compute_performance(res=res_lfmm,
                                     map=data$map,
                                     GTs=data$GTs)
  }else{
    perf_lfmm <- NULL
  }


  if(!is.null(res_emx)){
    perf_emx <- compute_performance(res=res_emx,
                                    map=data$map,
                                    GTs=data$GTs)
  }else{
    perf_emx <- NULL
  }



  out <- list(lfmm=list(res=res_lfmm,perf=perf_lfmm),emx=list(res=res_emx,perf=perf_emx),r2_F_ld=r2_F_ld,param_grid=param_grid)

  saveRDS(out,paste0(out_folder,basename(file)))
}

if(FALSE){
  perf <- perf_lfmm

  perf_full <- merge(
    param_grid,
    perf,
    by = c("qt", "alpha", "rho", "rho_ld"),
    all.x = TRUE
  )

  baseline <- perf_full[
    rho == "none",
    .(PR0 = mean(PR)),
    by = .(alpha, rho_ld)
  ]

  tmp <- merge(
    perf_full[rho != "none"],
    baseline,
    by = c("alpha","rho_ld")
  )

  # C_scores <- tmp[, .(
  #   C = uniqueN(param_id) / nrow(param_grid)
  # ), by = marker]


  tmp[, delta_PR := PR - PR0]

  tmp2 <- merge(
    tmp,
    r2_F_ld,
    by = c("rho","qt")
  )

  #tmp2[,hist(delta_PR)]

  tmp2[,E:=R2_F_ld*(n_loci/nrow(map))^2]
  summary(lm(delta_PR ~ E, data = tmp2))

  ggplot(tmp2,aes(E,delta_PR)) +
    geom_point() +
    geom_vline(xintercept = median(tmp2$E))

  hm_red <- tmp2[alpha==0.05,
                 .(mean_delta_PR = mean(delta_PR, na.rm = TRUE),
                   mean_n_loci = mean(n_loci, na.rm = TRUE),
                   mean_cor = mean(R2_F_ld, na.rm = TRUE),
                   mean_precision = mean(precision, na.rm = TRUE),
                   mean_recall = mean(recall, na.rm = TRUE),
                   mean_E = mean(E, na.rm = TRUE)),

                 by = .(rho, qt)
  ]

  hm_red[, rho := factor(rho, levels = sort(unique(rho)))]
  hm_red[, qt := factor(qt, levels = sort(unique(qt)))]
  hm_red[, better_than_full := mean_delta_PR > 0]

  hm_full <- tmp2[alpha==0.05,
                  .(mean_delta_PR = mean(delta_PR, na.rm = TRUE),
                    mean_n_loci = mean(n_loci, na.rm = TRUE),
                    mean_cor = mean(R2_F_ld, na.rm = TRUE),
                    mean_precision = mean(precision, na.rm = TRUE),
                    mean_recall = mean(recall, na.rm = TRUE),
                    mean_E = mean(E, na.rm = TRUE)),

                  by = .(rho, qt)
  ]

  hm_full[, rho := factor(rho, levels = sort(unique(rho)))]
  hm_full[, qt := factor(qt, levels = sort(unique(qt)))]
  hm_full[, better_than_full := mean_delta_PR > 0]

  p1 <- ggplot(hm_full, aes(rho, qt, fill = mean_delta_PR)) +
    geom_tile(color="white",fill="white") +
    geom_tile(data=hm_red,color="white") +
    geom_tile(
      data = hm_red[mean_delta_PR > 0],
      fill = NA,
      color = "black",
      linewidth = 0.8
    ) +
    scale_fill_viridis_c(option="turbo") +
    theme_bw() +
    ggtitle("E>median(E)")

  p1

}
