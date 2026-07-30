library(data.table)
library(SNPRelate)
# Helpers -----------------------------------------------------------------

default_cluster_colours <- function() {
  c(
    "#B2DF8A", "#FFD92F", "firebrick", "#33A02C", "#7FC97F", "#CAB2D6",
    "#FB8072", "grey30", "#E6AB02", "#FDC086", "steelblue", "#1F78B4",
    "#FB9A99", "#1B9E77", "#BC80BD", "#E31A1C", "#7570B3", "#A6761D",
    "#A6CEE3", "salmon", "#FFFF33", "forestgreen", "#FDCDAC", "#BF5B17",
    "#A6761D", "#FBB4AE", "#4DAF4A", "#B3E2CD", "#FDDAEC", "#BEBADA",
    "#FFF2AE", "#1F78B4", "#66C2A5", "#F0027F", "#E6AB02", "#E78AC3",
    "#FF7F00", "#8DA0CB", "#6A3D9A", "#B15928", "#E41A1C"
  )
}


detect_gt_input <- function(x, tol = 1e-8) {
  vals <- as.numeric(x)
  vals <- vals[is.finite(vals)]

  if (!length(vals)) return("hard")

  is_hard <- all(abs(vals - round(vals)) < tol) &&
    all(vals %in% c(0, 1, 2))

  if (is_hard) "hard" else "dosage"
}


polarize_genotypes <- function(x, cor_threshold = 0) {
  x <- as.matrix(x)

  if (ncol(x) <= 1L) return(x)

  ref_idx <- which.max(colSums(!is.na(x)))
  ref <- x[, ref_idx]

  flip <- vapply(seq_len(ncol(x)), function(j) {
    r <- suppressWarnings(
      stats::cor(ref, x[, j], use = "pairwise.complete.obs")
    )
    is.finite(r) && r < cor_threshold
  }, logical(1))

  if (any(flip)) {
    x[, flip] <- 2 - x[, flip]
  }

  x
}

expected_gt_hard <- function(x) {
  x <- polarize_genotypes(x)

  n <- rowSums(!is.na(x))
  c1 <- rowSums(x == 1, na.rm = TRUE)
  c2 <- rowSums(x == 2, na.rm = TRUE)

  E <- (c1 + 2 * c2) / n
  E[n == 0] <- NA_real_

  list(
    expected = E,
    n = n
  )
}

#x <- GTs[,cls[[2]]]
select_representative_snp <- function(x) {
  x <- polarize_genotypes(x)

  eg <- expected_gt_hard(x)
  E <- eg$expected
  n <- eg$n

  # Consensus genotype per individual
  consensus <- round(E)

  # Replace missing genotypes with consensus
  na_idx <- is.na(x)
  if (any(na_idx)) {
    row_ids <- row(x)[na_idx]
    x[na_idx] <- consensus[row_ids]
  }

  # Correlation to cluster-wide expected genotype
  cors <- apply(x, 2, function(g) {
    ok <- !is.na(E)
    if (sum(ok) < 3) return(NA_real_)
    cor(g[ok], E[ok])^2
  })

  best <- which.max(cors)

  # Representative SNP genotype with imputed missing values
  rep_gt <- x[, best]

  cluster_gt <- round(E)
  #cluster_gt[n < min_support] <- NA_integer_

  cor_E_hard <- cor(E, cluster_gt, use = "complete.obs")^2

  list(
    snp = colnames(x)[best],
    index = unname(best),
    MLG = cluster_gt,
    cor_best = cors[best],
    cor_E_hard = cor_E_hard,
    representative = rep_gt,
    eMLG = E,
    n = n
  )
}


score_eMLG <- function(x) {
  r2 <- suppressWarnings(
    stats::cor(round(x), x, use = "pairwise.complete.obs")^2
  )

  if (!is.finite(r2)) NA_real_ else r2
}


weighted_row_mean <- function(x, w) {
  x <- as.matrix(x)
  w <- as.numeric(w)

  ok <- !is.na(x)
  num <- rowSums(sweep(x, 2, w, `*`), na.rm = TRUE)
  den <- rowSums(sweep(ok, 2, w, `*`), na.rm = TRUE)

  y <- num / den
  y[den == 0] <- NA_real_

  y
}


cluster_level_map <- function(map_snp, id_col = "CL_id") {
  data.table::setDT(map_snp)

  map_snp[
    !is.na(get(id_col)),
    {
      pos <- as.numeric(Pos)

      .(
        Chr = Chr[1],
        Pos_min = min(pos, na.rm = TRUE),
        Pos_max = max(pos, na.rm = TRUE),
        Pos_mid = median(pos, na.rm = TRUE),
        span_bp = max(pos, na.rm = TRUE) - min(pos, na.rm = TRUE),
        n_loci = .N,
        rep_snp = rep_snp[1],
        r2_eMLG = if ("r2_eMLG" %in% names(.SD)) r2_eMLG[1] else NA_real_
      )
    },
    by = id_col
  ]
}

add_snp_positions <- function(map_snp, map_ref) {
  map_ref <- data.table::as.data.table(map_ref)
  map_snp <- data.table::as.data.table(map_snp)

  pos_cols <- intersect(c("marker", "Chr", "Pos"), names(map_ref))

  merge(
    unique(map_ref[, ..pos_cols]),
    map_snp,
    by = "marker",
    all.y = TRUE,
    sort = FALSE
  )
}

# Main functions -----------------------------------------------------------------
LD_clustering <- function(ld_decay,
                          map,
                          SNPs   = NULL,
                          rho_ld = 0.95,
                          rho_d  = 0.95,
                          ld_th  = NULL,
                          d_th   = NULL,
                          col_vector = NULL,
                          l_min  = 10,
                          cores  = 1) {
  data.table::setDT(map)

  #ch = "Chr1"
  out <- data.table::rbindlist(parallel::mclapply(names(ld_decay$by_chr), function(ch) {
    message(ch)
    map <- data.table::copy(map[Chr == ch])
    chr_obj <- ld_decay$by_chr[[ch]]
    if(is.data.frame(chr_obj$el)){
      el <- chr_obj$el
    }else{
      el <- data.table::fread(chr_obj$el, showProgress = FALSE)
    }



    ld_th_chr <- ld_th
    d_th_chr <- d_th

    if (is.null(ld_th_chr)) {
      b_chr <- ld_decay$decay_sum[Chr == ch, b]
      c_chr <- ld_decay$decay_sum[Chr == ch, c_pred]
      ld_th_chr <- ld_from_rho(b_chr, c_chr, rho = rho_ld)
    }

    if (is.null(d_th_chr)) {
      a_chr <- ld_decay$decay_sum[Chr == ch, a_pred]
      d_th_chr <- d_from_rho(a_chr, rho = rho_d)
    }

    if(!is.null(SNPs)){
      ed <- el[
        r2 > ld_th_chr &
          d < d_th_chr &
          (SNP1 %in% SNPs | SNP2 %in% SNPs),
        .(SNP1, SNP2)
      ]
    }else{
      ed <- el[
        r2 > ld_th_chr &
          d < d_th_chr,
        (SNP1 %in% SNPs | SNP2 %in% SNPs),
        .(SNP1, SNP2)
      ]
    }


    map[, `:=`(CL_id = NA_character_, CL_col = "grey")]

    if (nrow(ed) == 0L) {
      return(map)
    }

    g <- igraph::graph_from_data_frame(ed, directed = FALSE)
    comps <- igraph::components(g)
    cls <- split(names(comps$membership), comps$membership)
    cls <- cls[vapply(cls, length, integer(1)) >= l_min]

    if (!length(cls)) {
      return(map)
    }

    #cl_cols <- rep(col_vector, length.out = length(cls))

    cl_dt <- data.table::data.table(
      CL_id = rep(paste0(ch, "_", seq_along(cls)), lengths(cls)),
      marker = unlist(cls, use.names = FALSE),
      #CL_col = rep(cl_cols, lengths(cls)),
      n_loci = rep(lengths(cls),lengths(cls))
    )

    map[cl_dt, `:=`(CL_id = i.CL_id, n_loci=i.n_loci), on = "marker"]

  },mc.cores=cores))
  out[is.na(n_loci),CL_id:=marker]
  out[is.na(n_loci),n_loci:=1]
}

#' Construct expected multi-locus genotypes from LD clusters
#'
#' Converts SNP clusters into expected multi-locus genotypes (eMLGs) and
#' hard-coded consensus multi-locus genotypes (MLGs). For each cluster,
#' genotypes are polarized, summarized as a cluster-wide expected dosage,
#' and represented by the SNP most strongly correlated with this expected
#' genotype profile.
#'
#' @param cls Named list of SNP clusters. Each element contains marker IDs
#'   belonging to one LD cluster.
#' @param geno Genotype matrix with individuals in rows and SNPs in columns.
#'   Genotypes should be coded as allele dosages `0`, `1`, or `2`, with
#'   missing values allowed.
#' @param map_cl Marker map used to add chromosome and position information.
#' @param cor_th LD threshold used to define clusters.
#' @param cores Number of cores used for parallel processing.
#'
#' @return A list containing:
#' \describe{
#'   \item{eMLG}{Matrix of expected cluster genotypes.}
#'   \item{MLG}{Matrix of hard-coded consensus cluster genotypes.}
#'   \item{imp_rep_gt}{Matrix of representative SNP genotypes with missing
#'     values imputed from the cluster consensus.}
#'   \item{map_SNP}{Data table linking SNPs to clusters and positions.}
#'   \item{map_eMLG}{Data table describing each eMLG cluster, including genomic
#'     range, number of loci, representative SNP, and QC metrics.}
#'   \item{weight_mat}{Matrix of per-individual SNP support for each cluster.}
#'   \item{clusters}{Input SNP clusters.}
#'   \item{cor_th}{LD threshold used to generate clusters.}
#' }
#'
#' @details
#' For each SNP cluster, an expected genotype dosage is calculated for each
#' individual as the mean allele dosage across non-missing SNPs in the cluster.
#' A hard-coded MLG is then obtained by rounding this expected dosage to the
#' nearest genotype class. Missing representative SNP genotypes are imputed
#' using the rounded cluster consensus genotype.
#'
#' The representative SNP is selected as the marker with the highest Pearson
#' correlation to the cluster-wide expected genotype profile. This avoids
#' computing the full pairwise SNP correlation matrix and provides an efficient
#' approximation to choosing the SNP with the highest mean correlation to all
#' other SNPs in the cluster.
#'
#' The output includes QC metrics such as the correlation between the
#' representative SNP and the expected genotype, and the correlation between
#' the expected genotype and its rounded hard-call consensus.
#'
#' @export
make_eMLGs <- function(GTs,
                       map_cl,
                       cor_th = 0.8,
                       l_min = 0,
                       cores = 1) {



  map_cl <- data.table::copy(map_cl)
  data.table::setDT(map_cl)

  map_cl <- map_cl[!is.na(CL_id) & n_loci >= l_min]

  cls <- split(map_cl$marker, map_cl$CL_id)
  cls <- cls[lengths(cls) >= max(2,l_min)]

  if (!length(cls)) {
    stop("No LD clusters passed l_min.")
  }

  # #input = "rSNP"
  # expected_fun <- switch(
  #   input,
  #   hard = expected_gt_hard,
  #   dosage = expected_gt_dosage,
  #   rSNP = select_representative_snp
  # )

  process_cluster <- function(markers) {
    select_representative_snp(GTs[, markers, drop = FALSE])
  }


  if(length(cls)>1000){
    batch_size <- 1000

    idx_chunks <- split(seq_along(cls),ceiling(seq_along(cls) / batch_size))

    n_batches <- length(idx_chunks)

    #t1 <- Sys.time()
    res <- lapply(seq_along(idx_chunks), function(i) {
      idx <- idx_chunks[[i]]

      message(sprintf(
        "Starting batch %d/%d: entries %d-%d",
        i, n_batches, min(idx), max(idx)
      ))


      out <- parallel::mclapply(
        cls[idx],
        process_cluster,
        mc.cores = cores
      )
    })
    res <- unlist(res, recursive = FALSE)
  }else{
    message(sprintf(
      "Starting batch 1: entries %d-%d",1, length(cls)
    ))

    res <- parallel::mclapply(
      cls,
      process_cluster,
      mc.cores = cores
    )
  }

  # t2 <- Sys.time()



  #res$Chr1_1$
  # Representative SNP names
  rep_snps <- vapply(res, `[[`, character(1), "snp")

  # Representative SNP column indices
  rep_indices <- vapply(res, `[[`, integer(1), "index")

  # Correlation of representative SNP with expected genotype
  cor_best <- vapply(res, `[[`, numeric(1), "cor_best")
  cor_E_hard <- vapply(res, `[[`, numeric(1), "cor_E_hard")

  # Expected genotype matrix
  # rows = individuals, columns = clusters
  eMLG <- do.call(cbind, lapply(res, `[[`, "eMLG"))
  colnames(eMLG) <- names(res)

  MLG <- do.call(cbind, lapply(res, `[[`, "MLG"))
  colnames(MLG) <- names(res)

  # Representative SNP genotype matrix (NA replaced by consensus)
  rep_gt_mat <- do.call(cbind, lapply(res, `[[`, "representative"))
  colnames(rep_gt_mat) <- rep_snps

  # Weight matrix: number of SNPs contributing per individual
  # rows = individuals, columns = clusters
  weight_mat <- do.call(cbind, lapply(res, `[[`, "n"))
  colnames(weight_mat) <- names(res)


  map_snp <- data.table::data.table(
    CL_id = rep(names(cls), lengths(cls)),
    marker = unlist(cls, use.names = FALSE)
  )

  map_snp <- add_snp_positions(map_snp, map_cl)

  map_eMLG <- data.table(CL_id = names(cls),
                         Chr = map_snp[match(rep_snps,marker),Chr],
                         n_loci = lengths(cls),
                         rep_snp = rep_snps,
                         cor_best = cor_best,
                         cor_E_hard = cor_E_hard)


  positions <- map_snp[,.(Chr=Chr[1],min_pos=min(Pos),max_pos=max(Pos)),by=CL_id]
  positions[,range_pos:=max_pos-min_pos]

  map_eMLG <- positions[map_eMLG,on="CL_id"]

  list(
    eMLG = eMLG,
    MLG = MLG,
    imp_rep_gt = rep_gt_mat,
    map_SNP = map_snp,
    map_eMLG = map_eMLG,
    weight_mat = weight_mat,
    clusters = cls,
    cor_th = cor_th
  )
}


#' Collapse nearby correlated expected multi-locus genotypes
#'
#' Collapses expected multi-locus genotype clusters within chromosomes when
#' they are both physically close and highly correlated. This provides a
#' second-level reduction of correlated haplotype-like genotype blocks.
#'
#' @param MLGs Output list from `make_eMLGs()`.
#' @param distance_threshold Maximum physical distance, in base pairs, allowed
#'   between neighboring eMLGs for them to be considered part of the same
#'   candidate run. Default is `5e5`.
#' @param r2_threshold Minimum squared correlation required to collapse eMLGs.
#'   Default is `0.8`.
#' @param method Linkage method passed to [stats::hclust()]. Default is
#'   `"complete"`.
#' @param prefix Prefix for collapsed eMLG identifiers. Default is `"C_eMLG"`.
#'
#' @return A list containing:
#' \describe{
#'   \item{eMLG}{Matrix of collapsed hard-coded consensus eMLGs.}
#'   \item{map_eMLG}{Data table describing collapsed eMLGs, including member
#'     counts, number of loci, LD summaries, and QC metrics.}
#'   \item{map_SNP}{Data table linking SNPs to collapsed eMLGs.}
#'   \item{lookup}{Data table linking original eMLGs to collapsed eMLGs.}
#'   \item{clusters}{List of original eMLG IDs belonging to each collapsed
#'     eMLG.}
#'   \item{params}{List of parameters used for the collapse step.}
#' }
#'
#' @details
#' eMLGs are first ordered by chromosome and genomic position. Consecutive
#' eMLGs separated by more than `distance_threshold` base pairs are assigned to
#' different candidate runs. Within each run, pairwise squared correlations are
#' calculated between eMLG genotype vectors, and hierarchical clustering is
#' performed on the distance matrix `1 - r^2`.
#'
#' Groups are cut at `1 - r2_threshold`. With the default complete-linkage
#' method, all members of a collapsed group must satisfy the specified
#' correlation threshold to the group under the complete-linkage criterion.
#'
#' For collapsed groups containing multiple eMLGs, genotypes are summarized by
#' a weighted mean across member eMLGs, using the number of SNP loci in each
#' eMLG as weights. The weighted expected dosage is then rounded to produce a
#' hard-coded collapsed genotype.
#'
#' This step reduces pseudo-replication among nearby correlated haplotype-like
#' blocks while retaining QC metrics such as mean and minimum pairwise
#' \(r^2\), number of member eMLGs, total number of SNP loci, and hard-call
#' conformity.
#'
#' @export
collapse_eMLGs_by_chr_window <- function(MLGs,
                                         distance_threshold = 5e5,
                                         r2_threshold = 0.8,
                                         method = "complete",
                                         prefix = "C_eMLG") {
  eMLG <- as.matrix(MLGs$eMLG)

  map <- data.table::copy(MLGs$map_eMLG)
  data.table::setDT(map)

  # Standardize position column names
  if ("min_pos" %in% names(map)) data.table::setnames(map, "min_pos", "Pos_min")
  if ("max_pos" %in% names(map)) data.table::setnames(map, "max_pos", "Pos_max")

  data.table::setorderv(map, c("Chr", "Pos_min"))

  map[, gap_bp := Pos_min - data.table::shift(Pos_max), by = Chr]
  map[, new_run := is.na(gap_bp) | gap_bp > distance_threshold, by = Chr]
  map[, run := cumsum(new_run), by = Chr]

  offdiag_stats <- function(r2mat) {
    if (ncol(r2mat) <= 1L) {
      return(list(mean = 1, min = 1))
    }

    vals <- r2mat[upper.tri(r2mat)]
    vals <- vals[is.finite(vals)]

    list(
      mean = mean(vals, na.rm = TRUE),
      min = min(vals, na.rm = TRUE)
    )
  }

  collapse_run <- function(gmap) {
    ids <- intersect(gmap$CL_id, colnames(eMLG))

    if (length(ids) == 0L) return(NULL)

    if (length(ids) == 1L) {
      y <- eMLG[, ids]

      return(list(
        C_eMLG = matrix(y, ncol = 1, dimnames = list(rownames(eMLG), ids)),
        map = data.table::data.table(
          old_C_eMLG_id = ids,
          n_eMLGs = 1L,
          n_loci = gmap[match(ids, CL_id), n_loci],
          mean_r2_to_members = 1,
          min_r2_to_members = 1,
          r2_eMLG = score_eMLG(y),
          members = list(ids)
        )
      ))
    }

    x <- eMLG[, ids, drop = FALSE]

    r2 <- suppressWarnings(stats::cor(x, use = "pairwise.complete.obs")^2)
    r2[!is.finite(r2)] <- 0
    diag(r2) <- 1

    hc <- stats::hclust(
      stats::as.dist(1 - r2),
      method = method
    )

    groups <- split(ids, stats::cutree(hc, h = 1 - r2_threshold))

    score_group <- function(members) {
      z <- x[, members, drop = FALSE]

      if (length(members) == 1L) {
        y <- as.numeric(z[, 1])
        r2_members <- matrix(1, 1, 1)
      } else {
        w <- map[match(members, CL_id), n_loci]
        w[is.na(w)] <- 1

        E <- weighted_row_mean(
          polarize_genotypes(z),
          w
        )

        y <- round(E)
        y[is.na(E)] <- NA_real_

        r2_members <- suppressWarnings(
          stats::cor(z, use = "pairwise.complete.obs")^2
        )
        r2_members[!is.finite(r2_members)] <- NA_real_
        diag(r2_members) <- 1
      }

      od <- offdiag_stats(r2_members)

      list(
        y = y,
        n_loci = sum(map[match(members, CL_id), n_loci], na.rm = TRUE),
        mean_r2_to_members = od$mean,
        min_r2_to_members = od$min,
        r2_eMLG = score_eMLG(y)
      )
    }

    scores <- lapply(groups, score_group)

    C_eMLG <- do.call(cbind, lapply(scores, `[[`, "y"))

    run_prefix <- paste0(gmap$Chr[1], "_run", gmap$run[1], "_")
    colnames(C_eMLG) <- paste0(run_prefix, seq_along(groups))
    rownames(C_eMLG) <- rownames(eMLG)

    map_out <- data.table::data.table(
      old_C_eMLG_id = colnames(C_eMLG),
      n_eMLGs = lengths(groups),
      n_loci = vapply(scores, `[[`, numeric(1), "n_loci"),
      mean_r2_to_members = vapply(scores, `[[`, numeric(1), "mean_r2_to_members"),
      min_r2_to_members = vapply(scores, `[[`, numeric(1), "min_r2_to_members"),
      r2_eMLG = vapply(scores, `[[`, numeric(1), "r2_eMLG"),
      members = I(groups)
    )

    list(C_eMLG = C_eMLG, map = map_out)
  }

  runs <- split(map, by = c("Chr", "run"), keep.by = TRUE)

  runs <- split(map, by = c("Chr", "run"), keep.by = TRUE)

  chromosomes <- unique(map$Chr)
  n_chr <- length(chromosomes)

  res <- list()

  for (i in seq_along(chromosomes)) {
    chr <- chromosomes[i]

    message(sprintf(
      "Starting chromosome %s (%d/%d)",
      chr, i, n_chr
    ))

    chr_runs <- runs[vapply(runs, function(z) z$Chr[1] == chr, logical(1))]

    res_chr <- lapply(chr_runs, collapse_run)

    res <- c(res, res_chr)
  }

  res <- Filter(Negate(is.null), res)

  C_eMLG <- do.call(cbind, lapply(res, `[[`, "C_eMLG"))
  map_out <- data.table::rbindlist(lapply(res, `[[`, "map"), fill = TRUE)

  new_ids <- paste0(prefix, "_", seq_len(ncol(C_eMLG)))

  colnames(C_eMLG) <- new_ids
  map_out[, C_eMLG_id := new_ids]

  lookup <- data.table::data.table(
    C_eMLG_id = rep(map_out$C_eMLG_id, map_out$n_eMLGs),
    CL_id = unlist(map_out$members, use.names = FALSE)
  )

  map_snp_out <- merge(
    MLGs$map_SNP,
    lookup,
    by = "CL_id",
    all.x = FALSE,
    all.y = TRUE,
    sort = FALSE
  )

  list(
    eMLG = C_eMLG,
    map_eMLG = map_out,
    map_SNP = map_snp_out,
    lookup = lookup,
    clusters = split(lookup$CL_id, lookup$C_eMLG_id),
    params = list(
      stage = "collapse_eMLGs_by_chr_window",
      distance_threshold = distance_threshold,
      r2_threshold = r2_threshold,
      method = method
    )
  )
}

collapse_eMLGs_global <- function(MLGs,
                                  r2_threshold = 0.8,
                                  method = "complete",
                                  prefix = "G_eMLG") {
  eMLG <- as.matrix(MLGs$eMLG)

  r2 <- suppressWarnings(
    stats::cor(eMLG, use = "pairwise.complete.obs")^2
  )

  r2[!is.finite(r2)] <- 0
  diag(r2) <- 1

  hc <- stats::hclust(
    stats::as.dist(1 - r2),
    method = method
  )

  groups <- split(
    colnames(eMLG),
    stats::cutree(hc, h = 1 - r2_threshold)
  )

  members <- groups[[1]]
  MLGs$map_eMLG[match(members, C_eMLG_id)]

  score_group <- function(membrs) {
    x <- eMLG[, membrs, drop = FALSE]

    n_loci_group <- sum(
      MLGs$map_eMLG[
        match(membrs, C_eMLG_id),
        n_loci
      ],
      na.rm = TRUE
    )

    if (length(membrs) == 1L) {
      y <- as.numeric(x[, 1])
      r2_members <- 1
    } else {

      w <- MLGs$map_eMLG[
        match(membrs, C_eMLG_id),
        n_loci
      ]

      w[is.na(w)] <- 1

      y <- weighted_row_mean(
        polarize_genotypes(x),
        w
      )

      r2_members <- suppressWarnings(
        stats::cor(x, use = "pairwise.complete.obs")^2
      )
    }

    list(
      y = y,
      n_loci = n_loci_group,
      mean_r2_to_members = mean(r2_members, na.rm = TRUE),
      min_r2_to_members = min(r2_members, na.rm = TRUE),
      r2_eMLG = score_eMLG(y)
    )
  }

  scores <- lapply(groups, score_group)

  G_eMLG <- do.call(cbind, lapply(scores, `[[`, "y"))
  G_ids <- paste0(prefix, "_", seq_along(groups))

  colnames(G_eMLG) <- G_ids
  rownames(G_eMLG) <- rownames(eMLG)

  lookup <- data.table::data.table(
    G_eMLG_id = rep(G_ids, lengths(groups)),
    C_eMLG_id = unlist(groups, use.names = FALSE)
  )

  map_eMLG <- data.table::data.table(
    G_eMLG_id = G_ids,
    n_C_eMLGs = lengths(groups),
    n_loci = vapply(scores, `[[`, numeric(1), "n_loci"),
    mean_r2_to_members = vapply(scores, `[[`, numeric(1), "mean_r2_to_members"),
    min_r2_to_members = vapply(scores, `[[`, numeric(1), "min_r2_to_members"),
    r2_eMLG = vapply(scores, `[[`, numeric(1), "r2_eMLG"),
    members = I(groups)
  )

  map_snp <- merge(
    MLGs$map_SNP,
    lookup,
    by = "C_eMLG_id",
    all.x = FALSE,
    all.y = TRUE,
    sort = FALSE
  )

  data.table::setcolorder(
    map_snp,
    intersect(
      c(
        "marker", "Chr", "Pos",
        "CL_id", "C_eMLG_id", "G_eMLG_id",
        "CL_col", "n_loci", "r2_eMLG"
      ),
      names(map_snp)
    )
  )

  list(
    eMLG = G_eMLG,
    map_eMLG = map_eMLG,
    map_SNP = map_snp,
    clusters = split(lookup$C_eMLG_id, lookup$G_eMLG_id),
    params = list(
      stage = "collapse_eMLGs_global",
      r2_threshold = r2_threshold,
      method = method
    )
  )
}

annotate_C_pair_distance <- function(map_SNP) {
  C_map <- map_SNP[
    ,
    .(
      Chr = Chr[1],
      Pos_min = as.numeric(min(Pos, na.rm = TRUE)),
      Pos_max = as.numeric(max(Pos, na.rm = TRUE)),
      Pos_mid = as.numeric(stats::median(Pos, na.rm = TRUE))
    ),
    by = C_eMLG_id
  ]

  pairs <- data.table::CJ(
    C1 = C_map$C_eMLG_id,
    C2 = C_map$C_eMLG_id
  )[C1 < C2]

  pairs[C_map, `:=`(
    Chr1 = i.Chr,
    Pos_min1 = i.Pos_min,
    Pos_max1 = i.Pos_max
  ), on = .(C1 = C_eMLG_id)]

  pairs[C_map, `:=`(
    Chr2 = i.Chr,
    Pos_min2 = i.Pos_min,
    Pos_max2 = i.Pos_max
  ), on = .(C2 = C_eMLG_id)]

  pairs[
    ,
    same_chr := Chr1 == Chr2
  ]

  pairs[
    ,
    gap_bp := data.table::fifelse(
      same_chr,
      pmax(0, pmax(Pos_min1, Pos_min2) - pmin(Pos_max1, Pos_max2)),
      Inf
    )
  ]

  pairs[
    ,
    overlap := same_chr & gap_bp == 0
  ]

  pairs[]
}

# read in data ----------------------------------------------------------------
folder <- "./single_SNP_results/"
#out_folder <- "./OR_performance2/"
files <- list.files(folder,full.names = TRUE)
#done <- list.files(out_folder,full.names = TRUE)
#files <- files[!basename(files) %in% basename(done)]

file <- files[1]
#for(file in files){
  message("Processing ",basename(file))

  data <- readRDS(file)

  env <- data$env
  ld_decay_chr <- data$ld_decay_chr_9sp[[repl]]
  map_chr <- data$map[Chr_9sp==repl]
  GTs_chr <- data$GTs[,data$map[,Chr_9sp==repl]]


  #ld_decay_chr

  map_cl <- LD_clustering(ld_decay = ld_decay_chr,
                             map=map_chr,
                             ld_th = 0.8,
                             d_th=1e6,
                             l_min=1,
                             cores = 4)

  # map_cl[is.na(n_loci),n_loci:=1]
  # map_cl[n_loci==1,CL_id:=marker]
  #
  # map_cl[type=="QTN",]


  ## keep only hybrids without Sielva
  #keep_inds <- sample_info[,which(Species=="hybrid" & Population != "Sielva")]

  eMLGs <- make_eMLGs(GTs=GTs_chr,
                      map_cl = map_cl,
                      cor_th = 0.8,
                      l_min=2,
                      cores = 8)

  #eMLGs$map_eMLG[n_]


  #saveRDS(eMLGs_1mb,"eMLGs_1mb.rds")
  #eMLGs_1mb_hyb$map_eMLG
  #eMLGs$map_eMLG[,plot(cor_best, cor_E_hard)]
  #abline(0,1)
  #This approach prioritizes robust LD-supported haplotype signals and may reduce sensitivity to true isolated causal SNPs.
  C_eMLGs <- collapse_eMLGs_by_chr_window(MLGs = eMLGs,
                                                  distance_threshold = 5e5,
                                                  r2_threshold = 0.8,
                                                  method = "complete",
                                                  prefix = "C_eMLG")


  gds_path = tempfile(fileext = ".gds")

  gds <- create_gds_from_geno(geno = GTs_chr,map=map_chr,gds_path)
  on.exit({ snpgdsClose(gds); unlink(gds_path) }, add = TRUE)

  non_clustered_SNPs <- map_chr[!marker %in% eMLGs$map_SNP$marker,marker]

  gts_reduced <- cbind(eMLGs$eMLG,GTs_chr[,non_clustered_SNPs])

  GRM <- snpgdsGRM(gds,method="GCTA",autosome.only = FALSE)$grm
  emx_eMLG <- emmax(env,X = gts_reduced,K=GRM)
  length(emx_eMLG$pval)

  map_chr[,plot(-log10(p.adjust(emx_p,"fdr")))]
  plot(-log10(p.adjust(emx_eMLG$pval,"fdr")))


  tmp <- apply(data$ld_ws,2,function(x){
    cor(x,data$map[,max_LD_with_QTN],use="pair")^2
  })

  tmp2 <- apply(data$ld_ws,2,function(x){
    cor(x,data$map[,lfmm_F],use="pair")^2
  })
  ld_w <- data$ld_ws[,which.max(tmp2)]

  ld_w
  qt = seq(0,0.95,by=0.05)
  q = 0.9
  tmp <- sapply(qt,function(q){
    keep <- which(ld_w>quantile(ld_w,q,na.rm=TRUE))
    cor(ld_w[keep],data$map[keep,max_LD_with_QTN],use="pair")^2
  })

  tmp2 <- sapply(qt,function(q){
    keep <- which(ld_w>quantile(ld_w,q,na.rm=TRUE))
    cor(ld_w[keep],data$map[keep,lfmm_F],use="pair")^2
  })

  plot(tmp2,tmp)
  qt_opt <- qt[which.max(tmp2)]

  data$map[,indx:=.I]
  keep <- which(ld_w>quantile(ld_w,0.5,na.rm=TRUE))
  par(mfcol=c(2,1))
  data$map[keep,plot(indx,-log10(p.adjust(lfmm_p,"fdr")),pch=ifelse(type=="QTN",3,20))]
  data$map[,plot(indx,-log10(p.adjust(lfmm_p,"fdr")),pch=ifelse(type=="QTN",3,20))]

  ## close at exit


