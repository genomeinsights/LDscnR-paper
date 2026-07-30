#' Compute Per-SNP Local LD Support (ld_w)
#'
#' For each SNP and each requested relative LD threshold \eqn{\rho}, computes
#' local LD support as the median \eqn{r^2} with neighbouring SNPs within the
#' physical window \eqn{d_\rho} implied by that chromosome's fitted LD-decay
#' curve (see \code{compute_LD_decay()} / \code{d_from_rho()}).
#'
#' (Renamed from \code{precalculate_ld_w()} -- that name reflected an older
#' implementation and is no longer accurate.)
#'
#' @param rho Numeric vector of relative LD thresholds (0 < rho < 1).
#' @param ld_decay Object of class \code{"ld_decay"} produced by
#'   \code{compute_LD_decay()}, with \code{keep_el = TRUE} (or a valid
#'   \code{el_data_folder}) so that edge lists are available.
#' @param map Optional data.table with columns \code{marker} and \code{MAF},
#'   used only if \code{min_maf_ldw} is also supplied. MAF is matched to SNP
#'   ids in the edge list via \code{marker} (which should match \code{snp.id}
#'   in the GDS the edge lists were built from).
#' @param min_maf_ldw Optional MAF threshold. If supplied together with
#'   \code{map}, a SNP pair only contributes to the local-LD median if BOTH
#'   members have \code{MAF > min_maf_ldw} -- this both (a) excludes low-MAF
#'   neighbours from deflating a well-genotyped focal SNP's ld_w (low-MAF
#'   pairs mechanically cap r^2 regardless of true linkage), and (b) leaves
#'   \code{ld_w = NA} for any focal SNP whose own MAF falls below the
#'   threshold, since none of its pairs survive the filter. Default
#'   \code{NULL} applies no MAF filtering (matches prior behaviour).
#'
#' @return A numeric matrix of ld_w values, one row per SNP (explicitly
#'   \code{rownames()}'d by marker/snp id) and one column per requested
#'   \code{rho} (\code{colnames()} set to the rho values).
#'
#' @details
#' This filtering is independent of, and separate from, the MAF filter used
#' when fitting the decay curve itself (\code{min_maf_decay} in
#' \code{compute_LD_decay()}) -- the edge lists this function reads always
#' contain all SNPs at all MAFs; filtering happens only here, at the point of
#' summarizing local LD support.
#'
#' @export
compute_ld_w <- function(rho, ld_decay, map = NULL, min_maf_ldw = NULL) {

  message("Computing ld_w")

  if (!is.null(min_maf_ldw) && is.null(map)) {
    stop("`map` (with columns `marker` and `MAF`) must be supplied when `min_maf_ldw` is set.")
  }

  maf_vec <- NULL
  if (!is.null(map)) {
    if (!all(c("marker", "MAF") %in% names(map))) {
      stop("`map` must contain columns `marker` and `MAF`.")
    }
    maf_vec <- stats::setNames(map$MAF, map$marker)
  }

  if (length(rho) > 1) {
    pb <- txtProgressBar(min = 0, max = length(rho) - 1, style = 3)
    setTxtProgressBar(pb, 0)
  }

  ld_ws <- do.call(rbind, lapply(ld_decay$by_chr, function(chr_obj) {

    a <- ld_decay$decay_sum[Chr == chr_obj$decay_sum$Chr, a_pred]
    b <- ld_decay$decay_sum[Chr == chr_obj$decay_sum$Chr, b]
    d_window <- d_from_rho(a, rho)

    if (is.null(chr_obj$el)) stop("No edge list present")
    if (is.character(chr_obj$el)) chr_obj$el <- data.table::fread(chr_obj$el, showProgress = FALSE)

    el_sym <- data.table::rbindlist(list(
      chr_obj$el[, .(SNP = SNP1, pos = pos1, pos_other = pos2, SNP_other = SNP2, r2, d)],
      chr_obj$el[, .(SNP = SNP2, pos = pos2, pos_other = pos1, SNP_other = SNP1, r2, d)]
    ))

    if (!is.null(min_maf_ldw)) {
      maf_ok <- maf_vec[el_sym$SNP] > min_maf_ldw & maf_vec[el_sym$SNP_other] > min_maf_ldw
      maf_ok[is.na(maf_ok)] <- FALSE
      el_sym <- el_sym[maf_ok]
    }

    ld_w <- do.call(cbind, lapply(d_window, function(win) {
      ld_w_win <- el_sym[d < win, .(r2_median = median(r2)), by = SNP]
      ld_w_win[match(chr_obj$snp_ids, ld_w_win$SNP)]$r2_median
    }))

    ## explicit rownames -- removes the implicit dependency on ld_ws row
    ## order matching map/GTs column order elsewhere in the pipeline
    rownames(ld_w) <- chr_obj$snp_ids

    if (length(rho) > 1) {
      setTxtProgressBar(pb, which(names(ld_decay$by_chr) == chr_obj$decay_sum$Chr))
    }

    return(ld_w)
  }))

  if (length(rho) > 1) close(pb)
  colnames(ld_ws) <- rho

  return(ld_ws)
}
