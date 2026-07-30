######################################################
## Outlier-region (OR) analysis for the pooled
## simulated genome, reusing the LD-filtering machinery
## from ld_w_filtering_3sp.R (precompute_LD_edges,
## LD_igraph_components, run_one_grid, get_potential_outliers).
##
## Source ld_w_filtering_3sp.R BEFORE this script -- it relies
## on those functions being defined.
##
## Adds:
##   - build_shared_el(): precomputes one LD edge list ("el")
##                        up front, covering candidate markers
##                        from EITHER FDR mode, so run_one_grid
##                        doesn't rebuild it per grid point
##   - run_grid_sim():    runs the (rho x th_ldw x r2_th x l_min)
##                        grid with FDR done either "per_sim"
##                        (within each original 2-chromosome
##                        replicate, matching how emx_p/lfmm_p
##                        were actually computed) or "pooled"
##                        (across the whole concatenated genome)
##   - classify_ORs():    labels each called outlier region (OR)
##                        as true/false using map$true_pos, and
##                        computes OR-level recall / OR_FDR
######################################################

library(data.table)
library(igraph)
library(parallel)

#----------------------------------------------------------
# 0) Load one pooled/concatenated simulation genome
#----------------------------------------------------------
# conc_data <- readRDS("./parsed_sim_data_genomes/adapt_bgs_V2_c1_env1_genome.rds")
# map   <- conc_data$map
# GTs   <- conc_data$GTs      # NOTE: run_one_grid()'s on-the-fly `el`
#                             # builder (used only when el = NULL) looks
#                             # up a GLOBAL variable literally named
#                             # `GTs` -- inherited from how
#                             # ld_w_filtering_3sp.R was written. We
#                             # build `el` up front below, so this is
#                             # only a fallback path, but keep the name
#                             # if you ever hit it.
# ld_ws <- conc_data$ld_ws    # rows named by marker (from the
#                             # concat_genome_replicates.R script)

# map[, indx    := .I]
# map[, sim_id  := ceiling(Chr_9sp / 2)]     # groups the QTN+neutral
#                                             # chromosome pair back into
#                                             # one original 2-chromosome
#                                             # replicate (1..10)
# map[, true_pos := rho_d <= 0.99 & ld_rel > 0.25]
# map[chr_type == "ntrl", true_pos := FALSE]

#----------------------------------------------------------
# 1) Build one shared LD edge list for clustering, covering
#    every marker that could plausibly end up in an outlier
#    region under EITHER fdr_mode
#----------------------------------------------------------
build_shared_el <- function(map, GTs, ld_ws,
                            p_cols, th_ldw_grid,
                            alpha  = 0.05,
                            max_bp = 1e6,
                            cores  = 1) {

  ## candidates under a fully pooled FDR correction
  cand_pooled <- get_potential_outliers(
    map = map, ld_ws = ld_ws,
    th_ldw_grid = th_ldw_grid, p_cols = p_cols, alpha = alpha
  )

  ## candidates under a per-replicate FDR correction (more liberal,
  ## since m is ~10x smaller per group -- catches anything the
  ## pooled scan alone might miss)
  cand_per_sim <- unique(unlist(lapply(sort(unique(map$sim_id)), function(s) {
    get_potential_outliers(
      map = map[sim_id == s], ld_ws = ld_ws,
      th_ldw_grid = th_ldw_grid, p_cols = p_cols, alpha = alpha
    )
  })))

  all_candidates <- unique(c(cand_pooled, cand_per_sim))
  message(length(all_candidates), " candidate outlier markers found across both FDR modes")

  precompute_LD_edges(
    GTs    = GTs[, all_candidates, drop = FALSE],
    map    = map[marker %in% all_candidates],
    r2_min = 0.1,
    max_bp = max_bp,
    cores  = cores
  )
}

# el <- build_shared_el(map, GTs, ld_ws,
#                       p_cols = c(EMX = "emx_p"),
#                       th_ldw_grid = 1 - 10^seq(log10(1), log10(0.01), length.out = 20))

#----------------------------------------------------------
# 2) Run the (rho x th_ldw x r2_th x l_min) grid, with FDR
#    done either "per_sim" or "pooled"
#----------------------------------------------------------
run_grid_sim <- function(map, ld_ws, el,
                         fdr_mode = c("per_sim", "pooled"),
                         p_cols   = c(EMX = "emx_p"),
                         rho_grid, th_ldw_grid,
                         r2_grid, lmin_grid,
                         alpha = 0.05, bp_th = Inf, cores = 1) {

  fdr_mode   <- match.arg(fdr_mode)
  param_grid <- CJ(rho = rho_grid, th_ldw = th_ldw_grid)

  run_for_subset <- function(map_sub) {
    rbindlist(lapply(seq_len(nrow(param_grid)), function(i) {
      pars <- param_grid[i]
      run_one_grid(
        map       = map_sub,
        el        = el,
        ld_ws     = ld_ws,
        rho       = as.character(pars$rho),
        th_ldw    = pars$th_ldw,
        p_cols    = p_cols,
        alpha     = alpha,
        r2_grid   = r2_grid,
        lmin_grid = lmin_grid,
        bp_th     = bp_th,
        cores     = cores
      )
    }), fill = TRUE)
  }

  if (fdr_mode == "pooled") {
    out <- run_for_subset(map)
    out[, sim_id := "pooled"]
  } else {
    out <- rbindlist(lapply(sort(unique(map$sim_id)), function(s) {
      cbind(sim_id = s, run_for_subset(map[sim_id == s]))
    }), fill = TRUE)
  }

  out[, fdr_mode := fdr_mode]
  out[]
}

#----------------------------------------------------------
# 3) Classify each called outlier region as true/false and
#    compute OR-level recall / OR_FDR
#----------------------------------------------------------
#' @param grid_out output of run_grid_sim()
#' @param map      the (annotated) map, with true_pos and focal_QTN
#' @param method   name of the p_cols entry to score, e.g. "EMX"
classify_ORs <- function(grid_out, map, method = "EMX") {

  true_pos_vec <- setNames(map$true_pos,  map$marker)
  focal_vec    <- setNames(map$focal_QTN, map$marker)

  ## recall denominator: distinct causal loci (by focal_QTN) that
  ## have at least one true_pos-linked marker anywhere in the data
  n_true_loci_total <- length(unique(focal_vec[which(true_pos_vec)]))

  clusters_col <- grid_out[[method]]  # one element per grid row;
  # each element is a list of
  # marker-vectors (one per OR)

  res <- rbindlist(lapply(clusters_col, function(ors) {
    if (length(ors) == 0) {
      return(data.table(n_OR = 0L, n_OR_true = 0L, n_OR_false = 0L,
                        true_loci_recovered = 0L))
    }
    or_is_true <- vapply(ors, function(mk) any(true_pos_vec[mk], na.rm = TRUE), logical(1))
    true_hits  <- unlist(ors[or_is_true])
    true_hits  <- true_hits[true_pos_vec[true_hits]]
    n_loci_rec <- length(unique(focal_vec[true_hits]))

    data.table(
      n_OR      = length(ors),
      n_OR_true = sum(or_is_true),
      n_OR_false = sum(!or_is_true),
      true_loci_recovered = n_loci_rec
    )
  }))

  id_cols <- intersect(c("sim_id", "fdr_mode", "rho", "th_ldw", "r2_th", "l_min", "n_loci"),names(grid_out))
  out <- cbind(grid_out[, ..id_cols], res)
  out[, OR_FDR  := ifelse(n_OR > 0, n_OR_false / n_OR, NA_real_)]
  out[, recall  := true_loci_recovered / n_true_loci_total]
  out[]
}

#----------------------------------------------------------
# Example use
#----------------------------------------------------------
r2_grid     <- seq(0.6, 0.9, by = 0.1)
lmin_grid   <- c(1, 5, 10, 20)
rho_grid    <- colnames(ld_ws)          # as saved by precalculate_ld_w
th_ldw_grid <- 1 - 10^seq(log10(1), log10(0.01), length.out = 10)

grid_per_sim <- run_grid_sim(map, ld_ws, el, fdr_mode = "per_sim",
                              rho_grid = rho_grid, th_ldw_grid = th_ldw_grid,
                              r2_grid = r2_grid, lmin_grid = lmin_grid)

res_pooled  <- classify_ORs(grid_pooled,  map)
res_per_sim <- classify_ORs(grid_per_sim, map)

## e.g. recall vs OR_FDR across th_ldw at a fixed rho/r2_th/l_min:
res_per_sim[rho == "0.75" & r2_th == 0.7 & l_min == 5,
            .(th_ldw, recall, OR_FDR, n_OR)]

## or compare per_sim vs pooled FDR directly:
rbind(res_pooled, res_per_sim)[rho == "0.75" & r2_th == 0.7 & l_min == 5,
                               .(fdr_mode, sim_id, th_ldw, recall, OR_FDR)]
