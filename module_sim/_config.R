## module_sim/_config.R
## Shared configuration + the pooling helper for the simulated-data TP/FP
## validation of ld_outlier_clusters().  source() this at the top of every 0x
## script.  Run everything from LDscnR-paper/.
##
## Data (Nemo forward sims, one file per replicate chromosome-pair):
##   list(GTs[160 x ~30k], map, env, LD_decay, ld_ws) with PRECOMPUTED
##   emx_p/emx_F (EMMAX, GRM on LD-pruned) + lfmm_p/lfmm_F (LFMM K=5, GC) and
##   ground truth (type, Va, MAF, max_LD_with_QTN, ...).
##   V = INVERSE selection intensity (V1 strongest, V2 collapses);
##   c = dispersal (c1 some gene flow, c2 very limited -> structure-dominated).
## Two layouts are supported transparently:
##   - external full grid  /Volumes/Nemo/Nemo_sim/parsed_sim_data2
##     files adapt_bgs_chr{1..10}_V{}_c{}_env{}.rds ; ld_ws cols "rho_0.95"
##   - local V2 subset      parsed_sim_data
##     files adapt_bgs_chr{1..10}_V2_c{}_env{}.rds ; ld_ws cols "0.95"

suppressMessages({ library(LDscnR); library(data.table); library(parallel) })

mod      <- "/Users/petrikem/gitlab/LDscnR-paper/module_sim"
SIM_DATA <- Sys.getenv("LDSCNR_SIM_DATA",
                       unset = "/Volumes/Nemo/Nemo_sim/parsed_sim_data2")
if (!dir.exists(SIM_DATA)) SIM_DATA <- "parsed_sim_data"   # fall back to local V2 subset

## TP/FP scoring engine (LD+distance match of a cluster to a Va-qualified QTN)
source("R/define_ORs_functions.R")
source("R/Outlier_regions_simulation.R")

## Pool the 10 replicates that share (V, c, env) into one 20-chromosome genome:
## relabel Chr/marker R{i}_*, cbind genotypes, concatenate ld_w and decay rows.
## Same env => the replicate chromosomes are ~independent, which stabilises RMSC
## and yields 10 QTN-bearing + 10 neutral chromosomes as a TP/FP bed.
pool_group <- function(pattern, data_dir = SIM_DATA) {
  files <- list.files(data_dir, pattern = pattern, full.names = TRUE)
  if (!length(files)) stop("no files match ", pattern, " in ", data_dir)
  reps  <- as.integer(sub(".*_chr([0-9]+)_.*", "\\1", basename(files)))
  files <- files[order(reps)]
  maps <- gts <- decs <- vector("list", length(files)); ldw <- c()
  for (i in seq_along(files)) {
    d <- readRDS(files[i]); m <- as.data.table(d$map)
    m[, Chr := paste0("R", i, "_", Chr)][, marker := paste0("R", i, "_", marker)]
    G <- d$GTs; colnames(G) <- paste0("R", i, "_", colnames(G))
    lwc <- if ("rho_0.95" %in% colnames(d$ld_ws)) "rho_0.95" else "0.95"
    lw <- d$ld_ws[, lwc]; names(lw) <- m$marker
    ds <- as.data.table(d$LD_decay$decay_sum); ds[, Chr := paste0("R", i, "_", Chr)]
    maps[[i]] <- m; gts[[i]] <- G; ldw <- c(ldw, lw); decs[[i]] <- ds
  }
  list(map   = rbindlist(maps, fill = TRUE),
       GTs   = do.call(cbind, gts),
       ldw   = ldw,
       decay = list(decay_sum = rbindlist(decs, fill = TRUE)))
}

## per-QTN LD/distance tolerances for TP matching, from the pooled decay fit
score_thresholds <- function(decay_sum, rho_r2 = 0.75, rho_d = 0.95) {
  ds <- as.data.table(decay_sum)
  list(r2min = ld_from_rho(median(ds$b), median(ds$c), rho_r2),
       dmax  = d_from_rho(median(ds$a), rho_d))
}

## Method-AGNOSTIC second-level clustering: partition an arbitrary marker set
## (e.g. the UNION of every method's outlier SNPs) into regions with the SAME
## distance-restricted single-linkage rule ld_outlier_clusters() uses internally
## -- single-linkage on 1 - r^2 within chromosome (join at r2_link), then split
## each block wherever consecutive members are more than `dcap` apart. Returns a
## list of member-marker vectors (one per region). Because it is built once from
## the pooled outlier set, a single-SNP hit and an ld_w hit at the same locus
## land in the SAME region and are scored as one detection.
.split_gap <- function(cl, pos, dcap) {
  if (!is.finite(dcap)) return(cl)
  out <- integer(length(cl)); nxt <- 0L
  for (g in unique(cl)) {
    gi <- which(cl == g); o <- order(pos[gi])
    run <- integer(length(gi)); run[o] <- cumsum(c(TRUE, diff(pos[gi][o]) > dcap))
    out[gi] <- nxt + run; nxt <- nxt + max(run)
  }
  out
}
cluster_regions <- function(markers, map, GTs, r2_link = 0.5, dcap = Inf) {
  m <- as.data.table(map)[marker %in% markers]
  out <- list()
  for (ch in unique(m$Chr)) {
    mk <- m[Chr == ch, marker]; pos <- m[Chr == ch, Pos]
    if (length(mk) == 1L) { out <- c(out, list(mk)); next }
    R <- suppressWarnings(stats::cor(GTs[, mk], use = "pairwise.complete.obs")^2)
    R[!is.finite(R)] <- 0
    cl <- stats::cutree(stats::hclust(stats::as.dist(1 - R), method = "single"), h = 1 - r2_link)
    cl <- .split_gap(cl, pos, dcap)
    out <- c(out, split(mk, cl))
  }
  unname(out)
}
