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
dir_data <- file.path(mod, "data")        # caches / scored results (git-ignored)
dir_fig  <- file.path(mod, "figures")     # plots (git-ignored)
for (d in c(dir_data, dir_fig)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
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

## per-QTN LD/distance tolerances for TP matching, from the pooled decay fit.
## rho_r2=0.75 / rho_d=0.95 are decay-relative; dmax_cap bounds the distance (in
## low gene flow the decay is so slow d_from_rho blows up to Mb-scale).
score_thresholds <- function(decay_sum, rho_r2 = 0.75, rho_d = 0.95, dmax_cap = Inf) {
  ds <- as.data.table(decay_sum)
  list(r2min = ld_from_rho(median(ds$b), median(ds$c), rho_r2),
       dmax  = min(d_from_rho(median(ds$a), rho_d), dmax_cap))
}

## Precision-recall AUC (trapezoidal) for an ordered-threshold sweep: given the
## (recall, precision) operating points a single ordered knob traces (tau_C, or
## alpha), integrate precision over recall from 0. Ties in recall collapse to max
## precision (upper envelope); handles the mild non-monotonicity from re-clustering
## at each threshold. This is the STANDARD PR-AUC that the C-score enables (one
## ordered knob) -- replaces the random-search AUC-PR* used for the raw unordered grid.
pr_auc <- function(recall, precision) {
  ok <- is.finite(recall) & is.finite(precision)
  if (sum(ok) < 1L) return(NA_real_)
  d <- data.table::data.table(r = recall[ok], p = precision[ok])[, .(p = max(p)), by = r]
  data.table::setorder(d, r)
  r <- c(0, d$r); p <- c(d$p[1], d$p)                 # integrate from recall 0
  sum(diff(r) * (utils::head(p, -1) + utils::tail(p, -1)) / 2)
}

## Fast EMMAX for permutation/surrogate nulls: eigendecompose the (normalised)
## kinship and rotate the genotypes ONCE (fast_emmax_setup), then a whitened per-SNP
## F-test per phenotype (fast_emmax_p). Numerically identical to emmax() (cor 1.0,
## max |diff -log10 p| ~1e-9); ~25x faster per phenotype. Reuses the package's REML
## internals. Works for sim (per-file) and empirical (one genome) alike.
fast_emmax_setup <- function(GTs, K) {
  n <- nrow(GTs); one_n <- matrix(1, n, n) / n
  Kn <- (n - 1) / sum((diag(n) - one_n) * K) * K; Xo <- matrix(1, n, 1)
  ev <- eigen(Kn, symmetric = TRUE)
  list(Kn = Kn, Xo = Xo, eigR = LDscnR:::emma.eigen.R.wo.Z(Kn, Xo), V = ev$vectors, lam = ev$values,
       Xt = crossprod(ev$vectors, GTs), xot = as.numeric(crossprod(ev$vectors, Xo)), df2 = n - 2L)
}
fast_emmax_p <- function(pre, y) {
  re <- LDscnR:::emma.REMLE(y, pre$Xo, pre$Kn, eig.R = pre$eigR)
  wv <- 1 / sqrt(re$vg * pre$lam + re$ve); yt <- as.numeric(crossprod(pre$V, y)) * wv
  xo <- pre$xot * wv; Xtw <- pre$Xt * wv; a2 <- sum(xo * xo)
  ra <- yt - xo * (sum(xo * yt) / a2); RSS0 <- sum(ra * ra)
  Rb <- Xtw - outer(xo, as.numeric(crossprod(xo, Xtw)) / a2)
  RSSf <- RSS0 - as.numeric(crossprod(ra, Rb))^2 / colSums(Rb * Rb)
  stats::pf((RSS0 / RSSf - 1) * pre$df2, 1, pre$df2, lower.tail = FALSE)
}

## Dedup-NEUTRAL scoring: like evaluate_ORs() but a region matching an ALREADY-
## claimed true-positive QTN (a dedup loser -- an "extra" region for the same QTN)
## is counted as NEITHER TP nor FP, it is dropped. This removes the clustering-
## parameter sensitivity of performance: over-fragmenting a true region into N
## pieces still yields 1 TP and 0 extra FP, so there is no incentive to tune
## r2_link / distance to game the metric. Genuine FP (regions matching no true-
## positive QTN) still count; over-merging is still penalised via recall (one
## region can claim only one QTN). Uses diagnose_OR_classification() to identify
## the dedup losers. Returns the same fields as evaluate_ORs().
evaluate_ORs_qtn <- function(clusters, map, qtn_ld_table, r2_min_focal, d_max_focal) {
  n_true <- as.data.table(map)[true_pos_QTN == TRUE, .N]
  if (!length(clusters))
    return(list(TP = 0L, FP = 0L, FN = n_true, extras = 0L,
                Precision = NA_real_, Recall = 0, PR = 0))
  dg <- diagnose_OR_classification(clusters, map, qtn_ld_table, r2_min_focal, d_max_focal)
  extra <- dg$dropped_by_dedup == TRUE & dg$candidate_qtn_is_true_positive %in% TRUE
  TP <- sum(dg$is_TP)                       # distinct true-positive QTN (one kept per QTN)
  FP <- sum(!dg$is_TP & !extra)             # regions matching no true-positive QTN
  FN <- n_true - TP
  Precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  Recall    <- if (n_true > 0)   TP / n_true else NA_real_
  PR        <- if (!is.na(Precision) && !is.na(Recall)) Precision * Recall else NA_real_
  list(TP = TP, FP = FP, FN = FN, extras = sum(extra),
       Precision = Precision, Recall = Recall, PR = PR)
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
## Thresholds are decay-relative (rho) and self-calibrate per chromosome when
## `decay_sum` + `rho_ld`/`rho_d` are supplied: r2_link = ld_from_rho(b, c, rho_ld),
## distance = min(d_from_rho(a_pred, rho_d), dcap_max). The dcap_max cap is
## essential in low gene flow (slow decay -> d_from_rho ~Mb -> would merge a whole
## chromosome into one cluster). Falls back to the scalar r2_link/dcap otherwise.
cluster_regions <- function(markers, map, GTs, r2_link = 0.5, dcap = Inf,
                            decay_sum = NULL, rho_ld = NULL, rho_d = NULL, dcap_max = Inf) {
  m <- as.data.table(map)[marker %in% markers]
  r2_by <- d_by <- NULL
  if (!is.null(decay_sum) && (!is.null(rho_ld) || !is.null(rho_d))) {
    ds <- as.data.table(decay_sum); cc <- if ("c" %in% names(ds)) ds$c else rep(1, nrow(ds))
    ac <- if ("a_pred" %in% names(ds)) ds$a_pred else ds$a
    if (!is.null(rho_ld)) r2_by <- stats::setNames(ld_from_rho(ds$b, cc, rho_ld), ds$Chr)
    if (!is.null(rho_d))  d_by  <- stats::setNames(pmin(d_from_rho(ac, rho_d), dcap_max), ds$Chr)
  }
  out <- list()
  for (ch in unique(m$Chr)) {
    mk <- m[Chr == ch, marker]; pos <- m[Chr == ch, Pos]
    r2c <- if (!is.null(r2_by) && !is.na(r2_by[ch])) r2_by[[ch]] else r2_link
    dc  <- if (!is.null(d_by)  && !is.na(d_by[ch]))  d_by[[ch]]  else dcap
    if (length(mk) == 1L) { out <- c(out, list(mk)); next }
    R <- suppressWarnings(stats::cor(GTs[, mk], use = "pairwise.complete.obs")^2)
    R[!is.finite(R)] <- 0
    cl <- stats::cutree(stats::hclust(stats::as.dist(1 - R), method = "single"), h = 1 - r2c)
    cl <- .split_gap(cl, pos, dc)
    out <- c(out, split(mk, cl))
  }
  unname(out)
}
