## =====================================================================
## module_C2 / R/00_helpers.R
##
## Shared setup + primitives for the operating-grid stability exploration.
## Nothing here writes to the manuscript pipeline; every path is under module_C2/.
##
## Terminology used throughout (deliberately NOT "second-tier C-score"):
##   * cell        : one (tau_C, l_min) operating point on the grid G
##   * usable cell : a cell yielding >= 1 region with q_R < FDR
##   * anchor locus: a PRESPECIFIED region (member-marker vector) fixed at the
##                   primary operating point, independent of the grid search
##   * detection   : an anchor locus is called as a region in a cell
##   * significance: that called region has q_R < FDR in that cell
## =====================================================================
suppressMessages({ library(data.table); library(LDscnR) })

C2 <- new.env(parent = emptyenv())

C2$ROOT   <- "module_C2"
C2$RES    <- file.path(C2$ROOT, "results")
C2$FIG    <- file.path(C2$ROOT, "figures")
C2$CACHE  <- file.path(C2$ROOT, "cache")
C2$BUNDLE <- "module_sticklebacks_LDscnR/data/3sp_LDscnR_data.rds"
C2$NULLF  <- "module_sticklebacks_LDscnR/results/null_popperm_3sp.rds"
for (p in c(C2$RES, C2$FIG, C2$CACHE)) if (!dir.exists(p)) dir.create(p, recursive = TRUE)

## ---- fixed analysis constants (match the prototype except B) ---------
C2$TAUS    <- seq(0.02, 0.50, by = 0.02)
C2$LMINS   <- c(1L, 2L, 3L, 5L, 10L, 15L, 20L)
C2$RHO_LD  <- 0.60
C2$DCAP    <- 5e5
C2$FDR     <- 0.05
C2$OP_TAU  <- 0.05
C2$OP_LMIN <- 3L
C2$GAP     <- 1e4          # prototype's coordinate-merge tolerance
C2$ZISSOU  <- c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00")

## ---- data + null -----------------------------------------------------
## B is DETECTED, never capped: the prototype used BCAP = 100 while the bundle
## carries more surrogates. Using all of them lowers the p-value floor.
c2_load <- function(nullf = C2$NULLF, B = NULL) {
  d    <- readRDS(C2$BUNDLE)
  map  <- as.data.table(d$map)
  null <- readRDS(nullf)
  Bav  <- length(null$C_surr)
  Buse <- if (is.null(B)) Bav else min(B, Bav)
  edges <- ld_edges(null$universe, d$GTs, map[, .(marker, Chr, Pos)],
                    as.data.table(d$LD_decay$decay_sum),
                    rho_ld = C2$RHO_LD, dcap = C2$DCAP)
  list(map = map, edges = edges, C_obs = null$C_obs,
       surrs = null$C_surr[seq_len(Buse)], B = Buse, B_available = Bav,
       basis = if (is.null(null$basis)) "unknown" else null$basis,
       universe = null$universe,
       mpos = stats::setNames(map$Pos, map$marker),
       mchr = stats::setNames(as.integer(gsub("Chr", "", map$Chr)), map$marker))
}

## ---- cluster a C-vector at tau into a scored region table ------------
## `keep_markers = TRUE` also returns the member-marker vectors (needed for the
## anchor matching in Q3; the prototype kept coordinates only).
c2_cluster <- function(C, tau, D, keep_markers = FALSE) {
  mk <- names(C)[C >= tau]
  empty <- data.table(size = integer(), chr = integer(), lo = numeric(),
                      hi = numeric(), score = numeric(), maxC = numeric())
  if (!length(mk)) return(if (keep_markers) list(tab = empty, mk = list()) else empty)
  r <- ld_regions(mk, D$edges)
  tab <- rbindlist(lapply(r, function(x) data.table(
    size = length(x), chr = unname(D$mchr[x[1]]),
    lo = min(D$mpos[x]), hi = max(D$mpos[x]),
    score = sum(C[x]), maxC = max(C[x]))))
  if (keep_markers) list(tab = tab, mk = r) else tab
}

## ---- location-matched empirical p + BH q for one cell ----------------
## O    : observed region table for the cell (already size-filtered to l_min)
## SDT  : ONE data.table of all surrogate regions at this tau, with a `b` column
##        (surrogate index), already size-filtered to l_min
## B    : number of surrogates (so surrogates contributing no region still count)
## Returns emp_p and BH-adjusted q_R over the regions PRESENT IN THIS CELL --
## this is the per-cell hypothesis set whose changing size Q4 interrogates.
c2_emp_pq <- function(O, SDT, B) {
  if (!nrow(O)) return(data.table(emp_p = numeric(), q_R = numeric(),
                                  overlap_freq = numeric()))
  ep <- of <- numeric(nrow(O))
  for (i in seq_len(nrow(O))) {
    if (nrow(SDT)) {
      h <- SDT[chr == O$chr[i] & lo <= O$hi[i] & hi >= O$lo[i]]
      if (nrow(h)) {
        bb <- h[, .(m = max(score)), by = b]
        ep[i] <- (1 + sum(bb$m >= O$score[i])) / (1 + B)
        of[i] <- nrow(bb) / B
        next
      }
    }
    ep[i] <- 1 / (1 + B); of[i] <- 0
  }
  data.table(emp_p = ep, q_R = stats::p.adjust(ep, "BH"), overlap_freq = of)
}

## ---- anchor matching rules ------------------------------------------
## Given an anchor marker set `A` and a called region's markers `R`:
##   any     : |A n R| >= 1                      (most permissive)
##   overlap : |A n R| / min(|A|,|R|)   >= thr   (overlap coefficient)
##   recover : |A n R| / |A|            >= thr   (fraction of anchor recovered)
## Physical-span overlap is provided only as a COMPARATOR (it is what the
## prototype and the region-level p-value use), not as a recommended rule.
c2_match_score <- function(A, R, rule = c("any", "overlap", "recover")) {
  rule <- match.arg(rule)
  n <- length(intersect(A, R))
  if (!n) return(0)
  switch(rule,
         any     = 1,
         overlap = n / min(length(A), length(R)),
         recover = n / length(A))
}

## ---- small utilities -------------------------------------------------
c2_msg <- function(...) { cat(sprintf(...)); flush.console() }

## Spearman + top-k overlap between two named score vectors
c2_agree <- function(a, b, k = 5L) {
  ids <- intersect(names(a), names(b))
  ta <- names(sort(a[ids], decreasing = TRUE))[seq_len(min(k, length(ids)))]
  tb <- names(sort(b[ids], decreasing = TRUE))[seq_len(min(k, length(ids)))]
  list(spearman = suppressWarnings(stats::cor(a[ids], b[ids], method = "spearman")),
       top_k_overlap = length(intersect(ta, tb)) / min(k, length(ids)))
}

## ---- anchor loci: a PRESPECIFIED region set, fixed independently of the grid --
## The 17 EMMAX regions at the primary operating point (tau_C = 0.05, l_min = 3,
## rho_ld = 0.60), retained as MEMBER-MARKER VECTORS rather than coordinate spans.
c2_anchors <- function(core, tau = C2$OP_TAU, lmin = C2$OP_LMIN) {
  ## NB the primary operating point tau_C = 0.05 is NOT on the prototype's grid
  ## (seq(0.02, 0.50, 0.02) steps over it), so the anchors are clustered directly
  ## from C_obs + edges rather than looked up in the cached grid.
  ct <- core$by_tau[[as.character(tau)]]
  if (is.null(ct)) ct <- { x <- c2_cluster(core$D$C_obs, tau, core$D, keep_markers = TRUE)
                           list(obs = x$tab, mk = x$mk) }
  keep <- which(ct$obs$size >= lmin)
  mk  <- ct$mk[keep]
  ## data.tables restored from RDS lose their over-allocation, so `:=` on a
  ## subset errors on the internal selfref -- reallocate before modifying.
  tab <- data.table::setalloccol(data.table::copy(ct$obs[keep]))
  tab[, anchor := seq_along(keep)]
  tab[, label := sprintf("A%02d_Chr%d:%.2f-%.2fMb", anchor, chr, lo / 1e6, hi / 1e6)]
  tab[chr == 1 & lo < 21.93e6 & hi > 21.40e6, label := paste0(label, "*inv")]
  tab[chr == 4 & lo < 12.83e6 & hi > 12.79e6, label := paste0(label, "*Eda")]
  list(tab = tab, mk = stats::setNames(mk, tab$label))
}

## Add missing tau values to a cached core (e.g. the operating point 0.05, which
## the prototype grid steps over). Re-saves the cache.
c2_augment_core <- function(core, taus, path = file.path(C2$CACHE, "grid_core.rds")) {
  miss <- setdiff(as.character(taus), names(core$by_tau))
  if (!length(miss)) return(core)
  for (tc in miss) {
    tau <- as.numeric(tc)
    oc <- c2_cluster(core$D$C_obs, tau, core$D, keep_markers = TRUE)
    S  <- data.table::rbindlist(lapply(seq_along(core$D$surrs), function(b) {
            s <- c2_cluster(core$D$surrs[[b]], tau, core$D)
            if (nrow(s)) s[, b := b] else NULL }), fill = TRUE)
    if (!nrow(S)) S <- data.table(size = integer(), chr = integer(), lo = numeric(),
                                  hi = numeric(), score = numeric(), maxC = numeric(),
                                  b = integer())
    core$by_tau[[tc]] <- list(obs = oc$tab, mk = oc$mk, surr = S)
    c2_msg("[aug] added tau = %s (%d obs regions, %d surrogate regions)\n", tc, nrow(oc$tab), nrow(S))
  }
  saveRDS(core, path)
  core
}

## ---- match every anchor to the regions called in every cell ------------------
## For each (cell, anchor) records whether the anchor is DETECTED (some called
## region matches it under `rule`) and, if so, the matched region's statistics.
## `rule`: "any" | "overlap" | "recover" | "span" (span = coordinate overlap, the
## prototype's implicit rule, included only as a comparator).
## When several called regions match, the one with the largest match score is
## taken (ties broken by C-mass), and `n_match` records the ambiguity.
c2_anchor_grid <- function(core, anchors, B, rule = "any", thr = 0.5,
                           taus = C2$TAUS, lmins = C2$LMINS) {
  A <- anchors$mk; alab <- names(A)
  out <- list()
  for (tau in taus) {
    ct <- core$by_tau[[as.character(tau)]]
    Sb <- ct$surr[b <= B]
    ## match anchors to this tau's called regions ONCE (l_min only filters later)
    nO <- nrow(ct$obs)
    ms <- matrix(0, nrow = length(A), ncol = max(nO, 1L))
    if (nO) for (j in seq_len(nO)) {
      Rj <- ct$mk[[j]]
      for (i in seq_along(A)) {
        ms[i, j] <- if (rule == "span") {
          as.numeric(ct$obs$chr[j] == anchors$tab$chr[i] &&
                     ct$obs$lo[j] <= anchors$tab$hi[i] &&
                     ct$obs$hi[j] >= anchors$tab$lo[i])
        } else c2_match_score(A[[i]], Rj, rule = if (rule == "any") "any" else rule)
      }
    }
    hit <- if (rule %in% c("any", "span")) ms > 0 else ms >= thr
    for (lm in lmins) {
      keep <- which(ct$obs$size >= lm)
      if (!length(keep)) {
        out[[length(out) + 1L]] <- data.table(
          tau = tau, lmin = lm, label = alab, detected = FALSE, n_match = 0L,
          match_score = 0, reg = NA_integer_, size = NA_integer_, score = NA_real_,
          maxC = NA_real_, emp_p = NA_real_, q_R = NA_real_, n_tested = 0L)
        next
      }
      O  <- ct$obs[keep]
      pq <- c2_emp_pq(O, Sb[size >= lm], B)
      sub <- ms[, keep, drop = FALSE]; sh <- hit[, keep, drop = FALSE]
      best <- vapply(seq_along(A), function(i) {
        w <- which(sh[i, ])
        if (!length(w)) return(NA_integer_)
        w[which.max(sub[i, w] + 1e-9 * O$score[w])]
      }, integer(1))
      out[[length(out) + 1L]] <- data.table(
        tau = tau, lmin = lm, label = alab,
        detected = !is.na(best), n_match = rowSums(sh),
        match_score = vapply(seq_along(A), function(i)
          if (is.na(best[i])) 0 else sub[i, best[i]], numeric(1)),
        reg = keep[best],
        size = O$size[best], score = O$score[best], maxC = O$maxC[best],
        emp_p = pq$emp_p[best], q_R = pq$q_R[best], n_tested = nrow(O))
    }
  }
  res <- rbindlist(out)
  res[, sig := detected & !is.na(q_R) & q_R < C2$FDR]
  res[]
}

## =====================================================================
## Null-admissible operating grid (second iteration)
## =====================================================================

## The GRID is seq(0.02, 0.50, 0.02) x LMINS = 175 cells. tau_C = 0.05 is a
## REFERENCE coordinate appended to the cache for anchor construction -- it is NOT
## a grid member and must never enter an admissibility denominator.
C2$REF_TAU  <- 0.05
C2$is_grid_tau <- function(x) x %in% C2$TAUS

## ---- per-cell null behaviour ----------------------------------------
## Absent surrogate indices mean "this surrogate produced no region", so every
## statistic is taken over all B surrogates, not over the b values present.
c2_null_cell <- function(core, tau, lmin, B) {
  S <- core$by_tau[[as.character(tau)]]$surr
  S <- S[size >= lmin]
  nreg <- integer(B); nmk <- integer(B)
  if (nrow(S)) {
    a <- S[, .(n = .N, m = sum(size)), by = b]
    nreg[a$b] <- a$n; nmk[a$b] <- a$m
  }
  k <- sum(nreg > 0)
  ci <- stats::binom.test(k, B)$conf.int          # Clopper-Pearson (exact)
  data.table(tau = tau, lmin = lmin,
             p_null_any = k / B, p_null_lo = ci[1], p_null_hi = ci[2], k_null = k,
             mean_null_regions = mean(nreg),
             q50_null_regions = stats::quantile(nreg, 0.50, names = FALSE),
             q90_null_regions = stats::quantile(nreg, 0.90, names = FALSE),
             q95_null_regions = stats::quantile(nreg, 0.95, names = FALSE),
             q99_null_regions = stats::quantile(nreg, 0.99, names = FALSE),
             max_null_regions = max(nreg),
             mean_null_markers = mean(nmk),
             mean_null_coverage = mean(nmk) / length(core$D$universe))
}

## ---- markers covered by observed regions in a cell -------------------
c2_obs_markers <- function(core, tau, lmin) {
  ct <- core$by_tau[[as.character(tau)]]
  keep <- which(ct$obs$size >= lmin)
  if (!length(keep)) return(character(0))
  unique(unlist(ct$mk[keep], use.names = FALSE))
}

## ---- anchor detection over an arbitrary cell set ---------------------
## Detection ONLY (significance is handled separately, and is tested for
## redundancy rather than assumed). `thr` is the anchor-marker retention
## fraction |A n R| / |A|; `recip` additionally demands |A n R| / |R| >= thr.
c2_detect_grid <- function(core, anchors, cells, thr = 0.5, recip = FALSE) {
  A <- anchors$mk; alab <- names(A)
  out <- list()
  for (tau in unique(cells$tau)) {
    ct <- core$by_tau[[as.character(tau)]]
    lms <- cells[tau == ..tau]$lmin
    nO <- nrow(ct$obs)
    rec <- rcp <- matrix(0, length(A), max(nO, 1L))
    if (nO) for (j in seq_len(nO)) {
      Rj <- ct$mk[[j]]
      for (i in seq_along(A)) {
        n <- length(intersect(A[[i]], Rj))
        if (n) { rec[i, j] <- n / length(A[[i]]); rcp[i, j] <- n / length(Rj) }
      }
    }
    for (lm in lms) {
      keep <- which(ct$obs$size >= lm)
      if (!length(keep)) {
        out[[length(out) + 1L]] <- data.table(tau = tau, lmin = lm, label = alab,
                                              detected = FALSE, n_match = 0L, best_rec = 0)
        next
      }
      sr <- rec[, keep, drop = FALSE]; sp <- rcp[, keep, drop = FALSE]
      hit <- sr >= thr & (if (recip) sp >= thr else TRUE)
      out[[length(out) + 1L]] <- data.table(
        tau = tau, lmin = lm, label = alab,
        detected = rowSums(hit) > 0, n_match = rowSums(hit),
        best_rec = vapply(seq_along(A), function(i) max(sr[i, ]), numeric(1)))
    }
  }
  rbindlist(out)
}
