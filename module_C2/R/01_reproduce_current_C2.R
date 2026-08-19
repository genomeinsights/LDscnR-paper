## =====================================================================
## module_C2 / R/01_reproduce_current_C2.R    [Question 1]
##
## (a) Build ONE shared core cache over the (tau_C, l_min) grid, keeping the
##     things the prototype discarded: observed region MEMBER MARKERS, per-region
##     emp_p / q_R in every cell, and the per-cell hypothesis count.
## (b) Reproduce null_sig_landscape.R's score EXACTLY (B = 100 cap, coordinate
##     consensus merge) and check it against the committed prototype output.
## (c) Re-run the same definition at the FULL available B and report the shift.
##
## Rscript module_C2/R/01_reproduce_current_C2.R
## =====================================================================
source("module_C2/R/00_helpers.R")

CORE <- file.path(C2$CACHE, "grid_core.rds")

## ---------------------------------------------------------------------
## 1. core grid pass (cached)
## ---------------------------------------------------------------------
if (file.exists(CORE)) {
  core <- readRDS(CORE); D <- core$D
  c2_msg("[1] loaded cached grid core (B=%d)\n", core$B)
} else {
  D <- c2_load()                                   # B detected, never capped
  c2_msg("[1] %s null: B = %d of %d available ; universe = %d ; %d cells\n",
         D$basis, D$B, D$B_available, length(D$universe),
         length(C2$TAUS) * length(C2$LMINS))
  by_tau <- list()
  for (tau in C2$TAUS) {
    oc <- c2_cluster(D$C_obs, tau, D, keep_markers = TRUE)
    S  <- rbindlist(lapply(seq_along(D$surrs), function(b) {
            s <- c2_cluster(D$surrs[[b]], tau, D)
            if (nrow(s)) s[, b := b] else NULL }), fill = TRUE)
    if (!nrow(S)) S <- data.table(size = integer(), chr = integer(), lo = numeric(),
                                  hi = numeric(), score = numeric(), maxC = numeric(),
                                  b = integer())
    by_tau[[as.character(tau)]] <- list(obs = oc$tab, mk = oc$mk, surr = S)
    c2_msg("    tau=%.2f : %d obs regions, %d surrogate regions\n", tau, nrow(oc$tab), nrow(S))
  }
  core <- list(D = D, by_tau = by_tau, B = D$B)
  saveRDS(core, CORE)
  c2_msg("[1] wrote %s\n", CORE)
}
B_FULL <- core$B

## ---------------------------------------------------------------------
## 2. per-cell scoring at an arbitrary B  ->  long per-region table
## ---------------------------------------------------------------------
## Returns one row per (cell, observed region), carrying everything Q4 needs.
score_grid <- function(core, B) {
  out <- list()
  for (tau in C2$TAUS) {
    ct <- core$by_tau[[as.character(tau)]]
    Sb <- ct$surr[b <= B]
    for (lm in C2$LMINS) {
      keep <- which(ct$obs$size >= lm)
      if (!length(keep)) {
        out[[length(out) + 1L]] <- data.table(
          tau = tau, lmin = lm, reg = integer(), size = integer(), chr = integer(),
          lo = numeric(), hi = numeric(), score = numeric(), maxC = numeric(),
          emp_p = numeric(), q_R = numeric(), overlap_freq = numeric(), n_tested = 0L)
        next
      }
      O  <- ct$obs[keep]
      pq <- c2_emp_pq(O, Sb[size >= lm], B)
      out[[length(out) + 1L]] <- data.table(
        tau = tau, lmin = lm, reg = keep, O[, .(size, chr, lo, hi, score, maxC)],
        pq, n_tested = nrow(O))
    }
  }
  rbindlist(out)
}

## full |G|: cells producing NO observed region must still appear (with n_obs = 0),
## otherwise the grid silently shrinks to the cells that happen to be productive.
cells_from <- function(L) {
  full <- CJ(tau = C2$TAUS, lmin = C2$LMINS)
  got  <- L[!is.na(emp_p), .(n_obs = .N, n_sig = sum(q_R < C2$FDR)), by = .(tau, lmin)]
  out  <- merge(full, got, by = c("tau", "lmin"), all.x = TRUE)
  out[is.na(n_obs), `:=`(n_obs = 0L, n_sig = 0L)][]
}

## ---- 2a. the prototype's exact configuration: B = 100 ---------------
L100 <- score_grid(core, 100L)
res100 <- cells_from(L100)
NCELL <- length(C2$TAUS) * length(C2$LMINS)
stopifnot(nrow(res100) == NCELL)
NUSE100 <- res100[n_sig > 0, .N]
c2_msg("\n[2a] B=100 (prototype config): %d cells, %d usable (>=1 significant region)\n",
       NCELL, NUSE100)

## ---- 2b. full B ------------------------------------------------------
Lfull <- score_grid(core, B_FULL)
resF  <- cells_from(Lfull)
NUSEF <- resF[n_sig > 0, .N]
c2_msg("[2b] B=%d (all available)  : %d cells, %d usable\n", B_FULL, NCELL, NUSEF)
NPROD <- resF[n_obs > 0, .N]
c2_msg("[2c] cells producing >=1 OBSERVED region: %d ; of those, usable (>=1 significant): %d\n",
       NPROD, NUSEF)
c2_msg("     -> every productive cell is usable: %s. The 'usable' restriction therefore\n",
       NPROD == NUSEF)
c2_msg("        removes only the %d cells where the method finds NOTHING, i.e. |U| carries\n", NCELL - NUSEF)
c2_msg("        no information about null danger on this dataset.\n")
saveRDS(list(L100 = L100, Lfull = Lfull, res100 = res100, resF = resF, B_full = B_FULL),
        file.path(C2$CACHE, "grid_scored.rds"))

## ---------------------------------------------------------------------
## 3. the prototype's consensus merge -- reproduced verbatim
## ---------------------------------------------------------------------
## NB the prototype computes a coordinate-consensus cluster id `cl` and then
## aggregates `by = chr`, so `cl` is never used: its "loci" are CHROMOSOMES.
## Both aggregations are reproduced here so the discrepancy is explicit.
consensus <- function(L, NUSE, by_cl) {
  sig <- L[q_R < C2$FDR]
  if (!nrow(sig)) return(data.table())
  setorder(sig, chr, lo)
  sig[, cl := { o <- order(lo); s <- lo[o]; e <- hi[o]
                g <- c(TRUE, s[-1] > cummax(e)[-length(e)] + C2$GAP)
                cumsum(g)[order(o)] }, by = chr]
  key <- if (by_cl) c("chr", "cl") else "chr"
  loci <- sig[, .(lo = min(lo), hi = max(hi), n_cells_sig = uniqueN(paste(tau, lmin)),
                  best_q = min(q_R), max_score = round(max(score), 2),
                  lmin_max = max(lmin)), by = key]
  loci[, stability := round(n_cells_sig / NUSE, 3)]
  loci[, locus := sprintf("Chr%d:%.2f-%.2f", chr, lo / 1e6, hi / 1e6)]
  setorder(loci, -stability, best_q)
  loci[]
}

proto_chr <- consensus(L100, NUSE100, by_cl = FALSE)   # what the script actually does
proto_cl  <- consensus(L100, NUSE100, by_cl = TRUE)    # what its comment says it does
full_chr  <- consensus(Lfull, NUSEF,  by_cl = FALSE)

## ---- check against the committed prototype output -------------------
ref_f <- "module_sticklebacks_LDscnR/results/region_stability2.csv"
if (file.exists(ref_f)) {
  ref <- fread(ref_f)
  cmp <- merge(ref[, .(chr, ref_ncells = n_cells_sig, ref_stab = stability)],
               proto_chr[, .(chr, new_ncells = n_cells_sig, new_stab = stability)],
               by = "chr", all = TRUE)
  cmp[, d_ncells := new_ncells - ref_ncells]
  c2_msg("\n[3] reproduction vs committed region_stability2.csv (B=100, by=chr):\n")
  c2_msg("    loci: ref %d / new %d ; usable cells: ref 134 / new %d\n",
         nrow(ref), nrow(proto_chr), NUSE100)
  c2_msg("    max |delta n_cells_sig| = %d ; exact match on all loci: %s\n",
         max(abs(cmp$d_ncells), na.rm = TRUE), all(cmp$d_ncells == 0, na.rm = TRUE))
  if (any(cmp$d_ncells != 0, na.rm = TRUE)) print(cmp[d_ncells != 0])
}

c2_msg("\n[3] prototype aggregation `by = chr` -> %d 'loci' (i.e. chromosomes)\n", nrow(proto_chr))
c2_msg("    intended  `by = chr, cl` (10 kb consensus merge) -> %d loci\n", nrow(proto_cl))
c2_msg("    => the committed score is a per-CHROMOSOME statistic, not a per-locus one.\n")
print(proto_chr[, .(locus, stability, n_cells_sig, best_q = round(best_q, 4), max_score)])

## ---------------------------------------------------------------------
## 4. required reporting
## ---------------------------------------------------------------------
c2_msg("\n[4] Monte Carlo resolution: B=%d -> p floor = 1/(B+1) = %.5f ; B=100 -> %.5f\n",
       B_FULL, 1 / (B_FULL + 1), 1 / 101)
c2_msg("    regions with emp_p at the floor (full B, over all cells): %d of %d rows\n",
       Lfull[emp_p <= 1 / (B_FULL + 1) + 1e-12, .N], nrow(Lfull[!is.na(emp_p)]))

## can a single-cell locus score high when the usable grid is small?
one <- proto_chr[n_cells_sig == 1]
c2_msg("    single-usable-cell check: min attainable non-zero stability = 1/%d = %.4f\n",
       NUSE100, 1 / NUSE100)
if (nrow(one)) c2_msg("    loci significant in exactly 1 cell: %s (stability %.3f)\n",
                      paste(one$locus, collapse = ", "), one$stability[1])

## relationship of the score to size / strength / best q
ann <- L100[q_R < C2$FDR, .(max_size = max(size), max_maxC = max(maxC),
                            max_score = max(score), best_q = min(q_R)), by = chr]
rep1 <- merge(proto_chr[, .(chr, locus, stability, n_cells_sig)], ann, by = "chr")
setorder(rep1, -stability)
fwrite(rep1, file.path(C2$RES, "current_C2_reproduction.csv"))
cr <- function(a, b) round(stats::cor(a, b, method = "spearman"), 3)
c2_msg("\n[4] Spearman(stability, .) : span=%s  maxC=%s  C-mass=%s  -log10(best_q)=%s\n",
       cr(rep1$stability, rep1$max_size), cr(rep1$stability, rep1$max_maxC),
       cr(rep1$stability, rep1$max_score), cr(rep1$stability, -log10(rep1$best_q)))

fwrite(res100[, .(tau, lmin, n_obs, n_sig)], file.path(C2$RES, "grid_cells_B100.csv"))
fwrite(resF[,   .(tau, lmin, n_obs, n_sig)], file.path(C2$RES, "grid_cells_Bfull.csv"))
c2_msg("\n[4] wrote results/current_C2_reproduction.csv + grid_cells_B{100,full}.csv\n")
