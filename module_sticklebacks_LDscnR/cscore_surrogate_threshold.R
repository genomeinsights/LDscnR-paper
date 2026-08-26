## =====================================================================
## module_sticklebacks_LDscnR / cscore_surrogate_threshold.R
##
## EXPERIMENT: calibrate each C-score cell's rejection threshold from the
## SURROGATES rather than from the observed data's own BH.
##
## CURRENT (LDscnR::ld_cscore). For cell (rho, q*): candidates are the markers
## whose ld_w reaches the q* quantile; BH runs on the OBSERVED p-values among
## them; hits are q < alpha. Each surrogate is scored the same way but
## SEPARATELY, on its own p-values. Two consequences:
##   * BH assumes null p-values are uniform. Under residual structure they are
##     not, so the effective cutoff is off by an unknown amount.
##   * observed and surrogates are each scored at their OWN cutoff, so the
##     comparison is not like-for-like even before the region test.
##
## PROPOSED (this script). For cell (rho, q*) pick ONE threshold t from the
## pooled surrogate p-values among that cell's candidates, and apply it to
## observed and every surrogate alike. t comes from an EMPIRICAL-NULL FDR --
## BH's criterion with the null estimated rather than assumed:
##
##     O(t) = #{observed candidates with p <= t}
##     S(t) = (1/B) * #{pooled surrogate candidate p-values <= t}
##     t*   = max{ t : S(t) <= alpha * O(t) },  evaluated at the observed p_(i)
##
## The plain alpha-quantile of the surrogate p-values would instead make every
## surrogate marker a hit in ~alpha of cells BY CONSTRUCTION, so C_surr would sit
## near alpha everywhere and the sparsity that tau_C and the gate depend on would
## be gone. The FDR form keeps stringency.
##
## RESULT (2026-08-25, 3sp EMMAX regional permutation, B = 100): REJECTED.
## The criterion is satisfiable in only 4 of 400 cells, so C_obs is capped at
## 4/400 = 0.01 and NOTHING reaches tau_C = 0.05 -- zero regions called, against
## 426 markers at C >= 0.05 under the current scheme.
##
## It is not an implementation bug. The regional-permutation surrogates are
## severely non-uniform in the deep tail, where BH assumes uniformity:
##
##     threshold   observed n   surrogate n/B   inflation vs uniform
##     1e-5        80           25.7            3.3x
##     1e-6        13            5.3            6.7x
##     1e-7         2            1.2           15.6x
##     1e-8         0            0.4           49x
##
## At 1e-7 the observed data has 2 markers and the surrogates average 1.2; at
## 1e-8 observed has 0 and surrogates 0.4. GENOME-WIDE, the observed extreme tail
## is not more extreme than a structure-preserving permutation of it, so a correct
## empirical-null FDR rejects nothing -- and at i = 1 the criterion fails in 25 of
## 25 sampled cells.
##
## What it demonstrates is why the LOCATION-MATCHED test is necessary. Pooling
## surrogate p-values genome-wide reintroduces, one stage earlier, exactly the
## pooled count-FDR that framework section 4 removed. Both hold of the same null:
## pooled genome-wide, nothing survives; location-matched, the same regional
## permutation puts all 17 EMMAX regions at q_R <= 0.0149. Excess null peaks are
## scattered and rarely land on any particular observed locus.
##
## Untested caveat: this is specific to the regional-permutation basis, which
## preserves structure by construction. The genetic MVN null is far quieter
## (median 0 markers with C > 0 per surrogate vs 5) and might behave differently.
## Its surrogate p-values were not kept, only C-scores, so it could not be tried.
##
## PK settled the wider question on 2026-08-26 (HANDOFF_gc_decision.md): do not
## build an empirical null from pooled surrogate statistics, and do not replace
## the within-cell BH. This file stays a PROTOTYPE and must not be wired into the
## pipeline. It is kept because the measurement above is worth not repeating.
##
## This changes the DEFINITION of the C-score, which is Fang et al. (2021)'s
## (doi:10.1093/molbev/msab144) -- a more substantive departure than the
## continuous-vs-menu integration the package already documents. An experiment
## until shown to help.
##
## Run from the LDscnR-paper root:
##   Rscript module_sticklebacks_LDscnR/cscore_surrogate_threshold.R [B]
## =====================================================================
suppressMessages({ library(data.table); library(LDscnR) })

BND  <- "module_sticklebacks_LDscnR/data/3sp_LDscnR_data.rds"
PERM <- "emmax_perm_reginal.rds"
OUT  <- "/private/tmp/claude-539526166/-Users-petrikem-gitlab-LDscnR/f5c2953d-f16b-4266-bda5-08c843e9b161/scratchpad"
a <- commandArgs(trailingOnly = TRUE)
B <- if (length(a) >= 1) as.integer(a[1]) else 100L
ALPHA <- 0.05; QSTAR <- seq(0, 0.95, by = 0.05)
GCAP <- 1e-3     # bound on the search; the chosen t is checked against it

d <- readRDS(BND); map <- as.data.table(d$map); ld_ws <- d$ld_ws
p_obs <- map$emx_p
perm <- readRDS(PERM); perm <- perm[seq_len(min(B, length(perm)))]; B <- length(perm)
stopifnot(length(p_obs) == nrow(ld_ws), all(lengths(perm) == nrow(ld_ws)))
nmk <- nrow(ld_ws); rho_cols <- colnames(ld_ws); ncell <- length(rho_cols) * length(QSTAR)
cat(sprintf("[1] %d markers, %d rho x %d qstar = %d cells, B = %d\n",
            nmk, length(rho_cols), length(QSTAR), ncell, B)); flush.console()

## ---- sparse store of the only surrogate values that can be a threshold ------
t0 <- Sys.time()
mkl <- pl <- bl <- vector("list", B)
for (b in seq_len(B)) { w <- which(perm[[b]] < GCAP)
  mkl[[b]] <- w; pl[[b]] <- perm[[b]][w]; bl[[b]] <- rep.int(b, length(w)) }
S_mk <- unlist(mkl, use.names = FALSE); S_p <- unlist(pl, use.names = FALSE)
S_b  <- unlist(bl,  use.names = FALSE); rm(mkl, pl, bl); invisible(gc())
cat(sprintf("[2] surrogate p < %.0e: %d of %.3g (%.4f%%) in %.0f s\n", GCAP, length(S_p),
            as.numeric(nmk) * B, 100 * length(S_p) / (as.numeric(nmk) * B),
            as.numeric(Sys.time() - t0, units = "secs"))); flush.console()

## ---- proposed scheme -------------------------------------------------------
t0 <- Sys.time()
cnt_obs <- integer(nmk)
hit_mk <- hit_b <- vector("list", ncell)   # sparse (marker, surrogate) accumulation
k <- 0L; thr_used <- numeric(0); near_cap <- 0L; empty <- 0L
incand <- logical(nmk)

for (rc in rho_cols) {
  lw <- ld_ws[, rc]
  for (q in QSTAR) {
    k <- k + 1L
    cand <- which(lw >= stats::quantile(lw, q, na.rm = TRUE))
    if (!length(cand)) { empty <- empty + 1L; next }
    incand[] <- FALSE; incand[cand] <- TRUE

    pc <- p_obs[cand]; po <- sort(pc[!is.na(pc) & pc < GCAP])
    if (!length(po)) { empty <- empty + 1L; next }
    keep <- incand[S_mk]
    ps <- sort(S_p[keep])

    Scount <- findInterval(po, ps)                       # surrogate values <= p_(i)
    okk <- which(Scount <= ALPHA * B * seq_along(po))    # S(t)/B <= alpha * O(t)
    if (!length(okk)) { empty <- empty + 1L; next }
    t_star <- po[max(okk)]
    thr_used <- c(thr_used, t_star)
    if (t_star >= 0.5 * GCAP) near_cap <- near_cap + 1L

    h <- cand[which(pc <= t_star)]
    if (length(h)) cnt_obs[h] <- cnt_obs[h] + 1L
    sel <- which(keep & S_p <= t_star)
    if (length(sel)) { hit_mk[[k]] <- S_mk[sel]; hit_b[[k]] <- S_b[sel] }
  }
  cat(sprintf("   rho %s\n", rc)); flush.console()
}
C_obs_new <- cnt_obs / ncell

hm <- unlist(hit_mk, use.names = FALSE); hb <- unlist(hit_b, use.names = FALSE)
rm(hit_mk, hit_b); invisible(gc())
C_surr_new <- vector("list", B)
if (length(hm)) {
  key <- data.table(mk = hm, b = hb)[, .N, by = c("mk", "b")]
  ## split by surrogate index rather than filtering inside [ -- `b` is both a
  ## column name and the loop variable, which data.table's NSE cannot resolve
  sp <- split(seq_len(nrow(key)), key$b)
  for (bb in seq_len(B)) {
    ix <- sp[[as.character(bb)]]
    C_surr_new[[bb]] <- if (length(ix)) stats::setNames(key$N[ix] / ncell, map$marker[key$mk[ix]]) else numeric(0)
  }
} else for (b in seq_len(B)) C_surr_new[[b]] <- numeric(0)
names(C_obs_new) <- map$marker
cat(sprintf("[3] proposed scheme done in %.1f min | cells used %d/%d (%d empty) | t*: median %.3g range %.3g-%.3g | near GCAP %d\n",
            as.numeric(Sys.time() - t0, units = "mins"), ncell - empty, ncell, empty,
            stats::median(thr_used), min(thr_used), max(thr_used), near_cap)); flush.console()

## ---- current scheme, on the identical p-values ------------------------------
t0 <- Sys.time()
C_obs_cur <- ld_cscore(p_obs, ld_ws, alpha = ALPHA, qstar = QSTAR)
C_surr_cur <- lapply(seq_len(B), function(b) { C <- ld_cscore(perm[[b]], ld_ws, alpha = ALPHA, qstar = QSTAR); C[C > 0] })
cat(sprintf("[4] current scheme done in %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))

mk_null <- function(Co, Cs, basis) structure(list(
  C_obs = Co, C_surr = Cs,
  universe = unique(c(names(Co)[Co > 0], unlist(lapply(Cs, names), use.names = FALSE))),
  basis = basis, engine = "EMMAX", B = length(Cs)), class = "ld_null")
res <- list(current  = mk_null(C_obs_cur, C_surr_cur, "region_perm (BH on observed)"),
            proposed = mk_null(C_obs_new, C_surr_new, "region_perm (surrogate-calibrated t)"),
            thr_used = thr_used, ncell = ncell, alpha = ALPHA, B = B)
saveRDS(res, file.path(OUT, "cscore_surrthr.rds"))

cat("\n=== C-score comparison ===\n")
for (nm in c("current", "proposed")) { x <- res[[nm]]
  cs <- vapply(x$C_surr, function(s) sum(s > 0), numeric(1))
  cat(sprintf("  %-9s obs C>0 = %6d | obs C>=0.05 = %5d | surrogate C>0: median %6.0f, max %6.0f | universe %d\n",
              nm, sum(x$C_obs > 0), sum(x$C_obs >= 0.05), stats::median(cs), max(cs), length(x$universe))) }
cat("\n[done]\n")
