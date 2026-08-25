## =====================================================================
## module_sticklebacks_LDscnR / gate_background.R
##
## TASK 1 of the null rework (framework §4, "One gate, before any p-value is read").
##
## The gate: for each engine x basis, the MEDIAN PER-SURROGATE count of markers
## with C > 0 and of LD-regions of size >= l_min, against the observed counts. A
## basis whose surrogates approach the observed counts is mis-specified and its
## p-values are not interpretable.
##
## Previously the two engines were reported on different statistics (EMMAX as
## median peaks/surrogate, LFMM as POOLED n(C>0) across surrogates), so the
## numbers in the manuscript were not comparable. This puts every cell on one
## statistic, via LDscnR::ld_gate() -- the package function, not a local
## reimplementation, so the paper and the released pipeline cannot diverge.
##
## CANONICAL SET = the LEGACY EMMAX nulls (null_{uncapped,popperm,regionperm,
## spatial,latent}_3sp.rds). This is the set the manuscript's numbers come from:
## its regional-permutation null gives max q_R = 0.0149, which is the value in
## framework §4 and in the handoff. The newer engine-tagged set
## (null_emmax_*_3sp.rds) gives 0.0348 on the same statistic. Both are tabulated
## below, `canonical` marks which is which, and the gate VERDICTS agree across
## the two -- only the exact counts differ.
##
## CAVEAT to carry into the text: the legacy set does not share one B. Genetic is
## B = 100 and latent B = 50, while the rest are B = 200. The statistic is now
## uniform; the surrogate count behind it is not, and a median over 50 draws is a
## noisier estimate than one over 200.
##
## Run from the LDscnR-paper root:
##   Rscript module_sticklebacks_LDscnR/gate_background.R
## Writes results/gate_background.csv
## =====================================================================
suppressMessages({ library(data.table); library(LDscnR) })

RES <- "module_sticklebacks_LDscnR/results"
BND <- "module_sticklebacks_LDscnR/data/3sp_LDscnR_data.rds"
TAU <- 0.05; LMIN <- 3L; RHO_LD <- 0.60; DCAP <- 5e5

## ---- the (engine x basis) cells --------------------------------------------------
## `set` distinguishes the two EMMAX null generations; `canonical` marks the one
## the manuscript reports.
cells <- rbindlist(list(
  data.table(set="legacy",   canonical=TRUE,  engine="EMMAX", basis="genetic (MVN on K)",   file="null_uncapped_3sp.rds"),
  data.table(set="legacy",   canonical=TRUE,  engine="EMMAX", basis="global permutation",   file="null_popperm_3sp.rds"),
  data.table(set="legacy",   canonical=TRUE,  engine="EMMAX", basis="regional permutation", file="null_regionperm_3sp.rds"),
  data.table(set="legacy",   canonical=TRUE,  engine="EMMAX", basis="spatial (GP kernel)",  file="null_spatial_3sp.rds"),
  data.table(set="legacy",   canonical=TRUE,  engine="EMMAX", basis="latent (PC subspace)", file="null_latent_3sp.rds"),
  data.table(set="emmax_v2", canonical=FALSE, engine="EMMAX", basis="genetic (MVN on K)",   file="null_emmax_genetic_3sp.rds"),
  data.table(set="emmax_v2", canonical=FALSE, engine="EMMAX", basis="global permutation",   file="null_emmax_global_perm_3sp.rds"),
  data.table(set="emmax_v2", canonical=FALSE, engine="EMMAX", basis="regional permutation", file="null_emmax_region_perm_3sp.rds"),
  data.table(set="emmax_v2", canonical=FALSE, engine="EMMAX", basis="spatial (GP kernel)",  file="null_emmax_spatial_3sp.rds"),
  data.table(set="emmax_v2", canonical=FALSE, engine="EMMAX", basis="latent (PC subspace)", file="null_emmax_latent_3sp.rds"),
  data.table(set="lfmm",     canonical=TRUE,  engine="LFMM",  basis="regional permutation", file="null_lfmm_region_perm_3sp.rds")
))

d <- readRDS(BND); map <- as.data.table(d$map)

nulls <- lapply(seq_len(nrow(cells)), function(i) {
  x <- readRDS(file.path(RES, cells$file[i])); class(x) <- "ld_null"
  x$engine <- cells$engine[i]; x$basis <- cells$basis[i]; x
})
names(nulls) <- sprintf("%s | %s | %s", cells$engine, cells$basis, cells$set)

## ---- the LFMM four-null bundle is B = 1 per basis ---------------------------------
## Not a median in any meaningful sense: a single draw, reported so the cells are
## visible rather than silently absent. The B column in the output is the warning.
b1 <- readRDS(file.path(RES, "lfmm_b1_fournulls_Cs.rds"))
b1_map <- c(genetic = "genetic (MVN on K)", global_perm = "global permutation",
            regional_perm = "regional permutation", spatial = "spatial (GP kernel)")
for (s in names(b1_map)) {
  C <- b1[[s]]
  nm <- sprintf("LFMM | %s | lfmm_b1", b1_map[[s]])
  nulls[[nm]] <- structure(list(C_obs = b1$observed, C_surr = list(C[C > 0]),
                                universe = unique(c(names(b1$observed)[b1$observed > 0], names(C)[C > 0])),
                                basis = b1_map[[s]], engine = "LFMM", B = 1L),
                           class = "ld_null")
  cells <- rbind(cells, data.table(set = "lfmm_b1", canonical = FALSE, engine = "LFMM",
                                   basis = b1_map[[s]], file = "lfmm_b1_fournulls_Cs.rds"))
}

## ---- one shared edge graph -------------------------------------------------------
## Edge membership is a property of a marker pair, not of the marker set, and
## ld_regions() subsets to the markers it is handed -- so a graph over the union
## of every universe gives clustering identical to per-basis graphs, and makes
## the rows comparable by construction.
univ <- unique(unlist(lapply(nulls, `[[`, "universe"), use.names = FALSE))
cat(sprintf("[1] shared edge graph over %d markers\n", length(univ))); flush.console()
t0 <- Sys.time()
edges <- ld_edges(univ, d$GTs, map[, .(marker, Chr, Pos)],
                  as.data.table(d$LD_decay$decay_sum), rho_ld = RHO_LD, dcap = DCAP)
cat(sprintf("    built in %.1f s\n", as.numeric(Sys.time() - t0, units = "secs"))); flush.console()

## ---- the gate --------------------------------------------------------------------
g <- ld_gate(nulls, edges, tau = TAU, l_min = LMIN)      # warns on any basis that fails
out <- cbind(cells[, .(set, canonical)], as.data.table(g))
setnames(out, "basis", "cell")
setcolorder(out, c("canonical", "set", "engine", "cell", "B",
                   "obs_regions", "med_regions", "ratio", "pass"))
fwrite(out, file.path(RES, "gate_background.csv"))

cat(sprintf("\n=== GATE (tau_C = %.2f, l_min = %d, rho_ld = %.2f, d_cap = %.0e) ===\n",
            TAU, LMIN, RHO_LD, DCAP))
cat("\n--- CANONICAL (legacy EMMAX + LFMM) ---\n")
print(out[canonical == TRUE, .(engine, cell, B, obs_Cgt0, med_Cgt0, obs_regions,
                               med_regions, max_regions, frac_surr_any_region, ratio, pass)])
cat("\n--- comparison: the newer engine-tagged EMMAX set, and the LFMM B=1 bundle ---\n")
print(out[canonical == FALSE, .(set, engine, cell, B, obs_Cgt0, med_Cgt0, obs_regions,
                                med_regions, max_regions, ratio, pass)])
cat(sprintf("\n[2] wrote %s\n", file.path(RES, "gate_background.csv")))
