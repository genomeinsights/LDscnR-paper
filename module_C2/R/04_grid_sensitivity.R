## =====================================================================
## module_C2 / R/04_grid_sensitivity.R    [Question 5]
##
## Is the stability score invariant to grid design? Compare, on the fixed anchor
## loci, three prespecified grids:
##   G1 original    tau = seq(0.02,0.50,0.02) x l_min = 1,2,3,5,10,15,20   (175)
##   G2 coarse      tau = 0.05,0.10,0.20,0.30,0.50 x l_min = 1,2,3,5,10     (25)
##   G3 breakpoint  tau at distinct C-values where the observed region set actually
##                  CHANGES (so a long stretch of identical calls is one
##                  observation, not fifteen) x l_min = 1,2,3,5,10
## plus a jackknife over tau-axis halves, and a dense/sparse tau contrast.
##
## Rscript module_C2/R/04_grid_sensitivity.R
## =====================================================================
source("module_C2/R/00_helpers.R")
core <- readRDS(file.path(C2$CACHE, "grid_core.rds")); B <- core$B
anc  <- c2_anchors(core)
ag   <- readRDS(file.path(C2$CACHE, "anchor_grid.rds"))
gfull <- ag$G$rec50                                  # already computed on G1

score_on <- function(g, taus, lmins) {
  sub <- g[tau %in% taus & lmin %in% lmins]
  n <- length(taus) * length(lmins)
  sub[, .(S_G = sum(sig) / n, D = sum(detected) / n), by = label]
}

## ---- G1 original -----------------------------------------------------
G1 <- score_on(gfull, C2$TAUS, C2$LMINS)

## ---- G2 coarse -------------------------------------------------------
TAU2 <- c(0.05, 0.10, 0.20, 0.30, 0.50); LMIN2 <- c(1L, 2L, 3L, 5L, 10L)
G2 <- score_on(gfull, TAU2, LMIN2)

## ---- G3 breakpoint-oriented -----------------------------------------
## keep only tau values at which the observed region SET changes (region count or
## the multiset of sizes), so identical repeated calls are not counted as
## independent robustness observations.
sig_of_tau <- vapply(C2$TAUS, function(t) {
  o <- core$by_tau[[as.character(t)]]$obs
  paste(nrow(o), paste(sort(o$size), collapse = ","), sep = "|") }, character(1))
TAU3 <- C2$TAUS[!duplicated(sig_of_tau)]
c2_msg("[1] breakpoint tau axis: %d of %d tau values change the region set: %s\n",
       length(TAU3), length(C2$TAUS), paste(sprintf("%.2f", TAU3), collapse = " "))
G3 <- score_on(gfull, TAU3, LMIN2)

## ---- jackknife over the tau axis ------------------------------------
TL <- C2$TAUS[C2$TAUS <= 0.26]; TH <- C2$TAUS[C2$TAUS > 0.26]
GL <- score_on(gfull, TL, C2$LMINS); GH <- score_on(gfull, TH, C2$LMINS)
## dense vs sparse sampling of the SAME tau range
GD <- score_on(gfull, C2$TAUS, C2$LMINS)
GS <- score_on(gfull, C2$TAUS[seq(1, length(C2$TAUS), by = 3)], C2$LMINS)

GR <- list(G1_original = G1, G2_coarse = G2, G3_breakpoint = G3,
           G4_tau_low = GL, G5_tau_high = GH, G6_tau_sparse = GS)
tab <- Reduce(function(a, b) merge(a, b, by = "label"),
              lapply(names(GR), function(n)
                setnames(GR[[n]][, .(label, S_G)], "S_G", n)))
setorder(tab, -G1_original)
fwrite(tab, file.path(C2$RES, "grid_sensitivity.csv"))
c2_msg("\n[2] S_G under six grids (fixed anchors):\n")
print(tab[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])

## ---- agreement -------------------------------------------------------
vv <- function(n) stats::setNames(GR[[n]]$S_G, GR[[n]]$label)
c2_msg("\n[3] agreement with the original grid:\n")
agr <- rbindlist(lapply(setdiff(names(GR), "G1_original"), function(n) {
  a <- c2_agree(vv("G1_original"), vv(n), k = 5L)
  data.table(grid = n, spearman = round(a$spearman, 3), top5_overlap = round(a$top_k_overlap, 2),
             pearson = round(stats::cor(vv("G1_original"), vv(n)[GR$G1_original$label]), 3))
}))
print(agr)
fwrite(agr, file.path(C2$RES, "grid_sensitivity_agreement.csv"))

## ---- the named loci --------------------------------------------------
key <- tab[grepl("inv|Eda", label) | label %in% tab$label[anc$tab$size <= 4]]
c2_msg("\n[4] ranks of the named / small loci under each grid (1 = most stable):\n")
rk <- copy(tab); for (n in names(GR)) rk[[n]] <- frank(-tab[[n]], ties.method = "min")
sel <- rk[grepl("inv|Eda", label) | label %in% anc$tab$label[anc$tab$size <= 4]]
print(sel)
c2_msg("\n[4] (small = anchors with <= 4 SNPs: %d of %d)\n",
       sum(anc$tab$size <= 4), nrow(anc$tab))
c2_msg("[4] wrote results/grid_sensitivity.csv + grid_sensitivity_agreement.csv\n")
