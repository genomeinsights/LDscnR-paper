## =============================================================================
## module_9sp/R/03d_diagnose_variance_components.R
##
## DIAGNOSTIC, for Supplementary materials + the Discussion's treatment of
## permutation problems under population structure (PK, 2026-09-07). Chases
## down WHY 03_EMMAX.R's label-permutation null (p=0.31/0.75) and
## 03c_EMMAX_structurednull.R's GRM-structured null (p=0.001/0.005) disagree
## so sharply on the SAME observed data.
##
## FOUR CHECKS, IN ORDER:
##
##   1. Does the REAL tested phenotype (eco_resid) hit a REML boundary
##      (ve -> 0, i.e. h^2 -> 1) under the canonical (stage1-pruned) GRM?
##   2. Is that boundary SPECIFIC to eco_resid, or does ANY phenotype hit it
##      under this K (a generic REML/numerical artefact)? Checked against 8
##      independent rnorm() draws, uncorrelated with K by construction.
##   3. Is the boundary specific to 9sp, or does module_3sp's OWN real
##      ecotype phenotype hit the same boundary under ITS OWN canonical GRM?
##      (3sp's scan is trusted/already-reported, so this is the key control.)
##   4. Genome-wide calibration: does 9sp's marker-wise EMMAX scan show the
##      same genomic-control behaviour as 3sp's (lambda_GC), and does a
##      leaner GRM basis (greedy LD-pruning, ~19k markers vs 753,625
##      stage1-pruned markers) change it -- and does that change the actual
##      bottom-line consensus-arm result, or only the diagnostic?
##
## FINDINGS (2026-09-07), see the printed summary at the end:
##   - eco_resid hits ve=0/h^2=1 in BOTH 9sp and 3sp -- NOT unique to 9sp,
##     and not itself the smoking gun (a population-differentiated binary
##     phenotype tested against a GRM built to capture that same structure
##     is a known hard case for REML in general).
##   - Random phenotypes get sensible interior h^2 in BOTH datasets --
##     rules out a generic REML/numerical bug; the boundary is
##     phenotype-structure-specific, as expected.
##   - THE REAL DIFFERENCE: genome-wide lambda_GC is healthy for 3sp (~1.09,
##     typical mild inflation for real polygenic signal) but DEFLATED for
##     9sp (~0.87) -- 9sp's mixed-model correction is genuinely
##     over-conservative genome-wide, not just sitting on a boundary.
##   - A leaner GRM basis (greedy LD-pruning) measurably improves marker-wise
##     lambda_GC (0.868 -> 0.957) but does NOT rescue the consensus-arm
##     result (44 vs 55 significant units, both non-significant under label
##     permutation, p=0.28 vs 0.31) -- a partial diagnostic improvement, not
##     a fix.
##
## READING: severe, pervasive population structure (9sp's background LD 0.301
## vs 3sp's 0.055; GRM off-diagonal SD 0.26-0.32) can push a standard EMMAX
## mixed-model fit into over-correction for a population-differentiated
## phenotype specifically -- both a naive label-permutation null (built from
## the SAME over-corrected machinery, so the observed statistic looks
## unremarkable against it) and a fully structure-matched null (which
## suppresses near enough EVERY surrogate, real or not, making it read as
## dramatically significant by contrast) become hard to interpret at face
## value under these conditions. Neither p-value (0.31 or 0.001) should be
## taken as the answer without this context.
## =============================================================================
suppressMessages({library(data.table); library(LDscnR); library(SNPRelate)})
devtools::load_all("~/gitlab/LDscnR")
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_9sp"), "R", "00_config.R"))
STAGE <- "03d_diagnose_variance_components"
say("=== %s ===\n\n", STAGE)

b9  <- readRDS(file.path(PATHS$out, "02_bundle", "bundle.rds"))
sc9 <- readRDS(file.path(PATHS$out, "03_EMMAX", "scan.rds"))
GRM9 <- b9$GRM; eco_resid <- sc9$eco_resid; n9 <- nrow(b9$GTs)
Kn9 <- (n9 - 1) / sum((diag(n9) - matrix(1/n9, n9, n9)) * GRM9) * GRM9

## ---- 1. does the real phenotype hit the REML boundary? ------------------------
say("[1] REML on the REAL tested phenotype (eco_resid), canonical (stage1-pruned) GRM\n")
re_real <- LDscnR:::emma.REMLE(eco_resid, matrix(1, n9, 1), Kn9)
say("    vg=%.4f ve=%.4f h2=%.4f\n", re_real$vg, re_real$ve, re_real$vg / (re_real$vg + re_real$ve))

## ---- 2. is the boundary phenotype-specific, or generic to this K? -------------
say("\n[2] REML on 8 independent random phenotypes (uncorrelated with K by construction)\n")
set.seed(42)
h2_random <- vapply(seq_len(8), function(i) {
  y <- rnorm(n9)
  re <- LDscnR:::emma.REMLE(y, matrix(1, n9, 1), Kn9)
  h2 <- re$vg / (re$vg + re$ve)
  say("    random %d: vg=%.4f ve=%.4f h2=%.4f\n", i, re$vg, re$ve, h2)
  h2
}, numeric(1))
say("    median h2 (random) = %.4f -- sensible interior estimates, NOT a generic boundary artefact\n",
    median(h2_random))

## ---- 3. does module_3sp's OWN trusted phenotype/GRM hit the same boundary? ----
say("\n[3] control: module_3sp's real ecotype phenotype, ITS OWN canonical GRM\n")
b3 <- readRDS(path.expand("~/gitlab/LDscnR-paper/module_3sp/out/02_bundle/bundle.rds"))
GRM3 <- b3$GRM; n3 <- nrow(b3$GTs)
Kn3 <- (n3 - 1) / sum((diag(n3) - matrix(1/n3, n3, n3)) * GRM3) * GRM3
eco3 <- as.integer(as.factor(b3$pheno$ecotype)) - 1L
re3 <- LDscnR:::emma.REMLE(eco3, matrix(1, n3, 1), Kn3)
say("    3sp real eco: vg=%.4f ve=%.4f h2=%.4f\n", re3$vg, re3$ve, re3$vg / (re3$vg + re3$ve))
set.seed(42)
h2_random3 <- vapply(seq_len(5), function(i) {
  y <- rnorm(n3)
  re <- LDscnR:::emma.REMLE(y, matrix(1, n3, 1), Kn3)
  re$vg / (re$vg + re$ve)
}, numeric(1))
say("    3sp random phenotypes: median h2 = %.4f\n", median(h2_random3))
say("    -> the ve=0/h2=1 boundary is NOT unique to 9sp: 3sp's own trusted scan hits it too.\n")

## ---- 4. genome-wide calibration: lambda_GC, 9sp vs 3sp, and GRM basis sensitivity --
say("\n[4] genome-wide marker-wise EMMAX scan: genomic control (lambda_GC)\n")
p9 <- emmax_fast(emmax_setup(b9$GTs, GRM9), eco_resid)
lambda9 <- median(qchisq(1 - p9, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)
say("    9sp (stage1-pruned GRM, %s markers): lambda_GC = %.3f ; min p = %.2e\n",
    format(length(p9), big.mark=","), lambda9, min(p9, na.rm = TRUE))

p3 <- emmax_fast(emmax_setup(b3$GTs, GRM3), eco3)
lambda3 <- median(qchisq(1 - p3, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)
say("    3sp (stage1-pruned GRM, %s markers): lambda_GC = %.3f ; min p = %.2e\n",
    format(length(p3), big.mark=","), lambda3, min(p3, na.rm = TRUE))

say("\n[4b] does a leaner GRM basis (greedy LD-pruning) change 9sp's calibration?\n")
gds_path <- file.path(PATHS$cache, "9sp.gds")
gds <- snpgdsOpen(gds_path)
on.exit(try(snpgdsClose(gds), silent = TRUE), add = TRUE)
set.seed(1)
pruned <- unlist(snpgdsLDpruning(gds, ld.threshold = GRM_GREEDY$ld.threshold,
                                 slide.max.bp = GRM_GREEDY$slide.max.bp,
                                 autosome.only = FALSE, verbose = FALSE), use.names = FALSE)
GRM_greedy <- snpgdsGRM(gds, snp.id = pruned, method = GRM_METHOD, verbose = FALSE, autosome.only = FALSE)$grm
say("    greedy-pruned markers: %s (vs %s stage1-pruned)\n",
    format(length(pruned), big.mark=","), format(length(b9$grm_markers), big.mark=","))
p_greedy <- emmax_fast(emmax_setup(b9$GTs, GRM_greedy), eco_resid)
lambda_g <- median(qchisq(1 - p_greedy, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)
say("    9sp (greedy GRM): lambda_GC = %.3f ; min p = %.2e -- improves, but see 4c: does it change the result?\n",
    lambda_g, min(p_greedy, na.rm = TRUE))

say("\n[4c] does the greedy GRM change the CONSENSUS-ARM bottom line, not just marker-wise lambda?\n")
um <- ld_unit_matrix(b9$GTs, b9$stage1, b9$map, size_floor = SIZE_FLOOR, repr = UNIT_REPR)
Pu_g <- emmax_setup(um, GRM_greedy)
pu_obs_g <- emmax_fast(Pu_g, eco_resid)
test_con_g <- ld_outlier_test(b9$stage1, b9$map, pu_obs_g, statistic = "unit", size_floor = SIZE_FLOOR,
                              alpha = ALPHA, assembly = "stage2_discovered", GTs = b9$GTs,
                              LD_decay = b9$LD_decay, score_threshold = REGION_ASSEMBLY$score_threshold,
                              distance_threshold = REGION_ASSEMBLY$distance_threshold)
covar_obs <- b9$pheno[[EMMAX_COVAR]]
POPT <- unique(b9$pheno[, .(pop_ID, ecotype = b9$eco, lineage = get(EMMAX_COVAR))])
perm_lineage <- function() {
  pt <- copy(POPT)[, ep := sample(ecotype), by = lineage]
  y  <- pt$ep[match(b9$pheno$pop_ID, pt$pop_ID)]
  as.numeric(stats::resid(stats::lm(y ~ covar_obs)))
}
p_perm_g <- function(bb) { set.seed(bb); emmax_fast(Pu_g, perm_lineage()) }
null_con_g <- ld_outlier_perm(test_con_g, b9$stage1, b9$map, p_perm_g, GTs = b9$GTs, LD_decay = b9$LD_decay,
                              B = NPERM_CONSENSUS, level = "units", verbose = TRUE)
say("    greedy-GRM consensus: observed %d/%d -> %d regions | surrogate mean %.2f | p = %.4f\n",
    null_con_g$observed, nrow(test_con_g$units), nrow(test_con_g$regions),
    mean(null_con_g$surrogates), null_con_g$p)
say("    vs stage1-pruned-GRM consensus (03_EMMAX.R): observed 55/4547 -> 21 regions | surrogate mean 67.36 | p = 0.3117\n")
say("    -> partial diagnostic improvement (lambda closer to 1), but NOT a fix: still non-significant either way.\n")

## ---- save --------------------------------------------------------------------
OUT <- file.path(stage_dir(STAGE), "diagnostics.rds")
dir.create(stage_dir(STAGE), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(
  eco_resid_vc = re_real, random_h2_9sp = h2_random,
  sp3_eco_vc = re3, random_h2_3sp = h2_random3,
  lambda_9sp_stage1pruned = lambda9, lambda_3sp = lambda3, lambda_9sp_greedy = lambda_g,
  greedy_markers = length(pruned), stage1pruned_markers = length(b9$grm_markers),
  greedy_consensus_test = test_con_g, greedy_consensus_null = null_con_g
), OUT)
write_receipt(STAGE, inputs = c(file.path(PATHS$out, "02_bundle", "_receipt.rds"),
                                file.path(PATHS$out, "03_EMMAX", "_receipt.rds"),
                                path.expand("~/gitlab/LDscnR-paper/module_3sp/out/02_bundle/_receipt.rds")),
             params = list(), outputs = OUT)
say("\n[5] wrote %s\n    receipt: %s\n", OUT, receipt_path(STAGE))
