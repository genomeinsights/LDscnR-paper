## =====================================================================
## module_C2 / R/03_anchor_loci.R    [Question 3 + Question 4]
##
## Fix the loci INDEPENDENTLY of the grid search: the 17 EMMAX regions at the
## primary operating point (tau_C = 0.05, l_min = 3, rho_ld = 0.60), retained as
## member-marker vectors. Then match called regions to anchors in every cell by
## MARKER MEMBERSHIP (not physical span) under four rules, and record the full
## per-cell detail the BH diagnostics in Q4 need.
##
## Rscript module_C2/R/03_anchor_loci.R
## =====================================================================
source("module_C2/R/00_helpers.R")
core <- readRDS(file.path(C2$CACHE, "grid_core.rds")); B <- core$B
## the primary operating point is not on the prototype's tau grid -- add it so the
## coarse/breakpoint grids in 04 can actually include it.
core <- c2_augment_core(core, C2$OP_TAU)
NCELL <- length(C2$TAUS) * length(C2$LMINS)

## ---- 1. anchors ------------------------------------------------------
anc <- c2_anchors(core)
c2_msg("[1] anchor set: %d regions at tau=%.2f l_min=%d rho_ld=%.2f\n",
       nrow(anc$tab), C2$OP_TAU, C2$OP_LMIN, C2$RHO_LD)

## cross-check against the committed primary region table
ref_f <- "module_sticklebacks_LDscnR/results/regions_tau0.05_lmin3_rho0.60.csv"
if (file.exists(ref_f)) {
  ref <- fread(ref_f)[method == "EMMAX"]
  ref[, chr := as.integer(sub("Chr", "", Chr))]
  setorder(ref, chr, start)
  new <- anc$tab[order(chr, lo)]
  ok <- nrow(ref) == nrow(new) && all(ref$start == new$lo) && all(ref$end == new$hi) &&
        all(ref$n_snp == new$size)
  c2_msg("[1] matches committed EMMAX region set (%d rows): %s\n", nrow(ref), ok)
  if (!ok) { print(ref[, .(Chr, start, end, n_snp)]); print(new[, .(chr, lo, hi, size)]) }
}
c2_msg("    sizes: %s\n", paste(sort(anc$tab$size), collapse = " "))

## ---- 2. matching rules ----------------------------------------------
RULES <- list(
  any     = list(rule = "any",     thr = 0),
  ovl50   = list(rule = "overlap", thr = 0.50),
  rec50   = list(rule = "recover", thr = 0.50),
  span    = list(rule = "span",    thr = 0))       # comparator: coordinate overlap
G <- list()
for (nm in names(RULES)) {
  t0 <- Sys.time()
  G[[nm]] <- c2_anchor_grid(core, anc, B, rule = RULES[[nm]]$rule, thr = RULES[[nm]]$thr)
  c2_msg("[2] rule '%s': %.0f s ; detected cells = %d ; significant cells = %d\n",
         nm, as.numeric(Sys.time() - t0, "secs"), sum(G[[nm]]$detected), sum(G[[nm]]$sig))
}
saveRDS(list(anchors = anc, G = G, B = B), file.path(C2$CACHE, "anchor_grid.rds"))

## ---- 3. per-rule D / Q scores + rank agreement ----------------------
sc <- function(g) g[, .(D = sum(detected) / NCELL, Q = sum(sig) / NCELL,
                        mean_nmatch = mean(n_match[detected]),
                        ambig = mean(n_match[detected] > 1)), by = label]
S <- lapply(G, sc)
tabQ <- Reduce(function(a, b) merge(a, b, by = "label"),
               lapply(names(S), function(n) setnames(S[[n]][, .(label, Q)], "Q", paste0("Q_", n))))
tabD <- Reduce(function(a, b) merge(a, b, by = "label"),
               lapply(names(S), function(n) setnames(S[[n]][, .(label, D)], "D", paste0("D_", n))))
cmp <- merge(tabD, tabQ, by = "label")
setorder(cmp, -Q_any)
fwrite(cmp, file.path(C2$RES, "anchor_matching_rules.csv"))
c2_msg("\n[3] anchor scores by matching rule (D = detection, Q = significance):\n")
print(cmp)

c2_msg("\n[3] rank agreement between rules (Spearman on Q; top-5 overlap):\n")
for (a in names(S)) for (b in names(S)) if (a < b) {
  va <- stats::setNames(S[[a]]$Q, S[[a]]$label); vb <- stats::setNames(S[[b]]$Q, S[[b]]$label)
  ag <- c2_agree(va, vb, k = 5L)
  c2_msg("    %-6s vs %-6s : rho = %+.3f ; top5 overlap = %.2f\n", a, b, ag$spearman, ag$top_k_overlap)
}
c2_msg("\n[3] chaining/ambiguity (fraction of detected cells where >1 called region matches):\n")
print(S$any[order(-ambig), .(label, mean_nmatch = round(mean_nmatch, 2), ambig = round(ambig, 3))])

## ---- 3b. THE headline check: does significance add anything to detection? ----
c2_msg("\n[3b] detection vs significance, per rule (cells detected / cells significant):\n")
for (nm in names(G)) {
  gg <- G[[nm]]
  c2_msg("    %-6s : detected = %4d ; significant = %4d ; detected-but-NOT-significant = %d\n",
         nm, sum(gg$detected), sum(gg$sig), sum(gg$detected & !gg$sig))
}
c2_msg("    -> if the last column is 0, the location-matched significance layer is VACUOUS\n")
c2_msg("       on this dataset: the score collapses exactly onto detection stability.\n")

## ---- 4. Q4: what drives the score -- p, BH, or disappearance? -------
g <- G$rec50
det <- g[detected == TRUE]
det[, sig_p := emp_p < C2$FDR]                      # raw empirical significance
q4 <- det[, .(cells_detected = .N,
              frac_raw_p_sig = mean(sig_p),
              frac_BH_sig    = mean(q_R < C2$FDR),
              median_emp_p   = stats::median(emp_p),
              median_q       = stats::median(q_R),
              at_p_floor     = mean(emp_p <= 1 / (B + 1) + 1e-12),
              median_ntested = as.numeric(stats::median(n_tested))), by = label]
q4[, lost_to_BH := round(frac_raw_p_sig - frac_BH_sig, 3)]
setorder(q4, -frac_BH_sig)
fwrite(q4, file.path(C2$RES, "anchor_bh_behaviour.csv"))
c2_msg("\n[4] BH behaviour among DETECTED cells (rule = rec50):\n")
print(q4)

## per-cell: does BH get easier as the hypothesis set shrinks?
cellwise <- det[, .(n_tested = n_tested[1], n_sig = sum(q_R < C2$FDR),
                    n_rawsig = sum(emp_p < C2$FDR)), by = .(tau, lmin)]
c2_msg("\n[4] Spearman(n_tested, fraction of anchors BH-significant) = %.3f\n",
       stats::cor(cellwise$n_tested, cellwise$n_sig / pmax(cellwise$n_rawsig, 1), method = "spearman"))

## ---- 5. the required per-cell detail table --------------------------
detail <- G$rec50[, .(label, tau, lmin, detected, n_match, match_score,
                      size, score, maxC, emp_p, q_R, sig, n_tested)]
fwrite(detail, file.path(C2$RES, "anchor_locus_grid_detail.csv"))
c2_msg("\n[5] wrote results/anchor_locus_grid_detail.csv (%d rows), anchor_matching_rules.csv, anchor_bh_behaviour.csv\n",
       nrow(detail))
