## =============================================================================
## module_sim_LDscnR / sensitivity_grid_analyse.R
##
## Reduce the 80 panel files from sensitivity_grid.R (run on mini2) to the four
## things TODO_sensitivity.md asked for, reported on TWO AXES and never combined.
##
## 1. DOES THE HEADLINE HOLD ACROSS THE GRID? Precision over the 27 selection
##    cells runs 0.302-0.424 against the single-SNP baseline of 0.187, so the
##    advantage holds in EVERY cell at 1.62x-2.27x. Recall is flat at 0.20-0.25.
##
## 2. THE ld_w CONFOUND IS LARGELY DISMISSED. The concern was that V0.5_c1, the
##    cell where clustering wins, sits at an extreme ld_w operating point (5-6%
##    marker retention against 26-48% in V0.5_c2). Varying ld_w over the grid
##    moves V0.5_c1's precision by 1.09x (0.579/0.542/0.529) and V2_c1's by
##    1.04x, while the BETWEEN-CELL difference is 4-7x. So the cell effect is not
##    an ld_w artefact. V0.5_c2 IS ld_w-sensitive (2.1x), which is worth saying.
##
## 3. TEST AXIS. EMMAX and LFMM agree on precision (0.351 vs 0.337, r = 0.807)
##    but LFMM has markedly higher recall (0.325 vs 0.231, r = 0.882). THE ENGINE
##    MATTERS MORE FOR RECALL THAN ANY SELECTION PARAMETER DOES, which is the
##    clearest justification for reporting the two axes separately.
##
## 4. C, WITH THE PRE-REGISTERED READING AND NO MORE. C = the fraction of the 27
##    selection cells calling a 50 kb bin, at a fixed engine.
##
##      C            bins     % within dmax of a QTN
##      <= 0.25     36,133            2.1
##      0.25-0.5    24,035            2.5
##      0.5-0.75    48,619            1.8
##      0.75-<1     11,535            3.6
##      1.0         27,843            5.6
##
##    C = 1 is 4.0x enriched over bins called in a single setting. EVERYTHING
##    BELOW 0.75 IS FLAT, and 0.5-0.75 is LOWER than <=0.25, so C IS NOT A GRADED
##    SCORE -- it is an indicator that fires at 1. Spearman is 0.061. That is
##    exactly the reading TODO_sensitivity.md pre-registered ("regions close to 1
##    are probably safe; nothing is claimed about regions scoring low") and it
##    supports nothing stronger.
##
## A BIN-CONSTRUCTION TRAP, RECORDED BECAUSE IT NEARLY PRODUCED A NULL RESULT.
## `calls` expands each significant cluster over every bin it spans, so truth has
## to be on the same footing. Marking truth at the START BIN of a QTN-containing
## cluster instead diluted it by cluster span and gave 0.4% vs 0.9% with Spearman
## 0.021 -- a flat, uninformative C. Recomputing truth from QTN POSITIONS (within
## dmax) gives 1.4% vs 5.6%. Same data, opposite conclusion about whether C
## carries information.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/sensitivity_grid_analyse.R
## Env: DIR, OUT, SNP_PREC (single-SNP precision baseline from floor_sweep.R)
## =============================================================================
suppressMessages({library(data.table)})
DIR <- Sys.getenv("DIR", "module_sim_LDscnR/results/sensitivity_grid/panels")
OUT <- Sys.getenv("OUT", "module_sim_LDscnR/results/sensitivity_grid")
SNP_PREC <- as.numeric(Sys.getenv("SNP_PREC", "0.187"))
fs <- list.files(DIR, full.names = TRUE)

P <- rbindlist(lapply(fs, function(f) { z <- readRDS(f); if (is.null(z)) return(NULL)
  cbind(z$perf, cell = z$cell, tag = z$tag, env = z$env) }), fill = TRUE)
fwrite(P, file.path(OUT, "perf.csv"))

S <- P[engine == "emx", .(prec = mean(prec), rec = mean(rec)), by = .(ldw, rho, floor)]
cat(sprintf("1. precision over the 27 selection cells: %.3f - %.3f (baseline %.3f)\n",
            min(S$prec), max(S$prec), SNP_PREC))
cat(sprintf("   advantage %.2fx - %.2fx | holds in every cell: %s\n",
            min(S$prec)/SNP_PREC, max(S$prec)/SNP_PREC, all(S$prec > SNP_PREC)))
cat("\n2. precision by simulation cell x ld_w (the confound check)\n")
print(dcast(P[engine == "emx", .(prec = round(mean(prec), 3)), by = .(cell, ldw)],
            cell ~ ldw, value.var = "prec"))
E <- dcast(P[, .(prec = mean(prec), rec = mean(rec)), by = .(ldw, rho, floor, engine)],
           ldw + rho + floor ~ engine, value.var = c("prec", "rec"))
cat(sprintf("\n3. test axis: precision emx %.3f / lfm %.3f (r %.3f); recall emx %.3f / lfm %.3f (r %.3f)\n",
            mean(E$prec_emx), mean(E$prec_lfm), stats::cor(E$prec_emx, E$prec_lfm),
            mean(E$rec_emx), mean(E$rec_lfm), stats::cor(E$rec_emx, E$rec_lfm)))
fwrite(E, file.path(OUT, "two_axes.csv"))

## ---- C, on the selection axis only, scored against positional truth
CC <- rbindlist(lapply(fs, function(f) { z <- readRDS(f)
  if (is.null(z) || !nrow(z$calls)) return(NULL)
  cl <- z$calls[engine == "emx"]
  cl[, .(C = uniqueN(paste(ldw, rho, floor)) / 27), by = .(chrfile, Chr, bin)][
     , `:=`(cell = z$cell, tag = z$tag, env = z$env)][] }), fill = TRUE)
TB <- as.data.table(readRDS(file.path(OUT, "qtn_bins_true.rds")))[, true_bin := TRUE]
M  <- merge(CC, TB, by = c("cell","tag","env","chrfile","Chr","bin"), all.x = TRUE)
M[is.na(true_bin), true_bin := FALSE]
fwrite(M, file.path(OUT, "C_bins_scored.csv"))
M[, Cbin := cut(C, c(0, .25, .5, .75, .99, 1),
                labels = c("<=.25", ".25-.5", ".5-.75", ".75-<1", "1.0"), include.lowest = TRUE)]
cat("\n4. C against positional truth (within dmax of a QTN)\n")
print(M[, .(bins = .N, pct_true = round(100 * mean(true_bin), 1)), by = Cbin][order(Cbin)])
cat(sprintf("   Spearman %.3f | C=1 enrichment over single-setting bins %.1fx\n",
            stats::cor(M$C, as.numeric(M$true_bin), method = "spearman"),
            mean(M[C == 1]$true_bin) / mean(M[C <= 1/27 + 1e-9]$true_bin)))
