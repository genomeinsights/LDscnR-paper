## =====================================================================
## module_sticklebacks_LDscnR / kingman_cross_set_table.R
##
## Cross-set corroboration table for the 3sp outlier regions (manuscript Table S3,
## \label{tab:crossset}). Consensus loci across the four 3sp scans -- EMMAX and LFMM,
## each clustered at l_min = 3 and l_min = 10 -- with two Kingman *specific* EcoPeak
## sets (global, Pacific) as independent validation ticks.
##
## For each consensus locus (regions merged across sets by overlap, 10 kb join
## tolerance) it records which sets recover it; rows are the loci found by >= 2 sets,
## ordered by corroboration, and a final "unique" row counts regions found by only
## one set. All coordinates are gasAcu1; the Kingman peaks are the lifted BEDs from
## the kingman2021 module, read READ-ONLY (this script writes nothing there).
##
## Run from the LDscnR-paper root (needs the 3sp bundle + kingman2021/data/liftover):
##   Rscript module_sticklebacks_LDscnR/kingman_cross_set_table.R
## Writes module_sticklebacks_LDscnR/results/cross_set_kingman.csv.
## =====================================================================
suppressMessages({ library(data.table); library(LDscnR) })

BND  <- "module_sticklebacks_LDscnR/data/3sp_LDscnR_data.rds"
LIFT <- "kingman2021/data/liftover"                       # read-only (other module)
OUT  <- "module_sticklebacks_LDscnR/results/cross_set_kingman.csv"
ROMAN <- c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV",
           "XVI","XVII","XVIII","XIX","XX","XXI")
GAP <- 1e4                                                 # merge tolerance between sets
if (!dir.exists(dirname(OUT))) dir.create(dirname(OUT), recursive = TRUE)

## ---- 1. the four 3sp region sets (EMMAX/LFMM x l_min 3/10) --------------------------
d  <- readRDS(BND); m3 <- as.data.table(d$map)
C_emx  <- ld_cscore(emmax_fast(emmax_setup(d$GTs, d$GRM), d$eco), d$ld_ws, alpha = 0.05, qstar = seq(0, .95, .05))
C_lfmm <- ld_cscore(m3$lfmm_p, d$ld_ws, alpha = 0.05, qstar = seq(0, .95, .05))
names(C_emx) <- names(C_lfmm) <- m3$marker
uni   <- union(names(C_emx)[C_emx > 0], names(C_lfmm)[C_lfmm > 0])
edges <- ld_edges(uni, d$GTs, m3[, .(marker, Chr, Pos)], as.data.table(d$LD_decay$decay_sum),
                  rho_ld = 0.60, dcap = 1e5)
mpos <- stats::setNames(m3$Pos, m3$marker); mchr <- stats::setNames(as.character(m3$Chr), m3$marker)
regset <- function(C, lmin) { r <- ld_regions(names(C)[C >= 0.05], edges); r <- r[lengths(r) >= lmin]
  if (!length(r)) return(data.table(Chr = character(), start = numeric(), end = numeric()))
  rbindlist(lapply(r, function(x) data.table(Chr = unname(mchr[x[1]]), start = min(mpos[x]), end = max(mpos[x])))) }
sets <- list(E3 = regset(C_emx, 3L), E10 = regset(C_emx, 10L), L3 = regset(C_lfmm, 3L), L10 = regset(C_lfmm, 10L))
cat(sprintf("[1] set sizes: E3=%d E10=%d L3=%d L10=%d\n",
            nrow(sets$E3), nrow(sets$E10), nrow(sets$L3), nrow(sets$L10)))

## ---- 2. merge into consensus loci, record set membership ---------------------------
allr <- rbindlist(sets, idcol = "set"); allr[, chr_num := as.integer(gsub("Chr", "", Chr))]
allr[, cl := { o <- order(start); s <- start[o]; e <- end[o]
  g <- c(TRUE, s[-1] > cummax(e)[-length(e)] + GAP); cumsum(g)[order(o)] }, by = chr_num]
allr[, locus := paste(chr_num, cl, sep = "_")]
con  <- allr[, .(chr_num = chr_num[1], start = min(start), end = max(end)), by = locus]
pres <- dcast(unique(allr[, .(locus, set)]), locus ~ set, fun.aggregate = length, value.var = "set")
for (s in c("E3", "E10", "L3", "L10")) if (is.null(pres[[s]])) pres[[s]] <- 0L
con  <- merge(con, pres, by = "locus")
con[, n3sp := (E3 > 0) + (E10 > 0) + (L3 > 0) + (L10 > 0)]

## ---- 3. Kingman specific EcoPeaks (gasAcu1) -> validation ticks ---------------------
rd <- function(f) { b <- fread(file.path(LIFT, f), header = FALSE, col.names = c("chr", "start", "end", "pv"))
  b[, chr_num := match(sub("^chr", "", chr), ROMAN)]; b[!is.na(chr_num), .(chr_num, start, end)] }
tick <- function(P) { setkey(P, chr_num, start, end)
  !is.na(foverlaps(con[, .(chr_num, start, end, locus)], P, by.x = c("chr_num", "start", "end"),
                   type = "any", mult = "first", nomatch = NA)$start) }
con[, KGlob := tick(rd("pv_c155.specific.bed"))]
con[, KPac  := tick(rd("pv_c150.specific.bed"))]
con[, ntot := n3sp + (KGlob > 0) + (KPac > 0)]

## ---- 4. labels, shared table, unique row -------------------------------------------
con[, Mb := start / 1e6][, lab := sprintf("Chr%d:%.2f", chr_num, Mb)]
con[chr_num == 1 & start < 21.93e6 & end > 21.40e6, lab := paste0(lab, " (inversion)")]
con[chr_num == 4 & start < 12.83e6 & end > 12.79e6, lab := paste0(lab, " (Eda)")]
tk <- function(x) ifelse(x > 0 | x == TRUE, "Y", "")
shared <- con[ntot >= 2][order(-ntot, chr_num, start)]
tab <- shared[, .(Locus = lab, n = ntot, E3 = tk(E3), E10 = tk(E10), L3 = tk(L3), L10 = tk(L10),
                  KGlob = tk(KGlob), KPac = tk(KPac))]
uniq <- con[ntot == 1, .(E3 = sum(E3 > 0), E10 = sum(E10 > 0), L3 = sum(L3 > 0), L10 = sum(L10 > 0))]
cat(sprintf("[4] consensus loci=%d ; shared (>=2 sets)=%d ; LFMM l_min=3 singletons=%d\n",
            nrow(con), nrow(shared), uniq$L3))
cat("\n=== shared-loci table ===\n"); print(tab, nrow = 100)
cat("\n=== unique (found by only one of the four sets) ===\n"); print(uniq)
fwrite(tab, OUT)
fwrite(cbind(Locus = "UNIQUE (1 set only)", n = NA, uniq, KGlob = "", KPac = ""), OUT, append = TRUE)
cat(sprintf("[4] wrote %s\n", OUT))
