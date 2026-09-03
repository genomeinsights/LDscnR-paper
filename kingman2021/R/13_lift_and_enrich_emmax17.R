## Lift the 17 EMMAX l_min=3 regions gasAcu1 -> gasAcu1-4 with the UCSC liftOver
## binary (same tool R/08 used), VALIDATE the lift against the module's already-lifted
## LFMM regions, then recompute the two geography-matched positive controls against the
## 17-region set (for consistency with the main manuscript's structure-null analysis):
##   (b) rank-based enrichment of low Kingman p inside the regions   (per cohort)
##   (c) do the c151 N.-European suggestive peaks land in the regions more than chance
##
## PREREQ: `liftOver` on PATH (http://hgdownload.soe.ucsc.edu/admin/exe/) -- the same
## environment that ran R/08_liftover.sh. Also needs R/12_overlap_emmax17.R to have been
## run first (it writes data/regions_tau0.05_lmin3_rho0.60_emmax.csv).
##
## Writes data/liftover/emmax17_g14.bed, data/enrichment_emmax17.csv,
##        data/c151_peak_overlap_emmax17.csv
suppressMessages(library(data.table))
P    <- path.expand("~/gitlab/LDscnR-paper/kingman2021")
EP   <- path.expand("~/gitlab/LD-scaling-genome-scans/empirical_data/kingman2021/ecopeaks")
LIFT <- file.path(P, "data", "liftover")
CHAIN <- file.path(LIFT, "gasAcu1ToGasAcu1-4.chain")
ROMAN <- c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV",
           "XVI","XVII","XVIII","XIX","XX","XXI")
arab2rom <- function(chr) paste0("chr", ROMAN[as.integer(gsub("Chr","",chr))])

if (!nzchar(Sys.which("liftOver")))
  stop("liftOver not on PATH -- run this in the environment that ran R/08_liftover.sh")
lift_bed <- function(inbed, outbed) {
  ## -minMatch=0.5 to match R/08_liftover.sh, which produced the reference lfmm_g14.bed;
  ## liftOver's default (0.95) drops 3 of the 39 LFMM regions and fails the check below.
  system2("liftOver", c("-minMatch=0.5", inbed, CHAIN, outbed, tempfile()),
          stdout = FALSE, stderr = FALSE)
  fread(outbed, header = FALSE, col.names = c("Chr", "start", "end", "name"))
}

## ---- VALIDATE the lift: lfmm_g1 -> should reproduce lfmm_g14 ------------------------
myv <- lift_bed(file.path(LIFT, "lfmm_g1.bed"), tempfile(fileext = ".bed"))
g14 <- fread(file.path(LIFT, "lfmm_g14.bed"), header = FALSE, col.names = c("Chr", "start", "end", "name"))
v <- merge(myv, g14, by = "name", suffixes = c("_my", "_ref"))
cat(sprintf("[validate] %d/%d LFMM regions lifted; max |start diff| %d bp, max |end diff| %d bp\n",
            nrow(v), nrow(g14), max(abs(v$start_my - v$start_ref)), max(abs(v$end_my - v$end_ref))))
stopifnot(nrow(v) >= 39, max(abs(v$start_my - v$start_ref), abs(v$end_my - v$end_ref)) < 50)

## ---- lift the 17 EMMAX regions -----------------------------------------------------
E <- fread(file.path(P, "data", "regions_tau0.05_lmin3_rho0.60_emmax.csv"))
E[, chrR := arab2rom(Chr)]
fwrite(E[, .(chrR, start, end, name = paste0("EMMAX_", region))],
       file.path(LIFT, "emmax17_g1.bed"), sep = "\t", col.names = FALSE)
E14 <- lift_bed(file.path(LIFT, "emmax17_g1.bed"), file.path(LIFT, "emmax17_g14.bed"))
cat(sprintf("[lift] %d/%d EMMAX regions lifted to gasAcu1-4\n", nrow(E14), nrow(E)))

## ---- (b) rank-based enrichment of low Kingman p inside the regions -----------------
set.seed(1)
enr <- rbindlist(lapply(c("c151_nEur", "c155_global", "c150_pacNW"), function(coh) {
  d <- fread(file.path(EP, paste0(coh, ".snp_p.tsv.gz")), select = c("Chr", "Pos", "p"))
  RR <- E14[Chr %in% d$Chr]; RR[, span := end - start]
  rg <- d[, .(lo = min(Pos), hi = max(Pos)), by = Chr]; setkey(rg, Chr)
  mk <- function(X) { setkey(X, Chr, start, end)
    !is.na(foverlaps(d[, .(Chr, start = Pos, end = Pos)], X, by.x = c("Chr", "start", "end"),
                     type = "any", mult = "first", nomatch = NA)$name) }
  inreg <- mk(copy(RR)); r <- rg[J(RR$Chr)]
  rbindlist(lapply(c(0.001, 0.01), function(TOP) {
    thr <- quantile(d$p, TOP); obs <- mean(d$p[inreg] <= thr); nv <- numeric(500)
    for (b in 1:500) { st <- r$lo + floor(runif(nrow(RR)) * pmax(1, r$hi - r$lo - RR$span))
      nv[b] <- mean(d$p[mk(data.table(Chr = RR$Chr, start = st, end = st + RR$span, name = "x"))] <= thr) }
    data.table(cohort = coh, top = TOP, in_region_rate = obs, null_rate = mean(nv),
               fold = round(obs / mean(nv), 3), pval = (1 + sum(nv >= obs)) / 501) })) }))
fwrite(enr, file.path(P, "data", "enrichment_emmax17.csv"))
cat("\n=== enrichment of low Kingman p inside the 17 EMMAX regions ===\n"); print(enr[, .(cohort, top, fold, pval)])

## ---- (c) c151 suggestive-peak overlap ---------------------------------------------
pk <- fread(file.path(P, "data", "peaks", "c151_nEur_suggestivePeaks.bed"), header = FALSE,
            col.names = c("Chr", "start", "end", "n_top", "min_p"))
dc <- fread(file.path(EP, "c151_nEur.snp_p.tsv.gz"), select = c("Chr", "Pos"))
rng <- dc[, .(lo = min(Pos), hi = max(Pos)), by = Chr]; setkey(rng, Chr); setkey(E14, Chr, start, end)
hits <- function(X) nrow(foverlaps(X, E14, by.x = c("Chr", "start", "end"), type = "any", nomatch = NULL))
obs <- hits(pk); pk[, w := end - start]; set.seed(1); B <- 5000L; nl <- integer(B); r <- rng[J(pk$Chr)]
for (b in seq_len(B)) { st <- r$lo + floor(runif(nrow(pk)) * pmax(1, r$hi - r$lo - pk$w))
  nl[b] <- hits(data.table(Chr = pk$Chr, start = st, end = st + pk$w)) }
res <- data.table(n_peaks = nrow(pk), n_inside = obs, null_mean = round(mean(nl), 3),
                  fold = round(obs / mean(nl), 3), pval = (1 + sum(nl >= obs)) / (B + 1))
fwrite(res, file.path(P, "data", "c151_peak_overlap_emmax17.csv"))
cat("\n=== c151 suggestive peaks inside the 17 EMMAX regions ===\n"); print(res)
