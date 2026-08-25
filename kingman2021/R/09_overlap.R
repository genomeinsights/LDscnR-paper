## Overlap between the 3sp LD-aware outlier regions (module_sticklebacks_LDscnR,
## tau_C=0.05 / l_min=10 / rho_ld=0.60) and the Kingman EcoPeaks, in both directions:
##   (a) region-level overlap vs a within-chromosome shuffle null   [gasAcu1 space]
##   (b) rank-based enrichment of low Kingman p inside those regions [gasAcu1-4 space]
## Writes data/overlap_summary.csv, data/overlap_detail.csv, data/enrichment_summary.csv
suppressMessages(library(data.table))
P    <- path.expand("~/gitlab/LDscnR-paper/kingman2021")
EP   <- path.expand("~/gitlab/LD-scaling-genome-scans/empirical_data/kingman2021/ecopeaks")
LIFT <- file.path(P, "data", "liftover")
ROMAN <- c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV",
           "XVI","XVII","XVIII","XIX","XX","XXI")
rn <- function(x) match(sub("^chr", "", x), ROMAN)
B  <- 2000L; set.seed(1)

## ---------- (a) region-level overlap, gasAcu1 ----------------------------------------
L <- fread(file.path(P,"data","regions_tau0.05_lmin10_rho0.60.csv"))[method=="LFMM"]
L[, chr_num := as.integer(gsub("Chr","",Chr))]
d3 <- readRDS(path.expand("~/gitlab/LDscnR-paper/module_sticklebacks_LDscnR/data/3sp_LDscnR_data.rds"))
m3 <- as.data.table(d3$map); m3[, chr_num := as.integer(gsub("Chr","",Chr))]
rng <- m3[, .(lo=min(Pos), hi=max(Pos)), by=chr_num]; setkey(rng, chr_num); rm(d3, m3); gc()

rd <- function(s) { Pk <- fread(file.path(LIFT, paste0("pv_",s,".bed")), header=FALSE,
                                col.names=c("chr","start","end","pv"))
  Pk[, chr_num := rn(chr)]; Pk[, p_snp := as.numeric(tstrsplit(pv,"\\|")[[1]])]
  Pk <- Pk[!is.na(chr_num), .(chr_num,start,end,p_snp)]; setkey(Pk,chr_num,start,end); Pk }
ovl <- function(A, Pk) { A2 <- A[, .(chr_num,start,end,id=.I)]
  o <- foverlaps(A2, Pk, by.x=c("chr_num","start","end"), type="any", nomatch=NULL)
  if (!nrow(o)) return(list(nhit=0L, bp=0)); o[, w := pmin(i.end,end)-pmax(i.start,start)]
  list(nhit=uniqueN(o$id), bp=sum(o$w)) }

SETS <- c("c155.specific","c155.sensitive","c150.specific","c150.sensitive")
summ <- rbindlist(lapply(SETS, function(s) {
  Pk <- rd(s); o <- ovl(L, Pk)
  tot <- rng[, sum(hi-lo)]; pkbp <- Pk[chr_num %in% rng$chr_num, sum(end-start)]
  r <- rng[J(L$chr_num)]; nh <- integer(B); nb <- numeric(B)
  for (b in seq_len(B)) { st <- r$lo + floor(runif(nrow(L))*pmax(1, r$hi-r$lo-L$span))
    on <- ovl(data.table(chr_num=L$chr_num, start=st, end=st+L$span), Pk)
    nh[b] <- on$nhit; nb[b] <- on$bp }
  data.table(set=s, n_peaks=nrow(Pk), peak_Mb=pkbp/1e6, genome_pct=100*pkbp/tot,
             n_regions=nrow(L), regions_hit=o$nhit, regions_hit_null=mean(nh),
             fold_region=o$nhit/mean(nh), p_region=(1+sum(nh>=o$nhit))/(B+1),
             bp_Mb=o$bp/1e6, bp_pct=100*o$bp/sum(L$span), bp_pct_null=100*mean(nb)/sum(L$span),
             fold_bp=o$bp/mean(nb), p_bp=(1+sum(nb>=o$bp))/(B+1)) }))
fwrite(summ, file.path(P,"data","overlap_summary.csv")); print(summ[, 1:9])

## ---------- (a2) per-region detail with gene names -----------------------------------
G <- fread(file.path(P,"data","gname_g1.bed"), header=FALSE, col.names=c("chr","start","end","gene"))
G[, chr_num := rn(chr)]; G <- G[!is.na(chr_num)]; setkey(G, chr_num, start, end)
det <- rbindlist(lapply(c(GlobalSpec="c155.specific", PacSpec="c150.specific"), function(s) {
  Pk <- rd(s)
  o <- foverlaps(L[,.(chr_num,start,end,Chr,region)], Pk, by.x=c("chr_num","start","end"),
                 type="any", nomatch=NULL)
  if (!nrow(o)) return(NULL)
  o[, `:=`(ov_start=pmax(i.start,start), ov_end=pmin(i.end,end))][, ov_kb := (ov_end-ov_start)/1e3]
  o[, .(Chr, region, reg_start=i.start, reg_end=i.end, peak_start=start, peak_end=end,
        p_snp, ov_start, ov_end, ov_kb, chr_num)] }), idcol="set")
det[, ridx := .I]
gg <- foverlaps(det[,.(chr_num, start=ov_start, end=ov_end, ridx)], G,
                by.x=c("chr_num","start","end"), type="any", nomatch=NULL)
gg <- gg[!grepl("^(si:|zgc:|CABZ|BX|CU|CR|AL|FP|LO)", gene, ignore.case=TRUE)]
det <- merge(det, gg[, .(genes=paste(head(unique(gene),8), collapse=", ")), by=ridx],
             by="ridx", all.x=TRUE)
det[is.na(genes), genes := ""]; setorder(det, -ov_kb)
fwrite(det, file.path(P,"data","overlap_detail.csv"))
cat(sprintf("\noverlap_detail.csv: %d region x peak pairs\n", nrow(det)))

## ---------- (b) rank-based enrichment, gasAcu1-4 --------------------------------------
R14 <- fread(file.path(LIFT,"lfmm_g14.bed"), header=FALSE, col.names=c("Chr","start","end","name"))
enr <- rbindlist(lapply(c("c151_nEur","c155_global","c150_pacNW"), function(coh) {
  d <- fread(file.path(EP, paste0(coh,".snp_p.tsv.gz")), select=c("Chr","Pos","p"))
  RR <- R14[Chr %in% d$Chr]; RR[, span := end-start]
  rg <- d[, .(lo=min(Pos), hi=max(Pos)), by=Chr]; setkey(rg, Chr)
  mk <- function(X) { setkey(X, Chr, start, end)
    !is.na(foverlaps(d[, .(Chr, start=Pos, end=Pos)], X, by.x=c("Chr","start","end"),
                     type="any", mult="first", nomatch=NA)$name) }
  inreg <- mk(copy(RR)); r <- rg[J(RR$Chr)]
  rbindlist(lapply(c(0.001, 0.01), function(TOP) {
    thr <- quantile(d$p, TOP); obs <- mean(d$p[inreg] <= thr); nv <- numeric(500)
    for (b in 1:500) { st <- r$lo + floor(runif(nrow(RR))*pmax(1, r$hi-r$lo-RR$span))
      nv[b] <- mean(d$p[mk(data.table(Chr=RR$Chr, start=st, end=st+RR$span, name="x"))] <= thr) }
    data.table(cohort=coh, n_snp=nrow(d), top=TOP, p_thresh=thr, in_region_rate=obs,
               null_rate=mean(nv), fold=obs/mean(nv), pval=(1+sum(nv>=obs))/501) })) }))
fwrite(enr, file.path(P,"data","enrichment_summary.csv")); print(enr[, .(cohort, top, fold, pval)])
