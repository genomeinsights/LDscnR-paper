## How much do region sets agree with each other -- INCLUDING Kingman's own two peak sets
## with each other? That last pair is the ceiling: if the published Global-specific and
## Pacific-specific EcoPeaks only partly agree, then partial agreement between a 3sp region
## set and either of them is the expected result, not a failure.
##
## Everything is done in gasAcu1 space (Kingman peaks lifted by R/08) on the same
## SNP-covered chromosome ranges, with the SAME uniform-placement rotation null used in
## R/09 and R/12, so the numbers are commensurable with those sections.
##
## Reported per pair:
##   jaccard   bp intersection / bp union            -- but CAPPED at min/max span, so a
##             pair of very unequal sets can never score high; read overlap_coef instead
##   overlap_coef  bp intersection / bp of the SMALLER set (Szymkiewicz-Simpson) -- 1.0 means
##             the smaller set is entirely contained in the larger. This is the fair
##             cross-pair measure: it asks "how much of the smaller claim is corroborated?"
##   hit_A     fraction of A's features hitting B
##   fold      hit_A over the rotation null           -- span-sensitive, see ceiling
##   ceiling   the largest fold this pair could show (n_A / null_mean)
## Writes data/peakset_concordance.csv
suppressMessages(library(data.table))
P    <- path.expand("~/gitlab/LDscnR-paper/kingman2021")
LIFT <- file.path(P, "data", "liftover")
ROMAN <- c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV",
           "XVI","XVII","XVIII","XIX","XX","XXI")
B <- 2000L; set.seed(1)

d  <- readRDS(path.expand("~/gitlab/LDscnR-paper/module_sticklebacks_LDscnR/data/3sp_LDscnR_data.rds"))
m3 <- as.data.table(d$map); m3[, chr_num := as.integer(gsub("Chr","",Chr))]
rng <- m3[, .(lo=min(Pos), hi=max(Pos)), by=chr_num]; setkey(rng, chr_num); rm(d, m3); gc()

## --- the sets, all as (chr_num, start, end) in gasAcu1 -------------------------------
kpeak <- function(s) { x <- fread(file.path(LIFT, paste0("pv_",s,".bed")), header=FALSE,
                                  col.names=c("chr","start","end","pv"))
  x[, chr_num := match(sub("^chr","",chr), ROMAN)]; x[!is.na(chr_num), .(chr_num,start,end)] }
oreg <- function(f, meth) { x <- fread(file.path(P,"data",f))[method==meth]
  x[, .(chr_num = as.integer(gsub("Chr","",Chr)), start, end)] }

S <- list(
  `Kingman Global-spec`  = kpeak("c155.specific"),
  `Kingman Pacific-spec` = kpeak("c150.specific"),
  `Kingman Global-sens`  = kpeak("c155.sensitive"),
  `Kingman Pacific-sens` = kpeak("c150.sensitive"),
  `3sp EMMAX l_min=3`    = oreg("regions_tau0.05_lmin3_rho0.60_emmax.csv", "EMMAX"),
  `3sp LFMM l_min=10`    = oreg("regions_tau0.05_lmin10_rho0.60.csv",      "LFMM"),
  `3sp LFMM l_min=3`     = oreg("regions_tau0.05_lmin3_rho0.60_lfmm.csv",  "LFMM"))
S <- lapply(S, function(x) x[chr_num %in% rng$chr_num][order(chr_num, start)])

hit_frac <- function(A, Bk) { setkey(Bk, chr_num, start, end)
  o <- foverlaps(A[, .(chr_num,start,end,id=.I)], Bk, by.x=c("chr_num","start","end"),
                 type="any", nomatch=NULL)
  if (!nrow(o)) return(0); uniqueN(o$id)/nrow(A) }
bp_int <- function(A, Bk) { setkey(Bk, chr_num, start, end)
  o <- foverlaps(A[, .(chr_num,start,end)], Bk, by.x=c("chr_num","start","end"),
                 type="any", nomatch=NULL)
  if (!nrow(o)) return(0); o[, sum(pmin(i.end,end)-pmax(i.start,start))] }

pairs <- list(c(1,2), c(3,4), c(5,1), c(5,2), c(6,1), c(6,2), c(7,1), c(7,2), c(5,7), c(5,6))
res <- rbindlist(lapply(pairs, function(ij) {
  A <- copy(S[[ij[1]]]); Bk <- copy(S[[ij[2]]]); nm <- names(S)
  A[, w := end-start]
  obs <- hit_frac(A, Bk); inter <- bp_int(A, Bk)
  uni <- sum(A$w) + sum(Bk$end-Bk$start) - inter
  r <- rng[J(A$chr_num)]; nv <- numeric(B)
  for (b in seq_len(B)) { st <- r$lo + floor(runif(nrow(A))*pmax(1, r$hi-r$lo-A$w))
    nv[b] <- hit_frac(data.table(chr_num=A$chr_num, start=st, end=st+A$w), Bk) }
  data.table(A=nm[ij[1]], B=nm[ij[2]], n_A=nrow(A), n_B=nrow(Bk),
             span_A_Mb=round(sum(A$w)/1e6,2), span_B_Mb=round(sum(Bk$end-Bk$start)/1e6,2),
             inter_Mb=round(inter/1e6,3),
             overlap_coef=round(inter/min(sum(A$w), sum(Bk$end-Bk$start)),3),
             jaccard=round(inter/uni,4), hit_A=round(obs,3), null=round(mean(nv),3),
             fold=round(obs/mean(nv),2), ceiling=round(1/mean(nv),1),
             p=(1+sum(nv>=obs))/(B+1)) }))
fwrite(res, file.path(P,"data","peakset_concordance.csv"))
setorder(res, -overlap_coef)
print(res[, .(A, B, span_A_Mb, span_B_Mb, inter_Mb, overlap_coef, jaccard, fold, ceiling, p)])
