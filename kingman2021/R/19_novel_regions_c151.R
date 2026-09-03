## Are the 3sp regions that overlap NO published EcoPeak still real?
##
## The sensitivity reading of the LFMM surplus says its extra regions are genuine
## Atlantic marine-freshwater loci that Kingman's cohorts cannot see. That predicts they
## should carry signal in an INDEPENDENT, geography-matched cohort. c151 (Northern Europe,
## 9 marine / 18 freshwater) is exactly that: it is in the same VCF but was NOT used to
## call the published EcoPeaks (only c150 and c155 were), so it is independent of the
## truth set the regions are being scored against.
##
## Design: split each engine's regions into CORROBORATED (overlap a specific EcoPeak) and
## NOVEL (overlap none), then test each group for enrichment of low c151 p-values against a
## within-chromosome rotation null. The corroborated group is the INTERNAL POSITIVE CONTROL:
## c151 is underpowered (no FDR-significant SNP anywhere), so if the corroborated regions
## also show nothing, the test is uninformative rather than negative.
##
## CIRCULARITY: all 27 c151 samples are a SUBSET of the 84 c155 samples, so classifying
## regions by c155 peaks is circular in the deflationary direction -- a locus where the
## c151 samples carry signal is likelier to have been called a c155 peak and so labelled
## corroborated, depleting signal from the "novel" class by construction. c150 (Pacific)
## shares ZERO samples with c151, so PEAKSET="c150" gives a fully independent split and is
## the version to trust. PEAKSET="both" is reported for comparison only.
##
## Usage: PEAKSET=c150 Rscript kingman2021/R/19_novel_regions_c151.R [B]
## Writes data/novel_regions_c151.csv
suppressMessages(library(data.table))
a <- commandArgs(trailingOnly=TRUE); B <- if (length(a)>=1) as.integer(a[1]) else 500L
P    <- path.expand("~/gitlab/LDscnR-paper/kingman2021")
EP   <- path.expand("~/gitlab/LD-scaling-genome-scans/empirical_data/kingman2021/ecopeaks")
TR   <- path.expand("~/gitlab/LD-scaling-genome-scans/empirical_data/kingman2021/tracks")
LIFT <- file.path(P,"data","liftover")
CHAIN<- file.path(LIFT,"gasAcu1ToGasAcu1-4.chain")
ROMAN <- c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV",
           "XVI","XVII","XVIII","XIX","XX","XXI")
set.seed(1)
if (!nzchar(Sys.which("liftOver"))) stop("liftOver not on PATH")

## ---- lift both l_min=3 region sets into gasAcu1-4 ---------------------------------
lift_set <- function(csv, meth, tag) {
  R <- fread(file.path(P,"data",csv))[method==meth]
  R[, chrR := paste0("chr", ROMAN[as.integer(gsub("Chr","",Chr))])]
  fin <- file.path(LIFT, paste0(tag,"_g1.bed")); fout <- file.path(LIFT, paste0(tag,"_g14.bed"))
  fwrite(R[, .(chrR, start, end, name=paste0(tag,"_",region))], fin, sep="\t", col.names=FALSE)
  system2("liftOver", c("-minMatch=0.5", fin, CHAIN, fout, tempfile()), stdout=FALSE, stderr=FALSE)
  x <- fread(fout, header=FALSE, col.names=c("Chr","start","end","name"))
  cat(sprintf("%-14s %3d/%3d regions lifted, %.2f Mb\n", tag, nrow(x), nrow(R), sum(x$end-x$start)/1e6)); x
}
L <- lift_set("regions_tau0.05_lmin3_rho0.60_lfmm.csv",  "LFMM",  "lfmm136")
E <- lift_set("regions_tau0.05_lmin3_rho0.60_emmax.csv", "EMMAX", "emmax17")

## ---- classify: does the region overlap ANY specific EcoPeak? -----------------------
PEAKSET <- Sys.getenv("PEAKSET", "c150")
SETS <- if (PEAKSET == "c150") "c150.specific" else c("c155.specific","c150.specific")
cat(sprintf("\nclassifying regions by: %s\n", paste(SETS, collapse=" + ")))
pk <- rbindlist(lapply(SETS, function(s)
  fread(file.path(TR, sprintf("gasAcu1-4.%s.50kb.final.peaks.bed", s)), header=FALSE,
        select=1:3, col.names=c("Chr","start","end"))))
setkey(pk, Chr, start, end)
classify <- function(R) { R <- copy(R)
  o <- foverlaps(R[, .(Chr,start,end)], pk, by.x=c("Chr","start","end"), type="any", mult="first", nomatch=NA)
  R[, corroborated := !is.na(o$start)][] }
L <- classify(L); E <- classify(E)

## ---- c151 signal per group vs rotation null ----------------------------------------
d <- fread(file.path(EP,"c151_nEur.snp_p.tsv.gz"), select=c("Chr","Pos","p"))
rg <- d[, .(lo=min(Pos), hi=max(Pos)), by=Chr]; setkey(rg, Chr)
TOP <- 0.01; thr <- quantile(d$p, TOP)
mk <- function(X){ setkey(X, Chr, start, end)
  !is.na(foverlaps(d[, .(Chr, start=Pos, end=Pos)], X, by.x=c("Chr","start","end"),
                   type="any", mult="first", nomatch=NA)$name) }
score <- function(R, lab) {
  R <- R[Chr %in% d$Chr]; if (!nrow(R)) return(NULL)
  R[, span := end-start]
  inreg <- mk(copy(R)); obs <- mean(d$p[inreg] <= thr)
  r <- rg[J(R$Chr)]; nv <- numeric(B)
  for (b in seq_len(B)) { st <- r$lo + floor(runif(nrow(R))*pmax(1, r$hi-r$lo-R$span))
    nv[b] <- mean(d$p[mk(data.table(Chr=R$Chr, start=st, end=st+R$span, name="x"))] <= thr) }
  data.table(group=lab, n_regions=nrow(R), span_Mb=round(sum(R$span)/1e6,2), n_snp=sum(inreg),
             rate=round(obs,4), null=round(mean(nv),4), fold=round(obs/mean(nv),2),
             p=(1+sum(nv>=obs))/(B+1))
}
res <- rbindlist(list(
  score(E[corroborated==TRUE ], "EMMAX corroborated"),
  score(E[corroborated==FALSE], "EMMAX NOVEL"),
  score(L[corroborated==TRUE ], "LFMM corroborated"),
  score(L[corroborated==FALSE], "LFMM NOVEL")))
cat(sprintf("\nc151 top %.0f%% of p (p<=%.3g); rotation null B=%d\n", 100*TOP, thr, B))
print(res)
res[, peakset := PEAKSET]
fwrite(res, file.path(P,"data", sprintf("novel_regions_c151_%s.csv", PEAKSET)))
