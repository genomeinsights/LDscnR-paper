## c151 (Northern Europe) cannot reach genome-wide significance: with 9 marine / 18
## freshwater (median 7 marine actually called), the exact multivariate-hypergeometric
## test bottoms out at p = 3.2e-7, while BH at q<=0.01 over 1.7M tests would need ~54
## SNPs at that floor. Only 1 attains it. So there is no FDR-significant c151 peak set.
## This is very likely why the hub publishes EcoPeak BEDs for c150 and c155 only, even
## though Table S2 also defines c151/c153/c154 cohorts.
##
## What IS still usable is the *ranked* signal. This script emits a suggestive peak set
## from the top-q quantile of -log10 p, clearly flagged as NOT FDR-controlled.
## Usage: Rscript 06_c151_suggestive_peaks.R [top_frac] [merge_kb] [min_sig]
suppressMessages(library(data.table))
a <- commandArgs(trailingOnly=TRUE)
TOP   <- if (length(a)>=1) as.numeric(a[1]) else 0.001
MERGE <- (if (length(a)>=2) as.numeric(a[2]) else 50)*1000
MINSIG<- if (length(a)>=3) as.integer(a[3]) else 2L
DATA <- path.expand("~/gitlab/LD-scaling-genome-scans/empirical_data/kingman2021")
d <- fread(file.path(DATA,"ecopeaks","c151_nEur.snp_p.tsv.gz"))
thr <- quantile(d$p, TOP)
sig <- d[p <= thr]; setorder(sig, Chr, Pos)
cat(sprintf("c151: %d sites; top %.3g%% => p <= %.3g, %d SNPs\n", nrow(d), 100*TOP, thr, nrow(sig)))
sig[, newgrp := (Pos - shift(Pos, fill=-1000000000L) > MERGE) | (Chr != shift(Chr, fill=""))]
sig[, grp := cumsum(newgrp)]
pk <- sig[, .(start=min(Pos)-1L, end=max(Pos), n_top=.N, min_p=min(p)), by=.(Chr,grp)][n_top>=MINSIG]
pk[, width := end-start]
cat(sprintf("  -> %d suggestive peaks, %.2f Mb, median width %.1f kb  [NOT FDR-controlled]\n",
            nrow(pk), sum(pk$width)/1e6, median(pk$width)/1e3))
fwrite(pk[,.(Chr,start,end,n_top,min_p)],
       file.path(DATA,"ecopeaks","c151_nEur_suggestivePeaks.bed"), sep="\t", col.names=FALSE)
cat(sprintf("  wrote %s/ecopeaks/c151_nEur_suggestivePeaks.bed\n", DATA))
