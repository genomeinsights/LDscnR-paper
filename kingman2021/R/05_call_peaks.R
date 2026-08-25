## Call EcoPeaks from the SNP-based multivariate-hypergeometric p-values (04) and,
## when the cohort is c155_global, validate against the published peak set.
## Usage: Rscript 05_call_peaks.R <cohort> [fdr] [merge_kb] [min_sig]
suppressMessages(library(data.table))
a <- commandArgs(trailingOnly=TRUE)
COH <- if (length(a)>=1) a[1] else "c155_global"
FDR <- if (length(a)>=2) as.numeric(a[2]) else 0.01
MERGE <- (if (length(a)>=3) as.numeric(a[3]) else 50) * 1000
MINSIG <- if (length(a)>=4) as.integer(a[4]) else 2L
DATA <- "~/gitlab/LD-scaling-genome-scans/empirical_data/kingman2021"

d <- fread(file.path(path.expand(DATA), "ecopeaks", paste0(COH,".snp_p.tsv.gz")))
d[, q := p.adjust(p, "BH")]
sig <- d[q <= FDR]
cat(sprintf("%s: %d sites tested, %d significant at BH q<=%.3g (%.3f%%)\n",
            COH, nrow(d), nrow(sig), FDR, 100*nrow(sig)/nrow(d)))
if (!nrow(sig)) quit(save="no")

setorder(sig, Chr, Pos)
sig[, newgrp := (Pos - shift(Pos, fill=-1000000000L) > MERGE) | (Chr != shift(Chr, fill="")), by=NULL]
sig[, grp := cumsum(newgrp)]
pk <- sig[, .(start=min(Pos)-1L, end=max(Pos), n_sig=.N, min_p=min(p)), by=.(Chr, grp)][n_sig >= MINSIG]
pk[, width := end-start]
cat(sprintf("  -> %d peaks, %.2f Mb total, median width %.1f kb\n",
            nrow(pk), sum(pk$width)/1e6, median(pk$width)/1e3))
fwrite(pk[, .(Chr, start, end, n_sig, min_p)],
       file.path(path.expand(DATA), "ecopeaks", sprintf("%s_snpEcoPeaks_fdr%.3g.bed", COH, FDR)),
       sep="\t", col.names=FALSE)

## ---- validation against the published c155 peaks ----
if (COH == "c155_global") {
  TR <- file.path(path.expand(DATA), "tracks")
  setkey(pk, Chr, start, end)
  for (s in c("c155.specific","c155.sensitive")) {
    P <- fread(file.path(TR, sprintf("gasAcu1-4.%s.50kb.final.peaks.bed", s)),
               header=FALSE, select=1:3, col.names=c("Chr","start","end"))
    setkey(P, Chr, start, end)
    o <- foverlaps(P, pk, type="any", nomatch=NULL)
    rec <- uniqueN(o[, .(Chr, i.start)])
    o2 <- foverlaps(pk, P, type="any", nomatch=NULL)
    prec <- uniqueN(o2[, .(Chr, i.start)])
    cat(sprintf("  vs published %-15s : %d/%d published peaks recovered (%.0f%%); %d/%d of my peaks hit one (%.0f%%)\n",
                s, rec, nrow(P), 100*rec/nrow(P), prec, nrow(pk), 100*prec/nrow(pk)))
  }
}
