## Is the 3sp LFMM outlier-region set enriched for low Kingman marine-vs-freshwater
## p-values, and does the enrichment track the GEOGRAPHY of the Kingman cohort?
## Rank-based, so it works for c151 too, where no SNP can reach genome-wide significance.
## Usage: Rscript 07_signal_enrichment.R <cohort> [B]
suppressMessages(library(data.table))
a <- commandArgs(trailingOnly=TRUE)
COH <- if (length(a)>=1) a[1] else "c151_nEur"
B   <- if (length(a)>=2) as.integer(a[2]) else 1000L
DATA <- path.expand("~/gitlab/LD-scaling-genome-scans/empirical_data/kingman2021")
CH <- "/private/tmp/claude-539526166/-Users-petrikem-gitlab-LDscnR/e0f3b9b5-e5d3-4cab-8211-c3f2f61c09ff/scratchpad/chain"

d <- fread(file.path(DATA,"ecopeaks",paste0(COH,".snp_p.tsv.gz")), select=c("Chr","Pos","p"))
R <- fread(file.path(CH,"lfmm_g14.bed"), header=FALSE, col.names=c("Chr","start","end","name"))
R <- R[Chr %in% d$Chr]
setkey(d, Chr, Pos)
rng <- d[, .(lo=min(Pos), hi=max(Pos), n=.N), by=Chr]; setkey(rng, Chr)
R[, span := end-start]

mark <- function(RR) {
  setkey(RR, Chr, start, end)
  o <- foverlaps(d[, .(Chr, start=Pos, end=Pos)], RR, by.x=c("Chr","start","end"), type="any", mult="first", nomatch=NA)
  !is.na(o$name)
}
inreg <- mark(copy(R))
cat(sprintf("\n=== %s : %d SNPs tested, %d regions covering %d SNPs (%.2f%%)\n",
            COH, nrow(d), nrow(R), sum(inreg), 100*mean(inreg)))
for (TOP in c(0.001, 0.01)) {
  thr <- quantile(d$p, TOP)
  obs <- mean(d$p[inreg] <= thr)
  nullv <- numeric(B)
  r <- rng[J(R$Chr)]
  for (b in seq_len(B)) {
    st <- r$lo + floor(runif(nrow(R)) * pmax(1, r$hi - r$lo - R$span))
    m <- mark(data.table(Chr=R$Chr, start=st, end=st+R$span, name="x"))
    nullv[b] <- mean(d$p[m] <= thr)
  }
  cat(sprintf("  top %.1f%% of p (p<=%.2e): in-region rate %.4f  vs null %.4f  -> %.2fx   p=%.4f\n",
      100*TOP, thr, obs, mean(nullv), obs/mean(nullv), (1+sum(nullv>=obs))/(B+1)))
}
w <- suppressWarnings(wilcox.test(-log10(d$p[inreg]), -log10(d$p[!inreg]), alternative="greater"))
cat(sprintf("  median -log10(p): in-region %.3f  vs background %.3f   (Wilcoxon p=%.3g)\n",
    median(-log10(d$p[inreg])), median(-log10(d$p[!inreg])), w$p.value))
