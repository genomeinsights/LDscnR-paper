## Do the c151 (Northern Europe) suggestive peaks land inside the 3sp outlier regions
## more often than chance? Complements the rank-based test in 09_overlap.R by working on
## the thresholded peak set instead of the p-value ranks.
## Writes data/c151_peak_overlap.csv
suppressMessages(library(data.table))
P  <- path.expand("~/gitlab/LDscnR-paper/kingman2021")
EP <- path.expand("~/gitlab/LD-scaling-genome-scans/empirical_data/kingman2021/ecopeaks")
pk <- fread(file.path(P,"data","peaks","c151_nEur_suggestivePeaks.bed"), header=FALSE,
            col.names=c("Chr","start","end","n_top","min_p"))
R  <- fread(file.path(P,"data","liftover","lfmm_g14.bed"), header=FALSE,
            col.names=c("Chr","start","end","name"))
d  <- fread(file.path(EP,"c151_nEur.snp_p.tsv.gz"), select=c("Chr","Pos"))
rng <- d[, .(lo=min(Pos), hi=max(Pos)), by=Chr]; setkey(rng, Chr)
setkey(R, Chr, start, end)
hits <- function(X) nrow(foverlaps(X, R, by.x=c("Chr","start","end"), type="any", nomatch=NULL))
obs <- hits(pk); pk[, w := end-start]
set.seed(1); B <- 5000L; nl <- integer(B); r <- rng[J(pk$Chr)]
for (b in seq_len(B)) {
  st <- r$lo + floor(runif(nrow(pk))*pmax(1, r$hi-r$lo-pk$w))
  nl[b] <- hits(data.table(Chr=pk$Chr, start=st, end=st+pk$w))
}
res <- data.table(n_peaks=nrow(pk), n_inside=obs, null_mean=mean(nl),
                  fold=obs/mean(nl), pval=(1+sum(nl>=obs))/(B+1))
fwrite(res, file.path(P,"data","c151_peak_overlap.csv")); print(res)
