## =============================================================================
## stratified_paired_table.R -- the paired result, stratified, with magnitudes.
##
## Written in response to an external audit of doc/filter_then_test.pdf. Four
## defects it addresses, all valid:
##
##  1. WIN-RATE DENOMINATOR. winrate() conditioned on non-tied panels while the
##     caption said "percent of 80 panels". With 23% of panels carrying no
##     signal, ties are common and the two denominators differ materially.
##     Wins / losses / ties are now all reported.
##  2. MAGNITUDE. Signs alone treat a win by 1 and a win by 20 alike. The paired
##     distribution -- median change with IQR -- is reported alongside the sign
##     test, not instead of it.
##  3. POOLING. 80 panels span 4 parameter cells x 2 BGS arms. Stratified here.
##  4. ENDPOINT NAME. "true discovery" overstates bp proximity with causal
##     variants excluded. Called QTN-proximal throughout.
##
## Output feeds the document directly; nothing is transcribed by hand.
## =============================================================================
suppressMessages({library(data.table)})
R   <- Sys.getenv("RESULTS", "module_sim_LDscnR/results/filter_then_test")
OUT <- Sys.getenv("OUT", R)
IN  <- Sys.getenv("INFILE", "filter_then_test_clusters_emx.csv")
d   <- fread(file.path(R, IN))
key <- c("cell","tag","env","window_kb")
gw  <- d[method=="genome_wide", c(key,"n_tp","n_sig"), with=FALSE]
setnames(gw, c("n_tp","n_sig"), c("gw_tp","gw_sig"))
f   <- merge(d[method!="genome_wide"], gw, by=key)
f[, `:=`(dtp = n_tp - gw_tp, dsig = n_sig - gw_sig)]

## paired summary with ties made explicit
paired <- function(x) {
  x <- x[is.finite(x)]
  w <- sum(x > 0); l <- sum(x < 0); t <- sum(x == 0); nz <- w + l
  list(n = length(x), wins = w, losses = l, ties = t,
       win_pct_nonties = if (nz) 100*w/nz else NA_real_,
       win_pct_all = 100*w/max(length(x),1),
       med = median(x), q25 = quantile(x, .25), q75 = quantile(x, .75),
       p = if (nz) binom.test(w, nz)$p.value else NA_real_)
}

cat("=== A. POOLED, with ties explicit (50 kb window) ===\n")
a <- f[window_kb==50, paired(dtp), by=.(method,k)]
a[, `:=`(med = round(med,1), win_pct_nonties = round(win_pct_nonties), win_pct_all = round(win_pct_all),
         p = signif(p,3))]
print(a[method %in% c("ld_w","size"), .(method,k,n,wins,losses,ties,win_pct_nonties,win_pct_all,med,q25,q75,p)][order(method,k)])

cat("\n=== B. MAGNITUDE: median paired change in QTN-proximal discoveries [IQR], 50 kb ===\n")
b <- f[window_kb==50, .(med_dtp = median(dtp), q25 = quantile(dtp,.25), q75 = quantile(dtp,.75),
                        med_dsig = median(dsig)), by=.(method,k)]
print(dcast(b, k ~ method, value.var="med_dtp")[order(k)])
cat("\n  change in TOTAL discoveries (median), same strata:\n")
print(dcast(b, k ~ method, value.var="med_dsig")[order(k)])

cat("\n=== C. STRATIFIED BY BGS ARM (ld_w, 50 kb) ===\n")
cc <- f[window_kb==50 & method=="ld_w", paired(dtp), by=.(tag,k)]
print(cc[, .(tag,k,n,wins,losses,ties,win_pct_nonties=round(win_pct_nonties),med=round(med,1),p=signif(p,3))][order(tag,k)])

cat("\n=== D. STRATIFIED BY CELL x BGS (ld_w, 50 kb, k=5000) ===\n")
dd <- f[window_kb==50 & method=="ld_w" & k==5000, paired(dtp), by=.(cell,tag)]
print(dd[, .(cell,tag,n,wins,losses,ties,win_pct_nonties=round(win_pct_nonties),med=round(med,1),p=signif(p,3))][order(cell,tag)])

cat("\n=== E. SAME, cluster-size filter ===\n")
ee <- f[window_kb==50 & method=="size" & k==5000, paired(dtp), by=.(cell,tag)]
print(ee[, .(cell,tag,n,wins,losses,ties,win_pct_nonties=round(win_pct_nonties),med=round(med,1),p=signif(p,3))][order(cell,tag)])

cat("\n=== F. ld_w vs SIZE head to head, by cell (50 kb, all k) ===\n")
h <- merge(f[method=="ld_w",  c(key,"k","n_tp"), with=FALSE],
           f[method=="size", c(key,"k","n_tp"), with=FALSE],
           by=c(key,"k"), suffixes=c("_l","_s"))
hh <- h[window_kb==50, paired(n_tp_l - n_tp_s), by=.(cell,k)]
print(hh[, .(cell,k,wins,losses,ties,win_pct_nonties=round(win_pct_nonties),med=round(med,1),p=signif(p,3))][order(cell,k)])

fwrite(rbindlist(list(
  cbind(stratum="pooled",    a[, .(method,k,n,wins,losses,ties,win_pct_nonties,med,q25,q75,p)]),
  cbind(stratum="by_arm",    cc[, .(method="ld_w",k,n,wins,losses,ties,win_pct_nonties,med,q25,q75,p)]),
  cbind(stratum="cell_x_arm",dd[, .(method="ld_w",k=5000,n,wins,losses,ties,win_pct_nonties,med,q25,q75,p)])
), fill=TRUE), file.path(OUT, "stratified_paired.csv"))
cat(sprintf("\n  written: %s\n", file.path(OUT, "stratified_paired.csv")))
