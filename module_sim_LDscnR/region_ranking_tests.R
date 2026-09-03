## =====================================================================
## module_sim_LDscnR / region_ranking_tests.R
##
## Once regions are called, what should rank them? Tests s_R (summed C-mass),
## region size, mean C per marker, and C-squared against the truth.
##
## THE ANSWER IS THAT s_R IS SIZE. s_R sums C over member markers, so it is
## size-inflated by construction: Spearman(s_R, size) = +0.84, and once
## log(s_R) is residualised on log(size) the ranking AUC falls from 0.923 to
## 0.582, i.e. to near chance. At matched budget s_R and size are identical in
## 60 of 64 comparisons. Ranking by s_R is a size ordering and should be
## described as one, not presented as an independent statistic.
##
## mean_C (= s_R / size) IS independent of size and beats chance (AUC 0.804,
## 14 of 16 datasets) -- but ranks WORSE than size does, so it is informative
## rather than useful.
##
## The stratum table shows why: QTN content is almost entirely a size effect.
## 2 of 216 regions with 3-5 markers contain a QTN, against 42 of 52 with >50.
## Within a stratum s_R ranks well (0.81-0.84); the strata do the work.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/region_ranking_tests.R
## =====================================================================
## Does s_R carry information about QTN content BEYOND region size?
## s_R sums C over member markers, so it is size-inflated by construction.
## Three size-controlled views:
##   mean_C   = s_R / size, the per-marker version -- size divided out entirely
##   resid    = residual of log(s_R) on log(size), pooled within dataset
##   stratum  = rank by s_R only against regions of comparable size
suppressMessages(library(data.table))
D <- "module_sim_LDscnR/results/nulls_V2_c1"
fs <- list.files(D, pattern="^regions_.*[.]csv$", full.names=TRUE)
d <- rbindlist(lapply(fs, function(f) { r <- fread(f)
  r[, `:=`(env=as.integer(sub(".*_env([0-9]+)_.*","\\1",basename(f))),
           basis=sub(".*_emmax_(.+)_B[0-9]+\\.csv$","\\1",basename(f)))] }), fill=TRUE)
d <- d[!is.na(has_qtn) & size > 0 & s_R > 0]
d[, mean_C := s_R/size]
d[, resid := residuals(lm(log(s_R) ~ log(size))), by=.(env,basis)]
cat(sprintf("  %d regions, %d datasets, %d contain a QTN (%.1f%%)\n",
            nrow(d), uniqueN(d[,.(env,basis)]), sum(d$has_qtn), 100*mean(d$has_qtn)))
cat(sprintf("  Spearman(s_R, size) = %+.3f  <- how size-driven s_R is\n\n",
            cor(d$s_R, d$size, method="spearman")))

auc <- function(sc, pos) { ok <- is.finite(sc); sc <- sc[ok]; pos <- pos[ok]
  if (!any(pos) || all(pos)) return(NA_real_)
  r <- rank(sc); n1 <- sum(pos); (sum(r[pos]) - n1*(n1+1)/2)/(n1*sum(!pos)) }

cat("=== ranking AUC per dataset, then averaged (0.5 = chance) ===\n")
a <- d[, .(size=auc(size,has_qtn), s_R=auc(s_R,has_qtn),
           mean_C=auc(mean_C,has_qtn), resid=auc(resid,has_qtn)), by=.(env,basis)]
cat(sprintf("  size    %.3f  (n=%d usable)\n", mean(a$size,na.rm=TRUE), sum(!is.na(a$size))))
cat(sprintf("  s_R     %.3f\n", mean(a$s_R,na.rm=TRUE)))
cat(sprintf("  mean_C  %.3f   <- s_R with size divided out\n", mean(a$mean_C,na.rm=TRUE)))
cat(sprintf("  resid   %.3f   <- s_R residualised on size\n", mean(a$resid,na.rm=TRUE)))
cat(sprintf("\n  s_R beats size in %d of %d datasets; mean_C beats chance in %d\n",
            sum(a$s_R > a$size, na.rm=TRUE), sum(!is.na(a$s_R)), sum(a$mean_C > 0.5, na.rm=TRUE)))

cat("\n=== within size strata: does s_R rank inside a stratum? ===\n")
d[, stratum := cut(size, breaks=c(0,5,15,50,Inf), labels=c("3-5","6-15","16-50",">50"))]
st <- d[, .(n=.N, with_QTN=sum(has_qtn),
            auc_sR=round(auc(s_R,has_qtn),3)), by=stratum][order(stratum)]
print(st)
cat("\n=== matched budget: top N by each criterion ===\n")
topN <- function(x, sc, N) { o <- order(-sc); sum(x[o][seq_len(min(N,length(o)))]) }
b <- rbindlist(lapply(c(3,5,8,12), function(N)
  d[, .(budget=N, by_size=topN(has_qtn,size,N), by_s_R=topN(has_qtn,s_R,N),
        by_meanC=topN(has_qtn,mean_C,N)), by=.(env,basis)]))
print(b[, .(n=.N, size=round(mean(by_size),2), s_R=round(mean(by_s_R),2),
            mean_C=round(mean(by_meanC),2),
            sR_gt_size=sum(by_s_R>by_size), ties=sum(by_s_R==by_size),
            size_gt_sR=sum(by_size>by_s_R)), by=budget][order(budget)])
