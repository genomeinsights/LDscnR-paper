## =============================================================================
## module_sim_LDscnR / isolation_mechanism.R
##
## THE PAPER'S CENTRAL MECHANISM, TESTED DIRECTLY. The claim is that residual
## structure confounding produces ISOLATED false positives -- the kind of single
## unclustered hit a SNP scan reports -- while real adaptive loci produce
## CLUSTERED signal. If so, false positives should sit in smaller clusters.
## Neutral-chromosome discoveries are unambiguously false, QTN-tagging ones true,
## so the comparison needs no convention.
##
## CONFIRMED, on 2,547 discoveries over 40 panels:
##
##   route  class    n     median markers   % singleton   % >= 5 markers
##   cons   FALSE   310         14              15.5          69.0
##   cons   TRUE    205         28               4.4          87.8
##   best   FALSE  1057         26               8.5          78.2
##   best   TRUE    378         39               3.2          87.6
##   simes  FALSE   362         20              14.4          72.4
##   simes  TRUE    235         34               3.8          86.0
##
## Wilcoxon p = 4.5e-06, 1.3e-07, 1.3e-04. True positives sit in clusters about
## TWICE the size, and false positives are 3-4x more likely to be singletons.
##
## AND IT RECONCILES A RESULT THAT LOOKED LIKE A CONTRADICTION. The minimum
## cluster size floor is INERT above 2 (floor_sweep.R: p = 0.11, 0.40, 0.053).
## If false positives were isolated, a floor should remove them preferentially.
## It does not, because the distributions OVERLAP HEAVILY: only 15.5% of false
## positives are singletons and 69% still sit in clusters of five or more, while
## 12% of true positives fall below five. A threshold on size therefore trades
## one for the other and nets nothing.
##
## SO THE MECHANISM IS REAL BUT IT IS NOT A FILTER, and that distinction is the
## paper's point rather than a caveat on it. The value comes from CHANGING THE
## UNIT, not from filtering units: aggregation makes an isolated false positive
## contribute a single unit instead of surviving as a marker-level hit, and
## consolidates a clustered true positive into one. Nothing is thresholded.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/isolation_mechanism.R
## =============================================================================
suppressMessages(library(data.table))
D  <- "module_sim_LDscnR/results/structure_null"
CT <- as.data.table(readRDS(file.path(D, "cl_chrtype.rds")))
PC <- c(cons="p_cons", best="p_best", simes="p_simes")
out <- list()
for (f in list.files(D, pattern="spatial", full.names=TRUE)) {
  z <- readRDS(f); U <- z$units; k <- z$summary[1]
  ct <- CT[cell==k$cell & tag==k$tag & env==k$env, .(CL, chr_type)]
  U <- merge(U, ct, by="CL", all.x=TRUE)
  U[, tags_qtn := CL %in% unique(z$links$CL)]
  for (rt in names(PC)) {
    idx <- which(is.finite(U[[PC[rt]]]))
    sig <- idx[p.adjust(U[[PC[rt]]][idx], "BH") < 0.05]
    if (!length(sig)) next
    s <- U[sig]
    out[[length(out)+1]] <- data.table(route=rt, cell=k$cell, env=k$env, tag=k$tag,
      CL=s$CL, n_loci=s$n_loci,
      class=fifelse(s$tags_qtn, "TRUE (tags a QTN)",
              fifelse(s$chr_type=="ntrl", "FALSE (neutral chr)", "other")))
  }
}
R <- rbindlist(out)[class != "other"]
cat(sprintf("%d discoveries classified across %d panels\n\n", nrow(R), uniqueN(R[,.(cell,tag,env)])))
cat("== cluster size of TRUE vs FALSE discoveries\n")
print(R[, .(n=.N, median_markers=as.double(median(n_loci)), q25=as.double(quantile(n_loci,.25)),
            q75=as.double(quantile(n_loci,.75)), pct_singleton=round(100*mean(n_loci==1),1),
            pct_ge5=round(100*mean(n_loci>=5),1)), by=.(route, class)][order(route,class)])
for (rt in unique(R$route)) {
  a <- R[route==rt & class=="TRUE (tags a QTN)", n_loci]
  b <- R[route==rt & class=="FALSE (neutral chr)", n_loci]
  if (length(a) > 2 && length(b) > 2)
    cat(sprintf("\n  %-6s Wilcoxon TRUE vs FALSE: p = %.3g | median %g vs %g\n",
        rt, wilcox.test(a, b)$p.value, median(a), median(b)))
}
