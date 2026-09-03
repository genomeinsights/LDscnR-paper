## Precision, recall, and how the ranking depends on which you weight.
## F_beta = (1+b^2) P R / (b^2 P + R). b < 1 favours precision, b > 1 recall,
## b = 1 is the harmonic mean. The product P*R used elsewhere in this module is
## a third weighting again, and none of them is neutral -- the choice IS the
## claim, so the crossover is reported rather than a single number.
suppressMessages({library(data.table)})
d <- fread("module_sim_LDscnR/results/filter_then_test/snp_vs_cluster_dedup_allpanels.csv")
lab <- c(snp="single SNP", rep="representative", simes="Simes", emlg="eMLG")
d <- d[analysis %in% names(lab)]; d[, an := factor(lab[analysis], levels=lab)]
d[, `:=`(P = dedup_prec, R = dedup_rec)]
cat("=== precision and recall, medians ===\n")
print(dcast(d, an ~ tag, value.var=c("P","R"), fun.aggregate=function(z) round(median(z,na.rm=TRUE),3)))
cat("\n=== paired: is Simes' precision higher than single-SNP's? ===\n")
w <- dcast(d, cell+tag+env ~ an, value.var="P")
for (tg in c("nobgs","bgs")) {
  z <- w[tag==tg]; x <- z$Simes - z$`single SNP`; x <- x[is.finite(x)]; nz <- x[x!=0]
  cat(sprintf("  %-6s %2d W / %2d L / %2d T   median %+.3f   p %.3g\n", tg,
      sum(nz>0), sum(nz<0), sum(x==0), median(x), binom.test(sum(nz>0),length(nz))$p.value))
}
fb <- function(P, R, b) { v <- (1+b^2)*P*R / (b^2*P + R); v[!is.finite(v)] <- 0; v }
BET <- c(0.25, 0.5, 1, 1.5, 2, 3, 4)
cat("\n=== F_beta: Simes vs single-SNP, paired over 80 panels ===\n")
cat(sprintf("  %-6s %-22s %-22s %s\n", "beta", "median F (single SNP)", "median F (Simes)", "paired sign test"))
for (b in BET) {
  d[, f := fb(P, R, b)]
  v <- dcast(d, cell+tag+env ~ an, value.var="f")
  x <- v$Simes - v$`single SNP`; x <- x[is.finite(x)]; nz <- x[x!=0]
  cat(sprintf("  %-6.2f %-22.4f %-22.4f %2d W / %2d L  p %.4g%s\n", b,
      median(d[an=="single SNP"]$f, na.rm=TRUE), median(d[an=="Simes"]$f, na.rm=TRUE),
      sum(nz>0), sum(nz<0), binom.test(sum(nz>0),length(nz))$p.value,
      if (binom.test(sum(nz>0),length(nz))$p.value < 0.05) " *" else ""))
}
cat("\n  beta < 1 weights PRECISION more; beta > 1 weights RECALL more.\n")
cat("  The crossover is where the two analyses become indistinguishable.\n")
cat("\n=== same, nobgs only ===\n")
for (b in BET) {
  d[, f := fb(P, R, b)]
  v <- dcast(d[tag=="nobgs"], cell+env ~ an, value.var="f")
  x <- v$Simes - v$`single SNP`; x <- x[is.finite(x)]; nz <- x[x!=0]
  cat(sprintf("  beta %-5.2f  %2d W / %2d L   p %.4g\n", b, sum(nz>0), sum(nz<0),
      binom.test(sum(nz>0),length(nz))$p.value))
}
