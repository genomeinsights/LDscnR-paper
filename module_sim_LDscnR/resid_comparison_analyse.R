## =============================================================================
## module_sim_LDscnR / resid_comparison_analyse.R
##
## Residualising surrogates against the observed covariate: 3 bases x 2 settings
## x 20 panels at NSIM = 500. THIS REVERSES WHAT I CONCLUDED FROM SINGLE PANELS
## AT NSIM = 60 AND CIRCULATED TO THE PANEL SESSION.
##
##   basis    resid   E[V]     locking (x chance)   % draws with any discovery
##   vc       TRUE     0.12            0x                    5.9
##   vc       FALSE    0.13            3x                    6.2
##   pop      TRUE     0.86            1x                   18.9
##   pop      FALSE    0.89           90x  (SE 36)          19.7
##   spatial  TRUE   110.89          224x  (SE 69)          75.3
##   spatial  FALSE  243.32          681x  (SE 152)         94.1
##
## RESIDUALISING IS PROTECTIVE AGAINST REGION LOCKING, AND FOR THE `pop` BASIS IT
## IS THE DIFFERENCE BETWEEN A VALID NULL AND AN INVALID ONE -- 1x chance against
## 90x. I had reported the opposite ("orthogonal is not dissimilar"; that it cuts
## hit volume without reducing the overlap fraction), from a two-panel check at
## NSIM = 80 where `pop` produced 5-19 surrogate hits in total. The overlap
## fraction was estimated from single-digit counts and was noise.
##
## THE ORIGINAL CLAIM SURVIVES ONLY FOR `spatial`, where residualising helps 3x
## and leaves the null catastrophically locked either way. So the accurate
## statement is narrower: ORTHOGONALISING SUBSTANTIALLY REDUCES LOCKING BUT DOES
## NOT RESCUE A BASIS THAT IS BADLY LOCKED TO BEGIN WITH. Generalising the
## spatial result to all bases was the error.
##
## AND IT IS NEARLY FREE. Paired within panel, dropping residualisation raises
## E[V] by 3% for `pop` (0.861 -> 0.886, 15 of 20 panels, sign p = 0.041) and not
## detectably for `vc` (p = 0.82). So residualisation costs almost no power and
## buys two orders of magnitude of validity on the one basis where it matters.
## That cuts against the Joo prescription the panel session is following -- which
## avoids residualising precisely because it removes the confounded axis -- and
## the trade is measurable rather than a matter of principle.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/resid_comparison_analyse.R
## =============================================================================
D <- "module_sim_LDscnR/results/structure_null/resid_comparison"
S <- rbindlist(lapply(list.files(D, full.names=TRUE), function(f) {
  z <- readRDS(f); b <- basename(f)
  z$summary[, `:=`(basis = sub(".*simes_([a-z]+)_r[01]_.*","\\1", b),
                   resid = grepl("_r1_", b))][] }))
fwrite(S, "module_sim_LDscnR/results/structure_null/resid_comparison_summary.csv")
cat(sprintf("%d panels x %d cells | 20 panels expected per basis-resid\n\n",
            uniqueN(S[,.(cell,tag,env)]), uniqueN(S[,.(basis,resid,route,floor)])))
cat("== E[V] and locking, consensus route, floor 1, mean over 20 panels (SE)\n")
A <- S[route=="consensus" & floor==1, .(
  panels=.N, obs_R=round(mean(obs_R),1), EV=round(mean(EV),2),
  EV_se=round(sd(EV)/sqrt(.N),2),
  lock_x=round(mean(lock_overlap/pmax(lock_chance,1e-12), na.rm=TRUE)),
  lock_se=round(sd(lock_overlap/pmax(lock_chance,1e-12), na.rm=TRUE)/sqrt(.N)),
  pct_any=round(mean(pct_any),1)), by=.(basis,resid)]
print(A[order(basis,-resid)])
cat("\n== does dropping residualisation change E[V]? paired within panel\n")
W <- dcast(S[route=="consensus" & floor==1], basis+cell+tag+env ~ resid, value.var="EV")
setnames(W, c("FALSE","TRUE"), c("no_resid","resid"))
print(W[, .(panels=.N, mean_resid=round(mean(resid),3), mean_no_resid=round(mean(no_resid),3),
            n_higher_without=sum(no_resid>resid),
            sign_p=round(binom.test(sum(no_resid>resid), .N)$p.value,4)), by=basis])
