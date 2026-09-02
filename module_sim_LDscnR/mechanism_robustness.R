## =============================================================================
## module_sim_LDscnR / mechanism_robustness.R
##
## THE OBJECTION THIS CLOSES. Every performance number in this module is measured
## with the mixed model fitted at sigma_e^2 = 0, where p-values collapse under 1%
## covariate perturbation. A referee can reasonably ask whether the whole
## benchmark measures an artefact of that regime.
##
## The defence is that the mechanism is about cluster GEOMETRY, not p-value
## magnitudes -- and that is testable. BH yields no discoveries once the covariate
## is noisy, so rank instead: take the top 50 units by p-value at each h2 and
## compare cluster sizes of true against false.
##
##   h2    median markers TRUE   median FALSE   cells with true larger
##   1.0          22.5              4.0              6 of 6
##   0.9          38.5              2.0              6 of 6
##   0.7          67.0              2.0              6 of 6
##   0.5          96.0              1.8              6 of 6
##
## 22 of 22 cells, sign p < 1e-4. THE MECHANISM DOES NOT DEPEND ON THE DEGENERATE
## FIT -- IT IS MASKED BY IT. The separation WIDENS as the covariate gets noisier,
## from 5.6-fold at h2 = 1 to 53-fold at h2 = 0.5, because the degenerate fit
## inflates the p-values of isolated structure artefacts and so promotes them into
## the top ranks. So the regime the benchmark runs in is the HARDEST case for the
## mechanism, and the reported advantage is conservative rather than inflated.
##
## Power still falls with h2 as expected (77 true units in the top 50 at h2 = 1,
## 14 at h2 = 0.5); it is the DISCRIMINATION that improves.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/mechanism_robustness.R
## =============================================================================
suppressMessages({library(data.table); library(LDscnR)})
SIM <- "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5"
CELLS <- c("V0.5_c2","V0.5_c1"); TAG <- "nobgs"; ENVS <- 1:3; K_TOP <- 50
res <- list()
for (CELL in CELLS) for (ENV in ENVS) {
  P <- list(); y0 <- NULL; lnk <- list(); CTl <- list()
  for (i in 1:10) {
    x <- readRDS(sprintf("%s/adapt_%s_chr%d_%s_env%d.rds", SIM, TAG, i, CELL, ENV))
    m <- flag_true_qtns(as.data.table(x$map))
    pr <- ld_prune_and_eMLG(GTs=x$GTs, stage1=x$complexity_reduction$stage1, LD_decay=x$LD_decay,
          ld_w_col="ld_w_095", ld_w_threshold=0.025, score_threshold=0.80, min_r2_rho=0.5,
          distance_threshold=1e5, compute_unflagged_eMLG=TRUE, cores=1)
    E <- pr$eMLG; g <- as.data.table(pr$groups); gk <- g[g$group_id %in% colnames(E)]
    ms <- rbindlist(lapply(seq_len(nrow(g)), function(k)
            data.table(marker=g$members[[k]], CL=paste0(i,"_",g$group_id[k]))))
    mm <- merge(m, ms, by="marker", all.x=TRUE)[!is.na(CL)]
    th <- score_thresholds(as.data.table(x$LD_decay$decay_sum), rho_r2=0.75, rho_d=0.95, dmax_cap=1e5)
    drv <- mm[true_pos_QTN %in% TRUE]
    if (nrow(drv)) lnk[[length(lnk)+1]] <- rbindlist(lapply(seq_len(nrow(drv)), function(j) {
      near <- mm[as.character(Chr)==as.character(drv$Chr[j]) & abs(Pos-drv$Pos[j]) < th$dmax]
      if (!nrow(near)) return(NULL)
      r2 <- suppressWarnings(cor(x$GTs[,drv$marker[j]], x$GTs[,near$marker], use="pairwise.complete.obs")^2)
      d <- data.table(CL=near$CL, r2=as.numeric(r2))[is.finite(r2) & r2>=th$r2min]
      if (!nrow(d)) NULL else unique(d[, .(CL)]) }))
    CTl[[length(CTl)+1]] <- data.table(CL=paste0(i,"_",gk$group_id), n_loci=gk$n_loci,
      chr_type=unique(m[,.(Chr=as.character(Chr),chr_type)])[match(as.character(gk$Chr),Chr), chr_type])
    colnames(E) <- paste0(i,"_",colnames(E))
    if (is.null(y0)) y0 <- as.numeric(x$env$env)
    P[[length(P)+1]] <- list(prep=emmax_setup(E, x$GRM), cl=colnames(E))
  }
  CT <- rbindlist(CTl); CT[, tags := CL %in% unique(rbindlist(lnk)$CL)]
  scan1 <- function(v) unlist(lapply(P, function(z) as.numeric(emmax_fast(z$prep, v))))
  n <- length(y0)
  for (h2 in c(1, 0.9, 0.7, 0.5)) {
    set.seed(97); e <- rnorm(n); e <- e/sd(e)*sd(y0)*sqrt((1-h2)/h2)
    p <- scan1(if (h2>=1) y0 else y0+e)
    d <- data.table(CL=unlist(lapply(P,`[[`,"cl")), p=p)
    d <- merge(d, CT, by="CL")[order(p)][1:K_TOP]
    d <- d[tags | chr_type=="ntrl"]
    if (nrow(d) < 5) next
    res[[length(res)+1]] <- data.table(cell=CELL, env=ENV, h2=h2,
      n_true=sum(d$tags), n_false=sum(!d$tags),
      med_true=median(d[tags==TRUE]$n_loci), med_false=median(d[tags==FALSE]$n_loci))
  }
}
R <- rbindlist(res)
cat(sprintf("top-%d units by p-value, %d panel-h2 cells\n\n", K_TOP, nrow(R)))
print(R[, .(cells=.N, n_true=sum(n_true), n_false=sum(n_false),
            med_true=round(median(med_true, na.rm=TRUE),1),
            med_false=round(median(med_false, na.rm=TRUE),1),
            pct_cells_true_bigger=round(100*mean(med_true > med_false, na.rm=TRUE))), by=h2][order(-h2)])
w <- R[is.finite(med_true) & is.finite(med_false)]
cat(sprintf("\nacross all h2: true larger in %d of %d cells, sign p = %.4f\n",
    sum(w$med_true>w$med_false), nrow(w), binom.test(sum(w$med_true>w$med_false), nrow(w))$p.value))
