## Does selection elevate local LD -- and does background selection abolish it?
##
## The premise that selection elevates LD locally underlies a large class of
## LD-only outlier methods that use no association test at all. If it holds
## without background selection and fails with it, that bears on those methods
## directly.
##
## Test: the QTN enrichment of each cluster geometry class, computed SEPARATELY
## per BGS arm, under both denominators -- per unit (which units to examine) and
## per marker (is causal variation concentrated here). Classes are defined partly
## by density, so the marker denominator partly defines the effect away; both are
## reported because they answer different questions.
suppressMessages({library(data.table); library(LDscnR)})
SIM  <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
OUT  <- Sys.getenv("OUT", "module_sim_LDscnR/results/filter_then_test")
ENVS <- as.integer(strsplit(Sys.getenv("ENVS", "1,2,3,4,5"), ",")[[1]])
CELL <- Sys.getenv("CELL", "V0.5_c1")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
acc <- list()
for (tg in c("nobgs","bgs")) for (ev in ENVS) for (i in 1:10) {
  f <- sprintf("%s/adapt_%s_chr%d_%s_env%d.rds", SIM, tg, i, CELL, ev)
  if (!file.exists(f)) next
  x <- readRDS(f); m <- flag_true_qtns(as.data.table(x$map))
  pr <- ld_prune_and_eMLG(GTs=x$GTs, stage1=x$complexity_reduction$stage1, LD_decay=x$LD_decay,
        ld_w_col="ld_w_095", ld_w_threshold=0.025, score_threshold=0.80, min_r2_rho=0.5,
        distance_threshold=1e5, compute_unflagged_eMLG=FALSE, cores=1)
  g <- as.data.table(pr$groups)
  ms <- data.table(marker = unlist(g$members, use.names = FALSE),
            CL = paste0(tg,ev,i,"_",rep.int(g$group_id, lengths(g$members))))
  mm <- merge(m, ms, by="marker", all.x=TRUE)[!is.na(CL)]
  acc[[length(acc)+1]] <- mm[, .(n=.N, span=max(Pos)-min(Pos),
      ldw=median(ld_w_095,na.rm=TRUE), has_qtn=any(true_pos_QTN %in% TRUE)), by=CL][, arm := tg]
}
u <- rbindlist(acc)
u[, dens := fifelse(n==1, NA_real_, n/pmax(span/1e5,1e-9))]
u[, cls := fifelse(n==1,"singleton", fifelse(span<=1e5,
           fifelse(dens>=5,"short & dense","short & sparse"),
           fifelse(dens>=5,"long & dense","long & sparse")))]
cat(sprintf("  %d units, %d QTN-bearing (%d nobgs, %d bgs)\n\n", nrow(u), sum(u$has_qtn),
            sum(u$has_qtn & u$arm=="nobgs"), sum(u$has_qtn & u$arm=="bgs")))
res <- u[, {NQ <- sum(has_qtn); NM <- sum(n); NU <- .N
  .SD[, .(units=.N, markers=sum(n), qtn=sum(has_qtn),
          pct_units=100*.N/NU, pct_markers=100*sum(n)/NM, pct_qtn=100*sum(has_qtn)/NQ), by=cls]}, by=arm]
res[, `:=`(fold_unit = pct_qtn/pct_units, fold_marker = pct_qtn/pct_markers)]
cat("=== QTN enrichment by class, SEPARATELY per arm ===\n")
print(res[order(arm, -fold_unit), .(arm, cls, pct_units=round(pct_units,2),
      pct_markers=round(pct_markers,2), qtn, pct_qtn=round(pct_qtn,1),
      fold_unit=round(fold_unit,1), fold_marker=round(fold_marker,2))])
cat("\n=== the headline contrast: long & dense ===\n")
ld <- res[cls=="long & dense"]
print(ld[, .(arm, units, markers, qtn, pct_qtn=round(pct_qtn,1),
             fold_unit=round(fold_unit,1), fold_marker=round(fold_marker,2))])
## is the difference between arms significant? 2x2 on QTN in/out of the class
tab <- rbind(nobgs = c(ld[arm=="nobgs"]$qtn, sum(u$has_qtn & u$arm=="nobgs") - ld[arm=="nobgs"]$qtn),
             bgs   = c(ld[arm=="bgs"]$qtn,   sum(u$has_qtn & u$arm=="bgs")   - ld[arm=="bgs"]$qtn))
colnames(tab) <- c("in_long_dense","elsewhere")
cat("\n  QTN in long-and-dense vs elsewhere, by arm:\n"); print(tab)
ft <- fisher.test(tab)
cat(sprintf("  Fisher exact: OR %.2f [%.2f, %.2f], p = %.4g\n",
            ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))
cat("\n=== is the CLASS ITSELF smaller under BGS? ===\n")
print(u[, .(units=.N, pct_long_dense=round(100*mean(cls=="long & dense"),2),
            median_n_LD=as.numeric(median(n[cls=="long & dense"])),
            median_span_LD_kb=round(median(span[cls=="long & dense"])/1e3,1)), by=arm])
## An actual LD-only outlier scan takes the TOP TAIL of an LD statistic, not a
## whole geometric class. This is the closer proxy: top-k units by ld_w, per arm.
cat("\n=== LD-ONLY PROXY: top-k units by ld_w, per arm ===\n")
for (kk in c(1000, 5000, 20000)) {
  z <- u[, {ord <- order(-ldw); keep <- rep(FALSE, .N); keep[head(ord, kk)] <- TRUE
    .(k = kk, units = kk, pct_units = round(100*kk/.N, 2),
      qtn = sum(has_qtn[keep]), of = sum(has_qtn),
      pct_qtn = round(100*sum(has_qtn[keep])/sum(has_qtn), 1),
      fold = round((sum(has_qtn[keep])/kk) / (sum(has_qtn)/.N), 1))}, by = arm]
  print(z)
}
fwrite(res, file.path(OUT, sprintf("ld_elevation_by_arm_%s.csv", CELL)))
