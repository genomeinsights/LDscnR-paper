## Is a QTN CENTRED in its stage-2 cluster?
##
## Two accounts of why long-and-dense clusters hold 70% of QTN at 30-fold:
##   (a) selection extends LD around a causal site, so the QTN MANUFACTURES the
##       block -- the block should then be centred on the QTN;
##   (b) low-recombination regions form dense blocks anyway and QTN happen to be
##       in them -- such a block has no reason to be centred on anything.
##
## Suggested by the 3sp panel session as a test that does not depend on matching.
## Relative position = (QTN_pos - min) / (max - min) over the cluster's members,
## folded to |x - 0.5| so 0 = perfectly central, 0.5 = at an edge. Under (b) the
## folded statistic is Uniform(0, 0.5) with mean 0.25.
suppressMessages({library(data.table); library(LDscnR)})
SIM  <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
OUT  <- Sys.getenv("OUT", "module_sim_LDscnR/results/filter_then_test")
ENVS <- as.integer(strsplit(Sys.getenv("ENVS", "1,2,3,4,5"), ",")[[1]])
MINN <- as.integer(Sys.getenv("MINN", "5"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

acc <- list()
for (tg in c("nobgs","bgs")) for (ev in ENVS) for (i in 1:10) {
  f <- sprintf("%s/adapt_%s_chr%d_V0.5_c1_env%d.rds", SIM, tg, i, ev)
  if (!file.exists(f)) next
  x <- readRDS(f); m <- flag_true_qtns(as.data.table(x$map))
  pr <- ld_prune_and_eMLG(GTs=x$GTs, stage1=x$complexity_reduction$stage1, LD_decay=x$LD_decay,
        ld_w_col="ld_w_095", ld_w_threshold=0.025, score_threshold=0.80, min_r2_rho=0.5,
        distance_threshold=1e5, compute_unflagged_eMLG=FALSE, cores=1)
  g <- as.data.table(pr$groups)
  ms <- data.table(marker = unlist(g$members, use.names = FALSE),
            CL = paste0(tg,ev,i,"_",rep.int(g$group_id, lengths(g$members))))
  mm <- merge(m, ms, by="marker", all.x=TRUE)[!is.na(CL)]
  qcl <- mm[true_pos_QTN %in% TRUE, .(CL, qpos = Pos)]
  if (!nrow(qcl)) next
  gg <- mm[CL %in% qcl$CL, .(lo=min(Pos), hi=max(Pos), n=.N,
                             ldw=median(ld_w_095,na.rm=TRUE)), by=CL]
  z <- merge(qcl, gg, by="CL")[n >= MINN & hi > lo]
  if (nrow(z)) acc[[length(acc)+1]] <- z[, .(CL, arm=tg, env=ev, n, span=hi-lo, ldw,
                                             rel=(qpos-lo)/(hi-lo))]
}
d <- rbindlist(acc)
d[, folded := abs(rel - 0.5)]
cat(sprintf("  %d QTN-bearing clusters with >= %d markers\n\n", nrow(d), MINN))
cat(sprintf("  relative position: mean %.3f, median %.3f (0.5 = centre)\n", mean(d$rel), median(d$rel)))
cat(sprintf("  folded |rel-0.5|: mean %.3f, median %.3f\n", mean(d$folded), median(d$folded)))
cat("  under the no-centring account the folded statistic is Uniform(0,0.5), mean 0.25\n\n")
tt <- t.test(d$folded, mu = 0.25)
ks <- suppressWarnings(ks.test(d$rel, "punif"))
cat(sprintf("  t-test vs 0.25 : mean %.3f, 95%% CI [%.3f, %.3f], p = %.4g\n",
            tt$estimate, tt$conf.int[1], tt$conf.int[2], tt$p.value))
cat(sprintf("  KS vs Uniform  : D = %.3f, p = %.4g\n", ks$statistic, ks$p.value))
cat(sprintf("\n  in the central third (0.33-0.67): %d of %d = %.0f%% (expected 33%%)\n",
            sum(d$rel > 1/3 & d$rel < 2/3), nrow(d), 100*mean(d$rel > 1/3 & d$rel < 2/3)))
cat(sprintf("  in the outer 20%% (edges):        %d of %d = %.0f%% (expected 20%%)\n",
            sum(d$rel < .1 | d$rel > .9), nrow(d), 100*mean(d$rel < .1 | d$rel > .9)))
cat("\n  by arm:\n")
print(d[, .(clusters=.N, mean_rel=round(mean(rel),3), mean_folded=round(mean(folded),3),
            median_n=as.numeric(median(n)), median_span_kb=round(median(span)/1e3,1)), by=arm])
fwrite(d, file.path(OUT, "qtn_centrality.csv"))
