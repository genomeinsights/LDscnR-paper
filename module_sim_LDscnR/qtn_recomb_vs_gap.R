## =====================================================================
## module_sim_LDscnR / qtn_recomb_vs_gap.R
##
## Does the C-score's advantage over BH track the RECOMBINATION ENVIRONMENT of
## a genome's causal variants, rather than its power?
##
## Mechanism under test: ld_w can only mark a QTN that sits in a region with
## enough local LD to mark. Where recombination is high, ld_w is pinned at the
## chromosome background b and carries no information, so C discards markers BH
## would still test. Prediction: C - alpha should correlate with where the QTN
## sit (recombination / their own ld_w), and NOT with how much power there was.
##
## Truth uses flag_true_qtns(), i.e. the same detectable set the PR-AUC scoring
## uses -- not all QTN.
##
## Run from the LDscnR-paper root:
##   Rscript module_sim_LDscnR/qtn_recomb_vs_gap.R
## Env: SIM_DATA, RESULTS, CELLS (default: whatever result CSVs exist)
## =====================================================================
suppressMessages({ library(data.table); library(LDscnR) })
SIM  <- Sys.getenv("SIM_DATA", "/Volumes/Nemo/Nemo_sim/regen_sim_data_bgs5")
RES  <- Sys.getenv("RESULTS", "module_sim_LDscnR/results")

fs <- list.files(RES, pattern="^bgs_vs_nobgs_prauc[.]csv$", recursive=TRUE, full.names=TRUE)
if (!length(fs)) stop("no result CSVs under ", RES)
pr <- unique(rbindlist(lapply(fs, fread), fill=TRUE))
pr[, gap := PR_AUC_C - PR_AUC_alpha]
## one row per genome = (cell, tag); engine and l_min are averaged over
g <- pr[, .(gap = mean(gap, na.rm=TRUE), alpha = mean(PR_AUC_alpha, na.rm=TRUE),
            C = mean(PR_AUC_C, na.rm=TRUE), n_true = n_true[1]), by=.(cell, tag)]

pool_map <- function(tag, cell) {
  ff <- list.files(SIM, full.names=TRUE,
    pattern=sprintf("^adapt_%s_chr[0-9]+_%s[.]rds$", tag, gsub("\\.","[.]",cell)))
  if (!length(ff)) return(NULL)
  ff <- ff[order(as.integer(sub(".*_chr([0-9]+)_.*","\\1",basename(ff))))]
  rbindlist(lapply(seq_along(ff), function(i) {
    x <- readRDS(ff[i]); m <- as.data.table(x$map)
    b <- mean(as.data.table(x$LD_decay$decay_sum)$b, na.rm=TRUE)
    m[, `:=`(Chr=paste0("R",i,"_",Chr), marker=paste0("R",i,"_",marker), chr_b=b)]
    m
  }), fill=TRUE)
}

out <- rbindlist(lapply(seq_len(nrow(g)), function(i) {
  m <- pool_map(g$tag[i], g$cell[i]); if (is.null(m)) return(NULL)
  m <- flag_true_qtns(m)
  q <- m[true_pos_QTN %in% TRUE & is.finite(rec_rate)]
  if (!nrow(q)) return(NULL)
  data.table(cell=g$cell[i], tag=g$tag[i], gap=g$gap[i], alpha=g$alpha[i], C=g$C[i],
             n_qtn      = nrow(q),
             qtn_rec    = median(log1p(q$rec_rate)),              # where the QTN sit
             qtn_ldw    = median(q$ld_w_095, na.rm=TRUE),
             qtn_above_b= mean(q$ld_w_095 > q$chr_b, na.rm=TRUE), # reachable by ld_w at all
             gen_rec    = median(log1p(m$rec_rate), na.rm=TRUE))
}), fill=TRUE)
fwrite(out, file.path(RES, "qtn_recomb_vs_gap.csv"))

cat(sprintf("\n  %d genomes from %d cell(s)\n\n", nrow(out), uniqueN(out$cell)))
print(out[order(gap), .(cell, tag, n_qtn, qtn_rec=round(qtn_rec,2),
                        qtn_above_b=round(qtn_above_b,2), qtn_ldw=round(qtn_ldw,4),
                        alpha=round(alpha,3), gap=round(gap,3))])
ct <- function(x, y, lab) {
  if (sum(is.finite(x) & is.finite(y)) < 5) return(invisible())
  s <- suppressWarnings(cor.test(x, y, method="spearman"))
  cat(sprintf("  %-34s rho %+.3f   p %.3f\n", lab, s$estimate, s$p.value))
}
cat("\n=== does the gap track the QTN's recombination environment? ===\n")
ct(out$qtn_rec,     out$gap, "gap ~ median log1p(rec) at QTN")
ct(out$qtn_above_b, out$gap, "gap ~ frac QTN with ld_w > b")
ct(out$qtn_ldw,     out$gap, "gap ~ median ld_w at QTN")
cat("\n=== ... or does it just track power? ===\n")
ct(out$alpha,  out$gap, "gap ~ alpha PR-AUC (power proxy)")
ct(out$n_qtn,  out$gap, "gap ~ n detectable QTN")
