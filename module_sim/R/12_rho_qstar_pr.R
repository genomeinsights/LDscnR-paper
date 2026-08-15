## module_sim/R/12_rho_qstar_pr.R
## The sim (truth-aware) version of the poster heatmap: instead of TOTAL outlier
## regions (which conflate TP+FP and can be maximised by structure FP), colour the
## (rho, q*) grid by PR = Precision x Recall (dedup-neutral evaluate_ORs_qtn).
## PR peaks only where a cell finds many true QTN AND keeps FP low -> the real
## LD-filtering sweet spot. Precision and Recall saved too, as diagnostic panels.
## Run from LDscnR-paper/:  Rscript module_sim/R/12_rho_qstar_pr.R [V c env]
## Output (git-ignored): data/rhoqstar_pr_V*.rds, figures/rhoqstar_pr_V*.png

source("module_sim/R/_config.R")
suppressMessages({ library(ggplot2); library(patchwork) })
a   <- commandArgs(trailingOnly = TRUE)
V   <- if (length(a) >= 1) a[1] else "2"
cc  <- if (length(a) >= 2) a[2] else "1"
env <- if (length(a) >= 3) a[3] else "1"
QSTAR  <- seq(0, 0.95, by = 0.05)
R2LINK <- 0.4
FDR    <- 0.05

## pool GTs + map + full ld_ws (all rho columns)
files <- list.files(SIM_DATA, full.names = TRUE,
                    pattern = sprintf("^adapt_bgs_chr[0-9]+_V%s_c%s_env%s\\.rds$", V, cc, env))
reps <- as.integer(sub(".*_chr([0-9]+)_.*", "\\1", basename(files))); files <- files[order(reps)]
maps <- gts <- ldws <- decs <- vector("list", length(files))
for (i in seq_along(files)) {
  d <- readRDS(files[i]); m <- as.data.table(d$map)
  m[, `:=`(Chr = paste0("R", i, "_", Chr), marker = paste0("R", i, "_", marker))]
  G <- d$GTs; colnames(G) <- paste0("R", i, "_", colnames(G))
  lw <- d$ld_ws; rownames(lw) <- m$marker
  ds <- as.data.table(d$LD_decay$decay_sum); ds[, Chr := paste0("R", i, "_", Chr)]
  maps[[i]] <- m; gts[[i]] <- G; ldws[[i]] <- lw; decs[[i]] <- ds
}
map <- flag_true_positive_QTNs(rbindlist(maps, fill = TRUE))
GTs <- do.call(cbind, gts); LDW <- do.call(rbind, ldws)[map$marker, ]
th  <- score_thresholds(rbindlist(decs, fill = TRUE))
dcap <- d_from_rho(median(rbindlist(decs, fill = TRUE)$a_pred), 0.95)
RHO  <- colnames(LDW)
cat(sprintf("V%s_c%s_env%s: %d SNPs, %d chr, true_pos_QTN=%d | grid rho=%d x q*=%d\n",
    V, cc, env, nrow(map), uniqueN(map$Chr), map[true_pos_QTN==TRUE,.N], length(RHO), length(QSTAR)))

## pass 1: per cell -> significant markers + clusters (needs GTs, not qtab)
grid <- CJ(method = c("EMMAX", "LFMM"), rho = RHO, q = QSTAR, sorted = FALSE)
regl <- vector("list", nrow(grid)); allsig <- character(0)
for (r in seq_len(nrow(grid))) {
  pv <- map[[if (grid$method[r] == "EMMAX") "emx_p" else "lfmm_p"]]
  lw <- LDW[, grid$rho[r]]; thr <- stats::quantile(lw, grid$q[r], na.rm = TRUE)
  cand <- which(lw >= thr); qv <- stats::p.adjust(pv[cand], "BH"); sig <- map$marker[cand[qv < FDR]]
  regl[[r]] <- if (length(sig)) cluster_regions(sig, map, GTs, r2_link = R2LINK, dcap = dcap) else list()
  allsig <- c(allsig, sig)
}
qtab <- precompute_QTN_LD(GTs, map, unique(allsig), 2e6, cores = 4)

## pass 2: score each cell (dedup-neutral) -> TP/FP/precision/recall/PR
sc <- rbindlist(lapply(seq_len(nrow(grid)), function(r) {
  ev <- evaluate_ORs_qtn(regl[[r]], map, qtab, th$r2min, th$dmax)
  data.table(method = grid$method[r], rho = as.numeric(sub("rho_", "", grid$rho[r])), qstar = grid$q[r],
             TP = ev$TP, FP = ev$FP, Recall = ev$Recall, Precision = ev$Precision,
             PR = ifelse(is.na(ev$PR), 0, ev$PR))
}))
saveRDS(sc, file.path(dir_data, sprintf("rhoqstar_pr_V%s_c%s_env%s.rds", V, cc, env)))

tile <- function(fill, lab, opt = "turbo")
  ggplot(sc, aes(factor(rho), factor(qstar), fill = .data[[fill]])) +
    geom_tile() + facet_wrap(~method) + scale_fill_viridis_c(option = opt, name = lab) +
    labs(x = expression(LD~window~size~relative~to~LD~decay~(rho)), y = "local LD filtering threshold (q*)") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6), panel.grid = element_blank())
main <- tile("PR", "Precision\nx Recall") +
  labs(title = sprintf("V%s_c%s_env%s: PR across rho x q* (truth-aware sweet spot)", V, cc, env))
ggsave(file.path(dir_fig, sprintf("rhoqstar_pr_V%s_c%s_env%s.png", V, cc, env)), main, width = 12, height = 5, dpi = 150)
supp <- (tile("Precision", "Precision") + labs(title = "Precision")) /
        (tile("Recall", "Recall") + labs(title = "Recall"))
ggsave(file.path(dir_fig, sprintf("rhoqstar_pr_supp_V%s_c%s_env%s.png", V, cc, env)), supp, width = 12, height = 9, dpi = 150)

best <- sc[, .SD[which.max(PR)], by = method]
cat("\n=== best (rho, q*) by PR per method ===\n")
print(best[, .(method, rho, qstar, TP, FP, Precision = round(Precision,3), Recall = round(Recall,3), PR = round(PR,3))])
cat("wrote rhoqstar_pr figures\n")
