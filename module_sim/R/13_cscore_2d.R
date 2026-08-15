## module_sim/R/13_cscore_2d.R
## !! REPLICATE-AVERAGE: a SINGLE env is unreliable (env1 has repeatedly misled).
## Run across env1-5 and report mean +/- SE before making any claim. !!
## The poster's full consistency score: sweep BOTH ld_w window (rho) and filtering
## threshold (q*), and score each SNP by the proportion of (rho, q*) cells in which
## it is called an outlier (candidate & BH-FDR-significant within candidates). This
## generalises 10_cscore.R (q* only) to the whole (rho, q*) grid -- no single fragile
## (rho, q*) is trusted. Then gate on C to define outlier regions.
## Questions: (1) does 2D-C concentrate on QTN-linked SNPs even better than 1D?
##            (2) does C-gating beat single-SNP and the fixed (rho=0.95, RMSC q*)?
## Run from LDscnR-paper/:  Rscript module_sim/R/13_cscore_2d.R [V c env] [method]
##   method = emx_p (default) or lfmm_p
## Output (git-ignored): data/cscore2d_V*_<m>.rds, figures/cscore2d_*.png

source("module_sim/R/_config.R")
suppressMessages(library(ggplot2))
a    <- commandArgs(trailingOnly = TRUE)
V    <- if (length(a) >= 1) a[1] else "2"
cc   <- if (length(a) >= 2) a[2] else "1"
env  <- if (length(a) >= 3) a[3] else "1"
pcol <- if (length(a) >= 4) a[4] else "emx_p"
mtag <- sub("_p$", "", pcol)
QSTAR <- seq(0, 0.95, by = 0.05)
C_GRID <- seq(0.05, 1.0, by = 0.05)
R2LINK <- 0.4; FDR <- 0.05

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
RHO  <- colnames(LDW); pv <- map[[pcol]]; totQTN <- map[true_pos_QTN == TRUE, .N]
cat(sprintf("V%s_c%s_env%s [%s]: %d SNPs, true_pos_QTN=%d | grid rho=%d x q*=%d = %d cells\n",
    V, cc, env, mtag, nrow(map), totQTN, length(RHO), length(QSTAR), length(RHO) * length(QSTAR)))

## 2D C-score: fraction of (rho, q*) cells where a SNP is candidate & FDR-sig
ncell <- length(RHO) * length(QSTAR)
calls <- unlist(lapply(RHO, function(rc) { lw <- LDW[, rc]
  unlist(lapply(QSTAR, function(q) { thr <- stats::quantile(lw, q, na.rm = TRUE); cand <- which(lw >= thr)
    qv <- stats::p.adjust(pv[cand], "BH"); map$marker[cand[qv < FDR]] })) }))
tab <- table(calls)
Cs <- data.table(marker = names(tab), C = as.numeric(tab) / ncell)
map[, C := 0][Cs, on = "marker", C := i.C]

## (1) does 2D-C concentrate on QTN-linked SNPs?
map[, Cbin := cut(C, breaks = c(-.01, 0, .1, .3, .6, 1.01),
                  labels = c("0", "(0,.1]", "(.1,.3]", "(.3,.6]", "(.6,1]"))]
cat("\n=== 2D C-score vs QTN linkage ===\n")
print(map[, .(n = .N, med_maxLD_QTN = round(median(max_LD_with_QTN), 3),
              frac_r2gt.5 = round(mean(max_LD_with_QTN > 0.5), 3),
              on_neutral = sum(grepl("_Chr2$", Chr))), by = Cbin][order(Cbin)])

## (2) score: single-SNP, fixed (rho=0.95, RMSC q*), C-gated at each tau_C
mk_single <- map$marker[p.adjust(pv, "BH") < FDR]
rmsc95 <- rmsc_threshold(pv, LDW[, "rho_0.95"], grid = seq(0, 0.98, 0.01))
mk_fix <- { lw <- LDW[, "rho_0.95"]; thr <- quantile(lw, rmsc95$q_star); ci <- which(lw >= thr)
            qv <- p.adjust(pv[ci], "BH"); map$marker[ci[qv < FDR]] }
gated  <- lapply(C_GRID, function(ct) map[C >= ct, marker])
all_mk <- unique(c(mk_single, mk_fix, unlist(gated)))
qtab   <- precompute_QTN_LD(GTs, map, all_mk, 2e6, cores = 4)
score  <- function(mk, tag, tau = NA) { if (!length(mk)) return(NULL)
  reg <- cluster_regions(mk, map, GTs, r2_link = R2LINK, dcap = dcap)
  ev  <- evaluate_ORs_qtn(reg, map, qtab, th$r2min, th$dmax)
  data.table(set = tag, tau_C = tau, n_mk = length(mk), n_reg = length(reg), TP = ev$TP, FP = ev$FP,
             Recall = round(ev$Recall, 3), Precision = round(ev$Precision, 3), PR = round(ifelse(is.na(ev$PR),0,ev$PR),3)) }
res <- rbindlist(c(list(score(mk_single, "single-SNP"),
                        score(mk_fix, sprintf("fixed rho0.95 q*=%.2f", rmsc95$q_star))),
                   lapply(seq_along(C_GRID), function(i) score(gated[[i]], sprintf("C>=%.2f", C_GRID[i]), C_GRID[i]))),
                 fill = TRUE)
cat("\n=== scoring: single-SNP vs fixed (rho0.95,q*) vs 2D-C-gated (dedup-neutral) ===\n"); print(res)
best <- res[!is.na(tau_C)][which.max(PR)]
cat(sprintf("\nbest C-gate: %s -> TP=%d FP=%d Prec=%.3f Rec=%.3f PR=%.3f (vs single-SNP PR=%.3f, fixed PR=%.3f)\n",
    best$set, best$TP, best$FP, best$Precision, best$Recall, best$PR,
    res[set=="single-SNP", PR], res[grepl("fixed", set), PR]))
saveRDS(list(res = res, C = map[, .(marker, Chr, Pos, C, max_LD_with_QTN, true_pos_QTN)],
             params = list(V = V, c = cc, env = env, method = pcol, ncell = ncell)),
        file.path(dir_data, sprintf("cscore2d_V%s_c%s_env%s_%s.rds", V, cc, env, mtag)))

## PR path over tau_C + references
cg <- res[!is.na(tau_C)]; refs <- res[is.na(tau_C)]
p1 <- ggplot(cg, aes(Recall, Precision)) + geom_path(color = "grey60") +
  geom_point(aes(color = tau_C), size = 3) +
  geom_point(data = refs, aes(shape = set), size = 3.5, color = "black") +
  scale_color_viridis_c(name = "tau_C") + scale_shape_manual(values = c(15, 17), name = NULL) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = sprintf("V%s_c%s_env%s [%s]: 2D (rho,q*) C-gated vs single-SNP vs fixed rho0.95", V, cc, env, mtag),
       subtitle = "path = C gate tau_C; black = references") + theme_bw(base_size = 11)
ggsave(file.path(dir_fig, sprintf("cscore2d_pr_V%s_c%s_env%s_%s.png", V, cc, env, mtag)), p1, width = 7, height = 6, dpi = 150)

## C-score Manhattan: y = 2D consistency, colour = max r^2 with a QTN, triangle = true QTN
mm <- copy(map); mm[, chr_num := as.integer(factor(Chr, levels = unique(Chr)))]
clen <- mm[, .(len = max(Pos)), by = .(Chr, chr_num)][order(chr_num)]
spg <- 0.005 * sum(clen$len); clen[, start := c(0, head(cumsum(len + spg), -1))]
mm <- clen[, .(Chr, start)][mm, on = "Chr"][, x := Pos + start]
cmid <- clen[, .(Chr, mid = start + len / 2)]; shade <- clen[chr_num %% 2 == 0, .(xmin = start, xmax = start + len)]
qtnm <- mm[true_pos_QTN == TRUE]; bg <- mm[C == 0][sample(.N, min(.N, 15000))]
pp <- mm[C > 0]; setorder(pp, max_LD_with_QTN)
p2 <- ggplot() +
  geom_rect(data = shade, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = "grey93") +
  geom_point(data = bg, aes(x, C), color = "grey80", size = 0.3, alpha = 0.3) +
  geom_hline(yintercept = c(0.5, 0.9), linetype = 2, color = c("orange", "red"), linewidth = 0.3) +
  geom_point(data = pp, aes(x, C, color = max_LD_with_QTN), size = 1) +
  geom_point(data = qtnm, aes(x, y = -0.03), shape = 25, fill = "black", color = "black", size = 1.6) +
  scale_color_gradientn(colors = c("#3A3A98","#5DC8CD","#F2E205","#F29724","#D62828"),
                        limits = c(0, 1), name = expression(max~r^2~with~QTN)) +
  scale_x_continuous(breaks = cmid$mid, labels = gsub("R([0-9]+)_Chr", "\\1.", cmid$Chr), expand = c(0.01, 0.01)) +
  labs(x = "pooled chromosome", y = "C-score (2D rho x q* consistency)",
       title = sprintf("V%s_c%s_env%s [%s]: 2D (rho,q*) C-score Manhattan", V, cc, env, mtag),
       subtitle = "colour = max r^2 with a QTN; triangle = true QTN; dashed = tau_C 0.5 (orange) / 0.9 (red)") +
  theme_bw(base_size = 9) + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 6), legend.position = "top")
ggsave(file.path(dir_fig, sprintf("cscore2d_manhattan_V%s_c%s_env%s_%s.png", V, cc, env, mtag)), p2, width = 16, height = 5, dpi = 150)
cat("wrote cscore2d PR + Manhattan figures + rds\n")
