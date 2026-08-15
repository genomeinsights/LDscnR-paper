## module_sim/R/16_besttau_manhattan.R
## Manhattan comparison at the best operating points: C-score @ tau_C vs single-SNP
## @ best alpha, scored JOINTLY on one shared region frame so we can see which
## outlier regions (ORs) are shared between the two methods and which are unique.
## Four outlier sets = {C-score, single-SNP} x {EMMAX, LFMM}; union -> cluster once
## (fixed decay-relative rho_ld=0.75/rho_d=0.95<=500kb) -> ORs (>=l_min loci) ->
## TP/FP label -> colour by shared region identity (same colour across panels).
## Run from LDscnR-paper/:  Rscript module_sim/R/16_besttau_manhattan.R [V c env]
## Output (git-ignored): figures/besttau_manhattan_V*.png, data/besttau_V*.rds

source("module_sim/R/_config.R")
suppressMessages({ library(ggplot2); library(patchwork) })
a   <- commandArgs(trailingOnly = TRUE)
V   <- if (length(a) >= 1) a[1] else "2"
cc  <- if (length(a) >= 2) a[2] else "1"
env <- if (length(a) >= 3) a[3] else "1"
TAUC    <- 0.35                                # best C-score gate (both methods)
ALPHA   <- c(EMMAX = 0.05, LFMM = 0.02)        # best single-SNP alpha per method
QSTAR   <- seq(0, 0.95, by = 0.05); ALPHA_C <- c(0.001, 0.01, 0.05, 0.1)
LMIN    <- if (length(a) >= 4) as.integer(a[4]) else 2L   # l_min=1 -> no single-SNP filter
RHO_LD <- 0.75; RHO_D <- 0.95; DCAP <- 5e5

## pool GTs + map + full ld_ws + decay
files <- list.files(SIM_DATA, full.names = TRUE,
                    pattern = sprintf("^adapt_bgs_chr[0-9]+_V%s_c%s_env%s\\.rds$", V, cc, env))
reps <- as.integer(sub(".*_chr([0-9]+)_.*", "\\1", basename(files))); files <- files[order(reps)]
maps <- gts <- ldws <- decs <- vector("list", length(files))
for (i in seq_along(files)) {
  d <- readRDS(files[i]); m <- as.data.table(d$map)
  m[, `:=`(Chr = paste0("R", i, "_", Chr), marker = paste0("R", i, "_", marker))]
  G <- d$GTs; colnames(G) <- paste0("R", i, "_", colnames(G)); lw <- d$ld_ws; rownames(lw) <- m$marker
  ds <- as.data.table(d$LD_decay$decay_sum); ds[, Chr := paste0("R", i, "_", Chr)]
  maps[[i]] <- m; gts[[i]] <- G; ldws[[i]] <- lw; decs[[i]] <- ds
}
map <- flag_true_positive_QTNs(rbindlist(maps, fill = TRUE))
GTs <- do.call(cbind, gts); LDW <- do.call(rbind, ldws)[map$marker, ]
decay_sum <- rbindlist(decs, fill = TRUE)
th  <- score_thresholds(decay_sum, 0.75, 0.95, dmax_cap = DCAP)
RHO <- colnames(LDW); ncell <- length(RHO) * length(QSTAR) * length(ALPHA_C)

## outlier sets: C-score (gate tau_C) and single-SNP (best alpha), per method
Cgate <- function(pcol) { pv <- map[[pcol]]
  calls <- unlist(lapply(RHO, function(rc) { lw <- LDW[, rc]
    unlist(lapply(QSTAR, function(q) { thr <- stats::quantile(lw, q, na.rm = TRUE); cand <- which(lw >= thr)
      qv <- stats::p.adjust(pv[cand], "BH"); unlist(lapply(ALPHA_C, function(al) map$marker[cand[qv < al]])) })) }))
  tb <- table(calls); names(tb)[as.numeric(tb) / ncell >= TAUC] }
sets <- list("EMMAX C-score" = Cgate("emx_p"),
             "EMMAX single"  = map$marker[p.adjust(map$emx_p,  "BH") < ALPHA["EMMAX"]],
             "LFMM C-score"  = Cgate("lfmm_p"),
             "LFMM single"   = map$marker[p.adjust(map$lfmm_p, "BH") < ALPHA["LFMM"]])
cat("outlier SNPs:", paste(sprintf("%s=%d", names(sets), lengths(sets)), collapse = "  "), "\n")

## joint shared region frame: union -> cluster once (fixed rho) -> ORs >= l_min
union_mk <- unique(unlist(sets))
regions0 <- cluster_regions(union_mk, map, GTs, decay_sum = decay_sum, rho_ld = RHO_LD, rho_d = RHO_D, dcap_max = DCAP)
regions  <- regions0[lengths(regions0) >= LMIN]
qtab <- precompute_QTN_LD(GTs, map, union_mk, 2e6, cores = 4)
lab  <- classify_ORs(regions, map, qtab, th$r2min, th$dmax)
nsnp <- sapply(names(sets), function(nm) vapply(regions, function(mk) sum(mk %in% sets[[nm]]), integer(1)))
if (is.null(dim(nsnp))) nsnp <- matrix(nsnp, ncol = length(sets), dimnames = list(NULL, names(sets)))
support <- rowSums(nsnp > 0)
cat(sprintf("\nORs (>=%d loci): %d | TP: %s | FP: %s  (support = # of the 4 sets hitting)\n", LMIN, length(regions),
    paste(sprintf("%dx=%d", 1:4, tabulate(support[lab$is_TP], 4)), collapse = " "),
    paste(sprintf("%dx=%d", 1:4, tabulate(support[!lab$is_TP], 4)), collapse = " ")))
## how often do C-score and single-SNP hit the SAME OR (per method)?
for (m in c("EMMAX", "LFMM")) { hC <- nsnp[, paste(m,"C-score")] > 0; hS <- nsnp[, paste(m,"single")] > 0
  cat(sprintf("  %-5s: C-score ORs=%d, single ORs=%d, shared=%d (TP shared=%d)\n",
      m, sum(hC), sum(hS), sum(hC & hS), sum(hC & hS & lab$is_TP))) }

## Manhattan: 4 panels, y=-log10(q), colour = shared (>=2-set) region identity
map[, `:=`(x = NA_real_, nlq_emx = -log10(pmax(p.adjust(emx_p,"BH"),1e-30)), nlq_lfmm = -log10(pmax(p.adjust(lfmm_p,"BH"),1e-30)))]
map[, chr_num := as.integer(factor(Chr, levels = unique(Chr)))]
clen <- map[, .(len = max(Pos)), by = .(Chr, chr_num)][order(chr_num)]
spg <- 0.005 * sum(clen$len); clen[, start := c(0, head(cumsum(len + spg), -1))]
map <- clen[, .(Chr, start)][map, on = "Chr"][, x := Pos + start]
cmid <- clen[, .(Chr, mid = start + len/2)]; shade <- clen[chr_num %% 2 == 0, .(xmin = start, xmax = start + len)]
xkey <- setNames(map$x, map$marker); m2r <- integer(0); for (i in seq_along(regions)) m2r[regions[[i]]] <- i
reg_x <- vapply(regions, function(mk) stats::median(xkey[mk], na.rm=TRUE), numeric(1))
multi <- which(support >= 2L); multi <- multi[order(reg_x[multi])]
PAL <- rep(c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFD400","#A65628","#F781BF","#1B9E77","#D95F02","#66A61E","#E7298A"), 4)
reg_col <- rep("grey55", length(regions)); reg_col[multi] <- PAL[seq_along(multi)]
qtn <- map[true_pos_QTN == TRUE]
panel <- function(nm, ycol) { out <- map[marker %in% sets[[nm]]][, region := m2r[marker]]
  out <- out[!is.na(region)][, col := reg_col[region]]      # only SNPs that fall in an OR (>=l_min)
  ggplot() +
    geom_rect(data = shade, aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf), fill="grey93") +
    geom_point(data = map, aes(x, .data[[ycol]]), color="grey80", size=0.2, alpha=0.3) +
    geom_point(data = qtn, aes(x, y=-0.4), shape=25, fill="black", color="black", size=1.2) +
    geom_point(data = out, aes(x, .data[[ycol]], color=col), size=1.1) + scale_color_identity() +
    scale_x_continuous(breaks=cmid$mid, labels=gsub("R([0-9]+)_Chr","\\1.",cmid$Chr), expand=c(0.01,0.01)) +
    labs(x=NULL, y=expression(-log[10](q)), title=sprintf("%s  (ORs hit: %d)", nm, uniqueN(out$region))) +
    theme_bw(base_size=9) + theme(panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), axis.text.x=element_text(size=6)) }
P <- panel("EMMAX C-score","nlq_emx") / panel("EMMAX single","nlq_emx") /
     panel("LFMM C-score","nlq_lfmm") / panel("LFMM single","nlq_lfmm") +
  plot_annotation(title=sprintf("V%s_c%s_env%s: C-score(tau=%.2f) vs single-SNP(best alpha), joint OR frame", V, cc, env, TAUC),
                  subtitle="colour = shared (>=2-set) OR identity (same across panels); grey = set-unique OR; triangle = true QTN")
ggsave(file.path(dir_fig, sprintf("besttau_manhattan_lmin%d_V%s_c%s_env%s.png", LMIN, V, cc, env)), P, width=16, height=10, dpi=150)
saveRDS(list(sets=sets, regions=regions, lab=lab, nsnp=nsnp, support=support), file.path(dir_data, sprintf("besttau_lmin%d_V%s_c%s_env%s.rds", LMIN, V, cc, env)))
cat("wrote besttau_manhattan figure\n")
