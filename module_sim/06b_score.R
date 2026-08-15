## module_sim/06b_score.R
## CHEAP step (no genotypes). Reads the 06a cache and scores every method over
## the shared region frame under three views:
##   res        all outlier regions
##   res_snpN   single-SNP outliers removed: keep regions a method tags with
##              >= MIN_SNP SNPs  (MIN_SNP default 2; the recommended criterion)
##   res_multi  method-singletons removed: keep regions supported by >= 2 methods
## Also reports the two support diagnostics (are FPs shared? how many single-SNP?).
## Run from LDscnR-paper/:  Rscript module_sim/06b_score.R [V c env] [MIN_SNP]
## Output (git-ignored): module_sim/consensus_V{V}_c{c}_env{env}.rds

source("module_sim/_config.R")
a   <- commandArgs(trailingOnly = TRUE)
V   <- if (length(a) >= 1) a[1] else "2"
cc  <- if (length(a) >= 2) a[2] else "1"
env <- if (length(a) >= 3) a[3] else "1"
MIN_SNP <- if (length(a) >= 4) as.integer(a[4]) else 2L      # first-class >=k-SNP option

C <- readRDS(file.path(mod, sprintf("cache_V%s_c%s_env%s.rds", V, cc, env)))
map <- C$map; regions <- C$regions; sets <- C$sets; qtab <- C$qtab; th <- C$th
lab <- C$lab; totQTN <- C$total_true_QTN
cat(sprintf("V%s_c%s_env%s (cached): %d regions, %d true_pos_QTN | TP-match r2>%.2f dist<%.0fkb | MIN_SNP=%d\n",
    V, cc, env, length(regions), totQTN, th$r2min, th$dmax / 1e3, MIN_SNP))

## per-region: which methods hit it, and how many SNPs each contributes
nsnp <- sapply(names(sets), function(nm) vapply(regions, function(mk) sum(mk %in% sets[[nm]]), integer(1)))
if (is.null(dim(nsnp))) nsnp <- matrix(nsnp, ncol = length(sets), dimnames = list(NULL, names(sets)))
pos_of <- setNames(map$Pos, map$marker); chr_of <- setNames(map$Chr, map$marker)
region_dt <- data.table(region = seq_along(regions),
                        Chr = vapply(regions, function(mk) chr_of[[mk[1]]], character(1)),
                        center = vapply(regions, function(mk) as.numeric(stats::median(pos_of[mk])), numeric(1)),
                        n_loci = lengths(regions), is_TP = lab$is_TP, qtn = lab$qtn)
for (nm in names(sets)) {
  region_dt[[nm]] <- nsnp[, nm] > 0
  region_dt[[paste0("nsnp_", nm)]] <- nsnp[, nm]
}
region_dt[, support := rowSums(nsnp > 0)]
region_dt[, singleton := support == 1L]
region_dt[, neutral_chr := grepl("_Chr2$", Chr)]

## Q1: are FALSE positives shared across methods, or method-unique?
cat(sprintf("\nregion support by METHODS: TP: %s | FP: %s\n",
    paste(sprintf("%dx=%d", 1:4, tabulate(region_dt[is_TP == TRUE, support], 4)), collapse = " "),
    paste(sprintf("%dx=%d", 1:4, tabulate(region_dt[is_TP == FALSE, support], 4)), collapse = " ")))

## Q2: single-SNP outlier share per method
cat("single-SNP outliers (a region tagged by only 1 SNP for that method):\n")
for (nm in names(sets)) { h <- nsnp[, nm] > 0
  cat(sprintf("  %-13s %d/%d hit regions single-SNP (%.0f%%); of those %d are TP\n",
      nm, sum(nsnp[, nm] == 1), sum(h), 100 * mean(nsnp[h, nm] == 1),
      sum(nsnp[, nm] == 1 & region_dt$is_TP))) }

## dedup-neutral scoring: extras (duplicate regions for an already-claimed
## true-positive QTN) count as neither TP nor FP -> clustering-parameter robust.
score_mask <- function(mask_fun, tag) rbindlist(lapply(names(sets), function(nm) {
  hit <- mask_fun(nm); ev <- evaluate_ORs_qtn(regions[hit], map, qtab, th$r2min, th$dmax)
  data.table(filter = tag, method = nm, regions_hit = sum(hit), TP = ev$TP, FP = ev$FP,
             FN = ev$FN, extras = ev$extras,
             Precision = round(ev$Precision, 3), Recall = round(ev$Recall, 3),
             PR = round(ev$PR, 3), F1 = round(2 * ev$TP / (2 * ev$TP + ev$FP + ev$FN), 3)) }))

res       <- score_mask(function(nm) nsnp[, nm] >= 1L,                             "all")
res_snpN  <- score_mask(function(nm) nsnp[, nm] >= MIN_SNP,               sprintf(">=%d SNPs", MIN_SNP))
res_multi <- score_mask(function(nm) region_dt[[nm]] & region_dt$support >= 2L,    ">=2 methods")
cat("\n=== all outlier regions ===\n"); print(res)
cat(sprintf("\n=== single-SNP outliers removed (>=%d SNPs per method) ===\n", MIN_SNP)); print(res_snpN)
cat("\n=== method-singletons removed (>=2 methods) ===\n"); print(res_multi)

saveRDS(list(res = res, res_snpN = res_snpN, res_multi = res_multi, MIN_SNP = MIN_SNP,
             region_dt = region_dt, nsnp = nsnp, regions = regions, lab = lab, sets = sets,
             map = map, total_true_QTN = totQTN, params = C$params),
        file.path(mod, sprintf("consensus_V%s_c%s_env%s.rds", V, cc, env)))
