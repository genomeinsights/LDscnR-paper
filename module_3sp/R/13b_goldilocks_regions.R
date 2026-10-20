## =============================================================================
## module_3sp/R/13b_goldilocks_regions.R
##
## Companion to 13_goldilocks_units.R (Supplementary Table S7). S7 counts
## significant UNITS; several significant Stage-1 units can sit in one outlier
## region, so this recounts the same significant sets as REGIONS. Reads 13's saved
## per-unit p-values (out/13_goldilocks_units/goldilocks_units.rds) -- no test is
## rerun, so the significant sets are exactly S7's -- and writes only to
## out/13b_goldilocks_regions/.
##
## MERGING RULE: the main analysis's region assembly (ld_outlier_test(assembly =
## "stage2_discovered"), REGION_ASSEMBLY in 00_config.R), reproduced call-for-call:
## ld_prune_and_eMLG() over the significant units only, every unit flagged
## (ld_w_col = "n_loci", ld_w_threshold = 0, min_n_loci_flag = 1), min_r2 from
## ld_from_rho(b, c, rho = 0.50) per chromosome on the panel's own LD decay,
## distance_threshold = 1e5, score_threshold = 0.80. For Stage 2 the significant
## Stage-2 groups are fed in as the units to merge, under the same rule.
##
## Also writes out/13b_goldilocks_regions/table6_goldilocks_units.tex: 13's S7
## longtable (same rows, same formatting) with Regions (Bonf.) / Regions (BH)
## columns added. That file, not tables/table6_goldilocks_units.tex, is the one
## the supplement now uses.
## =============================================================================
suppressMessages({library(data.table); devtools::load_all("~/gitlab/LDscnR", quiet = TRUE)})
source(file.path(path.expand("~/gitlab/LDscnR-paper/module_3sp"), "R", "00_config.R"))
STAGE <- "13b_goldilocks_regions"
LDPP <- path.expand("~/gitlab/LD-pruning-paper")
BND  <- path.expand("~/gitlab/LDscnR-paper/3sp_data/3sp_LDscnR_data.rds")
ALPHA <- 0.05; INV <- list(Chr = "Chr1", from = 21.5e6, to = 21.9e6)
OUT_DIR <- file.path(PATHS$out, STAGE); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

G <- readRDS(file.path(PATHS$out, "13_goldilocks_units", "goldilocks_units.rds"))
K1 <- copy(G$K1); K2 <- copy(G$K2)

## members per unit, from the same stage-2 caches 13 read
fs <- sort(list.files(file.path(LDPP, "cache_3sp"), "^stage2_Chr[0-9]+_d1e\\+05\\.rds$", full.names = TRUE))
M1 <- list(); M2 <- list()
for (f in fs) { o <- readRDS(f); ch <- sub(".*stage2_(Chr[0-9]+)_.*", "\\1", f)
  a <- o$stage1$clusters; M1[[ch]] <- setNames(a$members, paste0(ch, ":", a$CL_id))
  b <- o$res$groups;      M2[[ch]] <- setNames(b$members, paste0(ch, ":", b$group_id)) }
M1 <- unlist(unname(M1), recursive = FALSE); M2 <- unlist(unname(M2), recursive = FALSE)
K1[, members := M1[uid]]; K2[, members := M2[uid]]
stopifnot(!any(vapply(K1$members, is.null, NA)), !any(vapply(K2$members, is.null, NA)),
          all(lengths(K1$members) == K1$n_loci), all(lengths(K2$members) == K2$n_loci))

say("[0] loading bundle (GTs, map, LD decay) ...\n")
d <- readRDS(BND); GTs <- d$GTs; map <- as.data.table(d$map); dec <- d$LD_decay
pos <- setNames(map$Pos, map$marker); rm(d); invisible(gc())

## the main-analysis assembly applied to an arbitrary set of significant units
assemble <- function(sig) {
  if (!nrow(sig)) return(data.table(Chr = character(), from = numeric(), to = numeric(),
                                    n_units = integer(), n_markers = integer(), units = character()))
  cl <- data.table(Chr = sig$Chr, CL_id = seq_len(nrow(sig)),
                   core_snp = vapply(sig$members, `[`, "", 1L), n_snps = lengths(sig$members),
                   members = sig$members)
  mk <- unlist(cl$members, use.names = FALSE)
  ms <- data.table(marker = mk, CL_id = rep.int(cl$CL_id, cl$n_snps),
                   n_loci = rep.int(cl$n_snps, cl$n_snps), Chr = rep.int(cl$Chr, cl$n_snps),
                   Pos = unname(pos[mk]))
  sub <- structure(list(map_snp = ms, clusters = cl, pruned = cl$core_snp),
                   class = "ld_complexity_reduction")
  pr <- suppressMessages(ld_prune_and_eMLG(
          GTs = GTs[, mk, drop = FALSE], stage1 = sub, ld_w_col = "n_loci", ld_w_threshold = 0,
          LD_decay = dec, min_r2_rho = CR_RHO, score_threshold = REGION_ASSEMBLY$score_threshold,
          distance_threshold = REGION_ASSEMBLY$distance_threshold,
          compute_unflagged_eMLG = FALSE, min_n_loci_eMLG = 1, min_n_loci_flag = 1, cores = 1))
  g <- as.data.table(pr$groups)
  rbindlist(lapply(seq_len(nrow(g)), function(i) { m <- g$members[[i]]
    hit <- vapply(sig$members, function(x) any(x %chin% m), NA)
    data.table(Chr = g$Chr[i], from = min(pos[m]), to = max(pos[m]), n_units = sum(hit),
               n_markers = length(m), units = paste(sig$uid[hit], collapse = ";")) }))[order(Chr, from)]
}

UN <- c("minp","simes","emlg","rand","ldw","centr")
REG <- list(); ROWS <- list()
for (v in UN) for (corr in c("bonferroni", "BH")) {
  res <- lapply(list(`Stage 1` = K1, `Stage 2` = K2), function(k) {
    p <- k[[paste0("p_", v)]]; ok <- !is.na(p)
    sig_flag <- if (corr == "bonferroni") ok & p <= ALPHA / sum(ok)
                else { q <- rep(NA_real_, length(p)); q[ok] <- p.adjust(p[ok], "BH"); ok & q <= ALPHA }
    sig <- k[sig_flag]
    list(sig = sig, reg = assemble(sig),
         inv = sig[Chr == INV$Chr & start <= INV$to & end >= INV$from, .N] > 0) })
  for (lab in names(res)) {
    REG[[paste(v, corr, lab, sep = "|")]] <- cbind(unit = v, correction = corr, partition = lab, res[[lab]]$reg)
  }
  ROWS[[paste(v, corr)]] <- data.table(unit = v, correction = corr,
    s1_units = nrow(res$`Stage 1`$sig), s1_regions = nrow(res$`Stage 1`$reg),
    s2_units = nrow(res$`Stage 2`$sig), s2_regions = nrow(res$`Stage 2`$reg),
    s1_inversion = res$`Stage 1`$inv, s2_inversion = res$`Stage 2`$inv)
  say("  %-6s %-10s S1 %2d units -> %2d regions ; S2 %2d units -> %2d regions\n", v, corr,
      nrow(res$`Stage 1`$sig), nrow(res$`Stage 1`$reg), nrow(res$`Stage 2`$sig), nrow(res$`Stage 2`$reg))
}
TAB <- rbindlist(ROWS); REGIONS <- rbindlist(REG, fill = TRUE)
print(TAB)
fwrite(TAB, file.path(OUT_DIR, "goldilocks_regions_summary.csv"))
fwrite(REGIONS, file.path(OUT_DIR, "goldilocks_regions_list.csv"))
saveRDS(list(summary = TAB, regions = REGIONS,
             rule = list(assembly = "stage2_discovered (ld_prune_and_eMLG over significant units)",
                         min_r2_rho = CR_RHO, score_threshold = REGION_ASSEMBLY$score_threshold,
                         distance_threshold = REGION_ASSEMBLY$distance_threshold,
                         inversion = INV)),
        file.path(OUT_DIR, "goldilocks_regions.rds"))

## ---- S7 with region columns (rows and formatting as in 13_goldilocks_units.R) ----
UN_LAB <- c(minp = "min-$p$ (uses phenotype)", simes = "Simes", emlg = "consensus",
           rand = "random member", ldw = "highest $\\ldw$", centr = "central member")
fmt_sci <- function(p) {
  e <- floor(log10(p)); m <- p / 10^e
  if (round(m, 2) >= 10) { m <- m / 10; e <- e + 1 }
  sprintf("$%.2f\\times10^{%d}$", m, e)
}
hdr <- paste("Partition & Unit & Tested & Min. $p$ & Sig. (Bonf.) & Sig. (BH) & Regions (Bonf.) &",
             "Regions (BH) & Median $p$ & Median rank/$n$ \\\\")
lines <- c(
  "% GENERATED by module_3sp/R/13_goldilocks_units.R and 13b_goldilocks_regions.R -- do not edit by hand.",
  "\\begin{longtable}{llrrrrrrrr}",
  paste0("\\caption{\\textbf{Unit representations compared at matched numbers of tests.} ",
         "Stage-1 units and the broader units produced by Stage 2 on all chromosomes were each ",
         "represented in six ways and tested genome-wide. In each partition, the minimum unit size ",
         "was raised until about 100 units remained, and units spanning more than 1 Mb were excluded. ",
         "Min-$p$ uses the phenotype to choose the member with the smallest p-value and is therefore ",
         "an optimistic reference, not a valid test; the random member is a chance baseline. ",
         "For the highest-$\\ldw$ and central members, Median rank/$n$ gives the median rank of the ",
         "chosen member's p-value within its unit (1 = smallest) and the median unit size. ",
         "The highest-$\\ldw$ and central-member Min. $p$ values are identical in the two partitions ",
         "because the same Chr5 region (2.48--2.50 Mb) was the top unit in both and the same member ",
         "was chosen; the random-member values agree only to three significant figures. ",
         "Regions: significant units assembled into outlier regions with the same Stage-2 rule as in ",
         "the main analysis.}",
         "\\label{tab:goldilocks-units}\\\\"),
  "\\toprule", hdr, "\\midrule", "\\endfirsthead",
  "\\multicolumn{10}{l}{\\tablename~\\thetable{} continued}\\\\",
  "\\toprule", hdr, "\\midrule", "\\endhead")
SUMMARY <- G$summary
for (lab in c("Stage 1", "Stage 2")) {
  rows <- SUMMARY[partition == lab]; rc <- if (lab == "Stage 1") "s1_regions" else "s2_regions"
  for (i in seq_len(nrow(rows))) {
    r <- rows[i]
    nreg <- function(corr) TAB[unit == r$unit & correction == corr][[rc]]
    rank_cell <- if (r$unit == "ldw") sprintf("%.0f / %.0f", r$median_rk_ldw, r$median_n)
                else if (r$unit == "centr") sprintf("%.0f / %.0f", r$median_rk_centr, r$median_n)
                else sprintf("-- / %.0f", r$median_n)
    lines <- c(lines, sprintf("%s & %s & %d & %s & %d & %d & %d & %d & %s & %s \\\\",
      if (i == 1) lab else "", UN_LAB[r$unit], r$n_with_p, fmt_sci(r$min_p), r$sig_bonf, r$sig_bh,
      nreg("bonferroni"), nreg("BH"), sprintf("%#.3g", r$median_p), rank_cell))
  }
  lines <- c(lines, "\\midrule")
}
lines[length(lines)] <- "\\bottomrule"
lines <- c(lines, "\\end{longtable}")
OUT_TEX <- file.path(OUT_DIR, "table6_goldilocks_units.tex")
writeLines(lines, OUT_TEX)
say("\nwrote %s\n", OUT_DIR)
