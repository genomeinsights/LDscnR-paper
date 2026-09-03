## kingman2021/R/02_build_rds.R
## Package one Kingman-2021 cohort into the object the 3sp scripts consume:
##   GTs  individuals x SNPs dosage matrix, colnames = marker. MEAN-IMPUTED and numeric,
##        not integer 0/1/2 -- see the missingness note below.
##   map  data.table(marker, Chr, chr_num, Pos, rec_rate, cM, af, f_missing)
##   eco  integer ecotype vector, 1 = Marine, 0 = Freshwater (matches 16_sim_machinery's
##        `eco <- as.integer(pheno$ecotype == "Marine")`)
##   meta the Table S2 rows for the cohort, in GTs row order
##
## rec_rate is the Rabbit Slough LDhelmet map in cM/Mb; cM is its running integral, so
## distance-restricted clustering can work in genetic rather than physical distance.
##
## MISSINGNESS. Mean sequencing coverage is 5.5x, so genotype missingness is pervasive:
## at MAF>=0.05 / F_MISSING<=0.20 roughly 12% of calls are missing, and tightening the
## filter is not a way out (F_MISSING<=0.05 keeps only ~4% of SNPs). We therefore keep
## the permissive site filter and mean-impute each SNP to its observed mean dosage,
## which is what gcta_grm()/fast_emmax_setup() in 16_sim_machinery_full_3sp.R require
## (they have no NA handling). Per-SNP `f_missing` is carried in `map` so any downstream
## step can filter on it. Caveat worth remembering: hard calls at 5.5x are themselves
## noisy, so r2/ld_w will be attenuated relative to a high-coverage dataset.
##
## Run from LDscnR-paper/ AFTER 01_extract_gts.sh for the same cohort:
##   Rscript kingman2021/R/02_build_rds.R c155_global
## Output (git-ignored): kingman2021/data/kingman2021_<cohort>.rds

suppressMessages({ library(data.table) })

args   <- commandArgs(trailingOnly = TRUE)
COHORT <- if (length(args) >= 1) args[1] else "c155_global"

DATA <- "/Users/petrikem/gitlab/LD-scaling-genome-scans/empirical_data/kingman2021"
MOD  <- "/Users/petrikem/gitlab/LDscnR-paper/kingman2021/data"
gtf  <- file.path(DATA, "gts", paste0(COHORT, ".gt.tsv.gz"))
smf  <- file.path(DATA, "gts", paste0(COHORT, ".samples.txt"))
stopifnot(file.exists(gtf), file.exists(smf))

samples <- readLines(smf)
cat(sprintf("cohort %s: %d samples\n", COHORT, length(samples)))

## ---- genotypes -----------------------------------------------------------------
## bcftools emits GT strings; map to dosage. Both phased (|) and unphased (/) appear
## in this VCF because GATK wrote PGT/PID physical phasing for some sites.
gt <- fread(gtf, header = FALSE, showProgress = TRUE)
setnames(gt, 1:4, c("Chr", "Pos", "REF", "ALT"))
stopifnot(ncol(gt) == 4L + length(samples))

dose_map <- c("0/0" = 0L, "0|0" = 0L,
              "0/1" = 1L, "0|1" = 1L, "1/0" = 1L, "1|0" = 1L,
              "1/1" = 2L, "1|1" = 2L)
gt_cols <- names(gt)[-(1:4)]
for (j in gt_cols) set(gt, j = j, value = unname(dose_map[gt[[j]]]))   # unmatched -> NA

GTs <- t(as.matrix(gt[, ..gt_cols]))
storage.mode(GTs) <- "integer"
rownames(GTs) <- samples

map <- gt[, .(Chr, Pos)]
map[, marker := paste0(Chr, ":", Pos)]
stopifnot(!anyDuplicated(map$marker))
colnames(GTs) <- map$marker
rm(gt); gc()

## roman -> 1..21 so downstream ordering/plotting has a numeric axis
ROMAN <- c("I","II","III","IV","V","VI","VII","VIII","IX","X",
           "XI","XII","XIII","XIV","XV","XVI","XVII","XVIII","XIX","XX","XXI")
map[, chr_num := match(sub("^chr", "", Chr), ROMAN)]
stopifnot(!anyNA(map$chr_num))
setorder(map, chr_num, Pos)
GTs <- GTs[, map$marker, drop = FALSE]

## ---- allele frequency, missingness, mean imputation --------------------------------
n_ind  <- nrow(GTs)
n_miss <- colSums(is.na(GTs))
map[, f_missing := n_miss / n_ind]
col_mean <- colMeans(GTs, na.rm = TRUE)
map[, af := col_mean / 2]

cat(sprintf("GTs: %d individuals x %d SNPs over %d chromosomes; %.2f%% missing calls\n",
            n_ind, ncol(GTs), uniqueN(map$Chr), 100 * sum(n_miss) / prod(dim(GTs))))
cat(sprintf("  per-SNP missingness: median %.3f, 90th pct %.3f, max %.3f\n",
            median(map$f_missing), quantile(map$f_missing, 0.9), max(map$f_missing)))

## a SNP with no calls at all would give NaN; the F_MISSING filter should exclude these
stopifnot(all(is.finite(col_mean)))
GTs <- matrix(as.numeric(GTs), nrow = n_ind,
              dimnames = list(rownames(GTs), colnames(GTs)))
na_idx <- which(is.na(GTs))
GTs[na_idx] <- col_mean[((na_idx - 1L) %/% n_ind) + 1L]     # column index of each NA
stopifnot(!anyNA(GTs))
cat(sprintf("  mean-imputed %d genotype cells; GTs is now numeric with no NA\n",
            length(na_idx)))

## ---- recombination: cM/Mb per SNP + cumulative cM --------------------------------
## Rabbit Slough rho-based map (Roberts Kingman et al. 2021), 2kb-ish intervals.
rec <- fread(file.path(DATA, "tracks", "gasAcu1-4.scaledRABSrecombRates.bedGraph"),
             header = FALSE, col.names = c("Chr", "start", "end", "rec_rate"))
rec <- rec[Chr %in% map$Chr]
setorder(rec, Chr, start)

## running cM along each chromosome: rate (cM/Mb) x interval length (Mb)
rec[, cM_end := cumsum((end - start) / 1e6 * rec_rate), by = Chr]
rec[, cM_start := cM_end - (end - start) / 1e6 * rec_rate]

## bedGraph is 0-based half-open; VCF Pos is 1-based. Pos p sits in [start, end) at p-1.
map[, pos0 := Pos - 1L]
setkey(rec, Chr, start, end)
hit <- foverlaps(map[, .(Chr, start = pos0, end = pos0 + 1L, marker)],
                 rec, type = "any", mult = "first", nomatch = NA)
map[, rec_rate := hit$rec_rate]
## linear interpolation of cM within the covering interval
map[, cM := hit$cM_start + (pos0 - hit$start) / 1e6 * hit$rec_rate]
map[, pos0 := NULL]

cat(sprintf("rec map: %.1f%% of SNPs covered; rec_rate median %.2f cM/Mb\n",
            100 * mean(!is.na(map$rec_rate)), median(map$rec_rate, na.rm = TRUE)))
if (anyNA(map$rec_rate))
  cat(sprintf("  NOTE: %d SNPs fall in gaps of the RABS map (rec_rate/cM = NA)\n",
              sum(is.na(map$rec_rate))))

## ---- ecotype ---------------------------------------------------------------------
meta <- fread(file.path(DATA, "meta", "tableS2_samples.tsv"))
meta <- meta[match(samples, seq_id)]
stopifnot(identical(meta$seq_id, samples))
eco <- as.integer(meta$ecotype == "M")          # 1 = Marine, 0 = Freshwater
stopifnot(!anyNA(eco))
cat(sprintf("ecotype: %d marine, %d freshwater\n", sum(eco == 1L), sum(eco == 0L)))

## the cohort flag is itself the ecotype; verify they agree rather than trusting one
if (COHORT %in% names(meta)) {
  flag <- suppressWarnings(as.integer(meta[[COHORT]]))
  stopifnot(identical(flag, eco))
  cat(sprintf("cohort flag %s agrees with Table S2 ecotype for all %d samples\n",
              COHORT, length(eco)))
}

## ---- save --------------------------------------------------------------------------
out <- file.path(MOD, sprintf("kingman2021_%s.rds", COHORT))
saveRDS(list(GTs    = GTs,
             map    = map[, .(marker, Chr, chr_num, Pos, rec_rate, cM, af, f_missing)],
             eco    = eco,
             meta   = meta,
             cohort = COHORT,
             source = "Roberts Kingman et al. 2021 Sci Adv 7:eabg5285; 227_genomes.final.filtered.vcf.gz (FigShare project 162634); gasAcu1-4"),
        out)
cat(sprintf("\nwrote %s (%.1f MB)\n", out, file.size(out) / 1e6))
