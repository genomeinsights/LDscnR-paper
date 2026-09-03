#!/usr/bin/env bash
## kingman2021/R/01_extract_gts.sh
## Slice the 227-genome joint VCF down to one analysis cohort and emit a compact
## genotype table for 02_build_rds.R. bcftools does the heavy filtering so that what
## reaches R is already the final SNP set.
##
## Usage:  bash 01_extract_gts.sh <cohort> [maf] [max_missing]
##   cohort      c155_global | c150_pacNW | all227
##   maf         within-cohort minor allele freq threshold   (default 0.05)
##   max_missing max fraction of missing genotypes per SNP   (default 0.20)
##
## Outputs (git-ignored, under the data root):
##   gts/<cohort>.samples.txt        sample order actually used
##   gts/<cohort>.gt.tsv.gz          CHROM POS REF ALT then one GT column per sample
##   gts/<cohort>.sites.tsv.gz       CHROM POS REF ALT AC AN F_MISSING (audit trail)
set -euo pipefail

COHORT="${1:?cohort required: c155_global | c150_pacNW | all227}"
MAF="${2:-0.05}"
MAXMISS="${3:-0.20}"

DATA=/Users/petrikem/gitlab/LD-scaling-genome-scans/empirical_data/kingman2021
VCF="$DATA/vcf/227_genomes.final.filtered.vcf.gz"
META="$DATA/meta/tableS2_samples.tsv"
OUT="$DATA/gts"; mkdir -p "$OUT"

command -v bcftools >/dev/null || { echo "bcftools not on PATH"; exit 1; }
[[ -s "$VCF" && -s "$VCF.tbi" ]] || { echo "missing VCF or .tbi"; exit 1; }

## --- sample list for the cohort -------------------------------------------------
## Table S2 cohort columns hold "1" (marine) / "0" (freshwater) / "" (not in cohort).
awk -F'\t' -v coh="$COHORT" '
  NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next }
  {
    if (coh=="all227") { print $(h["seq_id"]) }
    else { v = $(h[coh]); if (v=="0" || v=="1") print $(h["seq_id"]) }
  }' "$META" > "$OUT/$COHORT.samples.txt"

N=$(wc -l < "$OUT/$COHORT.samples.txt")
echo "cohort $COHORT: $N samples; maf>=$MAF, max missing $MAXMISS"
[[ "$N" -gt 0 ]] || { echo "no samples selected"; exit 1; }

## --- filter ---------------------------------------------------------------------
## Autosome-like scaffolds only: the 21 named chromosomes. chrM (mitochondrial),
## chrP (the Pitx1 BAC-derived scaffold) and chrUn (unplaced) are dropped -- they have
## no meaningful position axis for LD-decay or distance-restricted clustering.
CHRS=$(cut -f1 "$DATA/tracks/gasAcu1-4.chrom.sizes" | grep -vE '^chr(M|P|Un)$' | paste -sd, -)

## -S subsets samples; AC/AN are then recomputed within the cohort by `+fill-tags`,
## so the MAF filter is a within-cohort MAF, not the 227-genome one.
bcftools view -r "$CHRS" -S "$OUT/$COHORT.samples.txt" --force-samples -Ou "$VCF" \
  | bcftools view -m2 -M2 -v snps -Ou \
  | bcftools +fill-tags -Ou -- -t AF,AC,AN,F_MISSING \
  | bcftools view -i "MAF>=$MAF && F_MISSING<=$MAXMISS" -Oz -o "$OUT/$COHORT.filtered.vcf.gz"

bcftools index -t "$OUT/$COHORT.filtered.vcf.gz"

## --- emit genotype + site tables -------------------------------------------------
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
  "$OUT/$COHORT.filtered.vcf.gz" | gzip -c > "$OUT/$COHORT.gt.tsv.gz"

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\t%F_MISSING\n' \
  "$OUT/$COHORT.filtered.vcf.gz" | gzip -c > "$OUT/$COHORT.sites.tsv.gz"

echo "SNPs retained: $(bcftools index -n "$OUT/$COHORT.filtered.vcf.gz")"
echo "wrote $OUT/$COHORT.{samples.txt,gt.tsv.gz,sites.tsv.gz}"
