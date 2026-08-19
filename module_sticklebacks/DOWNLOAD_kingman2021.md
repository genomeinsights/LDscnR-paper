# Download instructions: Kingman et al. 2021 stickleback data (clear-peak contrast dataset)

**For:** a separate session tasked with acquiring and pre-processing the Kingman et al.
2021 stickleback data so we can run the LDscnR C-score outlier-region pipeline +
detectability gate on it. Paper: Kingman GAR et al. 2021, *Predicting future from past:
The genomic basis of recurrent and rapid stickleback evolution*, Sci Adv 7:eabg5285.

**Why:** it is the **clear-peak** end of the spectrum (recurrent marine→freshwater
adaptive loci standing sharply above a quiet background), the deliberate contrast to
our current **saturated** 3sp dataset (`module_sticklebacks/`, genome full of outlier
regions, no neutral floor). We want to confirm the detectability gate passes
emphatically and the negative-control-floor `tau_C` calibration works cleanly where a
real quiet background exists.

## What the LDscnR pipeline ultimately needs (the target)
A per-"genome" list saved as `.rds` mirroring the current 3sp format:
- `GTs` — individuals × SNPs genotype dosage matrix (0/1/2), genome-wide.
- `map` — data.table with `Chr`, `Pos`, `marker` (`Chr:Pos`), and ideally `rec_rate`.
- an **ecotype / environment vector** per individual (marine=0 / freshwater=1, or the
  freshwater-adaptation axis) — the phenotype for EMMAX/LFMM.
- a **recombination map** (cM or rate per position) for LD-decay + LD-aware clustering.
- (optional) the published **EcoPeak** coordinates as a near-truth peak set to score against.

## Data sources and the reality of each

### 1. Dryad — SNP genotyping array (CONVENIENT but TARGETED, not genome-wide)
DOI: https://doi.org/10.5061/dryad.pvmcvdnjm  (~471 MB total). Files:
- `SticklebackSNPGenotypingArrayE_AlaskanSamples.gasAcu1-4.{bed,bim,fam}` — PLINK, 1,643
  individuals, **only 363 SNPs** tagging 72 adaptive + 50 neutral regions (3 Alaskan pops).
- `gasAcu1-4_anc.fa` — imputed ancestral genome (FASTA).
- `Genotyping_of_marine_sticklebacks.README`.
> ⚠️ 363 SNPs is a validation panel, NOT a discovery scan. Download it (small, fast) and
> keep it — useful for checking the known peak/neutral regions and the ecotype labels in
> the `.fam` — but it CANNOT drive a genome-wide C-score outlier scan.

### 2. Dryad — reference genome
DOI: https://doi.org/10.5061/dryad.547d7wm6t — `gasAcu1-4` reference FASTA + liftOver
chains (to/from gasAcu1). Needed for alignment coordinates and to lift the EcoPeaks/recomb.

### 3. UCSC track hub — EcoPeak / TempoPeak / recombination (the peaks + recomb map)
Hub: `http://sbwdev.stanford.edu/kingsleyAssemblyHub/hub.txt`  (**HTTP only** — an
HTTPS fetch is refused; use `curl http://...` directly).
- Parse `hub.txt` → `genomes.txt` → the per-genome `trackDb.txt` to get the `bigDataUrl`s.
- Pull the **EcoPeak** and **TempoPeak** bigBed (the recurrent adaptive-loci calls = our
  near-truth peak set) and the **recombination-rate** bigWig.
- Convert with UCSC tools: `bigBedToBed`, `bigWigToBedGraph` (from
  http://hgdownload.soe.ucsc.edu/admin/exe/). These are in gasAcu1-4 coordinates.

### 4. NCBI SRA — WGS (the actual GENOME-WIDE discovery data; needs a calling pipeline)
- **PRJNA247503** — extant-population **individual** WGS (marine + freshwater). This is
  the one to scan genome-wide. Raw reads → requires alignment + variant calling.
- PRJNA671824 — contemporary **pool-seq** time series (pooled allele frequencies, no
  individual genotypes) — alternative / for the temporal axis.
- PRJNA671690 — RS 2009 BGI genome (assembly input, not needed here).

## Recommended acquisition path

**Step 0 — try to skip the heavy calling.** Before running a full WGS pipeline, check
whether a processed **genome-wide VCF / genotype matrix** already exists:
- the paper's Supplementary Materials at https://www.science.org/doi/10.1126/sciadv.abg5285
  (and PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC8213234/ — read its Data Availability),
- the Kingsley lab site / the UCSC hub directory (there may be a genotype/allele-freq track),
- email/repo of the corresponding author if nothing is posted.
If a processed VCF exists, go straight to Step 3.

**Step 1 — quick wins (do these regardless, small):**
```bash
# targeted array + ancestral ref + README
#   (download the 3 PLINK files + README from the Dryad pvmcvdnjm landing page)
# reference + liftover
#   (download gasAcu1-4 FASTA + chains from Dryad 547d7wm6t)
# peaks + recombination from the HTTP hub
curl -sO http://sbwdev.stanford.edu/kingsleyAssemblyHub/hub.txt
#   then follow genomes.txt -> trackDb.txt, curl the EcoPeak/TempoPeak bigBed + recomb bigWig,
#   and convert: bigBedToBed ecoPeak.bb ecoPeak.bed ; bigWigToBedGraph recomb.bw recomb.bedGraph
```

**Step 2 — genome-wide genotypes from WGS (if no processed VCF found):**
```bash
# tools: sra-toolkit, bwa-mem2 (or minimap2), samtools, bcftools (or GATK), plink2
prefetch PRJNA247503        # or per-run accessions; fasterq-dump each SRR
bwa-mem2 index gasAcu1-4.fa
# per sample: bwa-mem2 mem gasAcu1-4.fa R1.fq R2.fq | samtools sort -o s.bam ; samtools index s.bam
bcftools mpileup -f gasAcu1-4.fa *.bam | bcftools call -mv -Oz -o kingman.vcf.gz
bcftools view -m2 -M2 -v snps -q 0.05:minor kingman.vcf.gz -Oz -o kingman.biallelic.vcf.gz
plink2 --vcf kingman.biallelic.vcf.gz --recode A --out kingman   # -> .raw dosage matrix
```
Assign each sample its **ecotype** (marine vs freshwater) from the sample metadata in the
paper's Table S1 / the `.fam` population codes (Rabbit Slough / Resurrection Bay etc. are
marine ancestral; the freshwater lakes are the derived ecotype).

**Step 3 — package for the LDscnR pipeline:**
- Build `GTs` (indiv × SNP, 0/1/2) from the `.raw`; `map` (Chr/Pos/marker) from the `.bim`.
- Attach `rec_rate` per SNP by intersecting with the recombination bedGraph (Step 1).
- Build the ecotype vector.
- Save `list(GTs=, map=, env=, rec=...)` as `module_sticklebacks/kingman2021_3sp.rds`
  (or per-chromosome files), matching how `regen_stats.R` / the 3sp scripts consume data.
- Keep the EcoPeak bed as the truth-peak set for scoring.

## Handoff notes
- Prioritise finding a **processed genome-wide VCF** — it removes the entire Step-2
  bioinformatics burden. Only run the WGS calling if none is posted.
- Everything is in **gasAcu1-4** coordinates; keep that consistent across genotypes,
  recombination, and EcoPeaks (liftOver only if a source is in gasAcu1).
- The 363-SNP array alone is insufficient for the outlier scan — it is a sanity/validation
  panel, not the discovery data. Do not stop at Step 1.
- Report back: (a) whether a processed VCF was found, (b) #individuals and #SNPs
  genome-wide, (c) the ecotype breakdown, (d) paths to the packaged `.rds`, recomb map,
  and EcoPeak bed.
