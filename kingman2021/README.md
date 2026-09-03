# kingman2021

The **clear-peak** stickleback dataset — Roberts Kingman *et al.* 2021, *Predicting future
from past: the genomic basis of recurrent and rapid stickleback evolution*,
Sci Adv 7:eabg5285 ([10.1126/sciadv.abg5285](https://doi.org/10.1126/sciadv.abg5285)) —
acquired, packaged for the `LDscnR` pipeline, and compared against the 3sp outlier regions.

## The point

`module_sticklebacks_LDscnR/` is the **saturated** case: the whole 3sp genome carries real
outlier regions and there is no neutral floor. This module is the deliberate contrast —
recurrent marine→freshwater loci standing sharply above a **quiet background**, with an
independently published peak set to score against. It is what the detectability gate and
the negative-control-floor `tau_C` calibration need in order to be exercised somewhere a
floor actually exists.

Sanity check on that claim: chrI:1–3 Mb contains no EcoPeak, and **0 of 24,551 SNPs** there
reach p < 1e-6 in an uncorrected ecotype regression (λ ≈ 1.17). Inside the 39 Global-specific
EcoPeaks, extreme |ΔAF| is enriched **298×** over background.

## Headline findings

1. **No variant calling was needed.** A processed, GATK+VQSR joint VCF of all 227 genomes
   is mirrored on **FigShare project 162634** — reachable only via the hub description page,
   not from the paper's Data Availability statement. The 0.7 TB SRA path was superseded.
2. **Cohorts are one genome per population**, and Table S2's cohort columns hold the
   ecotype (`1`=marine, `0`=freshwater, blank = not in cohort). Use **`c155_global`**
   (84: 28 M / 56 F) as primary, paired with the 39 Global-specific EcoPeaks.
3. **Overlap with the 3sp outlier regions is weak** (1.1–2.5×, mostly ns) — **because of
   geography, not method failure.** Re-deriving the peak caller on the geographically
   matched Northern-European cohort roughly triples the agreement.
4. **c151 (N. Europe) cannot yield an FDR-significant peak set.** With 9 marine / 18
   freshwater the exact test bottoms out at p = 3.2e-7 while BH needs ≥54 SNPs at that
   floor. This is structural, and very likely why the hub ships peak BEDs for c150/c155 only.
5. **`l_min` is consequential.** The EMMAX `l_min=3` set — 17 regions over just 1.21 Mb —
   hits the Global-specific EcoPeaks **19.5×** over null (p=5e-4), far better than either
   LFMM set. And ***Eda* is a 4-SNP cluster**: recovered by EMMAX at `l_min=3`, discarded by
   `l_min=10`, and never found by LFMM at either floor. Any analysis demanding ≥10 linked
   SNPs throws away the textbook marine–freshwater locus.
6. **Kingman's own two datasets set a high ceiling.** Global-specific and Pacific-specific
   EcoPeaks are essentially *nested* — 95% of the Global-specific sequence lies inside the
   Pacific-specific set. Our best cross-dataset agreement is 0.44, less than half of that.
   So the 3sp regions are **not** reproducing Kingman at the achievable limit; the gap is
   real, and §"geographic gradient" argues it is largely geographic.
8. **The uncorroborated LFMM regions carry no signal in a matched cohort.** LFMM's 107
   regions overlapping no EcoPeak sit at **0.98×, p=0.27** in the geography-matched c151
   Northern-European cohort — over 21 Mb, invariant to how the split is cut — while both
   engines' corroborated regions light up (23.6× / 4.33×), so the cohort has the power.
   That closes the "novel Atlantic loci" reading of the LFMM surplus. It does **not** show
   they are false positives — real non-marine-freshwater structure looks the same.
7. **The EcoPeaks also calibrate, not just score.** Sweeping the τ_C × l_min grid against
   them: precision/fold are *not estimable* in most cells (a lone region that hits gives
   precision 1.0), and masking below 5 regions leaves **33 of 270 cells for EMMAX but 247
   for LFMM** — C discriminates for one engine and barely for the other. Caveat: LFMM has
   **no structure-aware null** here, so that comparison is confounded; see the report.

## Layout

```
kingman2021/
├── R/          scripts, numbered in run order
├── data/       small derived results (large .rds and .snp_p.tsv.gz are git-ignored)
├── figures/    report figures
└── doc/        LaTeX report + PDF, and the full provenance record
```

Bulk data (the 15 GB VCF, per-cohort genotype extracts, tracks, reference) lives outside
the repo at `~/gitlab/LD-scaling-genome-scans/empirical_data/kingman2021/`, which is
git-ignored there.

## Scripts

| Script | Does |
|---|---|
| `R/01_extract_gts.sh` | bcftools: subset the joint VCF to one cohort, biallelic SNPs, within-cohort MAF and missingness filters |
| `R/02_build_rds.R` | build `list(GTs, map, eco, meta, ...)`; **mean-imputes** `GTs` (see below) |
| `R/03_peaks_truth.R` | published EcoPeak/TempoPeak BEDs → truth table + per-SNP membership, row-aligned with `map` |
| `R/04_ecopeaks_snp_test.py` | re-implementation of Kingman's SNP-based multivariate-hypergeometric ecotype test |
| `R/05_call_peaks.R` | BH-FDR + 50 kb merge → peaks; validates against the published sets |
| `R/06_c151_suggestive_peaks.R` | rank-based (**not** FDR-controlled) c151 peak set |
| `R/07_signal_enrichment.R` | rank-based enrichment of low Kingman p inside the 3sp regions |
| `R/08_liftover.sh` | reconstruct the gasAcu1 ↔ gasAcu1-4 chains and lift everything |
| `R/09_overlap.R` | region-level overlap + per-region detail with gene names |
| `R/10_figures.R` | the four report figures |
| `R/11_c151_peak_overlap.R` | do the c151 suggestive peaks land in the 3sp regions more than chance? |
| `R/12_overlap_emmax17.R` | same overlap on the EMMAX `l_min=3` 17-region set |
| `R/13_lift_and_enrich_emmax17.R` | lift those 17 regions (validated against R/08) + enrichment + c151 overlap |
| `R/14_venn_emmax17.R` | Venn of a region set vs the two specific EcoPeak sets |
| `R/15_overlap_lfmm.R` | same overlap on the LFMM `l_min=3` 136-region set |
| `R/16_peakset_concordance.R` | pairwise concordance incl. Kingman-vs-Kingman as the ceiling |
| `R/17_overlap_sweep.R` | τ_C × l_min sweep vs the EcoPeaks, with region-floor masking (`REPLOT=1` re-renders from the saved grid) |
| `R/18_cscore_enrichment.R` | per-SNP C-score vs EcoPeak membership, EMMAX and LFMM |
| `R/19_novel_regions_c151.R` | do the regions overlapping no EcoPeak carry c151 signal? (`PEAKSET=both\|c150`) |
| `R/20_fig_novel_c151.R` | figure for the above |
| `R/run_all.sh` | waits for the VCF download, then builds both cohorts end to end |

## Run

```bash
Rscript kingman2021/R/02_build_rds.R c155_global
```

Full rebuild from the joint VCF (long — the VCF sweep dominates):

```bash
bash kingman2021/R/run_all.sh
```

## Two things to know before using the bundles

**`GTs` is mean-imputed and numeric, not integer 0/1/2.** Mean coverage is 5.5×, so ~12% of
calls are missing at `MAF>=0.05 / F_MISSING<=0.20`, and tightening does not help
(`F_MISSING<=0.05` keeps 3.8% of SNPs). `gcta_grm()`/`emmax_setup()` have no NA handling.
Per-SNP `f_missing` is in `map`. Expect r²/`ld_w` attenuation from low-coverage hard calls;
the VCF retains `PL` if genotype posteriors are ever needed.

**Coordinates differ between the datasets.** 3sp is **gasAcu1** (`Chr1..Chr21`, arabic);
everything here is **gasAcu1-4** (`chrI..chrXXI`, roman). `R/08_liftover.sh` reconstructs
the chain from the hub's bigChain tracks (the published chains are behind Dryad's anti-bot
gate) and validates it: peaks lift at ~100% with preserved span, and the *Eda* peak shifts
by 12 kb.

## Tooling

`R/08`, `R/13` and `R/19` need UCSC `liftOver` (and `08` also `bigBedToBed`, `chainSwap`).
These now live durably at
`~/gitlab/LD-scaling-genome-scans/empirical_data/kingman2021/tools/` — add that to `PATH`
before running them. They were previously only in a session scratch directory, which was
cleared; that was a reproducibility gap, now closed.

## Report

`doc/kingman2021_report.pdf` — build with:

```bash
cd kingman2021/doc && latexmk -pdf kingman2021_report.tex
```

`doc/00_PROVENANCE.md` carries the full acquisition record, including what was blocked and
why, and every validation number.
