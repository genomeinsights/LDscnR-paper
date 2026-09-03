# kingman2021 — liveness inventory

My half of the union needed before the consolidation deletes anything. Scope: everything
under `LDscnR-paper/kingman2021/` plus the bulk data at
`~/gitlab/LD-scaling-genome-scans/empirical_data/kingman2021/`.

**Verdict: every script in `R/` is live.** None is a candidate for deletion.

A note on the criterion, since it is what the consolidation got wrong: **"nothing imports
its output" does not mean dead.** These are analysis endpoints — they write a CSV or a
figure that a document cites. Having no downstream importer is the defining shape of the
class we are trying to keep, so call-graph reachability cannot distinguish live from dead
here. The test that works is: *does a document or another script cite its output?*

## Chain of production (each row consumed by the rows below it)

| Script | Writes | Consumed by |
|---|---|---|
| `01_extract_gts.sh` | `gts/<cohort>.{samples.txt,gt.tsv.gz,sites.tsv.gz}` | 02, 04 |
| `02_build_rds.R` | `data/kingman2021_<cohort>.rds` | 03 |
| `03_peaks_truth.R` | `data/kingman2021_<cohort>_truth.rds` | truth set for scoring |
| `04_ecopeaks_snp_test.py` | `ecopeaks/<cohort>.snp_p.tsv.gz` | 05, 06, 07, 11, 13, 19 |
| `05_call_peaks.R` | `ecopeaks/*_snpEcoPeaks_fdr*.bed` | validation vs published peaks (report §7) |
| `06_c151_suggestive_peaks.R` | `data/peaks/c151_nEur_suggestivePeaks.bed` | 11, 13, 19 |
| `07_signal_enrichment.R` | stdout | superseded in report by 09, kept as the standalone form |
| **`08_liftover.sh`** | **chains + `data/liftover/pv_*.bed`** | **09, 12, 13, 15, 16, 17, 18, 19 — and 11 of ldscnr-2c's 13 load-bearing scripts** |
| `09_overlap.R` | `overlap_summary.csv`, `overlap_detail.csv`, `enrichment_summary.csv` | report §5, §7 |
| `10_figures.R` | `fig1`–`fig4` | report |
| `11_c151_peak_overlap.R` | `c151_peak_overlap.csv` | report §7 |
| `12_overlap_emmax17.R` | `regions_…_emmax.csv`, `overlap_{summary,detail}_emmax17.csv` | report §8; feeds 13, 14, 16, 19 |
| `13_lift_and_enrich_emmax17.R` | `liftover/emmax17_g14.bed`, `enrichment_emmax17.csv`, `c151_peak_overlap_emmax17.csv` | report §8; feeds 19 |
| `14_venn_emmax17.R` | `fig5_venn_{emmax17,lfmm}.png` | report §8 **and the manuscript** |
| `15_overlap_lfmm.R` | `regions_…_lfmm.csv`, `overlap_{summary,detail}_lfmm.csv` | report §8; feeds 16, 19 |
| `16_peakset_concordance.R` | `peakset_concordance.csv` | report §9 (the concordance ceiling) |
| `17_overlap_sweep.R` | `overlap_sweep_*.csv`, `kingman_sweep_*.png` | report §10 **and the manuscript** |
| `18_cscore_enrichment.R` | `cscore_enrichment.csv`, `kingman_cscore_enrichment.png` | report §10 **and the manuscript** |
| `19_novel_regions_c151.R` | `novel_regions_c151_{both,c150}.csv` | latest result; report update pending |

## Read this before deleting anything

**`08_liftover.sh` is the load-bearing script nobody has listed.** It reconstructs the
gasAcu1↔gasAcu1-4 chains and produces `data/liftover/pv_c155.specific.bed` and
`pv_c150.specific.bed` — the files carrying **every external validation claim on the panel
side**, in both my inventory and ldscnr-2c's. It appeared on no deletion list only because
it was never traced: it writes its outputs as `"pv_$s.bed"` from a loop variable, so a grep
for `.specific.bed` finds no producer and the BEDs look like unregenerable supplied
artefacts. They are not. Provenance is `doc/00_PROVENANCE.md` §11.

**Three kingman2021 figures are in the manuscript right now** —
`fig5_venn_emmax17.png` (14), `kingman_cscore_enrichment.png` (18),
`kingman_sweep_EMMAX_c155specific.png` (17).

**Two files outside my folder should be preserved, not deleted:**
`module_sticklebacks_LDscnR/kingman_overlap_sweep.R` and `kingman_cscore_enrichment.R`.
Superseded by my 17 and 18 and they generate figures with known errors, so they must not be
re-run — but they are the provenance of figures that were in the manuscript until 19 Aug,
and report §10 records what was wrong with them. Move to a `superseded/` directory; do not
delete.

## Second-form dependencies — inputs whose PRODUCERS are gone

Applying ldscnr-2c's corrected rule (trace one level *past* each claim's producing script,
to what that script reads) to my own inventory found two I had missed. Both data files
survive; both **producers were deleted** in commit `3cb47ba`.

| Input | Read by | Producer | Status |
|---|---|---|---|
| `data/inputs/null_popperm_3sp.rds` | `R/12` → the **17 EMMAX regions**, the manuscript's headline set | `module_sticklebacks_LDscnR/permutation_null_3sp.R` | producer deleted; file **vendored** into `data/inputs/` and `R/12` repointed |
| `data/regions_tau0.05_lmin10_rho0.60.csv` | `R/08`, `R/09`, `R/16` | `module_sticklebacks_LDscnR/manhattan_regions.R` | producer deleted; CSV survives in `data/` |

Both producers are recoverable: `git show 3cb47ba^:module_sticklebacks_LDscnR/<file>`.
Neither input is regenerable from anything now on disk, so restore the producer first if
either ever needs re-deriving. This is the same blind spot as `08_liftover.sh`, one level
further out: a backwards trace from claims stops at the script producing the cited output
and never asks what that script *read*.

## External inputs this folder depends on

- `3sp_data/3sp_LDscnR_data.rds` (887 MB) — loaded by 09, 12, 15, 16, 17, 18. Repointed
  after the move out of `module_sticklebacks_LDscnR/data/`; verified by re-running 18.
- `~/gitlab/LD-scaling-genome-scans/empirical_data/kingman2021/` — the 15 GB VCF, tracks,
  reference, per-cohort extracts. Git-ignored in that repo. Not regenerable cheaply: the
  VCF is a 15 GB download and `04` is ~9 h of compute per cohort.
