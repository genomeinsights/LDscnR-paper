# Kingman et al. 2021 stickleback — acquired data (clear-peak contrast dataset)

Paper: Roberts Kingman GA et al. 2021, *Predicting future from past: The genomic basis of
recurrent and rapid stickleback evolution*, Sci Adv 7:eabg5285. DOI 10.1126/sciadv.abg5285

Data root: `/Users/petrikem/gitlab/LD-scaling-genome-scans/empirical_data/kingman2021/`
(sibling of the existing `3sp/` store; note `LD-scaling-genome-scans/.gitignore` ignores
`empirical_data/*sp/` only, so this directory is **not** git-ignored yet — see §5.)

## 1. Headline: no WGS calling needed

The handoff's Step 0 ("try to skip the heavy calling") **succeeded**. A processed,
genome-wide, GATK-called + VQSR-filtered joint VCF for all 227 genomes is published on
FigShare as part of the assembly-hub mirror. The entire Step-2 path (download 1758 SRA
runs / 0.7 TB, align, joint-call) is therefore **superseded and was not run**.

The source that made this possible is not in the paper's Data Availability statement:
the UCSC hub description page <https://web.stanford.edu/group/kingsley/hubDescription.html>
points to **FigShare project 162634**, which mirrors every file of the assembly hub.

## 2. What was acquired

| Path | Contents |
|---|---|
| `vcf/227_genomes.final.filtered.vcf.gz` (+`.tbi`) | 15.3 GB joint VCF, 227 genomes, gasAcu1-4 coords, GATK `ApplyVQSR` → `SelectVariants --exclude-filtered`. Pre-filter `AC>=5 & AC<=450`. |
| `meta/tableS2_samples.tsv` | Table S2, 227 rows. Sequence ID, population, GPS, **M/F ecotype**, water type, and the analysis-cohort flags c150/c151/c153/c154/c155. |
| `tracks/*.bb` + `*.bed` | EcoPeaks (truth set) and TempoPeaks, converted with `bigBedToBed`. |
| `tracks/gasAcu1-4.scaledRABSrecombRates.bw` + `.bedGraph` | Rabbit Slough recombination rate in **cM/Mb** (LDhelmet rho-based map, 818,839 intervals, covers 416.3/463 Mb, mean 4.28 cM/Mb). |
| `tracks/gasAcu1-4_LW_recomb.2kb.bw`, `..._PS_recomb.2kb.bw` | Shanfelter et al. 2019 Lake Washington / Puget Sound fine-scale maps (alternates). |
| `tracks/gasAcu1-4.chrom.sizes` | 21 chromosomes + chrM/chrP/chrUn. |
| `ref/gasAcu1-4.2bit` | Reference assembly (117,585,540 B). |

### Peak sets (region counts and span verified against the paper)

| BED | Track | Regions | Span |
|---|---|---|---|
| `gasAcu1-4.c155.specific.50kb.final.peaks.bed` | **Global Specific EcoPeaks** | 39 | 3.74 Mb |
| `gasAcu1-4.c155.sensitive.50kb.final.peaks.bed` | Global Sensitive EcoPeaks | 92 | 24.62 Mb |
| `gasAcu1-4.c150.specific.50kb.final.peaks.bed` | **N.E. Pacific Specific EcoPeaks** | 209 | 27.44 Mb |
| `gasAcu1-4.c150.sensitive.50kb.final.peaks.bed` | N.E. Pacific Sensitive EcoPeaks | 212 | 91.90 Mb |
| `CH_SC_LB.specific.final.bed` | Cheney-Scout-Loberg Specific TempoPeaks | 344 | 18.43 Mb |
| `CH_SC_LB.sensitive.final.bed` | C-S-L Sensitive TempoPeaks | 524 | 99.27 Mb |

Paper text: Global specific 39 regions/3.7 Mb; Pacific specific 209/27.4 Mb; Pacific
sensitive 212/91.9 Mb. **All match.** BED cols 4–5 carry the SNP-based and window-based
p-values.

## 3. Cohorts — this is the key design fact

The 227 genomes are **one genome per population** for the EcoPeak analyses. The Table S2
cohort flag columns encode the exact cohorts the EcoPeaks were called on, with the flag
value being the ecotype (`1`=marine, `0`=freshwater) and blank = not in that cohort:

* `c155_global` — **84 samples: 56 freshwater + 28 marine.** Exactly reproduces the hub's
  "Global selection: 56 freshwater and 28 marine populations". → truth set = the 39
  Global Specific EcoPeaks (0.8 % of the genome).
* `c150_pacNW` — **68 samples: 57 freshwater + 11 marine.** Hub says "57 freshwater and 12
  marine"; one marine is short — worth a sanity check, but the freshwater side matches
  exactly. → truth set = the 209 Pacific Specific EcoPeaks (5.9 % of the genome).

All 227 samples: 50 marine / 177 freshwater. Mean sequencing coverage 5.5× (paper).

**Sample-ID check: the 227 Table S2 `Sequence ID` values match the 227 VCF sample names
exactly — set equality, no leftovers on either side.**

## 4. Sources that were NOT usable

* **Dryad** (`10.5061/dryad.pvmcvdnjm` array panel, `10.5061/dryad.547d7wm6t` reference +
  liftOver chains) — the API download endpoint now returns
  `{"error":"Unauthorized, must have current bearer token"}` and the web download route is
  behind an Anubis anti-bot challenge. Not circumvented. The reference was obtained from
  FigShare instead. **Still missing: the 363-SNP array PLINK trio and the gasAcu1↔gasAcu1-4
  liftOver chains.** Neither is needed — the array is a validation panel, not discovery
  data, and everything here is already in gasAcu1-4 coordinates so no liftOver is required.
* **`sbwdev.stanford.edu`** (the hub URL in the paper's Data Availability statement) — the
  host is fully unreachable: no ICMP reply, TCP connect times out on both 80 and 443, and an
  off-network proxy also fails. It appears to be down, not merely blocked locally.
  FigShare project 162634 and the GitHub mirror `chenhijy/gasAcu1-4_assemblyHub` both serve
  the same files.
* **SRA PRJNA247503** — 1758 runs / 206 individuals / 0.7 TB, median 3.34 Gb per individual
  (≈7.2× on a 463 Mb genome). Not downloaded; the joint VCF supersedes it. (The 227 VCF
  samples = these 206 plus 21 earlier Jones et al. 2012 genomes.)

## 5. Built outputs

`run_all.sh` → `01_extract_gts.sh` (bcftools: 21 named chromosomes only, biallelic SNPs,
within-cohort `MAF>=0.05`, `F_MISSING<=0.20`) → `02_build_rds.R` → `03_peaks_truth.R`.

| File | Contents |
|---|---|
| `module_sticklebacks/kingman2021_c155_global.rds` (219 MB) | 84 indiv × **3,509,797** SNPs, 28 M / 56 F |
| `module_sticklebacks/kingman2021_c150_pacNW.rds` (223 MB) | 68 indiv × **4,073,192** SNPs, 11 M / 57 F |
| `module_sticklebacks/kingman2021_<cohort>_truth.rds` | peak table + per-SNP peak-membership flags, **row-aligned with `map`** |

Each `.rds` is `list(GTs, map, eco, meta, cohort, source)`; `map` is
`marker, Chr, chr_num, Pos, rec_rate, cM, af, f_missing`. `GTs` is numeric and
mean-imputed (see the missingness note in `02_build_rds.R`): 11.56 % of calls imputed for
c155, 7.96 % for c150. The cohort flag was checked against the Table S2 ecotype column for
every sample in both cohorts and agrees.

Recombination map covers **96.4 %** of SNPs; the remainder (126,913 c155 / 145,009 c150)
fall in gaps of the RABS map and carry `rec_rate = cM = NA`. Anything working in genetic
distance must handle those.

### Validation: the packaged data reproduces the published peaks

Naive |ΔAF| between ecotypes, specific EcoPeaks vs the rest of the genome:

| Cohort | median |ΔAF| in-peak | background | SNPs above background 99.9th pct | Enrichment |
|---|---|---|---|---|
| c155_global | 0.179 | 0.058 | 29.8 % in-peak vs 0.10 % out | **298×** |
| c150_pacNW | 0.196 | 0.088 | 1.61 % in-peak vs 0.10 % out | 16× |

c155 is the sharper contrast, as expected: a balanced 28/56 design and a truth set covering
only 0.8 % of the genome. **Use `c155_global` as the primary clear-peak dataset**; c150 is
the secondary (more truth regions, but only 11 marine individuals).

Background quietness, the property this dataset was acquired for: on chrI:1–3 Mb (no
EcoPeak), 0 of 24,551 SNPs reach p < 1e-6 and 8 reach p < 1e-4 in an uncorrected
ecotype regression, λ ≈ 1.17.

## 6. Housekeeping note

**Resolved:** `empirical_data/kingman2021/` has since been added to
`LD-scaling-genome-scans/.gitignore`, so the 15 GB VCF cannot be staged by accident.
This module moved from `module_sticklebacks/kingman2021/` to the top-level
`LDscnR-paper/kingman2021/` (R/ data/ figures/ doc/) — script paths were repointed and
re-verified after the move.

### superseded note

`empirical_data/kingman2021/` is currently **not** covered by
`LD-scaling-genome-scans/.gitignore` (which ignores `empirical_data/*sp/`). The 15 GB VCF
must not be committed — add `empirical_data/kingman2021/` to that .gitignore before any
`git add` in that repo.

## 7. Re-derived EcoPeaks (c151 Northern Europe) + 3sp overlap

### Method
`04_ecopeaks_snp_test.py` re-implements Kingman's **SNP-based** EcoPeak test exactly as
documented in the hub's `ecoPeakDescription.html`: holding the genotype-class totals
(n0,n1,n2) and marine group size m fixed, the marine ALT-allele count is scored against a
multivariate-hypergeometric null. The null depends only on (n0,n1,n2,m), so it is computed
once per distinct key and cached — that is what makes it tractable genome-wide.
`05_call_peaks.R` applies BH-FDR and merges significant SNPs within 50 kb.

The **window/CSS** method (Kingman's second method) is *not* implemented — the published
description does not give the CSS formula. So these peaks are the SNP-only component;
Kingman's "specific" set is the intersection of SNP-based **and** window-based calls.

### Validation (both published cohorts)
| Cohort | Setting | My peaks | Published "specific" recovered |
|---|---|---|---|
| c155_global | BH q≤0.01 | 286 peaks / 15.5 Mb | **36/39 (92 %)** |
| c155_global | BH q≤1e-4 | 46 peaks / 3.83 Mb | 34/39 (87 %), and 46/46 of mine hit a published sensitive peak |
| c150_pacNW | BH q≤1e-4 | 274 peaks / 35.2 Mb | **195/209 (93 %)** |

At q≤1e-4 the c155 output (46 peaks, 3.83 Mb) is near-identical in size to the published
specific set (39 peaks, 3.74 Mb). That is the calibrated "specific-like" operating point.

### c151 cannot yield an FDR-significant peak set — this is structural
c151 has 9 marine / 18 freshwater, and after missingness the **median marine n actually
called is 7** (range 2–9). The exact test therefore bottoms out at **p = 3.2e-7**, attained
by exactly 1 SNP. BH at q≤0.01 over 1.70 M tests would need ≥54 SNPs at that floor.
**Zero SNPs are significant at any FDR.** This is almost certainly why the hub publishes
EcoPeak BEDs for c150 and c155 only, although Table S2 also defines c151/c153/c154.

`06_c151_suggestive_peaks.R` therefore emits a **rank-based, NOT FDR-controlled** set:
top 0.1 % of p → 1,770 SNPs → **41 suggestive peaks, 1.31 Mb** (`c151_nEur_suggestivePeaks.bed`).

### The geographic gradient (`07_signal_enrichment.R`)
Rank-based enrichment of low Kingman p inside the 39 lifted 3sp LFMM regions (4.6–4.8 % of
tested SNPs), permutation null shuffling regions within chromosome:

| Kingman cohort | Geography vs 3sp | top 1 % of p | p | top 0.1 % | p |
|---|---|---|---|---|---|
| **c151 N. Europe** | **matches** | **4.64×** | 0.0040 | 1.07× | 0.28 |
| c155 Global | partly | 3.19× | 0.024 | 6.30× | 0.0080 |
| c150 N.E. Pacific | different basin | 1.47× | 0.22 | 0.71× | 0.39 |

And 11 of the 41 c151 suggestive peaks fall inside a 3sp LFMM region vs 2.8 expected —
**3.93×, p = 0.0002** (`R/11_c151_peak_overlap.R`).

Caveat: the top-0.1 % column is unreliable for c151, whose 1.70 M p-values take only 19,035
distinct values; that extreme tail is shaped by which genotype configurations can attain a
low p rather than by signal. The top-1 % column is the trustworthy comparison.


## 8. The l_min=3 region sets, and the Eda correction

The §7 comparison used the 40 LFMM `l_min=10` regions. The main manuscript's structure-null
analysis works from `l_min=3` sets, so the overlap was recomputed on those
(`R/12_overlap_emmax17.R`, `R/13_lift_and_enrich_emmax17.R`, `R/15_overlap_lfmm.R`).

| Set | l_min | Regions | Span | Global-specific hit | Fold | p |
|---|---|---|---|---|---|---|
| **EMMAX** | 3 | 17 | **1.21 Mb** | **5/17** | **19.5×** | **0.0005** |
| LFMM | 10 | 40 | 21.13 Mb | 4/40 | 1.38× | 0.33 |
| LFMM | 3 | 136 | 37.27 Mb | 9/136 | 1.87× | 0.036 |

The EMMAX `l_min=3` set claims 0.3% of the genome and hits 5 of the 39 Global-specific
EcoPeaks. The LFMM `l_min=3` set claims 30× more genome for a 1.87× enrichment — the
difference is precision, not coverage.

Rank-based enrichment inside the 17 EMMAX regions: c151 34.1× (p=0.0020), c155 25.9×
(p=0.0040), c150 17.6× (p=0.0020) at top 1%; c155 33.5×, c151 11.3×, c150 11.1× at top 0.1%.
c151 suggestive peaks inside the 17 regions: 5/41 vs 0.33 expected, 15.1×, p=0.0002.

**Do not compare folds across region sets** — fold scales inversely with claimed span, so
34× on 1.21 Mb and 4.64× on 19.92 Mb are not commensurable. The geographic gradient
reproduces only in the top-1% column; in the top-0.1% column c155 is highest in *both*
region sets, because c155 is the best-powered cohort and the extreme tail therefore ranks
power rather than geography.

### Eda — §7's "genuinely missed" was wrong

| Method | l_min | Eda recovered? | Why |
|---|---|---|---|
| EMMAX | 3 | **yes** | region 7 = Chr4:12,809,075–12,811,617, **4 SNPs** |
| EMMAX | 10 | no | same region, below the cluster-size floor |
| LFMM | 3 | no | Chr4 regions jump 12.284 → 14.034 Mb |
| LFMM | 10 | no | same gap |

The EMMAX region overlaps the `eda` gene body by **1,371 bp** and sits inside the
Global-specific EcoPeak with the most extreme SNP p-value in that set (7.9e-15; Pacific
2.2e-12). So it is a **method difference, not merely a threshold choice** — lowering l_min
does not rescue LFMM. This vindicates `module_sticklebacks_LDscnR/README.md`'s claim that
the C-score resolves Chr4/Eda; the earlier reading was off the wrong region set.

### A liftOver bug caught by validation

`R/13` originally called `liftOver` with the default `-minMatch` (0.95) while `R/08` — which
produced the reference `lfmm_g14.bed` it validates against — used `-minMatch=0.5`. Three of
39 regions failed to lift and the assertion halted the script; the 36 that lifted matched
the reference to **0 bp**. Fixed by adding `-minMatch=0.5`; validation now passes 39/39 at
0 bp. (An earlier attempt, preserved in `data/emmax17_enrich.log`, used a pure-R chain
lifter that disagreed with `liftOver` by up to 8.3 Mb and correctly aborted.)


## 9. How much agreement is achievable? (`R/16_peakset_concordance.R`)

Fold enrichment has a pair-specific ceiling and Jaccard is capped by the span ratio, so
cross-pair comparison uses the **overlap coefficient** = bp intersection / bp of the
*smaller* set (1.0 = the smaller claim is entirely corroborated).

| A | B | span A | span B | ∩ Mb | overlap coef | p |
|---|---|---|---|---|---|---|
| 3sp EMMAX l_min=3 | 3sp LFMM l_min=3 | 1.21 | 37.27 | 1.158 | **0.959** | 5e-4 |
| Kingman Global-spec | Kingman Pacific-spec | 3.52 | 24.96 | 3.340 | **0.950** | 5e-4 |
| Kingman Global-sens | Kingman Pacific-sens | 22.09 | 82.68 | 20.972 | 0.949 | 5e-4 |
| 3sp EMMAX l_min=3 | 3sp LFMM l_min=10 | 1.21 | 21.13 | 1.046 | 0.865 | 5e-4 |
| 3sp EMMAX l_min=3 | Kingman Pacific-spec | 1.21 | 24.96 | 0.536 | 0.444 | 0.037 |
| 3sp EMMAX l_min=3 | Kingman Global-spec | 1.21 | 3.52 | 0.497 | 0.411 | 5e-4 |
| 3sp LFMM l_min=3 | Kingman Global-spec | 37.27 | 3.52 | 0.972 | 0.277 | 0.049 |
| 3sp LFMM l_min=3 | Kingman Pacific-spec | 37.27 | 24.96 | 4.792 | 0.192 | 0.046 |
| 3sp LFMM l_min=10 | Kingman Pacific-spec | 21.13 | 24.96 | 3.173 | 0.150 | 0.028 |
| 3sp LFMM l_min=10 | Kingman Global-spec | 21.13 | 3.52 | 0.476 | 0.135 | 0.32 |

**Kingman's two published sets are essentially nested (0.950)** — the Global peaks are the
stringent subset of the Pacific peaks. That is the ceiling. Our best cross-dataset value is
0.444, less than half, so the 3sp sets are not reproducing Kingman at the achievable limit.

The only value matching the ceiling is our own two methods agreeing with each other (0.959),
which is a *within-dataset* comparison over a wider span ratio (1.21 in 37.27 vs 3.52 in
24.96) and so an easier bar. **No analysis of ours agrees with Kingman as well as Kingman's
two datasets agree with each other.**


## 10. Calibrating against the EcoPeaks (`R/17`, `R/18`)

Sweeping the τ_C × l_min grid and scoring each cell against the Global-specific EcoPeaks.

**Region-floor masking is essential.** Precision and rotation-null fold both divide by the
number of regions in a cell, so a cell holding one region that happens to hit scores
precision 1.0 and a fold limited only by the near-empty null (up to 200× on this grid).
Unmasked, these degenerate cells form a bright band across the high-τ_C / high-l_min
corner and invert the apparent optimum. `R/17` masks cells with `n_regions < 5` (a
pragmatic floor, not a derived one) and returns NA rather than dividing by a 1e-6 guard.

| Sweep | Estimable | F1-optimal | F1 max / at (0.05, 3) | Rank |
|---|---|---|---|---|
| EMMAX / global-spec | 33/270 | τ=0.10, l=1 | 0.197 / 0.189 | 3 |
| EMMAX / Pacific-spec | 33/270 | τ=0.05, l=1 | 0.092 / 0.048 | 2 |
| LFMM / global-spec | 247/270 | τ=0.10, l=6 | 0.135 / 0.102 | 15 |
| LFMM / Pacific-spec | 247/270 | τ=0.05, l=5 | 0.217 / 0.212 | 4 |

(0.05, 3) is not the F1 maximum in any of the four. Precision peaks in the low-τ_C corner
(0.40 at τ=0.05, l=6–9); fold does not — it favours the l_min=1 edge at high τ_C (46× at
τ=0.70), in cells at the floor. Recall is capped near 0.28 by biology, so F1 is a poor
objective here.

**33 vs 247 estimable cells is the headline.** Regions at l_min=3: EMMAX 17 → 8 → 2 → 1 → 0
across τ_C = 0.05/0.10/0.30/0.60/0.90; LFMM 136 → 133 → 107 → 84 → 29. Per-SNP, 782 LFMM
SNPs carry C > 0.5 against EMMAX's 13. And no threshold rescues LFMM: precision stays
0.029–0.068 and fold 1.8–6.4 across all 18 estimable τ_C, vs EMMAX 0.125–0.294 and
9.7–18.8 over its three.

**Per-SNP C-score.** The informative signal is the C=0 → C>0 step (EMMAX 24–53× over the
1.31% base; C=0 sits at 0.91×). Within C>0 neither method rises — EMMAX peaks in (0.05,0.1]
then falls to 23.5×, LFMM falls monotonically 21.7 → 3.8×. The Spearman ρ (0.21 / 0.06) is
the step restated, not a trend. Top bin: EMMAX 23.5× on 13 SNPs, LFMM 3.8× on 782.

**Open item: LFMM has no structure-aware null.** Its τ_C is transferred from EMMAX, so the
engine comparison is confounded. Note λ_GC(LFMM) = 1.008 vs EMMAX 1.093 — LFMM is *not*
inflated in the usual sense (genomic control forces the median); it has a heavy tail (207
SNPs at p<1e-6 vs 13; 952 at BH q<0.05 vs 3). The sweep shows LFMM's C does not *rank*
toward the EcoPeaks; a null would show whether its C is unusual *at all*. Those can
dissociate — LFMM could pass its negative control while failing the positive one, which
would mean it detects something real and structured that is not marine–freshwater
adaptation. `module_sticklebacks_LDscnR/lfmm_permutation_null_3sp.R` exists but has never
produced output; B ≈ 5–10 would settle the direction.


## 11. Provenance of `data/liftover/pv_*.specific.bed` (cite this section)

These four BEDs are the paper's external validation set. Three facts about them are asked
often enough to answer in one place.

**They are grep-invisible.** `R/08_liftover.sh` line 59 writes them as
`liftOver -minMatch=0.5 pin.bed gasAcu1-4ToGasAcu1.chain "pv_$s.bed" "pv_un_$s.bed"`
inside a loop over `c155.specific c155.sensitive c150.specific c150.sensitive`. The
filename is **constructed from a loop variable**, so grepping either repository for
`.specific.bed` returns nothing. That is not evidence of missing provenance.

**We did not call these peaks.** They are Roberts Kingman *et al.* 2021's published
EcoPeaks, downloaded as bigBed from the FigShare mirror of the authors' UCSC assembly hub
(project 162634; the `sbwdev.stanford.edu` host in the paper's Data Availability statement
is dead), converted with `bigBedToBed`, and lifted. So there is no threshold of ours to
report. **The authors' thresholds** (their §"EcoPeak identification"): two methods — a
CSS/genetic-distance test on 2500 bp windows and a per-base multivariate-hypergeometric
allele-count test — with *specific* = both agreeing at 1% FDR and *sensitive* = either at
5% FDR, then merged at 50 kb. Cols 4–5 of each BED carry the published SNP-based and
window-based p-values. Counts and spans were checked against the paper and match exactly
(§2).

**Chain and assembly.** Source `gasAcu1-4`; target `gasAcu1`, because the 3sp panel is in
gasAcu1 (`Chr1..Chr21`, arabic) while everything Kingman is gasAcu1-4 (`chrI..chrXXI`,
roman). The published chains sit behind Dryad's anti-bot gate, so
`gasAcu1-4ToGasAcu1.chain` was **reconstructed** from the hub's bigChain + bigChain.link
tracks (1,095 chains / 19,851 blocks) and validated: every peak set lifts at ~100% with
preserved span, chromosome assignment is preserved, and the *Eda* peak moves 12 kb. See
`data/liftover/README.md` and §4.

**Regeneration status, checked 2026-09-03.** The BEDs are all 13 tracked in git, so
their integrity does not depend on anyone's disk. They are **not, however, rebuildable
on this machine**: `08_liftover.sh` needs `bigBedToBed`, `liftOver` and `chainSwap` on
`PATH` and none of the three is installed anywhere under `~`, `/usr/local` or `/opt`.
The script would fail at its first pipeline stage. Install from
<http://hgdownload.soe.ucsc.edu/admin/exe/> before attempting any re-derivation, and
note the script also fetches two FigShare objects at run time, so it needs network.

This is recorded rather than fixed because the BEDs are tracked and no re-derivation is
pending. It matters only if the peaks are ever questioned: the answer to “can you
regenerate these?” is currently “yes, after installing three tools”, not “yes”.

**What the cohort labels denote** (§3, from Table S2): `c155` = the **global** cohort, 84
genomes, 28 marine / 56 freshwater, one genome per population. `c150` = the **N.E.
Pacific** cohort, 68 genomes, 11 marine / 57 freshwater. Both are one-genome-per-population
designs, which is why they are panels rather than population samples.

**Note for anyone summing peak widths across both cohort BEDs:** the two sets are largely
nested, not independent — §9 measures 95% of the Global-specific sequence lying inside the
Pacific-specific set. Summing unmerged widths across both therefore double-counts. Merge
within chromosome first.


## 12. Do the uncorroborated regions carry signal? (`R/19`, `R/20`)

Tests the reading that LFMM's surplus regions are real Atlantic loci the Kingman cohorts
cannot see. c151 (Northern Europe, 9 M / 18 F) is geography-matched to the 3sp panel and was
not used to call the published EcoPeaks.

| Classified by | Group | n | Span | Fold | p |
|---|---|---|---|---|---|
| c155+c150 (deflationary) | EMMAX corroborated | 7 | 0.87 Mb | 23.6× | 0.002 |
| | EMMAX novel | 10 | 0.33 Mb | 2.02× | 0.072 |
| | LFMM corroborated | 27 | 14.72 Mb | 4.33× | 0.016 |
| | **LFMM novel** | 107 | 20.98 Mb | **0.98×** | 0.27 |
| c150 only (inflationary) | EMMAX novel | 11 | 0.41 Mb | 62.5× | 0.002 |
| | **LFMM novel** | 107 | 20.98 Mb | **0.98×** | 0.27 |

**Circularity is unavoidable**: all 27 c151 samples are a subset of c155's 84, so
classifying by c155 is deflationary for the novel class and classifying by c150 alone is
inflationary. Both reported. The LFMM novel result is **identical under both**, which is
what makes it usable; the EMMAX novel figure is not, so 62.5× is not quoted (one region
crossing the boundary moves it from 2.02×).

The corroborated groups are the internal positive control and both are significant, so c151
has the power and a flat result is a real negative.

**Limits.** Zero c151 signal is equally consistent with false positives and with real
structure that is not marine–freshwater (ecotype-correlated demography). Only "carries no
detectable marine–freshwater signal" is established. Says nothing about whether those
regions clear their own structure null (§10, still open).

**Reconciliation with the sims.** LFMM buys recall at ~20% precision cost in simulation
(recall 0.325 vs 0.231; precision 0.291–0.337 vs 0.351–0.359) but its real-genome surplus
carries no signal at all. The gap between regimes is the structure-model-dependence
evidence: latent factors capture simple simulated structure, not continental IBD.

## 13. Tooling note

`R/08`, `R/13`, `R/19` require UCSC `liftOver` (plus `bigBedToBed`, `chainSwap` for `08`).
These were originally only in a session scratch directory under `/tmp`, which was cleared
between sessions — the scripts then failed with "liftOver not on PATH". Now installed
durably at `empirical_data/kingman2021/tools/`; add to `PATH` before running.
