# Instructions for the final LDscnR simulation re-analysis

## Objective

Rebuild the analyses in `~/gitlab/LDscnR-paper/module_sim_3sp53` as a small,
auditable manuscript pipeline. The analysis should answer one central question:

> Does phenotype-blind Stage-1 LD complexity reduction reduce false discoveries
> and improve precision relative to unrestricted marker-wise testing, while
> retaining useful recall across selection, dispersal and background-selection
> conditions?

The final pipeline should not attempt to preserve every exploratory analysis.
Keep the existing files and git history available, but make a clean distinction
between manuscript analyses and exploratory work.

Do not begin the full 1,400-run computation until the parser correction and the
small-grid validation gates below pass.

## Update after inspection of the Stage-2 Manhattan figures

The corrected parser and the full EMMAX/LFMM grids have now been run under
`final_analysis/`. Do not rerun those expensive association stages solely for
this update. Reuse their saved p-values, BH decisions and Stage-1 objects, then
rerun Stage-2 assembly, truth scoring, pooling and figures for all 1,400
combinations.

The present implementation is incomplete for the revised primary estimand:

- `R/05_score_truth.R` creates Stage-2 region arms for Simes and consensus but
  not for the unrestricted marker-wise comparators.
- Its current `.region_members()` assigns every assayed marker between a
  region's coordinate bounds to that region. Replace this with the actual
  constituent discovered Stage-1 cluster members so intervening markers cannot
  lend truth credit.
- `R/06_summarise.R` currently treats region-versus-unit results as a secondary
  contrast. Make comparable Stage-2 region results the primary performance
  table and retain marker/unit results as diagnostics.
- `figureS_simulation_manhattan_tpfp_stage2.R` currently prints counts of
  coloured region-member markers as “TP” and “FP”. Print TP/FP region counts
  prominently; label coloured-marker totals separately if retained.

Before pooling, validate the revised assembly/scoring on the illustrative
`nobgs/V0.5_c1/env3` data and print, for each engine and method: significant
hypotheses, Stage-2 regions, TP regions, FP regions and recovered QTN. Confirm
that every reported region maps to at least one discovered seed cluster and
that each seed cluster maps to exactly one reported region.

## Sources of truth

Use these in this order:

1. `~/gitlab/LDscnR-NEMO/make_prod.sh`, its production templates, and its current
   README for the generating design and NEMO locus convention.
2. Native NEMO `.map` and `.snp_geno` files plus the matching
   `params_3spC_53cM` reference maps and environmental surfaces.
3. The current LDscnR package implementation in `~/gitlab/LDscnR`.
4. `~/gitlab/LDscnR_manuscript/materials_and_methods.tex` for the intended
   inferential workflow.
5. Existing `module_sim_3sp53` code only as reusable implementation material,
   not as authoritative documentation.

Do not use old parsed bundles or summary RDS files as data sources. They were
created with a parser whose locus-index convention appears to be wrong.

## Preserve existing work

- Inspect `git status` before editing. There is currently an uncommitted change
  to `R/14_random_removal_control.R`; do not overwrite, move, restore or reformat
  that file.
- Do not delete old outputs. Write the rebuilt pipeline to a new versioned output
  root such as `out_final_v1/` and the parsed data to a correspondingly new root.
- Initially place the clean implementation in a clearly named subtree such as
  `final_analysis/`. After it has reproduced the complete result, propose a
  separate cleanup/promotion step; do not mix that with the statistical rebuild.
- Raw NEMO output is read-only.

## Analyses that belong in the final story

### Primary EMMAX comparison

Run all three methods on the same MAF-filtered genotypes and with the same
Stage-1-derived GCTA relationship matrix:

1. `emmax_snp`: unrestricted marker-wise EMMAX, with BH applied across all
   tested markers.
2. `emmax_simes`: marker-wise EMMAX p-values combined within phenotype-blind
   Stage-1 units by Simes, followed by BH across the tested Stage-1 units.
3. `emmax_consensus`: one consensus dosage per eligible Stage-1 unit, tested by
   EMMAX, followed by BH across units.

Stage 1 uses the manuscript defaults: `rho = 0.50`, MAF greater than 0.10, and a
minimum unit size of two markers. Stage 2 must not alter p-values, BH correction,
or the rejection set. It assembles discovered Stage-1 units into the candidate
regions that LDscnR ultimately reports. These Stage-2 regions are the primary
TP/FP scoring and reporting unit in the simulations. Marker- and Stage-1-unit-level
scores are retained only to show why region assembly is necessary.

### Testing unit versus reported-call unit

Use precise terminology throughout the code, figures and manuscript:

- BH is performed on the hypotheses actually tested: individual markers for the
  unrestricted scan, and Stage-1 units for Simes and consensus testing.
- `assembly = "stage2_discovered"` is applied only after those BH decisions.
- A resulting Stage-2 region is one reported candidate call and one primary
  TP/FP scoring unit. It is not a new hypothesis whose p-value was included in
  BH, because its membership depends on which Stage-1 units were significant.

Do not describe a post hoc Stage-2 region as “one tested hypothesis.” If the aim
were instead to make Stage 2 the literal inferential test unit, Stage 2 would have
to be constructed genome-wide without reference to significance before any
association test. That is a different method and is not what the supplied
Manhattan figures or `stage2_discovered` implementation do.

Apply a comparable post hoc assembly to every method before primary performance
comparison:

1. For Simes and consensus, seed Stage 2 with the BH-significant Stage-1 units.
2. For unrestricted marker-wise testing, mark every phenotype-blind Stage-1
   cluster, including singleton clusters, as discovered when it contains at least
   one BH-significant marker; then run the same Stage-2 assembly on those discovered
   clusters. This changes only reporting, never marker p-values or marker-wise BH.
3. Preserve the mapping from every reported region to its constituent discovered
   Stage-1 clusters and significant markers.

This prevents an unfair comparison in which one side counts every significant SNP
as a separate false call while the other side counts an assembled genomic region.

### Secondary association-engine check

Retain LFMM only as a portability analysis:

- unrestricted marker-wise LFMM2 plus BH;
- the same marker p-values combined within Stage-1 units by Simes plus BH;
- no LFMM consensus-dosage analysis.

Split LFMM into a separate stage so that the primary EMMAX analysis can finish
and be audited without rerunning LFMM. Before the full LFMM grid, verify the
choice `K = 5` using a structure diagnostic that is independent of association
truth. Do not choose K by maximizing QTN recovery. If K=5 remains a legacy fixed
choice rather than a newly justified one, label it as such in the receipt and
manuscript output.

### One compact BGS validation

Because BGS is a designed factor, retain one supplementary quality-control
analysis showing that it produced the intended genomic effect. Summarise paired
BGS versus no-BGS differences in nucleotide diversity overall and across
recombination-rate strata. One compact figure or table is sufficient. Do not
carry forward several overlapping BGS figures.

### Optional, truth-based size diagnostic

If useful after the main results are complete, report the proportion of
significant Stage-1 units linked to truth as a function of unit size, together
with the number of units in every bin. This is supplementary and descriptive.
Do not infer from it that a size threshold is self-calibrating.

## Analyses not required for the current manuscript

Do not rerun or place in the final outputs unless a later, explicit manuscript
question requires them:

- structured-null simulations (`R/12_*` and `R/13_*`);
- random-removal controls (`R/14_*`);
- size-floor permutation sweeps;
- `F_beta` curves;
- multiple overlapping BGS/recombination figures;
- range-expansion runs;
- a chromosome-2 “neutral chromosome” analysis.

The empirical three- and nine-spined panels already carry the paper's permutation
message. The simulations have known truth and should be used primarily to assess
recovery and false discoveries directly.

Chromosome 2 is not neutral by construction: the generating map contains one
potential QTN there. Do not call every chromosome-2 discovery a false positive.

## Phase 1: rebuild and validate the parser

Correct `R_parsing/01_parse_nemo.R` or replace it in the clean pipeline.

Nemo locus identifiers are one-based. In particular, `quant.101` exists when
`quanti_loci = 101`. Remove the current unconditional `+1` transformation before
joining loci to the type-specific reference-map rows.

Add hard assertions and a compact parser-QA table. At minimum verify:

- each NEMO genotype column has exactly one NEMO map record;
- every NEMO map record used in the genotype matrix has exactly one reference-map
  match within its locus type;
- type-specific indices are integers in `[1, n_loci_of_type]`;
- `quant.100` maps to chromosome 1 and `quant.101` maps to the single chromosome-2
  QTN in the unfiltered reference map;
- joined genotype columns and map rows remain in identical order after sorting
  and duplicate-position handling;
- chromosome, physical position, locus type and allelic effect come from the
  matched reference row;
- genotype values are valid diploid dosages and sample metadata remain aligned;
- the environmental value assigned to every individual matches its sampled
  population and the selected environmental surface;
- counts by chromosome and locus type are reported before and after removal of
  monomorphic loci, individual subsampling, duplicate positions and MAF filtering.

On at least three deliberately chosen native runs, compare the old and corrected
join. Print the mappings for the final few QTN indices so that the correction is
human-verifiable. Genotypes should be unchanged; coordinates, types and QTN truth
must be rebuilt.

Do not merely edit a comment. The receipt/fingerprint for all downstream stages
must change, and no old parsed or scored output may satisfy the new receipts.

## Phase 2: define the clean pipeline

Use a short, linear dependency graph. A suitable layout is:

```text
final_analysis/
  README.md
  R/
    00_config.R
    01_parse_nemo.R
    02_build_ld_units.R
    03_emmax.R
    04_lfmm.R
    05_score_truth.R
    06_summarise.R
    07_export_manuscript.R
  R_figures/
    figure_simulation_performance.R
    figureS_simulation_qc.R
    figureS_truth_sensitivity.R
  run_smoke_test.sh
  run_full_grid.sh
```

Centralise every result-affecting value in `00_config.R`. No stage script should
contain a second hard-coded copy. Include paths, cells, tags, maps, environments,
MAF, LD-decay settings, Stage-1 rho, size floor, alpha, GRM estimator, truth
thresholds, seeds, bootstrap settings, software versions and the LDscnR git SHA.

Parameterise all host-dependent paths. Do not assume that
`/Volumes/Large_storage`, `/Volumes/Nemo`, `~/LDscnR-NEMO` or
`~/gitlab/LDscnR-NEMO` exists on every machine. The driver should fail with a
clear message showing the missing configured path.

The full-grid driver must build a manifest of exactly
`7 cells x 2 BGS treatments x 10 map IDs x 10 environments = 1,400` combinations.
For the final run, missing inputs are failures, not silent skips. Each output is
keyed by all four identifiers. Preserve deterministic seeds and record status and
wall time per stage and combination.

Cache expensive phenotype-blind objects and association statistics separately.
Truth thresholds and summary choices should be rescorable without recomputing LD,
the GRM, EMMAX or LFMM.

## Phase 3: smoke-test gates

Run these gates in order:

1. Parser only on three runs spanning both BGS treatments, both chromosomes and
   at least two map IDs. Inspect the QTN mapping table manually.
2. End-to-end EMMAX on one combination. Confirm identical marker order through
   parsing, LD construction, GRM construction, association and truth scoring.
3. A 14-combination grid: all seven parameter cells, both BGS treatments, one
   map and one environment. Confirm that all expected methods and summary fields
   appear and that no partial-grid pooling is possible without an explicit
   `allow_partial = TRUE` development flag.
4. Only after gates 1--3 pass, run the full EMMAX grid.
5. Run LFMM separately after its K decision is documented.

Write a machine-readable validation report and stop immediately on a failed
invariant.

## Truth definition and scoring

For each analysed sample, calculate QTN additive-variance contribution as

```text
Va_j = 2 p_j (1 - p_j) a_j^2
```

where `p_j` is calculated from the analysed 160 individuals. The current primary
detectability definition is MAF greater than 0.10 and at least 5% of the total
QTN additive variance on that chromosome. Keep this as the manuscript operating
definition unless the rebuilt data expose a failure, and make it a named config
value rather than an inline constant.

First calculate marker- and Stage-1-unit-level truth linkage as diagnostics:

- A significant marker is truth-linked when that marker satisfies the specified
  LD and distance criteria to at least one detectable QTN.
- A significant Stage-1 unit is truth-linked when at least one of its member
  markers satisfies the same criteria.
- Diagnostic hypothesis-level precision is the number of significant hypotheses
  that are truth-linked divided by all significant hypotheses. Do not silently
  discard duplicate links to an already recovered QTN from this denominator.

The primary manuscript estimand is based on Stage-2 reported regions:

- A reported region is a TP when at least one marker belonging to one of its
  constituent discovered Stage-1 clusters satisfies the LD and distance criteria
  to a detectable QTN. Otherwise it is an FP.
- Use constituent cluster members for truth assignment. Do not use every assayed
  marker lying physically between the region bounds: an intervening, untested
  marker could otherwise lend truth credit to the region.
- Region precision is `TP regions / (TP regions + FP regions)`.
- Region recall is the number of unique detectable QTN recovered by at least one
  reported region divided by all detectable QTN.
- When several discovered Stage-1 units merge into one region, count the region
  once. A region containing both truth-linked and non-linked constituent units is
  one TP candidate region, not one TP plus one or more FP calls.
- A QTN that cannot be represented by any eligible Stage-1 unit remains in the
  primary recall denominator. Also report conditional recall among QTN covered by
  eligible units as a diagnostic. This separates coverage loss from test-power
  loss.
- Record TP regions, FP regions, FN QTN, number of inferential tests, number of
  significant inferential hypotheses, number of assembled regions, BH critical
  p-value, eligible-marker fraction, and detectable-QTN coverage for every method
  and combination.

Use the intended decay-relative matching thresholds (`rho_r2 = 0.75`,
`rho_d = 0.95`) and the manuscript distance cap from central config. Verify the
package helper's exact interpretation. Save the actual per-run `r2min` and `dmax`
values.

Keep hypothesis-level performance in a separate diagnostic table. Never combine
unique/deduplicated TP regions with undeduplicated false hypotheses and label that
hybrid ratio “precision.”

The Manhattan figure subtitles must report the number of TP and FP regions, not
the number of coloured member markers. Member-marker counts describe plotting
coverage and can be given separately, but they do not demonstrate the reduction
in false candidate calls.

For the diagnostic that asks how much is gained simply by removing unclustered
markers, do not reuse the current `*_snp_clustered` definition, which turns a
significant marker into discovery of its whole enclosing unit. Instead define an
optional `emmax_snp_nonsingleton` analysis that:

1. uses the unrestricted marker p-values and genome-wide marker BH correction;
2. excludes significant markers whose phenotype-blind Stage-1 cluster size is 1;
3. still scores every retained significant marker as a marker hypothesis.

This isolates singleton exclusion without granting a marker truth credit from
neighbouring markers.

## Pooling and uncertainty

Point estimates must be ratios of pooled counts, never means of per-run precision
or recall:

```text
precision = sum(TP) / (sum(TP) + sum(FP))
recall    = recovered unique QTN / detectable unique QTN
```

The ten environmental continuations for a map share its burn-in. They are not ten
independent burn-in replicates. Do not use the across-environment standard error
currently implemented in `R/05_pool.R`.

For the primary confidence interval, cluster-bootstrap the ten map/burn-in IDs and
retain all ten environmental continuations belonging to each sampled map. Keep
all methods paired within a resample. For BGS versus no-BGS contrasts, retain the
matched map, environment and seed identifiers in the same resample.

As a sensitivity analysis, use a crossed two-way bootstrap that independently
resamples the ten map IDs and ten environmental-surface IDs, then takes their
Cartesian product while retaining method and BGS pairing. This assesses whether
conclusions change when the ten environmental surfaces are also treated as a
sampled factor. Label the map-cluster bootstrap as conditional on the observed
environmental surfaces and the two-way result accordingly.

Use paired bootstrap contrasts as the main comparison between each reduced method
and its unrestricted marker-wise engine. Report differences in precision and
recall with percentile confidence intervals, not significance tests on 100
pseudo-independent map--environment cells.

## Sensitivity analyses worth retaining

Keep sensitivity work small and tied to a possible criticism:

1. Rescore truth using a small prespecified grid for the additive-variance-share
   threshold, including the primary 5% value. This is cheap because association
   statistics do not change.
2. Report conditional and unconditional recall with respect to Stage-1 coverage.
3. If time permits, compare primary map-cluster intervals with crossed map/environment
   intervals.

Do not add a large Stage-1 rho, unit-size, MAF or association-threshold grid unless
the main result depends on that choice. The empirical analysis already provides
the main Stage-1 parameter sensitivity check.

## Final outputs

Generate these from code, with no hand-copied numbers:

1. `results/simulation_performance.tsv` and `.rds`: one row per
   cell x BGS treatment x method, containing pooled Stage-2 region counts,
   region precision and recall, coverage, inferential test counts, and confidence
   intervals.
2. `results/simulation_method_contrasts.tsv`: paired differences from the
   unrestricted marker method.
3. `results/simulation_qc.tsv`: parser checks, input counts, and the compact BGS
   validation.
4. `results/simulation_truth_sensitivity.tsv`.
5. `figures/figure_simulation_performance.pdf`: the single main simulation figure.
   Prefer two aligned panels showing the paired change in Stage-2 region precision
   and recall from unrestricted marker-wise EMMAX for Stage-1 Simes and consensus
   across all seven cells, with BGS treatment distinguishable but not visually
   dominant. All methods must first be converted to comparable reported regions.
6. `figures/figureS_simulation_absolute_performance.pdf`: absolute precision and
   recall for every primary method and cell.
7. `figures/figureS_simulation_bgs_qc.pdf` or one compact supplementary table.
8. `figures/figureS_simulation_truth_sensitivity.pdf`, only if it adds information
   beyond its table.
9. `tables/tableS_simulation_performance.tex`: full exact results.
10. `values_simulation.tex`: manuscript macros for every value cited in prose.
11. `results/run_manifest.tsv` and `results/validation_report.txt`.

Do not select only high-dispersal cells for the main figure after seeing that the
other cells perform poorly. Show all designed cells. If some regimes fail, that
is part of the method's operating range and belongs in the result.

The main figure should make the method comparison primary. BGS is a robustness
factor, not the organising narrative. Avoid F-beta summaries; show precision and
recall directly.

## Acceptance criteria for the full run

- Correct one-based NEMO locus mapping is demonstrated by assertions and a saved
  QA table.
- No output produced by the old parser is reused.
- The manifest contains 1,400 expected combinations and 1,400 successful parse,
  LD, EMMAX and score records; LFMM completeness is reported separately.
- All methods within a combination use identical analysed individuals, marker
  order, phenotype and GRM basis.
- Stage 2 cannot change any tested-hypothesis p-value or BH decision; it changes
  only the number and bounds of reported candidate regions.
- Primary TP/FP counts are counts of comparable Stage-2 reported regions for every
  method, including the unrestricted marker-wise comparator.
- Region truth status is derived only from constituent discovered-cluster members,
  not unrelated markers that happen to fall inside the interval bounds.
- Pooled metrics are recomputed from counts.
- Uncertainty respects shared burn-ins and paired comparisons.
- No chromosome is labelled neutral merely from its chromosome number.
- Figures and tables read only the final summary objects and fail on incomplete
  data.
- A rerun with unchanged inputs and configuration is deterministic.
- The final README contains the exact execution order, expected run counts,
  approximate resource requirements, and the location of all manuscript outputs.

## Stop points requiring a decision from Petri

Stop and report rather than silently choosing if:

- the corrected one-based join does not place `quant.101` on chromosome 2;
- the native output/reference map is not one-to-one within locus type;
- the actual sample count is not 320 before or 160 after the intended subsampling;
- MAF 0.10 removes essentially all QTN in a substantial part of the design;
- the LFMM structure diagnostic materially contradicts K=5;
- fewer than 1,400 raw combinations are available;
- primary conclusions reverse between the map-cluster and crossed bootstrap;
- the definition of a post hoc reported region is needed to support a main-text
  performance claim.

At each stop point, provide the smallest diagnostic table needed for a decision;
do not launch a new exploratory branch automatically.
